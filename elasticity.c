//                        libCEED + PETSc Example: Elasticity
//
// This example demonstrates a simple usage of libCEED with PETSc to solve
//   elasticity problems.
//
// The code uses higher level communication protocols in DMPlex.
//
// Build with:
//
//     make elasticity [PETSC_DIR=</path/to/petsc>] [CEED_DIR=</path/to/libceed>]
//
// Sample runs:
//
//     ./elasticity -problem linElas -degree 2 -nu 0.3 -E 1 -forcing mms -boundary mms -mesh ./meshes/cylinder8_672e_4ss_us.exo
//     ./elasticity -problem hyperSS -ceed /cpu/self -degree 2 -nu 0.3 -E 1 -forcing mms -boundary mms -mesh ./meshes/cylinder8_672e_4ss_us.exo
//     ./elasticity -problem hyperFS -ceed /gpu/occa -degree 2 -nu 0.3 -E 1 -forcing mms -boundary mms -mesh ./meshes/cylinder8_672e_4ss_us.exo
//
//TESTARGS -ceed {ceed_resource} -test -degree 2 -nu 0.3 -E 1

/// @file
/// CEED elasticity example using PETSc with DMPlex
const char help[] = "Solve solid Problems with CEED and PETSc DMPlex\n";

#include <stdbool.h>
#include <string.h>
#include <petscksp.h>
#include <petscdmplex.h>
#include <ceed.h>
#include "setup.h"

int main(int argc, char **argv) {
  PetscInt    ierr;
  MPI_Comm    comm;
  AppCtx      appCtx;                 // contains problem options
  Physics     phys;                   // contains physical constants
  Units       units;                  // contains units scaling
  DM          dmOrig;
  DM          *levelDMs;
  PetscInt    ncompu = 3;             // 3 dofs in 3D
  Vec         *U, *Uloc;
  Vec         R, Rloc, F, Floc;       // loc: Local R:Residual
  PetscInt    *Ugsz, *Ulsz, *Ulocsz;  // g:global, sz: size
  UserMult    resCtx, jacobCtx;;
  Mat         mat;
  //Ceed constituents
  Ceed        ceed;
  CeedData    *ceeddata;
  CeedQFunction qf_restrict, qf_prolong;
  SNES        snes;


  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if (ierr)
    return ierr;

  // ---------------------------------------------------------------------------
  // Process command line options
  // ---------------------------------------------------------------------------
  comm = PETSC_COMM_WORLD;

  // -- Set mesh-file, polynomial degree, problem type
  ierr = processCommandLineOptions(comm, &appCtx); CHKERRQ(ierr);

  // -- Set Poison's ratio, Young's Modulus
  ierr = PetscMalloc1(1, &phys); CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &units); CHKERRQ(ierr);
  ierr = processPhysics(comm, phys, units); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Setup DM
  // ---------------------------------------------------------------------------
  // -- Create distributed DM from mesh file
  //      (interpolate if polynomial degree > 1)
  ierr = createDistributedDM(comm, &appCtx, &dmOrig); CHKERRQ(ierr);

  // -- Setup DM by polynomial degree
  ierr = PetscMalloc1(appCtx.numLevels, &levelDMs); CHKERRQ(ierr);
  for (int i = 0; i < appCtx.numLevels; i++) {
    ierr = DMClone(dmOrig, &levelDMs[i]); CHKERRQ(ierr);
    ierr = SetupDMByDegree(levelDMs[i], appCtx, appCtx.levelDegrees[i],
                           ncompu); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Setup solution and work vectors
  // ---------------------------------------------------------------------------
  // Allocate arrays
  ierr = PetscMalloc1(appCtx.numLevels, &U); CHKERRQ(ierr);
  ierr = PetscMalloc1(appCtx.numLevels, &Uloc); CHKERRQ(ierr);
  ierr = PetscMalloc1(appCtx.numLevels, &Ugsz); CHKERRQ(ierr);
  ierr = PetscMalloc1(appCtx.numLevels, &Ulsz); CHKERRQ(ierr);
  ierr = PetscMalloc1(appCtx.numLevels, &Ulocsz); CHKERRQ(ierr);

  // -- Setup solution vectors for each level
  for (int i = 0; i < appCtx.numLevels; i++) {
    // -- Create global unknown vector U
    ierr = DMCreateGlobalVector(levelDMs[i], &U[i]); CHKERRQ(ierr);
    ierr = VecGetSize(U[i], &Ugsz[i]); CHKERRQ(ierr);
    ierr = VecGetLocalSize(U[i], &Ulsz[i]); CHKERRQ(ierr); // For matShell

    // -- Create local unknown vector Uloc
    ierr = DMCreateLocalVector(levelDMs[i], &Uloc[i]); CHKERRQ(ierr);
    ierr = VecGetSize(Uloc[i], &Ulocsz[i]); CHKERRQ(ierr); // For libCeed
  }

  // -- Create residual and forcing vectors
  ierr = VecDuplicate(U[appCtx.numLevels-1], &R); CHKERRQ(ierr);
  ierr = VecDuplicate(U[appCtx.numLevels-1], &F); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc[appCtx.numLevels-1], &Rloc); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc[appCtx.numLevels-1], &Floc); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Set up libCEED
  // ---------------------------------------------------------------------------
  CeedInit(appCtx.ceedresource, &ceed);

  // -- Create libCEED local forcing vector
  CeedVector forceceed;
  CeedScalar *f;
  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecGetArray(Floc, &f); CHKERRQ(ierr);
    CeedVectorCreate(ceed, Ulocsz[appCtx.numLevels-1], &forceceed);
    CeedVectorSetArray(forceceed, CEED_MEM_HOST, CEED_USE_POINTER, f);
  }

  // -- Restriction and prolongation QFunction
  if (appCtx.multigridChoice != MULTIGRID_NONE) {
    CeedQFunctionCreateIdentity(ceed, ncompu, CEED_EVAL_NONE, CEED_EVAL_INTERP,
                                &qf_restrict);
    CeedQFunctionCreateIdentity(ceed, ncompu, CEED_EVAL_INTERP, CEED_EVAL_NONE,
                                &qf_prolong);
  }

  // -- Setup libCEED objects
  ierr = PetscMalloc1(appCtx.numLevels, &ceeddata); CHKERRQ(ierr);
  for (int i = 0; i < appCtx.numLevels; i++) {
    ierr = PetscCalloc1(1, &ceeddata[i]); CHKERRQ(ierr);
    ierr = SetupLibceedByDegree(levelDMs[i], ceed, appCtx, phys, ceeddata, i,
                                ncompu, Ugsz[i], Ulocsz[i], forceceed,
                                qf_restrict, qf_prolong);
    CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Setup global forcing vector
  // ---------------------------------------------------------------------------
  ierr = VecZeroEntries(F); CHKERRQ(ierr);

  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecRestoreArray(Floc, &f); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(levelDMs[appCtx.numLevels-1], Floc, ADD_VALUES,
                                F); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(levelDMs[appCtx.numLevels-1], Floc, ADD_VALUES,
                              F); CHKERRQ(ierr);
    CeedVectorDestroy(&forceceed);
  }

  // ---------------------------------------------------------------------------
  // Print problem summary
  // ---------------------------------------------------------------------------
  if (!appCtx.testMode) {
    const char *usedresource;
    CeedGetResource(ceed, &usedresource);

    ierr = PetscPrintf(comm,
                       "\n-- Elastisticy Example - libCEED + PETSc --\n"
                       "  libCEED:\n"
                       "    libCEED Backend                    : %s\n"
                       "  Problem:\n"
                       "    Problem Name                       : %s\n"
                       "    Forcing Function                   : %s\n"
                       "    Boundary Condition                 : %s\n"
                       "  Mesh:\n"
                       "    File                               : %s\n"
                       "    Number of 1D Basis Nodes (p)       : %d\n"
                       "    Number of 1D Quadrature Points (q) : %d\n"
                       "    Global nodes                       : %D\n"
                       "    Owned nodes                        : %D\n"
                       "    DoF per node                       : %D\n"
                       "  Multigrid:\n"
                       "    Type                               : %s\n"
                       "    Number of Levels                   : %d\n",
                       usedresource, problemTypesForDisp[appCtx.problemChoice],
                       forcingTypesForDisp[appCtx.forcingChoice],
                       boundaryTypesForDisp[appCtx.boundaryChoice],
                       appCtx.meshFile ? appCtx.meshFile : "Box Mesh",
                       appCtx.degree + 1, appCtx.degree + 1,
                       Ugsz[appCtx.numLevels-1]/ncompu,
                       Ulsz[appCtx.numLevels-1]/ncompu, ncompu,
                       multigridTypesForDisp[appCtx.multigridChoice],
                       appCtx.numLevels); CHKERRQ(ierr);

    if (appCtx.multigridChoice != MULTIGRID_NONE) {
      for (int i = 0; i < 2; i++) {  
        CeedInt level = i ? appCtx.numLevels - 1 : 0;
        ierr = PetscPrintf(comm,"    Level %D (%s):\n"
                           "      Number of 1D Basis Nodes (p)     : %d\n"
                           "      Global Nodes                     : %D\n"
                           "      Owned Nodes                      : %D\n",
                           level, i ? "fine" : "coarse",
                           appCtx.levelDegrees[level] + 1,
                           Ugsz[level]/ncompu, Ulsz[level]/ncompu);
                           CHKERRQ(ierr);
      }
    }
  }

  // Setup SNES
  ierr = SNESCreate(comm, &snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes, levelDMs[appCtx.numLevels-1]); CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &resCtx); CHKERRQ(ierr);
  // -- Jacobian context
  ierr = PetscMalloc1(1, &jacobCtx); CHKERRQ(ierr);
  ierr =  CreateMatrixFreeCtx(comm, levelDMs[appCtx.numLevels-1], Uloc[appCtx.numLevels-1], ceeddata[appCtx.numLevels-1], ceed, resCtx, jacobCtx);
  CHKERRQ(ierr);
  // -- Function that computes the residual
  ierr = SNESSetFunction(snes, R, FormResidual_Ceed, resCtx); CHKERRQ(ierr);
  // -- Form Action of Jacobian on delta_u
  ierr = MatCreateShell(comm, Ulsz[appCtx.numLevels-1], Ulsz[appCtx.numLevels-1], Ugsz[appCtx.numLevels-1], Ugsz[appCtx.numLevels-1], jacobCtx, &mat);
  CHKERRQ(ierr);
  ierr = MatShellSetOperation(mat, MATOP_MULT,
                              (void (*)(void))ApplyJacobian_Ceed); CHKERRQ(ierr);
  // pass NULL for Pmat: skips assembly when PCNONE, otherwise use coloring (-snes_fd_color) to
  // assemble a MATAIJ Jacobian for use with any real preconditioner
  ierr = SNESSetJacobian(snes, mat, NULL, FormJacobian, NULL); CHKERRQ(ierr);
  // -- Set inner KSP options
  {
    PC pc;
    KSP ksp;
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCNONE); CHKERRQ(ierr); //For Now No Preconditioner
    ierr = KSPSetFromOptions(ksp);
  }
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  // Set initial Guess
  ierr = VecSet(U[appCtx.numLevels-1], 1.0); CHKERRQ(ierr);

  // Solve SNES
  ierr = SNESSolve(snes, F, U[appCtx.numLevels-1]); CHKERRQ(ierr);

  // Compute error
  if (appCtx.forcingChoice == FORCE_MMS) {
    CeedScalar l2error = 1., l2Unorm = 1.;
    const CeedScalar *truearray;
    Vec errorVec, trueVec;
    ierr = VecDuplicate(U[appCtx.numLevels-1], &errorVec); CHKERRQ(ierr);
    ierr = VecSet(errorVec, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(U[appCtx.numLevels-1], &trueVec); CHKERRQ(ierr);
    ierr = VecSet(trueVec, 0.0); CHKERRQ(ierr);

    // Global true soltion vector
    CeedVectorGetArrayRead(ceeddata[appCtx.numLevels-1]->truesoln, CEED_MEM_HOST, &truearray);
    ierr = VecPlaceArray(resCtx->Yloc, truearray); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(resCtx->dm, resCtx->Yloc, INSERT_VALUES, trueVec);
    CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(resCtx->dm, resCtx->Yloc, INSERT_VALUES, trueVec);
    CHKERRQ(ierr);
    ierr = VecResetArray(resCtx->Yloc); CHKERRQ(ierr);
    CeedVectorRestoreArrayRead(ceeddata[appCtx.numLevels-1]->truesoln, &truearray);

    // Compute error
    ierr = VecWAXPY(errorVec, -1.0, U[appCtx.numLevels-1], trueVec); CHKERRQ(ierr);
    ierr = VecNorm(errorVec, NORM_2, &l2error); CHKERRQ(ierr);
    ierr = VecNorm(U[appCtx.numLevels-1], NORM_2, &l2Unorm); CHKERRQ(ierr);
    l2error /= l2Unorm;

    // Output
    if (!appCtx.testMode || l2error > 0.168) {
      ierr = PetscPrintf(comm, "  L2 Error: %f\n", l2error); CHKERRQ(ierr);
    }

    // Cleanup
    ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
    ierr = VecDestroy(&trueVec); CHKERRQ(ierr);
  }

  // Free objects
  for (int i=0; i < appCtx.numLevels; i++) {
    ierr = VecDestroy(&U[i]); CHKERRQ(ierr);
    ierr = VecDestroy(&Uloc[i]); CHKERRQ(ierr);
    ierr = DMDestroy(&levelDMs[i]); CHKERRQ(ierr);
    ierr = CeedDataDestroy(i, ceeddata[i]); CHKERRQ(ierr);
  }
  ierr = PetscFree(U); CHKERRQ(ierr);
  ierr = PetscFree(Uloc); CHKERRQ(ierr);
  ierr = PetscFree(Ugsz); CHKERRQ(ierr);
  ierr = PetscFree(Ulsz); CHKERRQ(ierr);
  ierr = PetscFree(Ulocsz); CHKERRQ(ierr);
  ierr = DMDestroy(&dmOrig); CHKERRQ(ierr);
  ierr = PetscFree(levelDMs); CHKERRQ(ierr);
  CeedQFunctionDestroy(&qf_restrict);
  CeedQFunctionDestroy(&qf_prolong);
  ierr = PetscFree(ceeddata); CHKERRQ(ierr);


  ierr = VecDestroy(&R); CHKERRQ(ierr);
  ierr = VecDestroy(&Rloc); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&Floc); CHKERRQ(ierr);
  ierr = VecDestroy(&resCtx->Yloc); CHKERRQ(ierr);
  ierr = MatDestroy(&mat); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = PetscFree(resCtx); CHKERRQ(ierr);
  ierr = PetscFree(jacobCtx); CHKERRQ(ierr);
  ierr = PetscFree(phys); CHKERRQ(ierr);
  ierr = PetscFree(units); CHKERRQ(ierr);
  ierr = PetscFree(appCtx.levelDegrees); CHKERRQ(ierr);
  CeedDestroy(&ceed);

  return PetscFinalize();
}
