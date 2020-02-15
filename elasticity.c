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

static PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat Jpre,void *ctx) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  // Nothing to do for J, which is of type MATSHELL and uses information from Newton.
  // Jpre might be AIJ (e.g., when using coloring), so we need to assemble it.
  ierr = MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  PetscInt    ierr;
  MPI_Comm    comm;
  AppCtx      appCtx; //contains polinomial basis degree, problem choice & mesh filename
  Physics     phys;   // contains physical constants
  Units       units;  // contains units scaling
  DM          dm;
  PetscInt    ncompu = 3;                 // 3 dofs in 3D
  Vec         U, Uloc, R, Rloc, F, Floc;  // loc: Local R:Residual
  PetscInt    Ugsz, Ulsz, Ulocsz;         // g:global, sz: size
  UserMult    resCtx, jacobCtx;;
  Mat         mat;
  //Ceed constituents
  Ceed        ceed;
  CeedData    ceeddata;
  SNES        snes;


  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if(ierr)
    return ierr;

  // Process command line options
  comm = PETSC_COMM_WORLD;
  // -- Set mesh-file, polynomial degree, problem type
  ierr = processCommandLineOptions(comm, &appCtx); CHKERRQ(ierr);
  // -- Set Poison's ratio, Young's Modulus
  ierr = PetscMalloc1(1, &phys); CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &units); CHKERRQ(ierr);
  ierr = processPhysics(comm, phys, units); CHKERRQ(ierr);

  // Setup DM
  // -- Create distributed DM from mesh file (interpolate if polynomial degree > 1)
  ierr = createDistributedDM(comm, &appCtx, &dm); CHKERRQ(ierr);
  // -- Setup DM by polinomial degree
  ierr = SetupDMByDegree(dm, &appCtx, ncompu);

  // Setup solution and work vectors
  // -- Create global unknown vector U, forcing vector F, and residual vector R
  ierr = DMCreateGlobalVector(dm, &U); CHKERRQ(ierr);
  ierr = VecGetSize(U, &Ugsz); CHKERRQ(ierr);
  ierr = VecGetLocalSize(U, &Ulsz); CHKERRQ(ierr); //For matShell
  ierr = VecDuplicate(U, &R); CHKERRQ(ierr);
  ierr = VecDuplicate(U, &F); CHKERRQ(ierr);
  // -- Create local Vectors
  ierr = DMCreateLocalVector(dm, &Uloc); CHKERRQ(ierr);
  ierr = VecGetSize(Uloc, &Ulocsz); CHKERRQ(ierr); //For libCeed
  ierr = VecZeroEntries(Uloc); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc, &Rloc); CHKERRQ(ierr);

  // Set up libCEED
  CeedInit(appCtx.ceedresource, &ceed);
  ierr = PetscCalloc1(1, &ceeddata); CHKERRQ(ierr);
  // -- Create local forcing vector
  CeedVector forceceed;
  CeedScalar *f;
  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecDuplicate(Uloc, &Floc); CHKERRQ(ierr);
    CeedVectorCreate(ceed, Ulocsz, &forceceed);
    ierr = VecGetArray(Floc, &f); CHKERRQ(ierr);
    CeedVectorSetArray(forceceed, CEED_MEM_HOST, CEED_USE_POINTER, f);
  }
  // -- libCEED objects setup
  ierr = SetupLibceedByDegree(dm, ceed, &appCtx, phys, ceeddata, ncompu, Ugsz,
                              Ulocsz, forceceed); CHKERRQ(ierr);
  // -- Setup global forcing vector
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecRestoreArray(Floc, &f); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(dm, Floc, ADD_VALUES, F); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(dm, Floc, ADD_VALUES, F); CHKERRQ(ierr);
    CeedVectorDestroy(&forceceed);
    ierr = VecDestroy(&Floc); CHKERRQ(ierr);
  }

  // Print problem summary
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
                       "    DoF per node                       : %D\n",
                       usedresource, problemTypesForDisp[appCtx.problemChoice],
                       forcingTypesForDisp[appCtx.forcingChoice],
                       boundaryTypesForDisp[appCtx.boundaryChoice],
                       appCtx.meshFile ? appCtx.meshFile : "Box Mesh",
                       appCtx.degree + 1, appCtx.degree + 1, Ugsz/ncompu,
                       Ulsz/ncompu, ncompu); CHKERRQ(ierr);
  }

  // Setup SNES
  ierr = SNESCreate(comm, &snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes, dm); CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &resCtx); CHKERRQ(ierr);
  // -- Jacobian context
  ierr = PetscMalloc1(1, &jacobCtx); CHKERRQ(ierr);
  ierr =  CreateMatrixFreeCtx(comm, dm, Uloc, ceeddata, ceed, resCtx, jacobCtx);
  CHKERRQ(ierr);
  // -- Function that computes the residual
  ierr = SNESSetFunction(snes, R, FormResidual_Ceed, resCtx); CHKERRQ(ierr);
  // -- Form Action of Jacobian on delta_u
  ierr = MatCreateShell(comm, Ulsz, Ulsz, Ugsz, Ugsz, jacobCtx, &mat);
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
  const CeedScalar *CeedTruSln_as_initial;
  Vec VecTruSln_as_initial;
  ierr = VecDuplicate(U, &VecTruSln_as_initial); CHKERRQ(ierr);
  //VecView(VecTruSln_as_initial,PETSC_VIEWER_STDOUT_WORLD);
  CeedVectorGetArrayRead(ceeddata->truesoln, CEED_MEM_HOST, &CeedTruSln_as_initial);
  ierr = VecPlaceArray(resCtx->Yloc, CeedTruSln_as_initial); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(resCtx->dm, resCtx->Yloc, INSERT_VALUES, VecTruSln_as_initial);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(resCtx->dm, resCtx->Yloc, INSERT_VALUES, VecTruSln_as_initial);CHKERRQ(ierr);
  //VecView(VecTruSln_as_initial,PETSC_VIEWER_STDOUT_WORLD);

  //ierr = VecSet(U, 1.0); CHKERRQ(ierr);
  // Solve SNES
  //ierr = SNESSolve(snes, F, U); CHKERRQ(ierr);
  ierr = VecDuplicate(VecTruSln_as_initial,&U); CHKERRQ(ierr);
  ierr = SNESSolve(snes, F, U); CHKERRQ(ierr);
  VecView(R,PETSC_VIEWER_STDOUT_WORLD);
  //clean up for initial guess as input
  ierr = VecResetArray(resCtx->Yloc); CHKERRQ(ierr);
  CeedVectorRestoreArrayRead(ceeddata->truesoln, &CeedTruSln_as_initial);
  ierr = VecDestroy(&VecTruSln_as_initial); CHKERRQ(ierr);



  // Compute error
  if (appCtx.forcingChoice == FORCE_MMS) {
    CeedScalar l2error, l2Unorm;
    const CeedScalar *truearray;
    Vec errorVec, trueVec;
    ierr = VecDuplicate(U, &errorVec); CHKERRQ(ierr);
    ierr = VecDuplicate(U, &trueVec); CHKERRQ(ierr);

    // Global true soltion vector
    CeedVectorGetArrayRead(ceeddata->truesoln, CEED_MEM_HOST, &truearray);
    ierr = VecPlaceArray(resCtx->Yloc, truearray); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(resCtx->dm, resCtx->Yloc, INSERT_VALUES, trueVec);
    CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(resCtx->dm, resCtx->Yloc, INSERT_VALUES, trueVec);
    CHKERRQ(ierr);
    ierr = VecResetArray(resCtx->Yloc); CHKERRQ(ierr);
    CeedVectorRestoreArrayRead(ceeddata->truesoln, &truearray);

    // Compute error
    ierr = VecWAXPY(errorVec, -1.0, U, trueVec); CHKERRQ(ierr);
    ierr = VecNorm(errorVec, NORM_2, &l2error); CHKERRQ(ierr);
    ierr = VecNorm(U, NORM_2, &l2Unorm); CHKERRQ(ierr);
    l2error /= l2Unorm;

    // Output
    if (!appCtx.testMode || l2error > 0.013) {
      ierr = PetscPrintf(comm, "  L2 Error: %f\n", l2error); CHKERRQ(ierr);
    }

    // Cleanup
    ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
    ierr = VecDestroy(&trueVec); CHKERRQ(ierr);
  }

  // Free objects
  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = VecDestroy(&Uloc); CHKERRQ(ierr);
  ierr = VecDestroy(&R); CHKERRQ(ierr);
  ierr = VecDestroy(&Rloc); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&resCtx->Yloc); CHKERRQ(ierr);
  ierr = MatDestroy(&mat); CHKERRQ(ierr);
  ierr = DMDestroy(&dm); CHKERRQ(ierr);
  ierr = PetscFree(resCtx); CHKERRQ(ierr);
  ierr = PetscFree(jacobCtx); CHKERRQ(ierr);
  ierr = PetscFree(phys); CHKERRQ(ierr);
  ierr = PetscFree(units); CHKERRQ(ierr);
  ierr = CeedDataDestroy(0, ceeddata); CHKERRQ(ierr);
  CeedDestroy(&ceed);
  SNESDestroy(&snes);

  return PetscFinalize();
}
