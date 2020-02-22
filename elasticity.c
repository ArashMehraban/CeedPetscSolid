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
  PetscInt       ierr;
  MPI_Comm       comm;
  AppCtx         appCtx;                 // contains problem options
  Physics        phys;                   // contains physical constants
  Units          units;                  // contains units scaling
  DM             dmOrig;                 // distributed DM to clone
  DM             *levelDMs;
  PetscInt       ncompu = 3;             // 3 dofs in 3D
  PetscInt       numLevels = 1, fineLevel = 0;
  Vec            U, *Ug, *Uloc;
  Vec            R, Rloc, F, Floc;       // loc: local, R: residual, F: forcing
  PetscInt       *Ugsz, *Ulsz, *Ulocsz;  // g: global, sz: size
  UserMult       resCtx, *jacobCtx;
  UserMultProlongRestr *prolongRestrCtx;
  FormJacobCtx   formJacobCtx;
  Mat            *jacobMat, *prolongRestrMat;
  PCMGCycleType  pcmgCycleType = PC_MG_CYCLE_V;
  //Ceed constituents
  Ceed           ceed;
  CeedData       *ceedData;
  CeedQFunction  qfRestrict, qfProlong;
  SNES           snes;


  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if (ierr)
    return ierr;

  // ---------------------------------------------------------------------------
  // Process command line options
  // ---------------------------------------------------------------------------
  comm = PETSC_COMM_WORLD;

  // -- Set mesh-file, polynomial degree, problem type
  ierr = processCommandLineOptions(comm, &appCtx); CHKERRQ(ierr);
  numLevels = appCtx.numLevels;
  fineLevel = numLevels - 1;

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
  ierr = PetscMalloc1(numLevels, &levelDMs); CHKERRQ(ierr);
  for (int level = 0; level < numLevels; level++) {
    ierr = DMClone(dmOrig, &levelDMs[level]); CHKERRQ(ierr);
    ierr = SetupDMByDegree(levelDMs[level], appCtx, appCtx.levelDegrees[level],
                           ncompu); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Setup solution and work vectors
  // ---------------------------------------------------------------------------
  // Allocate arrays
  ierr = PetscMalloc1(numLevels, &Ug); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &Uloc); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &Ugsz); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &Ulsz); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &Ulocsz); CHKERRQ(ierr);

  // -- Setup solution vectors for each level
  for (int level = 0; level < numLevels; level++) {
    // -- Create global unknown vector U
    ierr = DMCreateGlobalVector(levelDMs[level], &Ug[level]); CHKERRQ(ierr);
    ierr = VecGetSize(Ug[level], &Ugsz[level]); CHKERRQ(ierr);
    ierr = VecGetLocalSize(Ug[level], &Ulsz[level]); CHKERRQ(ierr); // For matShell

    // -- Create local unknown vector Uloc
    ierr = DMCreateLocalVector(levelDMs[level], &Uloc[level]); CHKERRQ(ierr);
    ierr = VecGetSize(Uloc[level], &Ulocsz[level]); CHKERRQ(ierr); // For libCeed
  }

  // -- Create residual and forcing vectors
  ierr = VecDuplicate(Ug[fineLevel], &U); CHKERRQ(ierr);
  ierr = VecDuplicate(Ug[fineLevel], &R); CHKERRQ(ierr);
  ierr = VecDuplicate(Ug[fineLevel], &F); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc[fineLevel], &Rloc); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc[fineLevel], &Floc); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Set up libCEED
  // ---------------------------------------------------------------------------
  CeedInit(appCtx.ceedresource, &ceed);

  // -- Create libCEED local forcing vector
  CeedVector forceCeed;
  CeedScalar *f;
  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecGetArray(Floc, &f); CHKERRQ(ierr);
    CeedVectorCreate(ceed, Ulocsz[fineLevel], &forceCeed);
    CeedVectorSetArray(forceCeed, CEED_MEM_HOST, CEED_USE_POINTER, f);
  }

  // -- Restriction and prolongation QFunction
  if (appCtx.multigridChoice != MULTIGRID_NONE) {
    CeedQFunctionCreateIdentity(ceed, ncompu, CEED_EVAL_NONE, CEED_EVAL_INTERP,
                                &qfRestrict);
    CeedQFunctionCreateIdentity(ceed, ncompu, CEED_EVAL_INTERP, CEED_EVAL_NONE,
                                &qfProlong);
  }

  // -- Setup libCEED objects
  ierr = PetscMalloc1(numLevels, &ceedData); CHKERRQ(ierr);
  // ---- Setup residual evaluator and geometric information
  ierr = PetscCalloc1(1, &ceedData[fineLevel]); CHKERRQ(ierr);
  ierr = SetupLibceedFineLevel(levelDMs[fineLevel], ceed, appCtx,
                               phys, ceedData, fineLevel, ncompu, 
                               Ugsz[fineLevel], Ulocsz[fineLevel], forceCeed,
                               qfRestrict, qfProlong);
  CHKERRQ(ierr);
  // ---- Setup Jacobian evaluator and prolongtion/restriction
  for (int level = 0; level < numLevels; level++) {
    if (level != fineLevel) {
      ierr = PetscCalloc1(1, &ceedData[level]); CHKERRQ(ierr);
    }
    ierr = SetupLibceedLevel(levelDMs[level], ceed, appCtx, phys,
                             ceedData,  level, ncompu, Ugsz[level],
                             Ulocsz[level], forceCeed, qfRestrict,
                             qfProlong); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Setup global forcing vector
  // ---------------------------------------------------------------------------
  ierr = VecZeroEntries(F); CHKERRQ(ierr);

  if (appCtx.forcingChoice != FORCE_NONE) {
    ierr = VecRestoreArray(Floc, &f); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(levelDMs[fineLevel], Floc, ADD_VALUES,
                                F); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(levelDMs[fineLevel], Floc, ADD_VALUES,
                              F); CHKERRQ(ierr);
    CeedVectorDestroy(&forceCeed);
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
                       Ugsz[fineLevel]/ncompu, Ulsz[fineLevel]/ncompu, ncompu,
                       multigridTypesForDisp[appCtx.multigridChoice],
                       numLevels); CHKERRQ(ierr);

    if (appCtx.multigridChoice != MULTIGRID_NONE) {
      for (int i = 0; i < 2; i++) {  
        CeedInt level = i ? fineLevel : 0;
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

  // ---------------------------------------------------------------------------
  // Setup SNES
  // ---------------------------------------------------------------------------
  ierr = SNESCreate(comm, &snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes, levelDMs[fineLevel]); CHKERRQ(ierr);

  // -- Jacobian evaluators
  ierr = PetscMalloc1(numLevels, &jacobCtx); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &jacobMat); CHKERRQ(ierr);
  for (int level = 0; level < numLevels; level++) {
    // -- Jacobian context for level
    ierr = PetscMalloc1(1, &jacobCtx[level]); CHKERRQ(ierr);
    ierr = SetupJacobianCtx(comm, levelDMs[level], Ug[level], Uloc[level],
                            ceedData[level], ceed, jacobCtx[level]);
    CHKERRQ(ierr);

    // -- Form Action of Jacobian on delta_u
    ierr = MatCreateShell(comm, Ulsz[level], Ulsz[level], Ugsz[level],
                          Ugsz[level], jacobCtx[level], &jacobMat[level]);
    CHKERRQ(ierr);
    ierr = MatShellSetOperation(jacobMat[level], MATOP_MULT,
                                (void (*)(void))ApplyJacobian_Ceed);
    CHKERRQ(ierr);
    ierr = MatShellSetOperation(jacobMat[level], MATOP_GET_DIAGONAL,
                                (void(*)(void))GetDiag_Ceed);

  }
  // FormJacobian restricts the gradient of the state vector and assembles
  //   the diagonal of each Jacobian
  ierr = PetscMalloc1(1, &formJacobCtx); CHKERRQ(ierr);
  formJacobCtx->jacobCtx = jacobCtx;
  formJacobCtx->numLevels = numLevels;
  ierr = SNESSetJacobian(snes, jacobMat[fineLevel], NULL,
                         FormJacobian, formJacobCtx); CHKERRQ(ierr);

  // -- Residual evaluation function
  ierr = PetscMalloc1(1, &resCtx); CHKERRQ(ierr);
  ierr = PetscMemcpy(resCtx, jacobCtx[fineLevel],
                     sizeof(*jacobCtx[fineLevel])); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, R, FormResidual_Ceed, resCtx); CHKERRQ(ierr);

  // -- Prolongation/Restriction evaluation
  ierr = PetscMalloc1(numLevels, &prolongRestrCtx); CHKERRQ(ierr);
  ierr = PetscMalloc1(numLevels, &prolongRestrMat); CHKERRQ(ierr);
  for (int level = 1; level < numLevels; level++) {
    // -- Prolongation/restriction context for level
    ierr = PetscMalloc1(1, &prolongRestrCtx[level]); CHKERRQ(ierr);
    ierr = SetupProlongRestrictCtx(comm, levelDMs[level-1], levelDMs[level],
                                   Ug[level], Uloc[level-1], Uloc[level],
                                   ceedData[level-1], ceedData[level], ceed,
                                   prolongRestrCtx[level]); CHKERRQ(ierr);

    // -- Form Action of Jacobian on delta_u
    ierr = MatCreateShell(comm, Ulsz[level], Ulsz[level-1], Ugsz[level],
                          Ugsz[level-1], prolongRestrCtx[level],
                          &prolongRestrMat[level]); CHKERRQ(ierr);
    // Note: In PCMG, restriction is the transpose of prolongation
    ierr = MatShellSetOperation(prolongRestrMat[level], MATOP_MULT,
                                (void (*)(void))Prolong_Ceed);
    ierr = MatShellSetOperation(prolongRestrMat[level], MATOP_MULT_TRANSPOSE,
                                (void (*)(void))Restrict_Ceed);
    CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Setup KSP
  // ---------------------------------------------------------------------------
  {
    PC pc;
    KSP ksp;

    // -- KSP
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetNormType(ksp, KSP_NORM_NATURAL); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT,
                            PETSC_DEFAULT); CHKERRQ(ierr);

    // -- Preconditioning
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);

    if (appCtx.multigridChoice == MULTIGRID_NONE) {
      // ---- No Multigrid
      ierr = PCSetType(pc, PCNONE); CHKERRQ(ierr);
    } else {
      // ---- PCMG
      ierr = PCSetType(pc, PCMG); CHKERRQ(ierr);

      // ------ PCMG levels
      ierr = PCMGSetLevels(pc, numLevels, NULL); CHKERRQ(ierr);
      for (int level = 0; level < numLevels; level++) {
        // -------- Smoother
        KSP smoother;
        PC smoother_pc;
        ierr = PCMGGetSmoother(pc, level, &smoother); CHKERRQ(ierr);
        ierr = KSPSetType(smoother, KSPCHEBYSHEV); CHKERRQ(ierr);
        ierr = KSPChebyshevEstEigSet(smoother, 0, 0.1, 0, 1.1); CHKERRQ(ierr);
        ierr = KSPChebyshevEstEigSetUseNoisy(smoother, PETSC_TRUE); CHKERRQ(ierr);
        ierr = KSPSetOperators(smoother, jacobMat[level], jacobMat[level]);
        CHKERRQ(ierr);
        ierr = KSPGetPC(smoother, &smoother_pc); CHKERRQ(ierr);
        ierr = PCSetType(smoother_pc, PCJACOBI); CHKERRQ(ierr);
        ierr = PCJacobiSetType(smoother_pc, PC_JACOBI_DIAGONAL); CHKERRQ(ierr);

      // -------- Work vector
      if (level != fineLevel) {
        ierr = PCMGSetX(pc, level, Ug[level]); CHKERRQ(ierr);
      }

      // -------- Level prolongation operator
      if (level > 0) {
        ierr = PCMGSetInterpolation(pc, level, prolongRestrMat[level]);
        CHKERRQ(ierr);
      }

      // -------- Coarse solve
      KSP coarse;
      PC coarse_pc;
      ierr = PCMGGetCoarseSolve(pc, &coarse); CHKERRQ(ierr);
      ierr = KSPSetType(coarse, KSPCG); CHKERRQ(ierr);
      ierr = KSPSetOperators(coarse, jacobMat[0], jacobMat[0]); CHKERRQ(ierr);
      ierr = KSPSetTolerances(coarse, 1e-10, 1e-10, PETSC_DEFAULT,
                              PETSC_DEFAULT); CHKERRQ(ierr);
      ierr = KSPGetPC(coarse, &coarse_pc); CHKERRQ(ierr);
      ierr = PCSetType(coarse_pc, PCJACOBI); CHKERRQ(ierr);
      ierr = PCJacobiSetType(coarse_pc, PC_JACOBI_DIAGONAL); CHKERRQ(ierr);
    }

    // -------- PCMG options
    ierr = PCMGSetType(pc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
    ierr = PCMGSetNumberSmooth(pc, 3); CHKERRQ(ierr);
    ierr = PCMGSetCycleType(pc, pcmgCycleType); CHKERRQ(ierr);
  }

    ierr = KSPSetFromOptions(ksp);
  }
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Set initial guess
  // ---------------------------------------------------------------------------
  ierr = VecSet(U, 1.0); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Solve SNES
  // ---------------------------------------------------------------------------
  ierr = SNESSolve(snes, F, U); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // Output summary
  // ---------------------------------------------------------------------------
  if (!appCtx.testMode) {
    // -- SNES
    SNESType snesType;
    SNESConvergedReason reason;
    PetscInt its;
    PetscReal rnorm;
    ierr = SNESGetType(snes, &snesType); CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes, &reason); CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes, &its); CHKERRQ(ierr);
    ierr = SNESGetFunctionNorm(snes, &rnorm); CHKERRQ(ierr);
    ierr = PetscPrintf(comm,
                       "  SNES:\n"
                       "    SNES Type                          : %s\n"
                       "    SNES Convergence                   : %s\n"
                       "    Total SNES Iterations              : %D\n"
                       "    Final rnorm                        : %e\n",
                       snesType, SNESConvergedReasons[reason], its,
                       (double)rnorm); CHKERRQ(ierr);
    // -- KSP
    KSP ksp;
    KSPType kspType;
    ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
    ierr = KSPGetType(ksp, &kspType); CHKERRQ(ierr);
    ierr = PetscPrintf(comm,
                       "  KSP:\n"
                       "    KSP Type                           : %s\n",
                       kspType); CHKERRQ(ierr);

    // -- PC
    if (appCtx.multigridChoice != MULTIGRID_NONE) {
      PC pc;
      PCMGType pcmgType;
      ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
      ierr = PCMGGetType(pc, &pcmgType); CHKERRQ(ierr);
      ierr = PetscPrintf(comm,
                         "  PCMG:\n"
                         "    PCMG Type                          : %s\n"
                         "    PCMG Cycle Type                    : %s\n",
                         PCMGTypes[pcmgType],
                         PCMGCycleTypes[pcmgCycleType]); CHKERRQ(ierr);
    }
  }

  // ---------------------------------------------------------------------------
  // Compute error
  // ---------------------------------------------------------------------------
  if (appCtx.forcingChoice == FORCE_MMS) {
    CeedScalar l2error = 1., l2Unorm = 1.;
    const CeedScalar *truearray;
    Vec errorVec, trueVec;

    // -- Work vectors
    ierr = VecDuplicate(U, &errorVec); CHKERRQ(ierr);
    ierr = VecSet(errorVec, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(U, &trueVec); CHKERRQ(ierr);
    ierr = VecSet(trueVec, 0.0); CHKERRQ(ierr);

    // -- Assebmle global true soltion vector
    CeedVectorGetArrayRead(ceedData[fineLevel]->truesoln, CEED_MEM_HOST,
                           &truearray);
    ierr = VecPlaceArray(resCtx->Yloc, truearray); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(resCtx->dm, resCtx->Yloc, INSERT_VALUES,
                                trueVec); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(resCtx->dm, resCtx->Yloc, INSERT_VALUES, trueVec);
    CHKERRQ(ierr);
    ierr = VecResetArray(resCtx->Yloc); CHKERRQ(ierr);
    CeedVectorRestoreArrayRead(ceedData[fineLevel]->truesoln, &truearray);

    // -- Compute l2 error
    ierr = VecWAXPY(errorVec, -1.0, U, trueVec); CHKERRQ(ierr);
    ierr = VecNorm(errorVec, NORM_2, &l2error); CHKERRQ(ierr);
    ierr = VecNorm(U, NORM_2, &l2Unorm); CHKERRQ(ierr);
    l2error /= l2Unorm;

    // -- Output
    if (!appCtx.testMode || l2error > 0.168) {
      ierr = PetscPrintf(comm,
                         "  Performance:\n"
                         "    L2 Error                         :  %f\n",
                         l2error); CHKERRQ(ierr);
    }

    // -- Cleanup
    ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
    ierr = VecDestroy(&trueVec); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // Free objects
  // ---------------------------------------------------------------------------
  // Data in arrays per level
  for (int level = 0; level < numLevels; level++) {
    // Vectors
    ierr = VecDestroy(&Ug[level]); CHKERRQ(ierr);
    ierr = VecDestroy(&Uloc[level]); CHKERRQ(ierr);

    // Jacobian matrix and data
    ierr = VecDestroy(&jacobCtx[level]->Yloc); CHKERRQ(ierr);
    ierr = VecDestroy(&jacobCtx[level]->diagVec); CHKERRQ(ierr);
    ierr = MatDestroy(&jacobMat[level]); CHKERRQ(ierr);
    ierr = PetscFree(jacobCtx[level]); CHKERRQ(ierr);

    // Prolongation/Restriction matrix and data
    if (level > 0) {
      ierr = VecDestroy(&prolongRestrCtx[level]->multVec); CHKERRQ(ierr);
      ierr = PetscFree(prolongRestrCtx[level]); CHKERRQ(ierr);
      ierr = MatDestroy(&prolongRestrMat[level]); CHKERRQ(ierr);
    }

    // DM
    ierr = DMDestroy(&levelDMs[level]); CHKERRQ(ierr);

    // libCEED objects
    ierr = CeedDataDestroy(level, ceedData[level]); CHKERRQ(ierr);
  }

  // Arrays
  ierr = PetscFree(Ug); CHKERRQ(ierr);
  ierr = PetscFree(Uloc); CHKERRQ(ierr);
  ierr = PetscFree(Ugsz); CHKERRQ(ierr);
  ierr = PetscFree(Ulsz); CHKERRQ(ierr);
  ierr = PetscFree(Ulocsz); CHKERRQ(ierr);
  ierr = PetscFree(jacobCtx); CHKERRQ(ierr);
  ierr = PetscFree(formJacobCtx); CHKERRQ(ierr);
  ierr = PetscFree(jacobMat); CHKERRQ(ierr);
  ierr = PetscFree(prolongRestrCtx); CHKERRQ(ierr);
  ierr = PetscFree(prolongRestrMat); CHKERRQ(ierr);
  ierr = PetscFree(appCtx.levelDegrees); CHKERRQ(ierr);
  ierr = PetscFree(ceedData); CHKERRQ(ierr);

  // libCEED QFunctions
  CeedQFunctionDestroy(&qfRestrict);
  CeedQFunctionDestroy(&qfProlong);
  CeedDestroy(&ceed);

  // PETSc Objects
  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = VecDestroy(&R); CHKERRQ(ierr);
  ierr = VecDestroy(&Rloc); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&Floc); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = DMDestroy(&dmOrig); CHKERRQ(ierr);
  ierr = PetscFree(levelDMs); CHKERRQ(ierr);

  // Structs
  ierr = PetscFree(resCtx); CHKERRQ(ierr);
  ierr = PetscFree(phys); CHKERRQ(ierr);
  ierr = PetscFree(units); CHKERRQ(ierr);

  return PetscFinalize();
}
