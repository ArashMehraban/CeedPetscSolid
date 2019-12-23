#ifndef setup_h
#define setup_h

#include <stdbool.h>
#include <string.h>
#include <petsc.h>
#include <petscdmplex.h>
#include <petscfe.h>
#include <ceed.h>
#include "qfunctions/solid/common.h"
#include "qfunctions/solid/linElas.h" //Linear Elasticity
#include "qfunctions/solid/hyperSS.h" //Hyperelasticity Small SmallStrain
#include "qfunctions/solid/hyperFS.h" //Hyperelasticity Small FiniteStrain

// Problem options
typedef enum {  //SmallStrain      FiniteStrain
  ELAS_LIN = 0, ELAS_HYPER_SS = 1, ELAS_HYPER_FS = 2
} problemType;
static const char *const problemTypes[] = {"linElas","hyperSS","hyperFS",
                                            "problemType","ELAS_",0};

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------

typedef struct{
  char          meshFile[PETSC_MAX_PATH_LEN]; // exodusII mesh file
  problemType   problemChoice;
  PetscInt      degree;
}AppCtx;

//Physics for Elasticty and Hyperelasticity probelms
typedef struct{
  PetscScalar   nu;      //Poisson's ratio
  PetscScalar   E;       //Young's Modulus
}Physics;

// Problem specific data
typedef struct {
  CeedInt ncompu, qdatasize;
  // CeedQFunctionUser setupgeo, setuprhs, apply, error;
  // const char *setupgeofname, *setuprhsfname, *applyfname, *errorfname;
  // CeedEvalMode inmode, outmode;
  // CeedQuadMode qmode;
  // PetscBool enforce_bc;
  PetscErrorCode (*bcs_func)(PetscInt, PetscReal, const PetscReal *,
                             PetscInt, PetscScalar *, void *);
}problemData;

problemData problemOptions[3] = {
  [ELAS_LIN] = {
      .ncompu = 3,
      // .qdatasize = 10,
      // .setupgeo = SetupDiffGeo,
      // .apply = Diff3,  //<---- linearElasticity
      // .error = Error3,
      // .setuprhs = SetupDiffRhs,
      // .setupgeofname = SetupDiffGeo_loc,
      // .setuprhsfname = SetupDiffRhs3_loc,
      // .applyfname = Diff_loc,
      // .errorfname = Error3_loc,
      // .inmode = CEED_EVAL_GRAD,
      // .outmode = CEED_EVAL_GRAD,
      // .qmode = CEED_GAUSS,
       .bcs_func = NULL//BcOne // Drichlet of all 1's
  },
  [ELAS_HYPER_SS] = {
      .ncompu = 3,
      // .qdatasize = 10,
      // .setupgeo = SetupDiffGeo,
      // .setuprhs = SetupDiffRhs3,
      // .apply = Diff3,  //<---- hyperelasticitySS
      // .error = Error3,
      // .setupgeofname = SetupDiffGeo_loc,
      // .setuprhsfname = SetupDiffRhs3_loc,
      // .applyfname = Diff_loc,
      // .errorfname = Error3_loc,
      // .inmode = CEED_EVAL_GRAD,
      // .outmode = CEED_EVAL_GRAD,
      // .qmode = CEED_GAUSS,
      .bcs_func = NULL //BCsDiff
  },
  [ELAS_HYPER_FS] = {
      .ncompu = 3,
      // .qdatasize = 10,
      // .setupgeo = SetupDiffGeo,
      // .setuprhs = SetupDiffRhs3,
      // .apply = Diff3,   //<---- hyperelasticityFS
      // .error = Error3,
      // .setupgeofname = SetupDiffGeo_loc,
      // .setuprhsfname = SetupDiffRhs3_loc,
      // .applyfname = Diff_loc,
      // .errorfname = Error3_loc,
      // .inmode = CEED_EVAL_GRAD,
      // .outmode = CEED_EVAL_GRAD,
      // .qmode = CEED_GAUSS,
      .bcs_func = NULL //BCsDiff
    }
};

// -----------------------------------------------------------------------------
// Helper Functions
// -----------------------------------------------------------------------------
static int processCommandLineOptions(MPI_Comm comm, AppCtx *appCtx){

  PetscErrorCode ierr;
  PetscBool meshFileFlag = PETSC_FALSE;
  PetscBool degreeFalg = PETSC_FALSE;
  appCtx->problemChoice = ELAS_LIN; //-problem = Linear Elasticity if not given
  appCtx->degree = 0;

  PetscFunctionBeginUser;

  ierr = PetscOptionsBegin(comm, NULL,
     "Elasticity / Hyperelasticity in PETSc with libCEED", NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-degree", "Polynomial degree of tensor product basis",
                         NULL, appCtx->degree, &appCtx->degree, &degreeFalg);
                         CHKERRQ(ierr);

  ierr = PetscOptionsString("-mesh", "Read mesh from file", NULL,
                           appCtx->meshFile, appCtx->meshFile,
                           sizeof(appCtx->meshFile), &meshFileFlag);
                           CHKERRQ(ierr);
  #if !defined(PETSC_HAVE_EXODUSII)
    SETERRQ(comm, PETSC_ERR_ARG_WRONG,
     "ExodusII support needed. Reconfigure your Arch with --download-exodusii");
  #endif

  ierr = PetscOptionsEnum("-problem",
                          "Solves Elasticity & Hyperelasticity Problems", NULL,
                          problemTypes,(PetscEnum)appCtx->problemChoice,
                          (PetscEnum *)&appCtx->problemChoice,
                           NULL); CHKERRQ(ierr);


  ierr = PetscOptionsEnd(); CHKERRQ(ierr);//End of setting AppCtx
  if(!degreeFalg){
      ierr = PetscPrintf(comm, "-degree option needed\n\n");CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
  }
  if(!meshFileFlag){
      ierr = PetscPrintf(comm, "-mesh option needed (file)\n\n");CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
  }

    PetscFunctionReturn(0);
}

static int processPhysics(MPI_Comm comm, Physics *phys){

    PetscErrorCode ierr;
    PetscBool nuFlag = PETSC_FALSE;
    PetscBool YoungFlag = PETSC_FALSE;
    phys->nu = 0;
    phys->E = 0;

    PetscFunctionBeginUser;

    ierr = PetscOptionsBegin(comm, NULL,
     "Elasticity / Hyperelasticity in PETSc with libCEED", NULL); CHKERRQ(ierr);
    ierr = PetscOptionsScalar("-nu", "Poisson's ratio", NULL, phys->nu,
                              &phys->nu, &nuFlag);CHKERRQ(ierr);
    ierr = PetscOptionsScalar("-E", "Young's Modulus", NULL, phys->E, &phys->E,
                               &YoungFlag);CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);//End of setting Physics
    if(!nuFlag){
        ierr = PetscPrintf(comm, "-nu option needed\n\n");CHKERRQ(ierr);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Poisson's ratio Error!");
    }
    if(!YoungFlag){
        ierr = PetscPrintf(comm, "-E option needed\n\n");CHKERRQ(ierr);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Young's Modulus Error!");
    }
    PetscFunctionReturn(0);
}

static int createDistributedDM(MPI_Comm comm, AppCtx *ctx, DM *dm){

    PetscErrorCode  ierr;
    const char      *filename = ctx->meshFile;
    PetscBool       interpolate = PETSC_FALSE;
    DM              distributedMesh = NULL;
    PetscPartitioner part;

    PetscFunctionBeginUser;
    if(ctx->degree >= 2)
            interpolate = PETSC_TRUE;

    ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm);CHKERRQ(ierr);
    ierr = DMPlexGetPartitioner(*dm, &part); CHKERRQ(ierr);
    ierr = PetscPartitionerSetFromOptions(part); CHKERRQ(ierr);
    ierr = DMPlexDistribute(*dm, 0, NULL, &distributedMesh);CHKERRQ(ierr);
    if (distributedMesh) {
      ierr = DMDestroy(dm);CHKERRQ(ierr);
      *dm  = distributedMesh;
    }
    PetscFunctionReturn(0);
}

// -----------------------------------------------------------------------------
// Problem Option Data
// -----------------------------------------------------------------------------
static int SetupDMByDegree(DM dm, AppCtx *appCtx, PetscInt ncompu){

 PetscErrorCode  ierr;
 PetscInt        dim;
 PetscInt        order = appCtx->degree;
 PetscSpace      P;
 PetscDualSpace  Q;
 DM              K;
 PetscFE         fe;
 PetscInt        quadPointsPerEdge;
 PetscQuadrature q;  //quadrature points
 PetscQuadrature fq; //face quadrature points
 // For Dirichlet (Essential) Boundary
 IS              faceSetIS;        //Index Set for Face Sets
 const char     *name="Face Sets"; //PETSc internal requirement
 PetscInt        numFaceSet;       //Number of FaceSets in faceSetIS
 const PetscInt *faceSetIds;       //id of each FaceSet

 PetscFunctionBeginUser;

  ierr = DMGetDimension(dm, &dim);
 
  //setup FE space (Space P) for tensor polynomials
  ierr = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &P); CHKERRQ(ierr);
  //ierr = PetscObjectSetOptionsPrefix((PetscObject) P, prefix); CHKERRQ(ierr);
  ierr = PetscSpacePolynomialSetTensor(P, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscSpaceSetFromOptions(P); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumComponents(P, ncompu); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumVariables(P, dim); CHKERRQ(ierr);
  ierr = PetscSpaceSetDegree(P, order, order); CHKERRQ(ierr);
  ierr = PetscSpaceSetUp(P); CHKERRQ(ierr);
  //setup FE dual space (Space Q) for tensor polynomials
  ierr = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &Q);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE); CHKERRQ(ierr);
  //ierr = PetscObjectSetOptionsPrefix((PetscObject) Q, prefix); CHKERRQ(ierr);
  ierr = PetscDualSpaceCreateReferenceCell(Q, dim, PETSC_FALSE, &K); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetDM(Q, K); CHKERRQ(ierr);
  ierr = DMDestroy(&K); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetNumComponents(Q, ncompu); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetOrder(Q, order); CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetTensor(Q, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetFromOptions(Q); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetUp(Q); CHKERRQ(ierr);
  /* Create element */
  ierr = PetscFECreate(PetscObjectComm((PetscObject) dm), &fe); CHKERRQ(ierr);
  //ierr = PetscObjectSetOptionsPrefix((PetscObject) fe, prefix); CHKERRQ(ierr);
  ierr = PetscFESetFromOptions(fe); CHKERRQ(ierr);
  ierr = PetscFESetBasisSpace(fe, P); CHKERRQ(ierr);
  ierr = PetscFESetDualSpace(fe, Q); CHKERRQ(ierr);
  ierr = PetscFESetNumComponents(fe, ncompu); CHKERRQ(ierr);
  ierr = PetscFESetUp(fe); CHKERRQ(ierr);
  ierr = PetscSpaceDestroy(&P); CHKERRQ(ierr);
  ierr = PetscDualSpaceDestroy(&Q); CHKERRQ(ierr);
  /* Create quadrature */
  quadPointsPerEdge = PetscMax(order + 1,1);
  ierr = PetscDTGaussTensorQuadrature(dim,   1, quadPointsPerEdge, -1.0, 1.0,
                                        &q); CHKERRQ(ierr);
  ierr = PetscDTGaussTensorQuadrature(dim-1, 1, quadPointsPerEdge, -1.0, 1.0,
                                        &fq); CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe, q); CHKERRQ(ierr);
  ierr = PetscFESetFaceQuadrature(fe, fq); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&q); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&fq); CHKERRQ(ierr);
  // Setup DM
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMAddField(dm, NULL, (PetscObject)fe); CHKERRQ(ierr);
  ierr = DMCreateDS(dm); CHKERRQ(ierr);
  // Add Dirichlet (Essential) boundray
  ierr = DMGetLabelIdIS(dm, name, &faceSetIS);CHKERRQ(ierr);
  ierr = ISGetLocalSize(faceSetIS,&numFaceSet);
  ierr = ISGetIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);
  ierr = DMAddBoundary(dm,DM_BC_ESSENTIAL,"wall","Face Sets",0,0,NULL,
                  (void(*)(void))problemOptions[appCtx->problemChoice].bcs_func,
                  numFaceSet,faceSetIds,NULL);CHKERRQ(ierr);
  ierr = ISRestoreIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);
  ierr = ISDestroy(&faceSetIS);CHKERRQ(ierr);
  ierr = DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}



#endif
