// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#ifndef setup_h
#define setup_h

#include <stdbool.h>
#include <string.h>
#include <petsc.h>
#include <petscdmplex.h>
#include <petscfe.h>
#include <ceed.h>
#include "qfunctions/solid/common.h"
#include "qfunctions/solid/linElas.h"
#include "qfunctions/solid/hyperSS.h" //Hyperelasticity Small SmallStrain
#include "qfunctions/solid/hyperFS.h" //Hyperelasticity Small FiniteStrain


// -----------------------------------------------------------------------------
// PETSc Operator Structs
// -----------------------------------------------------------------------------
//Physics Structs (No pointer to struct necessary as)
typedef struct{
  PetscScalar   nu;      //Poisson's ratio
  PetscScalar   E;       //Young's Modulus
}Physics;

// Data for PETSc Matshell
typedef struct UserO_ *UserO;
struct UserO_ {
  MPI_Comm comm;
  DM dm;
  Vec Xloc, Yloc, diag;
  CeedVector xceed, yceed;
  CeedOperator op;
  Ceed ceed;
};

// Data for PETSc Interp/Restrict Matshells
typedef struct UserIR_ *UserIR;
struct UserIR_ {
  MPI_Comm comm;
  DM dmc, dmf;
  Vec Xloc, Yloc, mult;
  CeedVector ceedvecc, ceedvecf;
  CeedOperator op;
  Ceed ceed;
};

// -----------------------------------------------------------------------------
// libCEED Data Struct
// -----------------------------------------------------------------------------

// libCEED data struct for level
typedef struct CeedData_ *CeedData;
struct CeedData_ {
  Ceed ceed;
  CeedBasis basisx, basisu, basisctof;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi, Erestrictui,
                      Erestrictqdi;
  CeedQFunction qf_apply;
  CeedOperator op_apply, op_restrict, op_interp;
  CeedVector qdata, xceed, yceed;
};

// -----------------------------------------------------------------------------
// Command Line Options
// -----------------------------------------------------------------------------

// Coarsening options
typedef enum {
  COARSEN_UNIFORM = 0, COARSEN_LOGARITHMIC = 1
} coarsenType;
static const char *const coarsenTypes [] = {"uniform","logarithmic",
                                            "coarsenType","COARSEN",0
                                           };

// -----------------------------------------------------------------------------
// Boundary Conditions
// -----------------------------------------------------------------------------
// test Boundary condition with all components being 1
PetscErrorCode BcOne(PetscInt dim, PetscReal time, const PetscReal x[],
                       PetscInt ncompu, PetscScalar *u, void *ctx){
    PetscFunctionBeginUser;

    for (PetscInt i = 0; i < ncompu; i++){
        PetscPrintf(PETSC_COMM_WORLD,"All Dirichlet boundary values = 1\n\n\n");
        u[i] = 1.0;
    }

    PetscFunctionReturn(0);
}


// Diff boundary condition function
PetscErrorCode BCsDiff(PetscInt dim, PetscReal time, const PetscReal x[],
                       PetscInt ncompu, PetscScalar *u, void *ctx) {
  // *INDENT-OFF*
  #ifndef M_PI
  #define M_PI    3.14159265358979323846
  #endif
  // *INDENT-ON*
  const CeedScalar c[3] = { 0, 1., 2. };
  const CeedScalar k[3] = { 1., 2., 3. };

  PetscFunctionBeginUser;

  for (PetscInt i = 0; i < ncompu; i++)
    u[i] = sin(M_PI*(c[0] + k[0]*x[0])) *
           sin(M_PI*(c[1] + k[1]*x[1])) *
           sin(M_PI*(c[2] + k[2]*x[2]));

  PetscFunctionReturn(0);
}

// Mass boundary condition function
PetscErrorCode BCsMass(PetscInt dim, PetscReal time, const PetscReal x[],
                       PetscInt ncompu, PetscScalar *u, void *ctx) {
  PetscFunctionBeginUser;

  for (PetscInt i = 0; i < ncompu; i++)
    u[i] = PetscSqrtScalar(PetscSqr(x[0]) + PetscSqr(x[1]) +
                           PetscSqr(x[2]));

  PetscFunctionReturn(0);
}

// -----------------------------------------------------------------------------
// Problem Option Data
// -----------------------------------------------------------------------------

// Problem options
typedef enum {  //SmallStrain      FiniteStrain
  ELAS_LIN = 0, ELAS_HYPER_SS = 1, ELAS_HYPER_FS = 2
} problemType;
static const char *const problemTypes[] = {"linElas","hyperSS","hyperFS",
                                            "problemType","ELAS_",0};

// Problem specific data
typedef struct {
  CeedInt ncompu, qdatasize;
  CeedQFunctionUser setupgeo, setuprhs, apply, error;
  const char *setupgeofname, *setuprhsfname, *applyfname, *errorfname;
  CeedEvalMode inmode, outmode;
  CeedQuadMode qmode;
  PetscBool enforce_bc;
  PetscErrorCode (*bcs_func)(PetscInt, PetscReal, const PetscReal *,
                             PetscInt, PetscScalar *, void *);
}problemData;

problemData problemOptions[3] = {
  [ELAS_LIN] = {
      .ncompu = 3,
      .qdatasize = 10,
      .setupgeo = SetupDiffGeo,
      .apply = Diff3,  //<---- linearElasticity
      .error = Error3,
      .setuprhs = SetupDiffRhs,
      .setupgeofname = SetupDiffGeo_loc,
      .setuprhsfname = SetupDiffRhs3_loc,
      .applyfname = Diff_loc,
      .errorfname = Error3_loc,
      .inmode = CEED_EVAL_GRAD,
      .outmode = CEED_EVAL_GRAD,
      .qmode = CEED_GAUSS,
      .bcs_func = BcOne // Drichlet of all 1's
  },
  [ELAS_HYPER_SS] = {
      .ncompu = 3,
      .qdatasize = 10,
      .setupgeo = SetupDiffGeo,
      .setuprhs = SetupDiffRhs3,
      .apply = Diff3,  //<---- hyperelasticitySS
      .error = Error3,
      .setupgeofname = SetupDiffGeo_loc,
      .setuprhsfname = SetupDiffRhs3_loc,
      .applyfname = Diff_loc,
      .errorfname = Error3_loc,
      .inmode = CEED_EVAL_GRAD,
      .outmode = CEED_EVAL_GRAD,
      .qmode = CEED_GAUSS,
      .bcs_func = NULL //BCsDiff
  },
  [ELAS_HYPER_FS] = {
      .ncompu = 3,
      .qdatasize = 10,
      .setupgeo = SetupDiffGeo,
      .setuprhs = SetupDiffRhs3,
      .apply = Diff3,   //<---- hyperelasticityFS
      .error = Error3,
      .setupgeofname = SetupDiffGeo_loc,
      .setuprhsfname = SetupDiffRhs3_loc,
      .applyfname = Diff_loc,
      .errorfname = Error3_loc,
      .inmode = CEED_EVAL_GRAD,
      .outmode = CEED_EVAL_GRAD,
      .qmode = CEED_GAUSS,
      .bcs_func = NULL //BCsDiff
    }
};

// -----------------------------------------------------------------------------
// PETSc FE Setup
// -----------------------------------------------------------------------------

// Create FE by degree
static int PetscFECreateByDegree(DM dm, PetscInt dim, PetscInt Nc,
                                 const char prefix[],PetscInt order,
                                 PetscFE *fem) {
  PetscQuadrature q, fq; //Quadrature , FaceQuadrature
  DM              K;
  PetscSpace      P;
  PetscDualSpace  Q;
  PetscInt        quadPointsPerEdge;
  PetscBool       tensor = PETSC_TRUE;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  /* Create space */
  ierr = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &P); CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) P, prefix); CHKERRQ(ierr);
  ierr = PetscSpacePolynomialSetTensor(P, tensor); CHKERRQ(ierr);
  ierr = PetscSpaceSetFromOptions(P); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumComponents(P, Nc); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumVariables(P, dim); CHKERRQ(ierr);
  ierr = PetscSpaceSetDegree(P, order, order); CHKERRQ(ierr);
  ierr = PetscSpaceSetUp(P); CHKERRQ(ierr);
  ierr = PetscSpacePolynomialGetTensor(P, &tensor); CHKERRQ(ierr);
  /* Create dual space */
  ierr = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &Q);
  CHKERRQ(ierr);
  ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE); CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) Q, prefix); CHKERRQ(ierr);
  ierr = PetscDualSpaceCreateReferenceCell(Q, dim, tensor, &K); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetDM(Q, K); CHKERRQ(ierr);
  ierr = DMDestroy(&K); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetNumComponents(Q, Nc); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetOrder(Q, order); CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetTensor(Q, tensor); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetFromOptions(Q); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetUp(Q); CHKERRQ(ierr);
  /* Create element */
  ierr = PetscFECreate(PetscObjectComm((PetscObject) dm), fem); CHKERRQ(ierr);
  ierr = PetscObjectSetOptionsPrefix((PetscObject) *fem, prefix); CHKERRQ(ierr);
  ierr = PetscFESetFromOptions(*fem); CHKERRQ(ierr);
  ierr = PetscFESetBasisSpace(*fem, P); CHKERRQ(ierr);
  ierr = PetscFESetDualSpace(*fem, Q); CHKERRQ(ierr);
  ierr = PetscFESetNumComponents(*fem, Nc); CHKERRQ(ierr);
  ierr = PetscFESetUp(*fem); CHKERRQ(ierr);
  ierr = PetscSpaceDestroy(&P); CHKERRQ(ierr);
  ierr = PetscDualSpaceDestroy(&Q); CHKERRQ(ierr);
  /* Create quadrature */
  quadPointsPerEdge = PetscMax(order + 1,1);
  ierr = PetscDTGaussTensorQuadrature(dim,   1, quadPointsPerEdge, -1.0, 1.0,
                                      &q); CHKERRQ(ierr);
  ierr = PetscDTGaussTensorQuadrature(dim-1, 1, quadPointsPerEdge, -1.0, 1.0,
                                      &fq); CHKERRQ(ierr);

  ierr = PetscFESetQuadrature(*fem, q); CHKERRQ(ierr);
  ierr = PetscFESetFaceQuadrature(*fem, fq); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&q); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&fq); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// -----------------------------------------------------------------------------
// PETSc Setup for Level
// -----------------------------------------------------------------------------

// This function sets up a DM for a given degree
static int SetupDMByDegree(DM dm, PetscInt degree, PetscInt ncompu,
                           problemType problemChoice) {
  PetscInt ierr, dim;
  PetscFE fe;

  PetscFunctionBeginUser;

  // Setup FE
  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ierr = PetscFECreateByDegree(dm, dim, ncompu, NULL, degree, &fe);
  CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMAddField(dm, NULL, (PetscObject)fe); CHKERRQ(ierr);

  // Setup DM
  ierr = DMCreateDS(dm); CHKERRQ(ierr);
  { // Adding Boundary
      IS              faceSetIS; //IS for Face Sets
      const char     *name="Face Sets"; //PETSc internal requirement
      PetscInt        numFaceSet;  //Number of Face Sets in faceSetIS
      const PetscInt *faceSetIds;  //id of each Face Set
      ierr = DMGetLabelIdIS(dm, name, &faceSetIS);CHKERRQ(ierr);
      ierr = ISGetLocalSize(faceSetIS,&numFaceSet);
      ierr = ISGetIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);
      ierr = DMAddBoundary(dm,DM_BC_ESSENTIAL,"wall","Face Sets",0,0,NULL,
                           (void(*)(void))problemOptions[problemChoice].bcs_func,
                           numFaceSet,faceSetIds,NULL);CHKERRQ(ierr);
      ierr = ISRestoreIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);
      ierr = ISDestroy(&faceSetIS);CHKERRQ(ierr);
  }
  ierr = DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// -----------------------------------------------------------------------------
// libCEED Setup for Level
// -----------------------------------------------------------------------------

// Destroy libCEED operator objects
static PetscErrorCode CeedDataDestroy(CeedInt i, CeedData data) {
  PetscInt ierr;

  CeedVectorDestroy(&data->qdata);
  CeedVectorDestroy(&data->xceed);
  CeedVectorDestroy(&data->yceed);
  CeedBasisDestroy(&data->basisx);
  CeedBasisDestroy(&data->basisu);
  CeedElemRestrictionDestroy(&data->Erestrictu);
  CeedElemRestrictionDestroy(&data->Erestrictx);
  CeedElemRestrictionDestroy(&data->Erestrictui);
  CeedElemRestrictionDestroy(&data->Erestrictxi);
  CeedElemRestrictionDestroy(&data->Erestrictqdi);
  CeedQFunctionDestroy(&data->qf_apply);
  CeedOperatorDestroy(&data->op_apply);
  if (i > 0) {
    CeedOperatorDestroy(&data->op_interp);
    CeedBasisDestroy(&data->basisctof);
    CeedOperatorDestroy(&data->op_restrict);
  }
  ierr = PetscFree(data); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Get CEED restriction data from DMPlex
static int CreateRestrictionPlex(Ceed ceed, CeedInt P, CeedInt ncomp,
                                 CeedElemRestriction *Erestrict, DM dm) {
  PetscInt ierr;
  PetscInt c, cStart, cEnd, nelem, nnodes, *erestrict, eoffset;
  PetscSection section;
  Vec Uloc;

  PetscFunctionBeginUser;

  // Get Nelem
  ierr = DMGetSection(dm, &section); CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart,& cEnd); CHKERRQ(ierr);
  nelem = cEnd - cStart;

  // Get indices
  ierr = PetscMalloc1(nelem*P*P*P, &erestrict); CHKERRQ(ierr);
  for (c=cStart, eoffset=0; c<cEnd; c++) {
    PetscInt numindices, *indices, i;
    ierr = DMPlexGetClosureIndices(dm, section, section, c, &numindices,
                                   &indices, NULL); CHKERRQ(ierr);
    for (i=0; i<numindices; i+=ncomp) {
      for (PetscInt j=0; j<ncomp; j++) {
        if (indices[i+j] != indices[i] + (PetscInt)(copysign(j, indices[i])))
          SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
                   "Cell %D closure indices not interlaced", c);
      }
      // Essential boundary conditions are encoded as -(loc+1)
      PetscInt loc = indices[i] >= 0 ? indices[i] : -(indices[i] + 1);
      erestrict[eoffset++] = loc/ncomp;
    }
    ierr = DMPlexRestoreClosureIndices(dm, section, section, c, &numindices,
                                       &indices, NULL); CHKERRQ(ierr);
  }

  // Setup CEED restriction
  ierr = DMGetLocalVector(dm, &Uloc); CHKERRQ(ierr);
  ierr = VecGetLocalSize(Uloc, &nnodes); CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(dm, &Uloc); CHKERRQ(ierr);
  CeedElemRestrictionCreate(ceed, nelem, P*P*P, nnodes/ncomp, ncomp,
                            CEED_MEM_HOST, CEED_COPY_VALUES, erestrict,
                            Erestrict);
  ierr = PetscFree(erestrict); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Set up libCEED for a given degree
static int SetupLibceedByDegree(DM dm, Ceed ceed, CeedInt degree,
                                PetscInt ncompu, PetscInt gsize,PetscInt xlsize,
                                problemType problemChoice, CeedData data,
                                PetscBool setup_rhs, CeedVector rhsceed,
                                CeedVector *target) {
  int ierr;
  DM dmcoord;
  PetscSection section;
  Vec coords;
  const PetscScalar *coordArray;
  CeedBasis basisx, basisu;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi,
                      Erestrictui, Erestrictqdi;
  CeedQFunction qf_setupgeo, qf_apply;
  CeedOperator op_setupgeo, op_apply;
  CeedVector xcoord, qdata, xceed, yceed;
  CeedInt dim;
  CeedInt qdatasize = problemOptions[problemChoice].qdatasize, ncompx, P, Q,
          cStart, cEnd, nelem;

  // CEED bases
  P = degree + 1;
  Q = P;
  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ncompx = dim;
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu, P, Q,
                                  problemOptions[problemChoice].qmode, &basisu);
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, Q,
                                  problemOptions[problemChoice].qmode, &basisx);

  // CEED restrictions
  ierr = DMGetCoordinateDM(dm, &dmcoord); CHKERRQ(ierr);
  ierr = DMPlexSetClosurePermutationTensor(dmcoord, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);

  CreateRestrictionPlex(ceed, 2, ncompx, &Erestrictx, dmcoord);
  CreateRestrictionPlex(ceed, P, ncompu, &Erestrictu, dm);

  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
  nelem = cEnd - cStart;

  CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, ncompu,
                                    &Erestrictui); CHKERRQ(ierr);
  CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q,
                                    qdatasize, &Erestrictqdi); CHKERRQ(ierr);
  CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, ncompx,
                                    &Erestrictxi); CHKERRQ(ierr);

  // Element coordinates
  ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
  ierr = VecGetArrayRead(coords, &coordArray); CHKERRQ(ierr);
  ierr = DMGetSection(dmcoord, &section); CHKERRQ(ierr);

  CeedElemRestrictionCreateVector(Erestrictx, &xcoord, NULL);
  CeedVectorSetArray(xcoord, CEED_MEM_HOST, CEED_COPY_VALUES,
                     (PetscScalar *)coordArray);
  ierr = VecRestoreArrayRead(coords, &coordArray); CHKERRQ(ierr);

  // Create the persistent vectors that will be needed in setup and apply
  CeedInt nqpts;
  CeedBasisGetNumQuadraturePoints(basisu, &nqpts);
  CeedVectorCreate(ceed, qdatasize*nelem*nqpts, &qdata);
  CeedVectorCreate(ceed, xlsize, &xceed);
  CeedVectorCreate(ceed, xlsize, &yceed);

  // Create the Q-function that builds the operator (i.e. computes its
  // quadrature data) and set its context data
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].setupgeo,
                              problemOptions[problemChoice].setupgeofname, &qf_setupgeo);
  CeedQFunctionAddInput(qf_setupgeo, "dx", ncompx*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_setupgeo, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(qf_setupgeo, "qdata", qdatasize, CEED_EVAL_NONE);

  // Set up PDE operator
  CeedInt inscale = problemOptions[problemChoice].inmode==CEED_EVAL_GRAD ? dim : 1;
  CeedInt outscale = problemOptions[problemChoice].outmode==CEED_EVAL_GRAD ? dim : 1;
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].apply,
                              problemOptions[problemChoice].applyfname, &qf_apply);
  CeedQFunctionAddInput(qf_apply, "u", ncompu*inscale,
                        problemOptions[problemChoice].inmode);
  CeedQFunctionAddInput(qf_apply, "qdata", qdatasize, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qf_apply, "v", ncompu*outscale,
                         problemOptions[problemChoice].outmode);

  // Create the operator that builds the quadrature data for the operator
  CeedOperatorCreate(ceed, qf_setupgeo, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &op_setupgeo);
  CeedOperatorSetField(op_setupgeo, "dx", Erestrictx, CEED_TRANSPOSE,
                       basisx, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setupgeo, "weight", Erestrictxi, CEED_NOTRANSPOSE,
                       basisx, CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setupgeo, "qdata", Erestrictqdi, CEED_NOTRANSPOSE,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

  // Create the operator
  CeedOperatorCreate(ceed, qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &op_apply);
  CeedOperatorSetField(op_apply, "u", Erestrictu, CEED_TRANSPOSE,
                       basisu, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_apply, "qdata", Erestrictqdi, CEED_NOTRANSPOSE,
                       CEED_BASIS_COLLOCATED, qdata);
  CeedOperatorSetField(op_apply, "v", Erestrictu, CEED_TRANSPOSE,
                       basisu, CEED_VECTOR_ACTIVE);

  // Setup qdata
  CeedOperatorApply(op_setupgeo, xcoord, qdata, CEED_REQUEST_IMMEDIATE);

  // Set up RHS if needed
  if (setup_rhs) {
    CeedQFunction qf_setuprhs;
    CeedOperator op_setuprhs;
    CeedVectorCreate(ceed, nelem*nqpts*ncompu, target);

    // Create the q-function that sets up the RHS and true solution
    CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].setuprhs,
                                problemOptions[problemChoice].setuprhsfname, &qf_setuprhs);
    CeedQFunctionAddInput(qf_setuprhs, "x", dim, CEED_EVAL_INTERP);
    CeedQFunctionAddInput(qf_setuprhs, "dx", ncompx*dim, CEED_EVAL_GRAD);
    CeedQFunctionAddInput(qf_setuprhs, "weight", 1, CEED_EVAL_WEIGHT);
    CeedQFunctionAddOutput(qf_setuprhs, "true_soln", ncompu, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf_setuprhs, "rhs", ncompu, CEED_EVAL_INTERP);

    // Create the operator that builds the RHS and true solution
    CeedOperatorCreate(ceed, qf_setuprhs, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_setuprhs);
    CeedOperatorSetField(op_setuprhs, "x", Erestrictx, CEED_TRANSPOSE,
                         basisx, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_setuprhs, "dx", Erestrictx, CEED_TRANSPOSE,
                         basisx, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_setuprhs, "weight", Erestrictxi, CEED_NOTRANSPOSE,
                         basisx, CEED_VECTOR_NONE);
    CeedOperatorSetField(op_setuprhs, "true_soln", Erestrictui, CEED_NOTRANSPOSE,
                         CEED_BASIS_COLLOCATED, *target);
    CeedOperatorSetField(op_setuprhs, "rhs", Erestrictu, CEED_TRANSPOSE,
                         basisu, CEED_VECTOR_ACTIVE);

    // Setup RHS and target
    CeedOperatorApply(op_setuprhs, xcoord, rhsceed, CEED_REQUEST_IMMEDIATE);
    CeedVectorSyncArray(rhsceed, CEED_MEM_HOST);

    // Cleanup
    CeedQFunctionDestroy(&qf_setuprhs);
    CeedOperatorDestroy(&op_setuprhs);
  }

  // Cleanup
  CeedQFunctionDestroy(&qf_setupgeo);
  CeedOperatorDestroy(&op_setupgeo);
  CeedVectorDestroy(&xcoord);

  // Save libCEED data required for level
  data->basisx = basisx; data->basisu = basisu;
  data->Erestrictx = Erestrictx;
  data->Erestrictu = Erestrictu;
  data->Erestrictxi = Erestrictxi;
  data->Erestrictui = Erestrictui;
  data->Erestrictqdi = Erestrictqdi;
  data->qf_apply = qf_apply;
  data->op_apply = op_apply;
  data->qdata = qdata;
  data->xceed = xceed;
  data->yceed = yceed;

  PetscFunctionReturn(0);
}

#ifdef multigrid
// Setup libCEED level transfer operator objects
static PetscErrorCode CeedLevelTransferSetup(Ceed ceed, CeedInt numlevels,
    CeedInt ncompu, problemType problemChoice, CeedData *data, CeedInt *leveldegrees,
    CeedQFunction qf_restrict, CeedQFunction qf_prolong) {
  // Return early if numlevels=1
  if (numlevels==1)
    PetscFunctionReturn(0);

  // Set up each level
  for (CeedInt i=1; i<numlevels; i++) {
    // P coarse and P fine
    CeedInt Pc = leveldegrees[i-1] + 1;
    CeedInt Pf = leveldegrees[i] + 1;

    // Restriction - Fine to corse
    CeedBasis basisctof;
    CeedOperator op_restrict;

    // Basis
    CeedBasisCreateTensorH1Lagrange(ceed, 3, ncompu, Pc, Pf,
                                    CEED_GAUSS_LOBATTO, &basisctof);

    // Create the restriction operator
    CeedOperatorCreate(ceed, qf_restrict, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_restrict);
    CeedOperatorSetField(op_restrict, "input", data[i]->Erestrictu,
                         CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED,
                         CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_restrict, "output", data[i-1]->Erestrictu,
                         CEED_TRANSPOSE, basisctof, CEED_VECTOR_ACTIVE);

    // Save libCEED data required for level
    data[i]->basisctof = basisctof;
    data[i]->op_restrict = op_restrict;

    // Interpolation - Corse to fine
    CeedOperator op_interp;

    // Create the prolongation operator
    CeedOperatorCreate(ceed, qf_prolong, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_interp);
    CeedOperatorSetField(op_interp, "input", data[i-1]->Erestrictu,
                         CEED_NOTRANSPOSE, basisctof, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_interp, "output", data[i]->Erestrictu,
                         CEED_TRANSPOSE, CEED_BASIS_COLLOCATED,
                         CEED_VECTOR_ACTIVE);

    // Save libCEED data required for level
    data[i]->op_interp = op_interp;
  }

  PetscFunctionReturn(0);
}
#endif

// -----------------------------------------------------------------------------
// Mat Shell Functions
// -----------------------------------------------------------------------------

#ifdef multigrid
// This function returns the computed diagonal of the operator
static PetscErrorCode MatGetDiag(Mat A, Vec D) {
  PetscErrorCode ierr;
  UserO user;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  ierr = VecCopy(user->diag, D); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#endif

// This function uses libCEED to compute the linear elasticity with
// Dirichlet boundary conditions
static PetscErrorCode FormResidual_Ceed(SNES snes, Mat A, Vec X, Vec Y, Physics phys) {
  PetscErrorCode ierr;
  UserO user;
  PetscScalar *x, *y;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Global-to-local
  ierr = DMGlobalToLocalBegin(user->dm, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->dm, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Yloc); CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  ierr = VecGetArray(user->Yloc, &y); CHKERRQ(ierr);
  CeedVectorSetArray(user->xceed, CEED_MEM_HOST, CEED_USE_POINTER, x);
  CeedVectorSetArray(user->yceed, CEED_MEM_HOST, CEED_USE_POINTER, y);

  // Apply CEED operator
  CeedOperatorApply(user->op, user->xceed, user->yceed, CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->yceed, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(user->dm, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(user->dm, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#ifdef multigrid
// This function uses libCEED to compute the action of the interp operator
static PetscErrorCode MatMult_Interp(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserIR user;
  PetscScalar *x, *y;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Global-to-local
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->dmc, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->dmc, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Yloc); CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  ierr = VecGetArray(user->Yloc, &y); CHKERRQ(ierr);
  CeedVectorSetArray(user->ceedvecc, CEED_MEM_HOST, CEED_USE_POINTER, x);
  CeedVectorSetArray(user->ceedvecf, CEED_MEM_HOST, CEED_USE_POINTER, y);

  // Apply CEED operator
  CeedOperatorApply(user->op, user->ceedvecc, user->ceedvecf,
                    CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->ceedvecf, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

  // Multiplicity
  ierr = VecPointwiseMult(user->Yloc, user->Yloc, user->mult);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(user->dmf, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(user->dmf, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// This function uses libCEED to compute the action of the restriction operator
static PetscErrorCode MatMult_Restrict(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserIR user;
  PetscScalar *x, *y;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Global-to-local
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->dmf, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->dmf, X, INSERT_VALUES, user->Xloc);
  CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Yloc); CHKERRQ(ierr);

  // Multiplicity
  ierr = VecPointwiseMult(user->Xloc, user->Xloc, user->mult); CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  ierr = VecGetArray(user->Yloc, &y); CHKERRQ(ierr);
  CeedVectorSetArray(user->ceedvecf, CEED_MEM_HOST, CEED_USE_POINTER, x);
  CeedVectorSetArray(user->ceedvecc, CEED_MEM_HOST, CEED_USE_POINTER, y);

  // Apply CEED operator
  CeedOperatorApply(user->op, user->ceedvecf, user->ceedvecc,
                    CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->ceedvecc, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(user->dmc, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(user->dmc, user->Yloc, ADD_VALUES, Y);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#endif

// This function calculates the error in the final solution
static PetscErrorCode ComputeErrorMax(UserO user, CeedOperator op_error,
                                      Vec X, CeedVector target,
                                      PetscReal *maxerror) {
  PetscErrorCode ierr;
  PetscScalar *x;
  CeedVector collocated_error;
  CeedInt length;

  PetscFunctionBeginUser;
  CeedVectorGetLength(target, &length);
  CeedVectorCreate(user->ceed, length, &collocated_error);

  // Global-to-local
  ierr = DMGlobalToLocal(user->dm, X, INSERT_VALUES, user->Xloc); CHKERRQ(ierr);

  // Setup CEED vector
  ierr = VecGetArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  CeedVectorSetArray(user->xceed, CEED_MEM_HOST, CEED_USE_POINTER, x);

  // Apply CEED operator
  CeedOperatorApply(op_error, user->xceed, collocated_error,
                    CEED_REQUEST_IMMEDIATE);

  // Restore PETSc vector
  VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);

  // Reduce max error
  *maxerror = 0;
  const CeedScalar *e;
  CeedVectorGetArrayRead(collocated_error, CEED_MEM_HOST, &e);
  for (CeedInt i=0; i<length; i++) {
    *maxerror = PetscMax(*maxerror, PetscAbsScalar(e[i]));
  }
  CeedVectorRestoreArrayRead(collocated_error, &e);
  ierr = MPI_Allreduce(MPI_IN_PLACE, maxerror,
                       1, MPIU_REAL, MPIU_MAX, user->comm); CHKERRQ(ierr);

  // Cleanup
  CeedVectorDestroy(&collocated_error);

  PetscFunctionReturn(0);
}
#endif
