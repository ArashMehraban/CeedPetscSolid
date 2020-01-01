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
static const char *const problemTypes[] = {"linElas","hyperSS","hyperFS", "problemType","ELAS_",0};

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------

typedef struct{
  char          meshFile[PETSC_MAX_PATH_LEN]; // exodusII mesh file
  problemType   problemChoice;
  PetscInt      degree;
}AppCtx;


//Physics struct for Each problem is moved to its .h file (libCEED necessity)
//linElas.h ,hyperSS.h, hyperFS.h

// Problem specific data
typedef struct {
  CeedInt           qdatasize;
  CeedQFunctionUser setupgeo, apply, error;
  const char        *setupgeofname, *applyfname, *errorfname;
  CeedQuadMode      qmode;
  PetscErrorCode    (*bcs_func)(PetscInt, PetscReal, const PetscReal *, PetscInt, PetscScalar *, void *);
}problemData;

problemData problemOptions[3] = {
  [ELAS_LIN] = {
      .qdatasize = 10, // For linear Elasticty we could do 6
      .setupgeo = SetupGeo,
      .apply = LinElas,
      .error = Error,
      .setupgeofname = SetupGeo_loc,
      .applyfname = LinElas_loc,
      .errorfname = Error_loc,
      .qmode = CEED_GAUSS,
      .bcs_func = NULL // Drichlet of all 1's
  },
  [ELAS_HYPER_SS] = {
     // .qdatasize = 10,
     // .setupgeo = SetupGeo,
     // .apply = HyperSS,
     // .error = Error,
     // .setupgeofname = SetupGeo_loc,
     // .applyfname = HyperSS_loc,
     // .errorfname = Error_loc,
     // .qmode = CEED_GAUSS,
     //  .bcs_func = NULL // Drichlet of all 1's
  },
  [ELAS_HYPER_FS] = {
     // .qdatasize = 10,
     // .setupgeo = SetupGeo,
     // .apply = HyperFS,
     // .error = Error,
     // .setupgeofname = SetupGeo_loc,
     // .applyfname = HyperFS_loc,
     // .errorfname = Error_loc,
     // .qmode = CEED_GAUSS,
     //  .bcs_func = NULL // Drichlet of all 1's
    }
};

// Data for PETSc Matshell
typedef struct UserMult_private *UserMult;
struct UserMult_private {
  MPI_Comm     comm;
  DM           dm;
  Vec          Xloc, Yloc, diag;
  CeedVector   xceed, yceed;
  CeedOperator op;
  Ceed         ceed;
};

// libCEED data struct for level
typedef struct CeedData_private *CeedData;
struct CeedData_private {
  Ceed                ceed;
  CeedBasis           basisx, basisu, basisctof;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi, Erestrictui, Erestrictqdi;
  CeedQFunction       qf_apply;
  CeedOperator        op_apply, op_restrict, op_interp;
  CeedVector          qdata, xceed, yceed;
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

  ierr = PetscOptionsBegin(comm, NULL, "Elasticity / Hyperelasticity in PETSc with libCEED", NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-degree", "Polynomial degree of tensor product basis", NULL, appCtx->degree, &appCtx->degree,
                         &degreeFalg); CHKERRQ(ierr);

  ierr = PetscOptionsString("-mesh", "Read mesh from file", NULL, appCtx->meshFile, appCtx->meshFile,
                            sizeof(appCtx->meshFile), &meshFileFlag); CHKERRQ(ierr);
  #if !defined(PETSC_HAVE_EXODUSII)
    SETERRQ(comm, PETSC_ERR_ARG_WRONG,
     "ExodusII support needed. Reconfigure your Arch with --download-exodusii");
  #endif

  ierr = PetscOptionsEnum("-problem", "Solves Elasticity & Hyperelasticity Problems", NULL, problemTypes,
                          (PetscEnum)appCtx->problemChoice,(PetscEnum *)&appCtx->problemChoice, NULL); CHKERRQ(ierr);


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

static int processPhysics(MPI_Comm comm, Physics phys){

    PetscErrorCode ierr;
    PetscBool nuFlag = PETSC_FALSE;
    PetscBool YoungFlag = PETSC_FALSE;
    phys->nu = 0;
    phys->E = 0;

    PetscFunctionBeginUser;

    ierr = PetscOptionsBegin(comm, NULL,"Elasticity / Hyperelasticity in PETSc with libCEED", NULL); CHKERRQ(ierr);
    ierr = PetscOptionsScalar("-nu", "Poisson's ratio", NULL, phys->nu, &phys->nu, &nuFlag);CHKERRQ(ierr);
    ierr = PetscOptionsScalar("-E", "Young's Modulus", NULL, phys->E, &phys->E, &YoungFlag);CHKERRQ(ierr);
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
  ierr = PetscSpacePolynomialSetTensor(P, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscSpaceSetFromOptions(P); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumComponents(P, ncompu); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumVariables(P, dim); CHKERRQ(ierr);
  ierr = PetscSpaceSetDegree(P, order, order); CHKERRQ(ierr);
  ierr = PetscSpaceSetUp(P); CHKERRQ(ierr);
  //setup FE dual space (Space Q) for tensor polynomials
  ierr = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &Q);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE); CHKERRQ(ierr);
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
  ierr = PetscFESetFromOptions(fe); CHKERRQ(ierr);
  ierr = PetscFESetBasisSpace(fe, P); CHKERRQ(ierr);
  ierr = PetscFESetDualSpace(fe, Q); CHKERRQ(ierr);
  ierr = PetscFESetNumComponents(fe, ncompu); CHKERRQ(ierr);
  ierr = PetscFESetUp(fe); CHKERRQ(ierr);
  ierr = PetscSpaceDestroy(&P); CHKERRQ(ierr);
  ierr = PetscDualSpaceDestroy(&Q); CHKERRQ(ierr);
  /* Create quadrature */
  quadPointsPerEdge = PetscMax(order + 1,1);
  ierr = PetscDTGaussTensorQuadrature(dim, 1, quadPointsPerEdge, -1.0, 1.0, &q); CHKERRQ(ierr);
  ierr = PetscDTGaussTensorQuadrature(dim-1, 1, quadPointsPerEdge, -1.0, 1.0, &fq); CHKERRQ(ierr);
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
                 (void(*)(void))problemOptions[appCtx->problemChoice].bcs_func,numFaceSet,faceSetIds,NULL);CHKERRQ(ierr);
  ierr = ISRestoreIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);
  ierr = ISDestroy(&faceSetIS);CHKERRQ(ierr);
  ierr = DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL); CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

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
static int CreateRestrictionPlex(Ceed ceed, CeedInt P, CeedInt ncomp, CeedElemRestriction *Erestrict, DM dm) {
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
    ierr = DMPlexGetClosureIndices(dm, section, section, c, &numindices, &indices, NULL); CHKERRQ(ierr);
    for (i=0; i<numindices; i+=ncomp) {
      for (PetscInt j=0; j<ncomp; j++) {
        if (indices[i+j] != indices[i] + (PetscInt)(copysign(j, indices[i])))
          SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell %D closure indices not interlaced", c);
      }
      // Essential boundary conditions are encoded as -(loc+1)
      PetscInt loc = indices[i] >= 0 ? indices[i] : -(indices[i] + 1);
      erestrict[eoffset++] = loc/ncomp;
    }
    ierr = DMPlexRestoreClosureIndices(dm, section, section, c, &numindices,  &indices, NULL); CHKERRQ(ierr);
  }

  // Setup CEED restriction
  ierr = DMGetLocalVector(dm, &Uloc); CHKERRQ(ierr);
  ierr = VecGetLocalSize(Uloc, &nnodes); CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(dm, &Uloc); CHKERRQ(ierr);
  CeedElemRestrictionCreate(ceed, nelem, P*P*P, nnodes/ncomp, ncomp,  CEED_MEM_HOST, CEED_COPY_VALUES, erestrict,Erestrict);
  ierr = PetscFree(erestrict); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

//Set up libCEED for a given degree
static int SetupLibceedByDegree(DM dm, Ceed ceed, AppCtx *appCtx, Physics phys, CeedData data, PetscInt ncompu,
                                PetscInt Ugsz,PetscInt Ulocsz){

int          ierr;
CeedInt      degree = appCtx->degree;
problemType  problemChoice = appCtx->problemChoice;
CeedInt      qdatasize = problemOptions[appCtx->problemChoice].qdatasize;
CeedInt      P, Q;
CeedInt      dim, ncompx;
CeedBasis    basisx, basisu;
DM           dmcoord;
CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi, Erestrictui, Erestrictqdi;
CeedInt      cStart, cEnd, nelem;
Vec          coords;
const PetscScalar *coordArray;
PetscSection section;
CeedVector   xcoord, qdata, xceed, yceed;
CeedQFunction qf_setupgeo, qf_apply;
CeedOperator op_setupgeo, op_apply;

PetscFunctionBeginUser;

ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
ncompx = dim;
// CEED bases (solid problems)
P = degree + 1;
Q = P;
CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu, P, Q,problemOptions[problemChoice].qmode, &basisu);
CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, Q,problemOptions[problemChoice].qmode, &basisx);

// CEED restrictions
ierr = DMGetCoordinateDM(dm, &dmcoord); CHKERRQ(ierr);
ierr = DMPlexSetClosurePermutationTensor(dmcoord, PETSC_DETERMINE, NULL);CHKERRQ(ierr);

CreateRestrictionPlex(ceed, 2, ncompx, &Erestrictx, dmcoord);
CreateRestrictionPlex(ceed, P, ncompu, &Erestrictu, dm);

ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
nelem = cEnd - cStart;
CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, ncompu, &Erestrictui); CHKERRQ(ierr);
CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, qdatasize, &Erestrictqdi); CHKERRQ(ierr);
CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, ncompx, &Erestrictxi); CHKERRQ(ierr);

// Element coordinates
ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
ierr = VecGetArrayRead(coords, &coordArray); CHKERRQ(ierr);
ierr = DMGetSection(dmcoord, &section); CHKERRQ(ierr);

CeedElemRestrictionCreateVector(Erestrictx, &xcoord, NULL);
CeedVectorSetArray(xcoord, CEED_MEM_HOST, CEED_COPY_VALUES, (PetscScalar *)coordArray);

// Create the persistent vectors that will be needed in setup and apply
ierr = VecRestoreArrayRead(coords, &coordArray); CHKERRQ(ierr);
CeedInt nqpts;
CeedBasisGetNumQuadraturePoints(basisu, &nqpts);
CeedVectorCreate(ceed, qdatasize*nelem*nqpts, &qdata);
CeedVectorCreate(ceed, Ulocsz, &xceed);
CeedVectorCreate(ceed, Ulocsz, &yceed);

// Create the Q-function that builds the operator (i.e. computes its
// quadrature data) and set its context data
// qdata returns dXdx_i,j and w * det.
CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].setupgeo,problemOptions[problemChoice].setupgeofname,
                            &qf_setupgeo);
CeedQFunctionAddInput(qf_setupgeo, "dx", ncompx*dim, CEED_EVAL_GRAD);
CeedQFunctionAddInput(qf_setupgeo, "weight", 1, CEED_EVAL_WEIGHT);
CeedQFunctionAddOutput(qf_setupgeo, "qdata", qdatasize, CEED_EVAL_NONE);
CeedOperatorCreate(ceed, qf_setupgeo, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_setupgeo);
CeedOperatorSetField(op_setupgeo, "dx", Erestrictx, CEED_NOTRANSPOSE, basisx, xcoord);
CeedOperatorSetField(op_setupgeo, "weight", Erestrictx, CEED_NOTRANSPOSE, basisx, CEED_VECTOR_NONE);
CeedOperatorSetField(op_setupgeo, "qdata", Erestrictqdi, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, qdata);

CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].apply,problemOptions[problemChoice].applyfname, &qf_apply);
CeedQFunctionAddInput(qf_apply, "du", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionAddInput(qf_apply, "qdata", qdatasize, CEED_EVAL_NONE);
CeedQFunctionAddOutput(qf_apply, "dv", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionSetContext(qf_apply, &phys, sizeof(phys));
CeedOperatorCreate(ceed, qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_apply);
CeedOperatorSetField(op_apply, "dv", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
CeedOperatorSetField(op_apply, "qdata", Erestrictqdi, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, qdata);
CeedOperatorSetField(op_apply, "du", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);

CeedOperatorApply(op_setupgeo, xcoord, qdata, CEED_REQUEST_IMMEDIATE);

PetscFunctionReturn(0);
}


#endif
