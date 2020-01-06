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
#include "qfunctions/solid/constantForce.h"     // Constant forcing function
#include "qfunctions/solid/manufacturedForce.h" // Manufactured solution

// Problem options
typedef enum {  //SmallStrain      FiniteStrain
  ELAS_LIN = 0, ELAS_HYPER_SS = 1, ELAS_HYPER_FS = 2
} problemType;
static const char *const problemTypes[] = {"linElas","hyperSS","hyperFS", "problemType","ELAS_",0};

// Forcing function options
typedef enum {
  FORCE_NONE = 0, FORCE_CONST = 1, FORCE_MMS = 2
} forcingType;
static const char *const forcingTypes[] = {"none","constant","manufactured","FORCE_",0};

typedef enum{
   BDRY_WALL = 0, BDRY_WALL_FORCE = 1 , BDRY_MMS = 2
}boundaryType;
static const char *const boundaryTypes[] = {"wall","force","mms", "boundaryType","BDRY_",0};

typedef PetscErrorCode BCFunc(PetscInt, PetscReal, const PetscReal*, PetscInt, PetscScalar*, void*);
BCFunc BCBend1_ss, BCBend2_ss, BCMMS;
BCFunc* boundaryOptions[] = {BCBend1_ss, BCBend2_ss, BCMMS};

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------

typedef struct{
  char          meshFile[PETSC_MAX_PATH_LEN]; // exodusII mesh file
  problemType   problemChoice;
  forcingType   forcingChoice;
  PetscInt      degree;
  boundaryType  boundaryChoice;
}AppCtx;


//Physics struct for Each problem is moved to its .h file (libCEED necessity)
//linElas.h ,hyperSS.h, hyperFS.h

// Problem specific data
typedef struct {
  CeedInt           qdatasize;
  CeedQFunctionUser setupgeo, apply, jacob, error;
  const char        *setupgeofname, *applyfname, *jacobfname, *errorfname;
  CeedQuadMode      qmode;
}problemData;

problemData problemOptions[3] = {
  [ELAS_LIN] = {
      .qdatasize = 10, // For linear Elasticty we could do 6
      .setupgeo = SetupGeo,
      .apply = LinElasF,
      .jacob = LinElasdF,
      .error = Error,
      .setupgeofname = SetupGeo_loc,
      .applyfname = LinElasF_loc,
      .jacobfname = LinElasdF_loc,
      .errorfname = Error_loc,
      .qmode = CEED_GAUSS
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
  },
  [ELAS_HYPER_FS] = {
     // .qdatasize = 10,
     // .setupgeo = SetupGeo,
     // .apply = HyperFS,
     // .error = Error,
     // .setupgeofname = SetupGeo_loc,
     // .applyfname = HyperFS_loc,
     // .errorfname = Error_loc,
     // .qmode = CEED_GAUSS
    }
};

// Forcing function data

typedef struct{
  CeedQFunctionUser setupforcing;
  const char        *setupforcingfname;
}forcingData;

forcingData forcingOptions[3] = {
  [FORCE_NONE] = {
      .setupforcing = NULL,
      .setupforcingfname = NULL
  },
  [FORCE_CONST] = {
      .setupforcing = SetupConstantForce,
      .setupforcingfname = SetupConstantForce_loc
  },
  [FORCE_MMS] = {
      .setupforcing = SetupMMSForce,
      .setupforcingfname = SetupMMSForce_loc
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
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi, Erestrictui, Erestrictqdi, ErestrictGradui;
  CeedQFunction       qf_apply, qf_jacob;
  CeedOperator        op_apply, op_jacob, op_restrict, op_interp;
  CeedVector          qdata, gradu, xceed, yceed, truesolution;
};

// -----------------------------------------------------------------------------
// Helper Functions
// -----------------------------------------------------------------------------
static int processCommandLineOptions(MPI_Comm comm, AppCtx *appCtx){

  PetscErrorCode ierr;
  PetscBool meshFileFlag = PETSC_FALSE;
  PetscBool degreeFalg = PETSC_FALSE;
  PetscBool boundaryFlag = PETSC_FALSE;
  appCtx->problemChoice = ELAS_LIN; //-problem = Linear Elasticity if not given
  appCtx->degree = 0;
  appCtx->boundaryChoice = BDRY_WALL;
  appCtx->forcingChoice = FORCE_NONE;

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

  ierr = PetscOptionsEnum("-forcing", "Set forcing function option", NULL, forcingTypes,
                          (PetscEnum)appCtx->forcingChoice,(PetscEnum *)&appCtx->forcingChoice, NULL); CHKERRQ(ierr);


 ierr = PetscOptionsEnum("-boundary", "Set Dirichlet (Essential) Boundary option", NULL, boundaryTypes,
                (PetscEnum)appCtx->boundaryChoice,(PetscEnum *)&appCtx->boundaryChoice, &boundaryFlag); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);//End of setting AppCtx
  if(!degreeFalg){
      ierr = PetscPrintf(comm, "-degree option needed\n\n");CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
  }
  if(!meshFileFlag){
      ierr = PetscPrintf(comm, "-mesh option needed (file)\n\n");CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
  }
  if(!boundaryFlag){
      ierr = PetscPrintf(comm, "-boundary option needed\n\n");CHKERRQ(ierr);
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
 PetscQuadrature fq; //face quadrature points (For future: Nuemman boundary)
 // For Dirichlet (Essential) Boundary
 IS              faceSetIS;         //Index Set for Face Sets
 const char      *name="Face Sets"; //PETSc internal requirement
 PetscInt        numFaceSets;        //Number of FaceSets in faceSetIS
 const PetscInt  *faceSetIds;       //id of each FaceSet

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
  //ierr = DMCreateDS(dm); CHKERRQ(ierr);

  // Add Dirichlet (Essential) boundray
  ierr = DMCreateDS(dm); CHKERRQ(ierr);
  ierr = DMGetLabelIdIS(dm, name, &faceSetIS);CHKERRQ(ierr);
  ierr = ISGetLocalSize(faceSetIS,&numFaceSets);
  ierr = ISGetIndices(faceSetIS, &faceSetIds);CHKERRQ(ierr);

  for (PetscInt i=0; i<numFaceSets; i++){
    ierr = DMAddBoundary(dm,DM_BC_ESSENTIAL,NULL,"Face Sets",0,0,NULL,
        (void(*)(void))boundaryOptions[appCtx->boundaryChoice],1,&faceSetIds[i],(void *)(&faceSetIds[i]));CHKERRQ(ierr);
   }
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
  CeedVectorDestroy(&data->gradu);
  CeedVectorDestroy(&data->xceed);
  CeedVectorDestroy(&data->yceed);
  CeedVectorDestroy(&data->truesolution);
  CeedBasisDestroy(&data->basisx);
  CeedBasisDestroy(&data->basisu);
  CeedElemRestrictionDestroy(&data->Erestrictu);
  CeedElemRestrictionDestroy(&data->Erestrictx);
  CeedElemRestrictionDestroy(&data->Erestrictui);
  CeedElemRestrictionDestroy(&data->ErestrictGradui);
  CeedElemRestrictionDestroy(&data->Erestrictxi);
  CeedElemRestrictionDestroy(&data->Erestrictqdi);
  CeedQFunctionDestroy(&data->qf_apply);
  CeedOperatorDestroy(&data->op_apply);
  CeedQFunctionDestroy(&data->qf_jacob);
  CeedOperatorDestroy(&data->op_jacob);
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
                                PetscInt Ugsz, PetscInt Ulocsz, CeedVector forceceed){

int          ierr;
CeedInt      degree = appCtx->degree;
problemType  problemChoice = appCtx->problemChoice;
forcingType  forcingChoice = appCtx->forcingChoice;
CeedInt      qdatasize = problemOptions[appCtx->problemChoice].qdatasize;
CeedInt      P, Q;
CeedInt      dim, ncompx;
CeedBasis    basisx, basisu;
DM           dmcoord;
CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi, Erestrictui, Erestrictqdi, ErestrictGradui;
CeedInt      cStart, cEnd, nelem;
Vec          coords;
const PetscScalar *coordArray;
PetscSection section;
CeedVector   xcoord, qdata, gradu, xceed, yceed;
CeedQFunction qf_setupgeo, qf_apply, qf_jacob;
CeedOperator op_setupgeo, op_apply, op_jacob;

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
CeedElemRestrictionCreateIdentity(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, dim* ncompu, &ErestrictGradui); CHKERRQ(ierr);

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
CeedVectorCreate(ceed, dim*ncompu*nelem*nqpts, &gradu);
CeedVectorCreate(ceed, Ulocsz, &xceed);
CeedVectorCreate(ceed, Ulocsz, &yceed);

// Create the QFunction and Operator that computes the quadrature data
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

// Create the QFunction and Operator that evaluates the residual
CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].apply,problemOptions[problemChoice].applyfname, &qf_apply);
CeedQFunctionAddInput(qf_apply, "du", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionAddInput(qf_apply, "qdata", qdatasize, CEED_EVAL_NONE);
CeedQFunctionAddOutput(qf_apply, "dv", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionAddOutput(qf_apply, "gradu", ncompu*dim, CEED_EVAL_NONE);
CeedQFunctionSetContext(qf_apply, &phys, sizeof(phys));
CeedOperatorCreate(ceed, qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_apply);
CeedOperatorSetField(op_apply, "du", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
CeedOperatorSetField(op_apply, "qdata", Erestrictqdi, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, qdata);
CeedOperatorSetField(op_apply, "dv", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
CeedOperatorSetField(op_apply, "gradu", ErestrictGradui, CEED_NOTRANSPOSE, basisu, gradu);

// Create the QFunction and Operator that calculates the Jacobian
CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].jacob,problemOptions[problemChoice].jacobfname, &qf_jacob);
CeedQFunctionAddInput(qf_jacob, "deltadu", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionAddInput(qf_jacob, "qdata", qdatasize, CEED_EVAL_NONE);
CeedQFunctionAddInput(qf_jacob, "gradu", ncompu*dim, CEED_EVAL_NONE);
CeedQFunctionAddOutput(qf_jacob, "deltadv", ncompu*dim, CEED_EVAL_GRAD);
CeedQFunctionSetContext(qf_jacob, &phys, sizeof(phys));
CeedOperatorCreate(ceed, qf_jacob, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_jacob);
CeedOperatorSetField(op_jacob, "deltadu", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
CeedOperatorSetField(op_jacob, "qdata", Erestrictqdi, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, qdata);
CeedOperatorSetField(op_jacob, "deltadv", Erestrictu, CEED_NOTRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
CeedOperatorSetField(op_jacob, "gradu", ErestrictGradui, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, gradu);

// Compute the quadrature data
CeedOperatorApply(op_setupgeo, xcoord, qdata, CEED_REQUEST_IMMEDIATE);

  // Set up forcing term, if needed
  if (forcingChoice != FORCE_NONE) {
    CeedQFunction qf_setupforce;
    CeedOperator op_setupforce;
    if (forcingChoice == FORCE_MMS)
      CeedVectorCreate(ceed, nelem*nqpts*ncompu, &(data->truesolution));

    // Create the q-function that sets up the forcing vector (and true solution for MMS)
    CeedQFunctionCreateInterior(ceed, 1, forcingOptions[forcingChoice].setupforcing, forcingOptions[forcingChoice].setupforcingfname, &qf_setupforce);
    CeedQFunctionAddInput(qf_setupforce, "dx", ncompx*dim, CEED_EVAL_GRAD);
    CeedQFunctionAddInput(qf_setupforce, "weight", 1, CEED_EVAL_WEIGHT);
    CeedQFunctionAddInput(qf_setupforce, "x", dim, CEED_EVAL_INTERP);
    CeedQFunctionAddOutput(qf_setupforce, "rhs", ncompu, CEED_EVAL_INTERP);
    if (forcingChoice == FORCE_MMS)
      CeedQFunctionAddOutput(qf_setupforce, "true_soln", ncompu, CEED_EVAL_NONE);
    CeedQFunctionSetContext(qf_setupforce, &phys, sizeof(phys));

    // Create the operator that builds the forcing vector (and true solution for MMS)
    CeedOperatorCreate(ceed, qf_setupforce, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &op_setupforce);
    CeedOperatorSetField(op_setupforce, "dx", Erestrictx, CEED_TRANSPOSE, basisx, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_setupforce, "weight", Erestrictxi, CEED_NOTRANSPOSE, basisx, CEED_VECTOR_NONE);
    CeedOperatorSetField(op_setupforce, "x", Erestrictx, CEED_TRANSPOSE, basisx, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_setupforce, "rhs", Erestrictu, CEED_TRANSPOSE, basisu, CEED_VECTOR_ACTIVE);
    if (forcingChoice == FORCE_MMS)
      CeedOperatorSetField(op_setupforce, "true_soln", Erestrictui, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, data->truesolution);

    // Setup forcing vector (and true solution, for MMS)
    CeedOperatorApply(op_setupforce, xcoord, forceceed, CEED_REQUEST_IMMEDIATE);
    CeedVectorSyncArray(forceceed, CEED_MEM_HOST);

    // Cleanup
    CeedQFunctionDestroy(&qf_setupforce);
    CeedOperatorDestroy(&op_setupforce);
  }

  // Cleanup
  CeedQFunctionDestroy(&qf_setupgeo);
  CeedOperatorDestroy(&op_setupgeo);
  CeedVectorDestroy(&xcoord);

  // Save libCEED data required for level
  data->basisx = basisx;
  data->basisu = basisu;
  data->Erestrictx = Erestrictx;
  data->Erestrictu = Erestrictu;
  data->Erestrictxi = Erestrictxi;
  data->Erestrictui = Erestrictui;
  data->Erestrictqdi = Erestrictqdi;
  data->ErestrictGradui = ErestrictGradui;
  data->qf_apply = qf_apply;
  data->op_apply = op_apply;
  data->qf_jacob = qf_jacob;
  data->op_jacob = op_jacob;
  data->qdata = qdata;
  data->gradu = gradu;
  data->xceed = xceed;
  data->yceed = yceed;

PetscFunctionReturn(0);
}

static PetscErrorCode ApplyLocalCeedOp(Vec X, Vec Y, UserMult user){

    PetscErrorCode ierr;
    PetscScalar *x, *y;

    PetscFunctionBeginUser;
    ierr = DMGlobalToLocalBegin(user->dm, X, INSERT_VALUES, user->Xloc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->dm, X, INSERT_VALUES, user->Xloc); CHKERRQ(ierr);
    ierr = DMPlexInsertBoundaryValues(user->dm, PETSC_TRUE, user->Xloc, 0,NULL,NULL,NULL); CHKERRQ(ierr);
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
    ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
    ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

    // Local-to-global
    ierr = VecZeroEntries(Y); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(user->dm, user->Yloc, ADD_VALUES, Y); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->dm, user->Yloc, ADD_VALUES, Y); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode FormResidual_Ceed(SNES snes, Vec X, Vec Y, void *ptr) {
  PetscErrorCode ierr;
  UserMult user = (UserMult)ptr;

  PetscFunctionBeginUser;
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);
  // ierr = VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  // PetscPrintf(PETSC_COMM_WORLD, "\n\n X: \n\n");
  // ierr = VecView(Y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode ApplyJacobian_Ceed(Mat A, Vec X, Vec Y){
    PetscErrorCode ierr;
    UserMult user;

    PetscFunctionBeginUser;
    ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);
    ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode CreateMatrixFreeCtx(MPI_Comm comm, DM dm, Vec vloc, CeedData ceeddata, Ceed ceed, UserMult residualCtx, UserMult jacobianCtx){

    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    residualCtx->comm = comm;
    residualCtx->dm = dm;
    residualCtx->Xloc = vloc;
    ierr = VecDuplicate(vloc, &residualCtx->Yloc);CHKERRQ(ierr);
    residualCtx->xceed = ceeddata->xceed;
    residualCtx->yceed = ceeddata->yceed;
    residualCtx->op = ceeddata->op_apply;
    residualCtx->ceed = ceed;
    ierr = PetscMemcpy(jacobianCtx, residualCtx, sizeof(*residualCtx)); CHKERRQ(ierr);
    jacobianCtx->op = ceeddata->op_jacob;
    PetscFunctionReturn(0);
}

// boundary Functions
//
/* BCMMS boundary function explanation
ss : (sideset)
MMS: Boundary corresponding to the Method of Manufactured Solutions
Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
Also check ..\meshes\cyl-hol.8.jou

 left: sideset 999
 right: sideset 998
 outer: sideset 997
 inner: sideset 996
  / \-------------------\              y
 /   \                   \             |
(  O  )                   )      x ____|
 \   /                   /              \
  \ /-------------------/                \ z

 values on all points of the mesh is set beased on given solution below
 on u[0], u[1], u[2]
*/
PetscErrorCode BCMMS(PetscInt dim, PetscReal time, const PetscReal coords[],
                       PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscScalar x = coords[0];
  PetscScalar y = coords[1];
  PetscScalar z = coords[2];

  PetscFunctionBeginUser;

  u[0]= exp(2*x)*sin(3*y)*cos(4*z);
  u[1]= exp(3*y)*sin(4*z)*cos(2*x);
  u[2]= exp(4*z)*sin(2*x)*cos(3*y);

  PetscFunctionReturn(0);
}
/*BCBend2_ss boundary function explanation
ss : (sideset)
2_ss : two sides of the geometry
Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
Also check ..\meshes\cyl-hol.8.jou

 left: sideset 999
 right: sideset 998
  / \-------------------\              y
 /   \                   \             |
(  O  )                   )      x ____|
 \   /                   /              \
  \ /-------------------/                \ z

  0 values on the left side of the cyl-hole (sideset 999)
 -1 values on y direction of the right side of the cyl-hole (sideset 999)
*/
PetscErrorCode BCBend2_ss(PetscInt dim, PetscReal time, const PetscReal coords[],
                       PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscInt* faceID = (PetscInt*)ctx;
  PetscFunctionBeginUser;

  switch (*faceID){
  case 999: // left side of the cyl-hol
    u[0]= 0;
    u[1]= 0;
    u[2]= 0;
    break;
   case 998: //right side of the cyl-hol
    u[0]= 0;
    u[1]= -1; //bend in the -y direction
    u[2]= 0;
   break;
}
  PetscFunctionReturn(0);
}
/*BCBend1_ss boundary function explanation
ss : (sideset)
1_ss : 1 side (left side) of the geometry
Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
Also check ..\meshes\cyl-hol.8.jou

 left: sideset 999

  / \-------------------\              y
 /   \                   \             |
(  O  )                   )      x ____|
 \   /                   /              \
  \ /-------------------/                \ z

  0 values on the left side of the cyl-hole (sideset 999)
*/
PetscErrorCode BCBend1_ss(PetscInt dim, PetscReal time, const PetscReal coords[],
                       PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscFunctionBeginUser;

  u[0]= 0;
  u[1]= 0;
  u[2]= 0;

  PetscFunctionReturn(0);
}


#endif setup_h
