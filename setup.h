#ifndef setup_h
#define setup_h

#include <stdbool.h>
#include <string.h>

#include <petsc.h>
#include <petscdmplex.h>
#include <petscfe.h>

#include <ceed.h>

#include "qfunctions/solid/common.h"            // Geometric factors
#include "qfunctions/solid/linElas.h"           // Linear elasticity
#include "qfunctions/solid/hyperSS.h"           // Hyperelasticity small strain
#include "qfunctions/solid/hyperFS.h"           // Hyperelasticity finite strain
#include "qfunctions/solid/constantForce.h"     // Constant forcing function
#include "qfunctions/solid/manufacturedForce.h" // Manufactured solution forcing
#include "qfunctions/solid/manufacturedTrue.h"  // Manufactured true solution

// -----------------------------------------------------------------------------
// Command Line Options
// -----------------------------------------------------------------------------
// Problem options
typedef enum {
  ELAS_LIN = 0, ELAS_HYPER_SS = 1, ELAS_HYPER_FS = 2
} problemType;
static const char *const problemTypes[] = {"linElas",
                                           "hyperSS",
                                           "hyperFS",
                                           "problemType","ELAS_",0};
static const char *const problemTypesForDisp[] = {"Linear elasticity",
                                                  "Hyper elasticity small strain",
                                                  "Hyper elasticity finite strain"};

// Forcing function options
typedef enum {
  FORCE_NONE = 0, FORCE_CONST = 1, FORCE_MMS = 2
} forcingType;
static const char *const forcingTypes[] = {"none",
                                           "constant",
                                           "mms",
                                           "forcingType","FORCE_",0};
static const char *const forcingTypesForDisp[] = {"None",
                                                  "Constant",
                                                  "Manufactured solution"};

// Boundary condition options
typedef enum {
  BDRY_WALL_NONE = 0, BDRY_WALL_WEIGHT = 1, BDRY_MMS = 2
} boundaryType;
static const char *const boundaryTypes[] = {"wall_none",
                                            "wall_weight",
                                            "mms",
                                            "boundaryType","BDRY_",0};
static const char *const boundaryTypesForDisp[] = {"Wall with free end",
                                                   "Wall with weighted end",
                                                   "Manufactured solution"};

typedef PetscErrorCode BCFunc(PetscInt, PetscReal, const PetscReal *, PetscInt,
                              PetscScalar *, void *);
BCFunc BCBend1_ss, BCBend2_ss, BCMMS;
BCFunc *boundaryOptions[] = {BCBend1_ss, BCBend2_ss, BCMMS};

// Multigrid options
typedef enum {
  MULTIGRID_LOGARITHMIC = 0, MULTIGRID_UNIFORM = 1, MULTIGRID_NONE = 2
} multigridType;
static const char *const multigridTypes [] = {"logarithmic",
                                              "uniform",
                                              "none",
                                              "multigridType","MULTIGRID",0
                                             };
static const char *const multigridTypesForDisp[] = {"P-multigrid, logarithmic coarsening",
                                                    "P-multigrind, uniform coarsening",
                                                    "No multigrid"};

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------
// Units
typedef struct Units_private *Units;
struct Units_private {
  // Fundamental units
  PetscScalar meter;
  PetscScalar kilogram;
  PetscScalar second;
  // Derived unit
  PetscScalar Pascal;
};

// Application context from user command line options
typedef struct {
  char          ceedresource[PETSC_MAX_PATH_LEN]; // libCEED backend
  char          meshFile[PETSC_MAX_PATH_LEN];     // exodusII mesh file
  PetscBool     testMode;
  problemType   problemChoice;
  forcingType   forcingChoice;
  boundaryType  boundaryChoice;
  multigridType multigridChoice;
  PetscInt      degree;
  PetscInt      numLevels;
  PetscInt      *levelDegrees;
  PetscInt      maxDiagState;
} AppCtx;

// Problem specific data
typedef struct {
  CeedInt           qdatasize;
  CeedQFunctionUser setupgeo, apply, jacob;
  const char        *setupgeofname, *applyfname, *jacobfname;
  CeedQuadMode      qmode;
} problemData;

// Data specific to each problem option
problemData problemOptions[3] = {
  [ELAS_LIN] = {
    .qdatasize = 10, // For linear elasticity, 6 would be sufficient
    .setupgeo = SetupGeo,
    .apply = LinElasF,
    .jacob = LinElasdF,
    .setupgeofname = SetupGeo_loc,
    .applyfname = LinElasF_loc,
    .jacobfname = LinElasdF_loc,
    .qmode = CEED_GAUSS
  },
  [ELAS_HYPER_SS] = {
    .qdatasize = 10,
    .setupgeo = SetupGeo,
    .apply = HyperSSF,
    .jacob = HyperSSdF,
    .setupgeofname = SetupGeo_loc,
    .applyfname = HyperSSF_loc,
    .jacobfname = HyperSSdF_loc,
    .qmode = CEED_GAUSS
  },
  [ELAS_HYPER_FS] = {
     .qdatasize = 10,
     .setupgeo = SetupGeo,
     .apply = HyperFSF,
     .jacob = HyperFSdF,
     .setupgeofname = SetupGeo_loc,
     .applyfname = HyperFSF_loc,
     .jacobfname = HyperFSdF_loc,
     .qmode = CEED_GAUSS
  }
};

// Forcing function data
typedef struct {
  CeedQFunctionUser setupforcing;
  const char        *setupforcingfname;
} forcingData;

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
  Vec          Xloc, Yloc, diagVec;
  PetscInt     diagState, maxDiagState;
  CeedVector   Xceed, Yceed;
  CeedOperator op;
  Ceed         ceed;
};

// Data for Jacobian setup routine
typedef struct FormJacobCtx_private *FormJacobCtx;
struct FormJacobCtx_private {
  UserMult     *jacobCtx;
  PetscInt     numLevels;
  SNES         snesCoarse;
  Mat          jacobMatMF, jacobMatCoarse;
  Vec          Ucoarse;
};

// Data for PETSc Prolongation/Restriction Matshell
typedef struct UserMultProlongRestr_private *UserMultProlongRestr;
struct UserMultProlongRestr_private {
  MPI_Comm     comm;
  DM           dmC, dmF;
  Vec          locVecC, locVecF, multVec;
  CeedVector   ceedVecC, ceedVecF;
  CeedOperator opProlong, opRestrict;
  Ceed         ceed;
};

// libCEED data struct for level
typedef struct CeedData_private *CeedData;
struct CeedData_private {
  Ceed                ceed;
  CeedBasis           basisx, basisu, basisCtoF;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictqdi, ErestrictGradui;
  CeedQFunction       qfApply, qfJacob;
  CeedOperator        opApply, opJacob, opRestrict, opProlong;
  CeedVector          qdata, gradu, xceed, yceed, truesoln;
};

// -----------------------------------------------------------------------------
// Process command line options
// -----------------------------------------------------------------------------
// Process general command line options
static int processCommandLineOptions(MPI_Comm comm, AppCtx *appCtx) {

  PetscErrorCode ierr;
  PetscBool meshFileFlag = PETSC_FALSE;
  PetscBool degreeFalg   = PETSC_FALSE;
  PetscBool boundaryFlag = PETSC_FALSE;
  PetscBool ceedFlag     = PETSC_FALSE;
  appCtx->problemChoice  = ELAS_LIN;       // Default - Linear Elasticity
  appCtx->degree         = 3;
  appCtx->boundaryChoice = BDRY_WALL_NONE; // Related to mesh choice
  appCtx->forcingChoice  = FORCE_NONE;     // Default - no forcing term
  appCtx->maxDiagState   = 1;              // Default - no diagonal reuse

  PetscFunctionBeginUser;

  ierr = PetscOptionsBegin(comm, NULL,
                           "Elasticity / Hyperelasticity in PETSc with libCEED",
                           NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-ceed", "CEED resource specifier",
                            NULL, appCtx->ceedresource, appCtx->ceedresource,
                            sizeof(appCtx->ceedresource), &ceedFlag);
  CHKERRQ(ierr);

  ierr = PetscOptionsInt("-degree", "Polynomial degree of tensor product basis",
                         NULL, appCtx->degree, &appCtx->degree,
                         &degreeFalg); CHKERRQ(ierr);

  ierr = PetscOptionsString("-mesh", "Read mesh from file", NULL,
                            appCtx->meshFile, appCtx->meshFile,
                            sizeof(appCtx->meshFile), &meshFileFlag);
  CHKERRQ(ierr);
  #if !defined(PETSC_HAVE_EXODUSII)
  SETERRQ(comm, PETSC_ERR_ARG_WRONG,
          "ExodusII support needed. Reconfigure your Arch with --download-exodusii");
  #endif

  ierr = PetscOptionsEnum("-problem",
                          "Solves Elasticity & Hyperelasticity Problems",
                          NULL, problemTypes, (PetscEnum)appCtx->problemChoice,
                          (PetscEnum *)&appCtx->problemChoice, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsEnum("-forcing", "Set forcing function option", NULL,
                          forcingTypes, (PetscEnum)appCtx->forcingChoice,
                          (PetscEnum *)&appCtx->forcingChoice, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsEnum("-boundary",
                          "Set Dirichlet (Essential) Boundary option", NULL,
                          boundaryTypes, (PetscEnum)appCtx->boundaryChoice,
                          (PetscEnum *)&appCtx->boundaryChoice, &boundaryFlag);
  CHKERRQ(ierr);

  appCtx->multigridChoice = MULTIGRID_LOGARITHMIC;
  ierr = PetscOptionsEnum("-multigrid", "Set multigrid type option", NULL,
                          multigridTypes, (PetscEnum)appCtx->multigridChoice,
                          (PetscEnum *)&appCtx->multigridChoice, NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsInt("-maxDiagState", "Set number of times to use Jacobian diagonal before recalculation",
                         NULL, appCtx->maxDiagState, &appCtx->maxDiagState,
                         NULL); CHKERRQ(ierr);

  appCtx->testMode = PETSC_FALSE;
  ierr = PetscOptionsBool("-test",
                          "Testing mode (do not print unless error is large)",
                          NULL, appCtx->testMode, &(appCtx->testMode), NULL);
  CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr); // End of setting AppCtx

  // Check for all required values set
  if (!appCtx->testMode) {
    if(!degreeFalg) {
      ierr = PetscPrintf(comm, "-degree option needed\n\n"); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
    }
    if(!meshFileFlag) {
      ierr = PetscPrintf(comm, "-mesh option needed (file)\n\n"); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
    }
    if(!boundaryFlag) {
      ierr = PetscPrintf(comm, "-boundary option needed\n\n"); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "AppCtx ERROR!");
    }
  } else {
    appCtx->boundaryChoice = BDRY_MMS;
    appCtx->forcingChoice = FORCE_MMS;
  }

  // Provide default ceed resource if not specified
  if (!ceedFlag) {
    const char* ceedresource = "/cpu/self";
    strncpy(appCtx->ceedresource, ceedresource, 10);
  }

  // Determine number of levels
  switch (appCtx->multigridChoice) {
  case MULTIGRID_LOGARITHMIC:
    appCtx->numLevels = ceil(log(appCtx->degree)/log(2)) + 1;
    break;
  case MULTIGRID_UNIFORM:
    appCtx->numLevels = appCtx->degree;
    break;
  case MULTIGRID_NONE:
    appCtx->numLevels = 1;
    break;
  }

  // Populate array of degrees for each level for multigrid
  ierr = PetscMalloc1(appCtx->numLevels, &(appCtx->levelDegrees));
  CHKERRQ(ierr);

  switch (appCtx->multigridChoice) {
  case MULTIGRID_LOGARITHMIC:
    for (int i=0; i<appCtx->numLevels-1; i++)
      appCtx->levelDegrees[i] = pow(2,i);
    appCtx->levelDegrees[appCtx->numLevels-1] = appCtx->degree;
    break;
  case MULTIGRID_UNIFORM:
    for (int i=0; i<appCtx->numLevels; i++)
      appCtx->levelDegrees[i] = i + 1;
    break;
  case MULTIGRID_NONE:
    appCtx->levelDegrees[0] = appCtx->degree;
    break;
  }

  PetscFunctionReturn(0);
};

// Process physics options
static int processPhysics(MPI_Comm comm, Physics phys, Units units) {

  PetscErrorCode ierr;
  PetscBool nuFlag = PETSC_FALSE;
  PetscBool YoungFlag = PETSC_FALSE;
  phys->nu = 0;
  phys->E = 0;
  units->meter     = 1;        // 1 meter in scaled length units
  units->second    = 1;        // 1 second in scaled time units
  units->kilogram  = 1;        // 1 kilogram in scaled mass units

  PetscFunctionBeginUser;

  ierr = PetscOptionsBegin(comm, NULL,
                           "Elasticity / Hyperelasticity in PETSc with libCEED",
                           NULL); CHKERRQ(ierr);

  ierr = PetscOptionsScalar("-nu", "Poisson's ratio", NULL, phys->nu, &phys->nu,
                            &nuFlag); CHKERRQ(ierr);

  ierr = PetscOptionsScalar("-E", "Young's Modulus", NULL, phys->E, &phys->E,
                            &YoungFlag); CHKERRQ(ierr);

  ierr = PetscOptionsScalar("-units_meter", "1 meter in scaled length units",
                            NULL, units->meter, &units->meter, NULL);
  CHKERRQ(ierr);
  units->meter = fabs(units->meter);

  ierr = PetscOptionsScalar("-units_second","1 second in scaled time units",
                            NULL, units->second, &units->second, NULL);
  CHKERRQ(ierr);
  units->second = fabs(units->second);

  ierr = PetscOptionsScalar("-units_kilogram","1 kilogram in scaled mass units",
                            NULL, units->kilogram, &units->kilogram, NULL);
  CHKERRQ(ierr);
  units->kilogram = fabs(units->kilogram);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr); // End of setting Physics

  // Check for all required options to be set
  if(!nuFlag) {
    ierr = PetscPrintf(comm, "-nu option needed\n\n"); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Poisson's ratio Error!");
  }
  if(!YoungFlag) {
    ierr = PetscPrintf(comm, "-E option needed\n\n"); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Young's Modulus Error!");
  }

  // Define derived units
  units->Pascal = units->kilogram / (units->meter * PetscSqr(units->second));

  // Scale E to GPa
  phys->E *= units->Pascal;

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Setup DM
// -----------------------------------------------------------------------------
static PetscErrorCode CreateBCLabel(DM dm, const char name[]) {
  int ierr;
  DMLabel label;

  PetscFunctionBeginUser;

  ierr = DMCreateLabel(dm, name); CHKERRQ(ierr);
  ierr = DMGetLabel(dm, name, &label); CHKERRQ(ierr);
  ierr = DMPlexMarkBoundaryFaces(dm, 1, label); CHKERRQ(ierr);
  ierr = DMPlexLabelComplete(dm, label); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// Read mesh and distribute DM in parallel
static int createDistributedDM(MPI_Comm comm, AppCtx *ctx, DM *dm) {

  PetscErrorCode  ierr;
  const char      *filename = ctx->meshFile;
  // Note: interpolate if polynomial degree > 1
  PetscBool       interpolate = PETSC_TRUE;
  DM              distributedMesh = NULL;
  PetscPartitioner part;

  PetscFunctionBeginUser;

  // Read mesh
  if (ctx->degree >= 2)
    interpolate = PETSC_TRUE;

  if (ctx->testMode) {
    PetscInt dim = 3, cells[3] = {3, 3, 3};
    ierr = DMPlexCreateBoxMesh(comm, dim, PETSC_FALSE, cells, NULL,
                               NULL, NULL, interpolate, dm); CHKERRQ(ierr);
  } else {
    ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm); CHKERRQ(ierr);
  }

  // Distribute DM in parallel
  ierr = DMPlexGetPartitioner(*dm, &part); CHKERRQ(ierr);
  ierr = PetscPartitionerSetFromOptions(part); CHKERRQ(ierr);
  ierr = DMPlexDistribute(*dm, 0, NULL, &distributedMesh); CHKERRQ(ierr);
  if (distributedMesh) {
    ierr = DMDestroy(dm); CHKERRQ(ierr);
    *dm  = distributedMesh;
  }

  PetscFunctionReturn(0);
};

// Setup DM with FE space of apropriate degree
static int SetupDMByDegree(DM dm, AppCtx appCtx, PetscInt order,
                           PetscInt ncompu) {
  PetscErrorCode  ierr;
  PetscInt        dim;
  PetscSpace      P;
  PetscDualSpace  Q;
  DM              K;
  PetscFE         fe;
  PetscInt        quadPointsPerEdge;
  PetscQuadrature q;  // quadrature points
  PetscQuadrature fq; // face quadrature points (For future: Nuemman boundary)
  // Variables for Dirichlet (Essential) Boundary
  IS              faceSetIS;           // Index Set for Face Sets
  const char      *name = "Face Sets"; // PETSc internal requirement
  PetscInt        numFaceSets;         // Number of FaceSets in faceSetIS
  const PetscInt  *faceSetIds;         // id of each FaceSet

  PetscFunctionBeginUser;

  ierr = DMGetDimension(dm, &dim);

  // Setup FE space (Space P) for tensor polynomials
  ierr = PetscSpaceCreate(PetscObjectComm((PetscObject) dm), &P); CHKERRQ(ierr);
  ierr = PetscSpacePolynomialSetTensor(P, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscSpaceSetFromOptions(P); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumComponents(P, ncompu); CHKERRQ(ierr);
  ierr = PetscSpaceSetNumVariables(P, dim); CHKERRQ(ierr);
  ierr = PetscSpaceSetDegree(P, order, order); CHKERRQ(ierr);
  ierr = PetscSpaceSetUp(P); CHKERRQ(ierr);
  // Setup FE dual space (Space Q) for tensor polynomials
  ierr = PetscDualSpaceCreate(PetscObjectComm((PetscObject) dm), &Q);
  CHKERRQ(ierr);
  ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE); CHKERRQ(ierr);
  ierr = PetscDualSpaceCreateReferenceCell(Q, dim, PETSC_FALSE, &K);
  CHKERRQ(ierr);
  ierr = PetscDualSpaceSetDM(Q, K); CHKERRQ(ierr);
  ierr = DMDestroy(&K); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetNumComponents(Q, ncompu); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetOrder(Q, order); CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetTensor(Q, PETSC_TRUE); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetFromOptions(Q); CHKERRQ(ierr);
  ierr = PetscDualSpaceSetUp(Q); CHKERRQ(ierr);
  // Create element
  ierr = PetscFECreate(PetscObjectComm((PetscObject) dm), &fe); CHKERRQ(ierr);
  ierr = PetscFESetFromOptions(fe); CHKERRQ(ierr);
  ierr = PetscFESetBasisSpace(fe, P); CHKERRQ(ierr);
  ierr = PetscFESetDualSpace(fe, Q); CHKERRQ(ierr);
  ierr = PetscFESetNumComponents(fe, ncompu); CHKERRQ(ierr);
  ierr = PetscFESetUp(fe); CHKERRQ(ierr);
  ierr = PetscSpaceDestroy(&P); CHKERRQ(ierr);
  ierr = PetscDualSpaceDestroy(&Q); CHKERRQ(ierr);
  // Create quadrature
  quadPointsPerEdge = PetscMax(order + 1,1);
  ierr = PetscDTGaussTensorQuadrature(dim, 1, quadPointsPerEdge, -1.0, 1.0, &q);
  CHKERRQ(ierr);
  ierr = PetscDTGaussTensorQuadrature(dim-1, 1, quadPointsPerEdge, -1.0, 1.0,
                                      &fq); CHKERRQ(ierr);
  ierr = PetscFESetQuadrature(fe, q); CHKERRQ(ierr);
  ierr = PetscFESetFaceQuadrature(fe, fq); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&q); CHKERRQ(ierr);
  ierr = PetscQuadratureDestroy(&fq); CHKERRQ(ierr);
  // Setup DM
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMAddField(dm, NULL, (PetscObject)fe); CHKERRQ(ierr);

  // Add Dirichlet (Essential) boundray
  ierr = DMCreateDS(dm); CHKERRQ(ierr);
  if (appCtx.testMode) {
    // Test mode - box mesh
    PetscBool hasLabel;
    PetscInt marker_ids[1] = {1};
    DMHasLabel(dm, "marker", &hasLabel);
    if (!hasLabel)
      CreateBCLabel(dm, "marker");
    ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL,
                         (void(*)(void))boundaryOptions[appCtx.boundaryChoice],
                         1, marker_ids, NULL);
    CHKERRQ(ierr);
  } else {
    // ExodusII mesh
    ierr = DMGetLabelIdIS(dm, name, &faceSetIS); CHKERRQ(ierr);
    ierr = ISGetLocalSize(faceSetIS,&numFaceSets);
    ierr = ISGetIndices(faceSetIS, &faceSetIds); CHKERRQ(ierr);

    for (PetscInt i=0; i<numFaceSets; i++) {
      ierr = DMAddBoundary(dm,DM_BC_ESSENTIAL,NULL,"Face Sets",0,0,NULL,
                           (void(*)(void))boundaryOptions[appCtx.boundaryChoice],
                           1, &faceSetIds[i], (void *)(&faceSetIds[i]));
      CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(faceSetIS, &faceSetIds); CHKERRQ(ierr);
    ierr = ISDestroy(&faceSetIS); CHKERRQ(ierr);
  }
  ierr = DMPlexSetClosurePermutationTensor(dm, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);
  ierr = PetscFEDestroy(&fe); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// libCEED Functions
// -----------------------------------------------------------------------------
// Destroy libCEED objects
static PetscErrorCode CeedDataDestroy(CeedInt level, CeedData data) {
  PetscInt ierr;

  // Vectors
  CeedVectorDestroy(&data->qdata);
  CeedVectorDestroy(&data->gradu);
  CeedVectorDestroy(&data->xceed);
  CeedVectorDestroy(&data->yceed);
  CeedVectorDestroy(&data->truesoln);

  // Restrictions
  CeedElemRestrictionDestroy(&data->Erestrictu);
  CeedElemRestrictionDestroy(&data->Erestrictx);
  CeedElemRestrictionDestroy(&data->ErestrictGradui);
  CeedElemRestrictionDestroy(&data->Erestrictqdi);

  // Bases
  CeedBasisDestroy(&data->basisx);
  CeedBasisDestroy(&data->basisu);

  // QFunctions
  CeedQFunctionDestroy(&data->qfJacob);
  CeedQFunctionDestroy(&data->qfApply);

  // Operators
  CeedOperatorDestroy(&data->opJacob);
  CeedOperatorDestroy(&data->opApply);

  // Restriction and Prolongation data
  CeedBasisDestroy(&data->basisCtoF);
  CeedOperatorDestroy(&data->opProlong);
  CeedOperatorDestroy(&data->opRestrict);

  ierr = PetscFree(data); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// Get libCEED restriction data from DMPlex
static int CreateRestrictionPlex(Ceed ceed, CeedInterlaceMode imode, CeedInt P,
                                 CeedInt ncomp, CeedElemRestriction *Erestrict,
                                 DM dm) {
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
  CeedElemRestrictionCreate(ceed, imode, nelem, P*P*P, nnodes/ncomp, ncomp,
                            CEED_MEM_HOST, CEED_COPY_VALUES, erestrict,
                            Erestrict);
  ierr = PetscFree(erestrict); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// Set up libCEED for a given degree
static int SetupLibceedFineLevel(DM dm, Ceed ceed, AppCtx appCtx, Physics phys,
                                 CeedData *data, PetscInt fineLevel,
                                 PetscInt ncompu, PetscInt Ugsz, PetscInt Ulocsz,
                                 CeedVector forceCeed, CeedQFunction qfRestrict,
                                 CeedQFunction qfProlong) {
  int           ierr;
  CeedInt       P = appCtx.levelDegrees[fineLevel] + 1;
  CeedInt       Q = appCtx.levelDegrees[fineLevel] + 1;
  CeedInt       dim, ncompx;
  CeedInt       nqpts;
  CeedInt       qdatasize = problemOptions[appCtx.problemChoice].qdatasize;
  problemType   problemChoice = appCtx.problemChoice;
  forcingType   forcingChoice = appCtx.forcingChoice;
  DM            dmcoord;
  PetscSection  section;
  Vec           coords;
  PetscInt      cStart, cEnd, nelem;
  const PetscScalar *coordArray;
  CeedVector    xcoord;
  CeedBasis     basisx;
  CeedQFunction qfSetupGeo, qfApply;
  CeedOperator  opSetupGeo, opApply;

  PetscFunctionBeginUser;

  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ncompx = dim;

  // ---------------------------------------------------------------------------
  // libCEED restrictions
  // ---------------------------------------------------------------------------
  ierr = DMGetCoordinateDM(dm, &dmcoord); CHKERRQ(ierr);
  ierr = DMPlexSetClosurePermutationTensor(dmcoord, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);

  // -- Coordinate restriction
  ierr = CreateRestrictionPlex(ceed, CEED_INTERLACED, 2, ncompx,
                               &(data[fineLevel]->Erestrictx), dmcoord);
  CHKERRQ(ierr);

  // -- Solution restriction
  ierr = CreateRestrictionPlex(ceed, CEED_INTERLACED, P, ncompu,
                               &data[fineLevel]->Erestrictu, dm); CHKERRQ(ierr);

  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
  nelem = cEnd - cStart;

  // -- Geometric data restriction
  CeedElemRestrictionCreateStrided(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, qdatasize,
                                   CEED_STRIDES_BACKEND,
                                   &data[fineLevel]->Erestrictqdi);
  // -- State vector gradient restriction
  if (problemChoice != ELAS_LIN)
    CeedElemRestrictionCreateStrided(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q,
                                     dim*ncompu, CEED_STRIDES_BACKEND,
                                     &data[fineLevel]->ErestrictGradui);

  // ---------------------------------------------------------------------------
  // Element coordinates
  // ---------------------------------------------------------------------------
  ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
  ierr = VecGetArrayRead(coords, &coordArray); CHKERRQ(ierr);
  ierr = DMGetSection(dmcoord, &section); CHKERRQ(ierr);

  CeedElemRestrictionCreateVector(data[fineLevel]->Erestrictx, &xcoord, NULL);
  CeedVectorSetArray(xcoord, CEED_MEM_HOST, CEED_COPY_VALUES,
                     (PetscScalar *)coordArray);
  ierr = VecRestoreArrayRead(coords, &coordArray); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // libCEED bases
  // ---------------------------------------------------------------------------
  // -- Solution basis
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu, P, Q,
                                  problemOptions[problemChoice].qmode,
                                  &data[fineLevel]->basisu);
  // -- Coordinate basis
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, Q,
                                  problemOptions[problemChoice].qmode,
                                  &data[fineLevel]->basisx);

  // ---------------------------------------------------------------------------
  // Persistent libCEED vectors
  // ---------------------------------------------------------------------------
  CeedBasisGetNumQuadraturePoints(data[fineLevel]->basisu, &nqpts);
  // -- Geometric data vector
  CeedVectorCreate(ceed, qdatasize*nelem*nqpts, &data[fineLevel]->qdata);
  // -- State gradient vector
  if (problemChoice != ELAS_LIN)
    CeedVectorCreate(ceed, dim*ncompu*nelem*nqpts, &data[fineLevel]->gradu);

  // ---------------------------------------------------------------------------
  // Geometric factor computation
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the quadrature data
  //   qdata returns dXdx_i,j and w * det.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].setupgeo,
                              problemOptions[problemChoice].setupgeofname,
                              &qfSetupGeo);
  CeedQFunctionAddInput(qfSetupGeo, "dx", ncompx*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qfSetupGeo, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(qfSetupGeo, "qdata", qdatasize, CEED_EVAL_NONE);

  // -- Operator
  CeedOperatorCreate(ceed, qfSetupGeo, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &opSetupGeo);
  CeedOperatorSetField(opSetupGeo, "dx", data[fineLevel]->Erestrictx,
                       data[fineLevel]->basisx, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(opSetupGeo, "weight", CEED_ELEMRESTRICTION_NONE,
                       data[fineLevel]->basisx, CEED_VECTOR_NONE);
  CeedOperatorSetField(opSetupGeo, "qdata", data[fineLevel]->Erestrictqdi,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

  // -- Compute the quadrature data
  CeedOperatorApply(opSetupGeo, xcoord, data[fineLevel]->qdata,
                    CEED_REQUEST_IMMEDIATE);

  // -- Cleanup
  CeedQFunctionDestroy(&qfSetupGeo);
  CeedOperatorDestroy(&opSetupGeo);

  // ---------------------------------------------------------------------------
  // Local residual evaluator
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the residual of the
  //   non-linear PDE.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].apply,
                              problemOptions[problemChoice].applyfname,
                              &qfApply);
  CeedQFunctionAddInput(qfApply, "du", ncompu*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qfApply, "qdata", qdatasize, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qfApply, "dv", ncompu*dim, CEED_EVAL_GRAD);
  if (problemChoice != ELAS_LIN)
    CeedQFunctionAddOutput(qfApply, "gradu", ncompu*dim, CEED_EVAL_NONE);
  CeedQFunctionSetContext(qfApply, phys, sizeof(phys));

  // -- Operator
  CeedOperatorCreate(ceed, qfApply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &opApply);
  CeedOperatorSetField(opApply, "du", data[fineLevel]->Erestrictu,
                       data[fineLevel]->basisu, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(opApply, "qdata", data[fineLevel]->Erestrictqdi,
                       CEED_BASIS_COLLOCATED, data[fineLevel]->qdata);
  CeedOperatorSetField(opApply, "dv", data[fineLevel]->Erestrictu,
                       data[fineLevel]->basisu, CEED_VECTOR_ACTIVE);
  if (problemChoice != ELAS_LIN)
    CeedOperatorSetField(opApply, "gradu", data[fineLevel]->ErestrictGradui,
                         data[fineLevel]->basisu, data[fineLevel]->gradu);
  // -- Save libCEED data
  data[fineLevel]->qfApply = qfApply;
  data[fineLevel]->opApply = opApply;

  // ---------------------------------------------------------------------------
  // Forcing term, if needed
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the forcing term (RHS)
  //   for the non-linear PDE.
  // ---------------------------------------------------------------------------
  if (forcingChoice != FORCE_NONE) {
    CeedQFunction qfSetupForce;
    CeedOperator opSetupForce;

    // -- QFunction
    CeedQFunctionCreateInterior(ceed, 1,
                                forcingOptions[forcingChoice].setupforcing,
                                forcingOptions[forcingChoice].setupforcingfname,
                                &qfSetupForce);
    CeedQFunctionAddInput(qfSetupForce, "x", ncompx, CEED_EVAL_INTERP);
    CeedQFunctionAddInput(qfSetupForce, "qdata", qdatasize, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qfSetupForce, "force", ncompu, CEED_EVAL_INTERP);
    CeedQFunctionSetContext(qfSetupForce, phys, sizeof(phys));

    // -- Operator
    CeedOperatorCreate(ceed, qfSetupForce, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &opSetupForce);
    CeedOperatorSetField(opSetupForce, "x", data[fineLevel]->Erestrictx,
                         data[fineLevel]->basisx, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(opSetupForce, "qdata", data[fineLevel]->Erestrictqdi,
                         CEED_BASIS_COLLOCATED, data[fineLevel]->qdata);
    CeedOperatorSetField(opSetupForce, "force", data[fineLevel]->Erestrictu,
                         data[fineLevel]->basisu, CEED_VECTOR_ACTIVE);

    // -- Compute forcing term
    CeedOperatorApply(opSetupForce, xcoord, forceCeed, CEED_REQUEST_IMMEDIATE);
    CeedVectorSyncArray(forceCeed, CEED_MEM_HOST);

    // -- Cleanup
    CeedQFunctionDestroy(&qfSetupForce);
    CeedOperatorDestroy(&opSetupForce);
  }

  // ---------------------------------------------------------------------------
  // True solution, for MMS
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the true solution at
  //   the mesh nodes for validation with the manufactured solution.
  // ---------------------------------------------------------------------------
  if (forcingChoice == FORCE_MMS) {
    CeedScalar *truearray;
    const CeedScalar *multarray;
    CeedVector multvec;
    CeedBasis basisxtrue;
    CeedQFunction qfTrue;
    CeedOperator opTrue;

    // -- Solution vector
    CeedVectorCreate(ceed, Ulocsz, &(data[fineLevel]->truesoln));

    // -- Basis
    CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, P, CEED_GAUSS_LOBATTO,
                                    &basisxtrue);

    // QFunction
    CeedQFunctionCreateInterior(ceed, 1, MMSTrueSoln, MMSTrueSoln_loc,
                                &qfTrue);
    CeedQFunctionAddInput(qfTrue, "x", ncompx, CEED_EVAL_INTERP);
    CeedQFunctionAddOutput(qfTrue, "true_soln", ncompu, CEED_EVAL_NONE);

    // Operator
    CeedOperatorCreate(ceed, qfTrue, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                       &opTrue);
    CeedOperatorSetField(opTrue, "x", data[fineLevel]->Erestrictx, basisxtrue,
                         CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(opTrue, "true_soln", data[fineLevel]->Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

    // -- Compute true solution
    CeedOperatorApply(opTrue, xcoord, data[fineLevel]->truesoln,
                      CEED_REQUEST_IMMEDIATE);

    // -- Multiplicity calculation
    CeedElemRestrictionCreateVector(data[fineLevel]->Erestrictu, &multvec, NULL);
    CeedVectorSetValue(multvec, 0.);
    CeedElemRestrictionGetMultiplicity(data[fineLevel]->Erestrictu, multvec);

    // -- Multiplicity correction
    CeedVectorGetArray(data[fineLevel]->truesoln, CEED_MEM_HOST, &truearray);
    CeedVectorGetArrayRead(multvec, CEED_MEM_HOST, &multarray);
    for (int i = 0; i < Ulocsz; i++)
      truearray[i] /= multarray[i];
    CeedVectorRestoreArray(data[fineLevel]->truesoln, &truearray);
    CeedVectorRestoreArrayRead(multvec, &multarray);

    // -- Cleanup
    CeedVectorDestroy(&multvec);
    CeedBasisDestroy(&basisxtrue);
    CeedQFunctionDestroy(&qfTrue);
    CeedOperatorDestroy(&opTrue);
  }

  // ---------------------------------------------------------------------------
  // Cleanup
  // ---------------------------------------------------------------------------
  CeedBasisDestroy(&basisx);
  CeedVectorDestroy(&xcoord);

  PetscFunctionReturn(0);
};

// Set up libCEED for a given degree
static int SetupLibceedLevel(DM dm, Ceed ceed, AppCtx appCtx, Physics phys,
                             CeedData *data, PetscInt level,
                             PetscInt ncompu, PetscInt Ugsz, PetscInt Ulocsz,
                             CeedVector forceCeed, CeedQFunction qfRestrict,
                             CeedQFunction qfProlong) {
  int           ierr;
  CeedInt       fineLevel = appCtx.numLevels - 1;
  CeedInt       P = appCtx.levelDegrees[level] + 1;
  CeedInt       Q = appCtx.levelDegrees[fineLevel] + 1;
  CeedInt       dim;
  CeedInt       qdatasize = problemOptions[appCtx.problemChoice].qdatasize;
  problemType   problemChoice = appCtx.problemChoice;
  CeedQFunction qfJacob;
  CeedOperator  opJacob, opProlong = NULL, opRestrict = NULL;

  PetscFunctionBeginUser;

  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // libCEED restrictions
  // ---------------------------------------------------------------------------
  if (level != fineLevel) {
    // -- Solution restriction
    ierr = CreateRestrictionPlex(ceed, CEED_INTERLACED, P, ncompu,
                                 &data[level]->Erestrictu, dm); CHKERRQ(ierr);
  }

  // ---------------------------------------------------------------------------
  // libCEED bases
  // ---------------------------------------------------------------------------
  // -- Solution basis
  if (level != fineLevel)
    CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu, P, Q,
                                    problemOptions[problemChoice].qmode,
                                    &data[level]->basisu);

  // -- Prolongation basis
  if (level != 0)
    CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu,
                                    appCtx.levelDegrees[level-1] + 1, P,
                                    CEED_GAUSS_LOBATTO,
                                    &data[level]->basisCtoF);

  // ---------------------------------------------------------------------------
  // Persistent libCEED vectors
  // ---------------------------------------------------------------------------
  CeedVectorCreate(ceed, Ulocsz, &data[level]->xceed);
  CeedVectorCreate(ceed, Ulocsz, &data[level]->yceed);

  // ---------------------------------------------------------------------------
  // Jacobian evaluator
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the action of the
  //   Jacobian for each linear solve.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].jacob,
                              problemOptions[problemChoice].jacobfname,
                              &qfJacob);
  CeedQFunctionAddInput(qfJacob, "deltadu", ncompu*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qfJacob, "qdata", qdatasize, CEED_EVAL_NONE);
  if (problemChoice != ELAS_LIN)
    CeedQFunctionAddInput(qfJacob, "gradu", ncompu*dim, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qfJacob, "deltadv", ncompu*dim, CEED_EVAL_GRAD);
  CeedQFunctionSetContext(qfJacob, phys, sizeof(phys));

  // -- Operator
  CeedOperatorCreate(ceed, qfJacob, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &opJacob);
  CeedOperatorSetField(opJacob, "deltadu", data[level]->Erestrictu,
                       data[level]->basisu, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(opJacob, "qdata", data[fineLevel]->Erestrictqdi,
                       CEED_BASIS_COLLOCATED, data[fineLevel]->qdata);
  CeedOperatorSetField(opJacob, "deltadv", data[level]->Erestrictu,
                       data[level]->basisu, CEED_VECTOR_ACTIVE);
  if (problemChoice != ELAS_LIN)
    CeedOperatorSetField(opJacob, "gradu", data[fineLevel]->ErestrictGradui,
                         CEED_BASIS_COLLOCATED, data[fineLevel]->gradu);

  // ---------------------------------------------------------------------------
  // Restriction and Prolongation
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the prolongation and
  //   restriction between the p-multigrid levels.
  // ---------------------------------------------------------------------------
  if ((level != 0) && appCtx.multigridChoice != MULTIGRID_NONE) {
    // -- Restriction
    CeedOperatorCreate(ceed, qfRestrict, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &opRestrict);
    CeedOperatorSetField(opRestrict, "input", data[level]->Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(opRestrict, "output", data[level-1]->Erestrictu,
                         data[level]->basisCtoF, CEED_VECTOR_ACTIVE);

    // -- Prolongation
    CeedOperatorCreate(ceed, qfProlong, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &opProlong);
    CeedOperatorSetField(opProlong, "input", data[level-1]->Erestrictu,
                         data[level]->basisCtoF, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(opProlong, "output", data[level]->Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  }

  // ---------------------------------------------------------------------------
  // Save libCEED data required for level
  // ---------------------------------------------------------------------------
  // -- QFunctions
  data[level]->qfJacob = qfJacob;

  // -- Operators
  data[level]->opJacob = opJacob;
  if (opProlong)
    data[level]->opProlong = opProlong;
  if (opRestrict)
    data[level]->opRestrict = opRestrict;

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Apply libCEED Operators
// -----------------------------------------------------------------------------
// This function returns the computed diagonal of the operator
static PetscErrorCode GetDiag_Ceed(Mat A, Vec D) {
  PetscErrorCode ierr;
  UserMult user;

  PetscFunctionBeginUser;

  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Recompute diagonal if it is too stale
  if (user->diagState > user->maxDiagState) {
    CeedVector ceedDiagVec;
    const CeedScalar *diagArray;

    // -- Compute Diagonal
    CeedOperatorAssembleLinearDiagonal(user->op, &ceedDiagVec,
                                       CEED_REQUEST_IMMEDIATE);

    // -- Place in PETSc vector
    CeedVectorGetArrayRead(ceedDiagVec, CEED_MEM_HOST, &diagArray);
    ierr = VecPlaceArray(user->Xloc, diagArray); CHKERRQ(ierr);

    // -- Local-to-Global
    ierr = VecZeroEntries(user->diagVec); CHKERRQ(ierr);
    ierr = DMLocalToGlobal(user->dm, user->Xloc, ADD_VALUES, user->diagVec);
    CHKERRQ(ierr);

    // -- Cleanup
    ierr = VecResetArray(user->Xloc); CHKERRQ(ierr);
    CeedVectorRestoreArrayRead(ceedDiagVec, &diagArray);
    CeedVectorDestroy(&ceedDiagVec);  

    // -- Update state
    user->diagState = 0;  
  }

  // Copy diagonal
  ierr = VecCopy(user->diagVec, D); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// This function uses libCEED to compute the local action of an operator
static PetscErrorCode ApplyLocalCeedOp(Vec X, Vec Y, UserMult user) {
  PetscErrorCode ierr;
  PetscScalar *x, *y;

  PetscFunctionBeginUser;

  // Global-to-local
  ierr = DMGlobalToLocal(user->dm, X, INSERT_VALUES, user->Xloc); CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Yloc); CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  ierr = VecGetArray(user->Yloc, &y); CHKERRQ(ierr);
  CeedVectorSetArray(user->Xceed, CEED_MEM_HOST, CEED_USE_POINTER, x);
  CeedVectorSetArray(user->Yceed, CEED_MEM_HOST, CEED_USE_POINTER, y);

  // Apply CEED operator
  CeedOperatorApply(user->op, user->Xceed, user->Yceed, CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->Yceed, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobal(user->dm, user->Yloc, ADD_VALUES, Y); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the non-linear residual
static PetscErrorCode FormResidual_Ceed(SNES snes, Vec X, Vec Y, void *ctx) {
  PetscErrorCode ierr;
  UserMult user = (UserMult)ctx;

  PetscFunctionBeginUser;

  // Use computed BCs
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);
  ierr = DMPlexInsertBoundaryValues(user->dm, PETSC_TRUE, user->Xloc, 0, NULL,
                                    NULL, NULL); CHKERRQ(ierr);

  // libCEED for local action of residual evaluator
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the non-linear residual
static PetscErrorCode ApplyJacobianCoarse_Ceed(SNES snes, Vec X, Vec Y,
                                               void *ctx) {
  PetscErrorCode ierr;
  UserMult user = (UserMult)ctx;

  PetscFunctionBeginUser;

  // Use computed BCs
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);

  // libCEED for local action of residual evaluator
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the action of the Jacobian
static PetscErrorCode ApplyJacobian_Ceed(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserMult user;

  PetscFunctionBeginUser;

  // Zero boundary values
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);

  // libCEED for local action of Jacobian
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the action of the prolongation operator
static PetscErrorCode Prolong_Ceed(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserMultProlongRestr user;
  PetscScalar *c, *f;

  PetscFunctionBeginUser;

  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Global-to-local
  ierr = VecZeroEntries(user->locVecC); CHKERRQ(ierr);
  ierr = DMGlobalToLocal(user->dmC, X, INSERT_VALUES, user->locVecC); 
  CHKERRQ(ierr);
  ierr = VecZeroEntries(user->locVecF); CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->locVecC, (const PetscScalar **)&c);
  CHKERRQ(ierr);
  ierr = VecGetArray(user->locVecF, &f); CHKERRQ(ierr);
  CeedVectorSetArray(user->ceedVecC, CEED_MEM_HOST, CEED_USE_POINTER, c);
  CeedVectorSetArray(user->ceedVecF, CEED_MEM_HOST, CEED_USE_POINTER, f);

  // Apply CEED operator
  CeedOperatorApply(user->opProlong, user->ceedVecC, user->ceedVecF,
                    CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->ceedVecF, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->locVecC, (const PetscScalar **)c);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->locVecF, &f); CHKERRQ(ierr);

  // Multiplicity
  ierr = VecPointwiseMult(user->locVecF, user->locVecF, user->multVec);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobal(user->dmF, user->locVecF, ADD_VALUES, Y);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// This function uses libCEED to compute the action of the restriction operator
static PetscErrorCode Restrict_Ceed(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserMultProlongRestr user;
  PetscScalar *c, *f;

  PetscFunctionBeginUser;

  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);

  // Global-to-local
  ierr = VecZeroEntries(user->locVecF); CHKERRQ(ierr);
  ierr = DMGlobalToLocal(user->dmF, X, INSERT_VALUES, user->locVecF);
  CHKERRQ(ierr);
  ierr = VecZeroEntries(user->locVecC); CHKERRQ(ierr);

  // Multiplicity
  ierr = VecPointwiseMult(user->locVecF, user->locVecF, user->multVec);
  CHKERRQ(ierr);

  // Setup CEED vectors
  ierr = VecGetArrayRead(user->locVecF, (const PetscScalar **)&f); CHKERRQ(ierr);
  ierr = VecGetArray(user->locVecC, &c); CHKERRQ(ierr);
  CeedVectorSetArray(user->ceedVecF, CEED_MEM_HOST, CEED_USE_POINTER, f);
  CeedVectorSetArray(user->ceedVecC, CEED_MEM_HOST, CEED_USE_POINTER, c);

  // Apply CEED operator
  CeedOperatorApply(user->opRestrict, user->ceedVecF, user->ceedVecC,
                    CEED_REQUEST_IMMEDIATE);
  CeedVectorSyncArray(user->ceedVecC, CEED_MEM_HOST);

  // Restore PETSc vectors
  ierr = VecRestoreArrayRead(user->locVecF, (const PetscScalar **)&f);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(user->locVecC, &c); CHKERRQ(ierr);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobal(user->dmC, user->locVecC, ADD_VALUES, Y);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Setup context data for Jacobian evaluation
static PetscErrorCode SetupJacobianCtx(MPI_Comm comm, AppCtx appCtx, DM dm,
    Vec V, Vec Vloc, CeedData ceedData, Ceed ceed, UserMult jacobianCtx) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // PETSc objects
  jacobianCtx->comm = comm;
  jacobianCtx->dm = dm;

  // Work vectors
  jacobianCtx->Xloc = Vloc;
  ierr = VecDuplicate(Vloc, &jacobianCtx->Yloc); CHKERRQ(ierr);
  jacobianCtx->Xceed = ceedData->xceed;
  jacobianCtx->Yceed = ceedData->yceed;

  // Diagonal vector
  ierr = VecDuplicate(V, &jacobianCtx->diagVec); CHKERRQ(ierr);
  jacobianCtx->diagState = appCtx.maxDiagState + 1; // Force initial assembly
  jacobianCtx->maxDiagState = appCtx.maxDiagState;

  // libCEED operator
  jacobianCtx->op = ceedData->opJacob;

  // Ceed
  jacobianCtx->ceed = ceed;

  PetscFunctionReturn(0);
};

// Setup context data for prolongation and restriction operators
static PetscErrorCode SetupProlongRestrictCtx(MPI_Comm comm, DM dmC, DM dmF,
    Vec VF, Vec VlocC, Vec VlocF, CeedData ceedDataC, CeedData ceedDataF,
    Ceed ceed, UserMultProlongRestr prolongRestrCtx) {
  PetscErrorCode ierr;
  PetscScalar *multArray;

  PetscFunctionBeginUser;

  // PETSc objects
  prolongRestrCtx->comm = comm;
  prolongRestrCtx->dmC = dmC;
  prolongRestrCtx->dmF = dmF;

  // Work vectors
  prolongRestrCtx->locVecC = VlocC;
  prolongRestrCtx->locVecF = VlocF;
  prolongRestrCtx->ceedVecC = ceedDataC->xceed;
  prolongRestrCtx->ceedVecF = ceedDataF->xceed;

  // libCEED operators
  prolongRestrCtx->opProlong = ceedDataF->opProlong;
  prolongRestrCtx->opRestrict = ceedDataF->opRestrict;

  // Ceed
  prolongRestrCtx->ceed = ceed;

  // Multiplicity vector
  // -- Set libCEED vector
  ierr = VecZeroEntries(VlocF);
  ierr = VecGetArray(VlocF, &multArray); CHKERRQ(ierr);
  CeedVectorSetArray(ceedDataF->xceed, CEED_MEM_HOST, CEED_USE_POINTER,
                     multArray);

  // -- Compute multiplicity
  CeedElemRestrictionGetMultiplicity(ceedDataF->Erestrictu, ceedDataF->xceed);

  // -- Restore PETSc vector
  ierr = VecRestoreArray(VlocF, &multArray); CHKERRQ(ierr);

  // -- Local-to-global
  ierr = VecZeroEntries(VF); CHKERRQ(ierr);
  ierr = DMLocalToGlobal(dmF, VlocF, ADD_VALUES, VF); CHKERRQ(ierr);

  // -- Global-to-local
  ierr = VecDuplicate(VlocF, &prolongRestrCtx->multVec); CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dmF, VF, INSERT_VALUES, prolongRestrCtx->multVec);
  CHKERRQ(ierr);

  // -- Reciprocal
  ierr = VecReciprocal(prolongRestrCtx->multVec); CHKERRQ(ierr);

  // -- Reset work arrays
  ierr = VecZeroEntries(VF); CHKERRQ(ierr);
  ierr = VecZeroEntries(VlocF); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Jacobian setup
// -----------------------------------------------------------------------------
static PetscErrorCode FormJacobian(SNES snes, Vec U, Mat J, Mat Jpre,
                                   void *ctx) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  // Context data
  FormJacobCtx  formJacobCtx = (FormJacobCtx)ctx;
  PetscInt      numLevels = formJacobCtx->numLevels;
  UserMult      *jacobCtx = formJacobCtx->jacobCtx;

  // Update diagonal state counter
  for (int level = 0; level < numLevels; level++)
    jacobCtx[level]->diagState++;

  // Form coarse assembled matrix
  ierr = VecZeroEntries(formJacobCtx->Ucoarse); CHKERRQ(ierr);
  ierr = SNESComputeJacobianDefaultColor(formJacobCtx->snesCoarse,
                                         formJacobCtx->Ucoarse,
                                         formJacobCtx->jacobMatMF,
                                         formJacobCtx->jacobMatCoarse, NULL);
  CHKERRQ(ierr);

  // Jpre might be AIJ (e.g., when using coloring), so we need to assemble it
  ierr = MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Boundary Functions
// -----------------------------------------------------------------------------
// BCMMS boundary function
// ss : (sideset)
// MMS: Boundary corresponding to the Method of Manufactured Solutions
// Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
// Also check ..\meshes\cyl-hol.8.jou
//
// left: sideset 999
// right: sideset 998
// outer: sideset 997
// inner: sideset 996
//
//   / \-------------------\              y
//  /   \                   \             |
// (  O  )                   )      x ____|
//  \   /                   /              \    Coordinate axis
//   \ /-------------------/                \ z
//
// Values on all points of the mesh is set based on given solution below
// for u[0], u[1], u[2]
PetscErrorCode BCMMS(PetscInt dim, PetscReal time, const PetscReal coords[],
                     PetscInt ncompu, PetscScalar *u, void *ctx) {
  PetscScalar x = coords[0];
  PetscScalar y = coords[1];
  PetscScalar z = coords[2];

  PetscFunctionBeginUser;

  u[0] = exp(2*x)*sin(3*y)*cos(4*z)/1e8;
  u[1] = exp(3*y)*sin(4*z)*cos(2*x)/1e8;
  u[2] = exp(4*z)*sin(2*x)*cos(3*y)/1e8;

  PetscFunctionReturn(0);
};

// BCBend2_ss boundary function
// ss : (sideset)
// 2_ss : two sides of the geometry
// Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
// Also check ..\meshes\cyl-hol.8.jou
//
// left: sideset 999
// right: sideset 998
//
//   / \-------------------\              y
//  /   \                   \             |
// (  O  )                   )      x ____|
//  \   /                   /              \    Coordinate axis
//   \ /-------------------/                \ z
//
//  0 values on the left side of the cyl-hole (sideset 999)
// -1 values on y direction of the right side of the cyl-hole (sideset 999)
PetscErrorCode BCBend2_ss(PetscInt dim, PetscReal time,
                          const PetscReal coords[],
                          PetscInt ncompu, PetscScalar *u, void *ctx) {
  PetscInt *faceID = (PetscInt *)ctx;

  PetscFunctionBeginUser;

  switch (*faceID) {
  case 999:    // left side of the cyl-hol
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
    break;
  case 998:    // right side of the cyl-hol
    u[0] = 0;
    u[1] = -1; // bend in the -y direction
    u[2] = 0;
    break;
  }
  PetscFunctionReturn(0);
};

// BCBend1_ss boundary function
// ss : (sideset)
// 1_ss : 1 side (left side) of the geometry
// Cylinder with a whole in the middle (see figure ..\meshes\surface999-9.png)
// Also check ..\meshes\cyl-hol.8.jou
//
// left: sideset 999
//
//   / \-------------------\              y
//  /   \                   \             |
// (  O  )                   )      x ____|
//  \   /                   /              \    Coordinate axis
//   \ /-------------------/                \ z
//
//  0 values on the left side of the cyl-hole (sideset 999)
PetscErrorCode BCBend1_ss(PetscInt dim, PetscReal time,
                          const PetscReal coords[],
                          PetscInt ncompu, PetscScalar *u, void *ctx) {
  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  PetscFunctionReturn(0);
};

#endif //setup_h
