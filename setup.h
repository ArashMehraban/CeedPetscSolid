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
  BDRY_WALL_NONE = 0, BDRY_WALL_WEIGHT = 1, BDRY_MMS = 2, BDRY_CUBE = 3
} boundaryType;
static const char *const boundaryTypes[] = {"wall_none",
                                            "wall_weight",
                                            "mms",
                                            "cube",
                                            "boundaryType","BDRY_",0};
static const char *const boundaryTypesForDisp[] = {"Wall with free end",
                                                   "Wall with weighted end",
                                                   "Manufactured solution",
                                                   "Cube"};

typedef PetscErrorCode BCFunc(PetscInt, PetscReal, const PetscReal *, PetscInt,
                              PetscScalar *, void *);
BCFunc BCBend1_ss, BCBend2_ss, BCMMS, BCCube;
BCFunc *boundaryOptions[] = {BCBend1_ss, BCBend2_ss, BCMMS, BCCube};

// Multigrid options
typedef enum {
  MULTIGRID_LOGARITHMIC = 0, MULTIGRID_UNIFORM = 1, MULTIGRID_NONE = 2
} multigridType;
static const char *const multigridTypes [] = {"logarithmic",
                                              "uniform",
                                              "none",
                                              "multigridType","MULTIGRID",0
                                             };
static const char *const multigridTypesForDisp[] = {"P-multigrid with logarithmic coarsening",
                                                    "P-multigrind with uniform coarsening",
                                                    "No multigrid"};

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------
// Units
typedef struct Units_private *Units;
struct Units_private {
  // fundamental units
  PetscScalar meter;
  PetscScalar kilogram;
  PetscScalar second;
  // derived unit
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
    .qdatasize = 10, // For linear Elasticty we could do 6
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
  Vec          Xloc, Yloc, diag;
  CeedVector   xceed, yceed;
  CeedOperator op;
  Ceed         ceed;
};

// libCEED data struct for level
typedef struct CeedData_private *CeedData;
struct CeedData_private {
  Ceed                ceed;
  CeedBasis           basisx, basisu, basisCtoF;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi,
                      Erestrictqdi, ErestrictGradui;
  CeedQFunction       qf_apply, qf_jacob;
  CeedOperator        op_apply, op_jacob, op_restrict, op_prolong;
  CeedVector          qdata, gradu, xceed, yceed, truesoln;
};

// -----------------------------------------------------------------------------
// Helper Functions
// -----------------------------------------------------------------------------
static int processCommandLineOptions(MPI_Comm comm, AppCtx *appCtx) {

  PetscErrorCode ierr;
  PetscBool meshFileFlag = PETSC_FALSE;
  PetscBool degreeFalg = PETSC_FALSE;
  PetscBool boundaryFlag = PETSC_FALSE;
  PetscBool ceedFlag = PETSC_FALSE;
  appCtx->problemChoice = ELAS_LIN; //-problem = Linear Elasticity if not given
  appCtx->degree = 3;
  appCtx->boundaryChoice = BDRY_WALL_NONE;
  appCtx->forcingChoice = FORCE_NONE;

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

static int createDistributedDM(MPI_Comm comm, AppCtx *ctx, DM *dm) {

  PetscErrorCode  ierr;
  const char      *filename = ctx->meshFile;
  PetscBool       interpolate = PETSC_TRUE;
  DM              distributedMesh = NULL;
  PetscPartitioner part;

  PetscFunctionBeginUser;
  if (ctx->degree >= 2)
    interpolate = PETSC_TRUE;

  if (ctx->testMode) {
    PetscInt dim = 3, cells[3] = {3, 3, 3};
    ierr = DMPlexCreateBoxMesh(comm, dim, PETSC_FALSE, cells, NULL,
                               NULL, NULL, interpolate, dm); CHKERRQ(ierr);
  } else {
    ierr = DMPlexCreateFromFile(comm, filename, interpolate, dm); CHKERRQ(ierr);
  }
  ierr = DMPlexGetPartitioner(*dm, &part); CHKERRQ(ierr);
  ierr = PetscPartitionerSetFromOptions(part); CHKERRQ(ierr);
  ierr = DMPlexDistribute(*dm, 0, NULL, &distributedMesh); CHKERRQ(ierr);
  if (distributedMesh) {
    ierr = DMDestroy(dm); CHKERRQ(ierr);
    *dm  = distributedMesh;
  }
  PetscFunctionReturn(0);
};

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
  // For Dirichlet (Essential) Boundary
  IS              faceSetIS;          // Index Set for Face Sets
  const char      *name="Face Sets";  // PETSc internal requirement
  PetscInt        numFaceSets;        // Number of FaceSets in faceSetIS
  const PetscInt  *faceSetIds;        // id of each FaceSet

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
// Destroy libCEED operator objects
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
  CeedElemRestrictionDestroy(&data->Erestrictxi);
  CeedElemRestrictionDestroy(&data->Erestrictqdi);

  // Bases
  CeedBasisDestroy(&data->basisx);
  CeedBasisDestroy(&data->basisu);

  // QFunctions
  CeedQFunctionDestroy(&data->qf_jacob);
  CeedQFunctionDestroy(&data->qf_apply);

  // Operators
  CeedOperatorDestroy(&data->op_jacob);
  CeedOperatorDestroy(&data->op_apply);

  // Restriction and Prolongation data
  if (level > 0) {
    CeedBasisDestroy(&data->basisCtoF);
    CeedOperatorDestroy(&data->op_prolong);
    CeedOperatorDestroy(&data->op_restrict);
  }

  ierr = PetscFree(data); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// Get CEED restriction data from DMPlex
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
static int SetupLibceedByDegree(DM dm, Ceed ceed, AppCtx appCtx, Physics phys,
                                CeedData *data, PetscInt level, PetscInt ncompu,
                                PetscInt Ugsz, PetscInt Ulocsz,
                                CeedVector forceceed, CeedQFunction qf_restrict,
                                CeedQFunction qf_prolong) {
  int          ierr;
  CeedInt      P = appCtx.levelDegrees[level] + 1;
  CeedInt      Q = P;
  CeedInt      dim, ncompx;
  CeedInt      nqpts;
  CeedInt      qdatasize = problemOptions[appCtx.problemChoice].qdatasize;
  problemType  problemChoice = appCtx.problemChoice;
  forcingType  forcingChoice = appCtx.forcingChoice;
  DM           dmcoord;
  PetscSection section;
  Vec          coords;
  PetscInt     cStart, cEnd, nelem;
  const PetscScalar *coordArray;
  CeedVector   xcoord, qdata, gradu, xceed, yceed;
  CeedElemRestriction Erestrictx, Erestrictu, Erestrictxi,
                      Erestrictqdi, ErestrictGradui;
  CeedBasis    basisx, basisu, basisCtoF;
  CeedQFunction qf_setupgeo, qf_apply, qf_jacob;
  CeedOperator op_setupgeo, op_apply, op_jacob;

  PetscFunctionBeginUser;

  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ncompx = dim;

  // ---------------------------------------------------------------------------
  // CEED restrictions
  // ---------------------------------------------------------------------------
  ierr = DMGetCoordinateDM(dm, &dmcoord); CHKERRQ(ierr);
  ierr = DMPlexSetClosurePermutationTensor(dmcoord, PETSC_DETERMINE, NULL);
  CHKERRQ(ierr);

  // -- Coordinate restriction
  ierr = CreateRestrictionPlex(ceed, CEED_INTERLACED, 2, ncompx, &Erestrictx,
                               dmcoord); CHKERRQ(ierr);
  // -- Solution restriction
  ierr = CreateRestrictionPlex(ceed, CEED_INTERLACED, P, ncompu, &Erestrictu,
                               dm); CHKERRQ(ierr);

  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd); CHKERRQ(ierr);
  nelem = cEnd - cStart;
  // -- Geometric data restriction
  CeedElemRestrictionCreateStrided(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, qdatasize,
                                   CEED_STRIDES_BACKEND, &Erestrictqdi);
  // -- Dummy weights restriction
  CeedElemRestrictionCreateStrided(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q, ncompx,
                                   CEED_STRIDES_BACKEND, &Erestrictxi);
  // -- State vector gradient restriction
  if (problemChoice != ELAS_LIN)
    CeedElemRestrictionCreateStrided(ceed, nelem, Q*Q*Q, nelem*Q*Q*Q,
                                     dim*ncompu, CEED_STRIDES_BACKEND,
                                     &ErestrictGradui);

  // ---------------------------------------------------------------------------
  // Element coordinates
  // ---------------------------------------------------------------------------
  ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
  ierr = VecGetArrayRead(coords, &coordArray); CHKERRQ(ierr);
  ierr = DMGetSection(dmcoord, &section); CHKERRQ(ierr);

  CeedElemRestrictionCreateVector(Erestrictx, &xcoord, NULL);
  CeedVectorSetArray(xcoord, CEED_MEM_HOST, CEED_COPY_VALUES,
                     (PetscScalar *)coordArray);
  ierr = VecRestoreArrayRead(coords, &coordArray); CHKERRQ(ierr);

  // ---------------------------------------------------------------------------
  // CEED bases
  // ---------------------------------------------------------------------------
  // -- Solution basis
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu, P, Q,
                                  problemOptions[problemChoice].qmode, &basisu);
  // -- Coordinate basis
  CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, Q,
                                  problemOptions[problemChoice].qmode, &basisx);
  // -- Prolongation basis
  if (level != 0)
    CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompu,
                                    appCtx.levelDegrees[level-1] + 1, P,
                                    CEED_GAUSS_LOBATTO, &basisCtoF);

  // ---------------------------------------------------------------------------
  // Persistent libCEED vectors
  // ---------------------------------------------------------------------------
  CeedBasisGetNumQuadraturePoints(basisu, &nqpts);
  // -- Geometric data vector
  CeedVectorCreate(ceed, qdatasize*nelem*nqpts, &qdata);
  // -- State gradient vector
  if (problemChoice != ELAS_LIN)
    CeedVectorCreate(ceed, dim*ncompu*nelem*nqpts, &gradu);
  // -- Work vectors
  CeedVectorCreate(ceed, Ulocsz, &xceed);
  CeedVectorCreate(ceed, Ulocsz, &yceed);

  // ---------------------------------------------------------------------------
  // Geometric factor computation
  // ---------------------------------------------------------------------------
  // Create the QFunction and Operator that computes the quadrature data
  //   qdata returns dXdx_i,j and w * det.
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].setupgeo,
                              problemOptions[problemChoice].setupgeofname,
                              &qf_setupgeo);
  CeedQFunctionAddInput(qf_setupgeo, "dx", ncompx*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_setupgeo, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(qf_setupgeo, "qdata", qdatasize, CEED_EVAL_NONE);

  // -- Operator
  CeedOperatorCreate(ceed, qf_setupgeo, CEED_QFUNCTION_NONE,
                     CEED_QFUNCTION_NONE, &op_setupgeo);
  CeedOperatorSetField(op_setupgeo, "dx", Erestrictx, basisx,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setupgeo, "weight", Erestrictxi, basisx,
                       CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setupgeo, "qdata", Erestrictqdi,
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

  // -- Compute the quadrature data
  CeedOperatorApply(op_setupgeo, xcoord, qdata, CEED_REQUEST_IMMEDIATE);

  // -- Cleanup
  CeedQFunctionDestroy(&qf_setupgeo);
  CeedOperatorDestroy(&op_setupgeo);

  // ---------------------------------------------------------------------------
  // Local residual evaluator
  // ---------------------------------------------------------------------------
  if (level == appCtx.numLevels - 1) {
    // -- QFunction
    CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].apply,
                                problemOptions[problemChoice].applyfname,
                                &qf_apply);
    CeedQFunctionAddInput(qf_apply, "du", ncompu*dim, CEED_EVAL_GRAD);
    CeedQFunctionAddInput(qf_apply, "qdata", qdatasize, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf_apply, "dv", ncompu*dim, CEED_EVAL_GRAD);
    if (problemChoice != ELAS_LIN)
      CeedQFunctionAddOutput(qf_apply, "gradu", ncompu*dim, CEED_EVAL_NONE);
    CeedQFunctionSetContext(qf_apply, phys, sizeof(phys));

    // -- Operator
    CeedOperatorCreate(ceed, qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                       &op_apply);
    CeedOperatorSetField(op_apply, "du", Erestrictu, basisu, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_apply, "qdata", Erestrictqdi, CEED_BASIS_COLLOCATED,
                         qdata);
    CeedOperatorSetField(op_apply, "dv", Erestrictu, basisu, CEED_VECTOR_ACTIVE);
    if (problemChoice != ELAS_LIN)
      CeedOperatorSetField(op_apply, "gradu", ErestrictGradui, basisu, gradu);
  }

  // ---------------------------------------------------------------------------
  // Jacobian evaluator
  // ---------------------------------------------------------------------------
  // -- QFunction
  CeedQFunctionCreateInterior(ceed, 1, problemOptions[problemChoice].jacob,
                              problemOptions[problemChoice].jacobfname,
                              &qf_jacob);
  CeedQFunctionAddInput(qf_jacob, "deltadu", ncompu*dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_jacob, "qdata", qdatasize, CEED_EVAL_NONE);
  if (problemChoice != ELAS_LIN)
    CeedQFunctionAddInput(qf_jacob, "gradu", ncompu*dim, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qf_jacob, "deltadv", ncompu*dim, CEED_EVAL_GRAD);
  CeedQFunctionSetContext(qf_jacob, phys, sizeof(phys));

  // -- Operator
  CeedOperatorCreate(ceed, qf_jacob, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                     &op_jacob);
  CeedOperatorSetField(op_jacob, "deltadu", Erestrictu, basisu,
                       CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_jacob, "qdata", Erestrictqdi, CEED_BASIS_COLLOCATED,
                       qdata);
  CeedOperatorSetField(op_jacob, "deltadv", Erestrictu, basisu,
                       CEED_VECTOR_ACTIVE);
  if (problemChoice != ELAS_LIN)
    CeedOperatorSetField(op_jacob, "gradu", ErestrictGradui,
                         CEED_BASIS_COLLOCATED, gradu);

  // ---------------------------------------------------------------------------
  // Forcing term, if needed
  // ---------------------------------------------------------------------------
  if ((level == appCtx.numLevels - 1) && (forcingChoice != FORCE_NONE)) {
    CeedQFunction qf_setupforce;
    CeedOperator op_setupforce;

    // -- QFunction
    CeedQFunctionCreateInterior(ceed, 1,
                                forcingOptions[forcingChoice].setupforcing,
                                forcingOptions[forcingChoice].setupforcingfname,
                                &qf_setupforce);
    CeedQFunctionAddInput(qf_setupforce, "x", ncompx, CEED_EVAL_INTERP);
    CeedQFunctionAddInput(qf_setupforce, "qdata", qdatasize, CEED_EVAL_NONE);
    CeedQFunctionAddOutput(qf_setupforce, "force", ncompu, CEED_EVAL_INTERP);
    CeedQFunctionSetContext(qf_setupforce, phys, sizeof(phys));

    // -- Operator
    CeedOperatorCreate(ceed, qf_setupforce, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_setupforce);
    CeedOperatorSetField(op_setupforce, "x", Erestrictx, basisx,
                         CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_setupforce, "qdata", Erestrictqdi,
                         CEED_BASIS_COLLOCATED, qdata);
    CeedOperatorSetField(op_setupforce, "force", Erestrictu, basisu,
                         CEED_VECTOR_ACTIVE);

    // -- Compute forcing term
    CeedOperatorApply(op_setupforce, xcoord, forceceed, CEED_REQUEST_IMMEDIATE);
    CeedVectorSyncArray(forceceed, CEED_MEM_HOST);

    // -- Cleanup
    CeedQFunctionDestroy(&qf_setupforce);
    CeedOperatorDestroy(&op_setupforce);
  }

  // ---------------------------------------------------------------------------
  // True solution, for MMS
  // ---------------------------------------------------------------------------
  if ((level == appCtx.numLevels - 1) && (forcingChoice == FORCE_MMS)) {
    CeedScalar *truearray;
    const CeedScalar *multarray;
    CeedVector multvec;
    CeedBasis basisxtrue;
    CeedQFunction qf_true;
    CeedOperator op_true;

    // -- Solution vector
    CeedVectorCreate(ceed, Ulocsz, &(data[level]->truesoln));

    // -- Basis
    CeedBasisCreateTensorH1Lagrange(ceed, dim, ncompx, 2, P, CEED_GAUSS_LOBATTO,
                                    &basisxtrue);

    // QFunction
    CeedQFunctionCreateInterior(ceed, 1, MMSTrueSoln, MMSTrueSoln_loc,
                                &qf_true);
    CeedQFunctionAddInput(qf_true, "x", ncompx, CEED_EVAL_INTERP);
    CeedQFunctionAddOutput(qf_true, "true_soln", ncompu, CEED_EVAL_NONE);

    // Operator
    CeedOperatorCreate(ceed, qf_true, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE,
                       &op_true);
    CeedOperatorSetField(op_true, "x", Erestrictx, basisxtrue,
                         CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_true, "true_soln", Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

    // -- Compute true solution
    CeedOperatorApply(op_true, xcoord, data[level]->truesoln,
                      CEED_REQUEST_IMMEDIATE);

    // -- Multiplicity calculation
    CeedElemRestrictionCreateVector(Erestrictu, &multvec, NULL);
    CeedVectorSetValue(multvec, 0.);
    CeedElemRestrictionGetMultiplicity(Erestrictu, multvec);

    // -- Multiplicity correction
    CeedVectorGetArray(data[level]->truesoln, CEED_MEM_HOST, &truearray);
    CeedVectorGetArrayRead(multvec, CEED_MEM_HOST, &multarray);
    for (int i = 0; i < Ulocsz; i++)
      truearray[i] /= multarray[i];
    CeedVectorRestoreArray(data[level]->truesoln, &truearray);
    CeedVectorRestoreArrayRead(multvec, &multarray);

    // -- Cleanup
    CeedVectorDestroy(&multvec);
    CeedBasisDestroy(&basisxtrue);
    CeedQFunctionDestroy(&qf_true);
    CeedOperatorDestroy(&op_true);
  }

  // ---------------------------------------------------------------------------
  // Restriction and Prolongation
  // ---------------------------------------------------------------------------
  if ((level != 0) && appCtx.multigridChoice != MULTIGRID_NONE) {
    CeedOperator op_restrict, op_prolong;

    // -- Restriction
    CeedOperatorCreate(ceed, qf_restrict, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_restrict);
    CeedOperatorSetField(op_restrict, "input", Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_restrict, "output", data[level-1]->Erestrictu,
                         basisCtoF, CEED_VECTOR_ACTIVE);

    // -- Prolongation
    CeedOperatorCreate(ceed, qf_prolong, CEED_QFUNCTION_NONE,
                       CEED_QFUNCTION_NONE, &op_prolong);
    CeedOperatorSetField(op_prolong, "input", data[level-1]->Erestrictu,
                         basisCtoF, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_prolong, "output", Erestrictu,
                         CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

    // -- Save libCEED data required for level
    data[level]->op_restrict = op_restrict;
    data[level]->op_prolong = op_prolong;
  }

  // ---------------------------------------------------------------------------
  // Cleanup
  // ---------------------------------------------------------------------------
  CeedVectorDestroy(&xcoord);

  // ---------------------------------------------------------------------------
  // Save libCEED data required for level
  // ---------------------------------------------------------------------------
  // -- Vectors
  data[level]->qdata = qdata;
  if (problemChoice != ELAS_LIN)
    data[level]->gradu = gradu;
  data[level]->xceed = xceed;
  data[level]->yceed = yceed;

  // -- Restrictions
  data[level]->Erestrictx = Erestrictx;
  data[level]->Erestrictu = Erestrictu;
  data[level]->Erestrictxi = Erestrictxi;
  data[level]->Erestrictqdi = Erestrictqdi;
  if (problemChoice != ELAS_LIN)
    data[level]->ErestrictGradui = ErestrictGradui;

  // -- Bases
  data[level]->basisx = basisx;
  data[level]->basisu = basisu;
  if (level > 0)
    data[level]->basisCtoF = basisCtoF;

  // -- QFunctions
  data[level]->qf_jacob = qf_jacob;
  if (level == appCtx.numLevels - 1)
    data[level]->qf_apply = qf_apply;

  // -- Operators
  data[level]->op_jacob = op_jacob;
  if (level == appCtx.numLevels - 1)
    data[level]->op_apply = op_apply;

  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Apply libCEED Operators
// -----------------------------------------------------------------------------
// This function uses libCEED to compute the local action of an operator
static PetscErrorCode ApplyLocalCeedOp(Vec X, Vec Y, UserMult user) {
  PetscErrorCode ierr;
  PetscScalar *x, *y;

  PetscFunctionBeginUser;
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
  ierr = VecRestoreArrayRead(user->Xloc, (const PetscScalar **)&x); CHKERRQ(ierr);
  ierr = VecRestoreArray(user->Yloc, &y); CHKERRQ(ierr);

  // Local-to-global
  ierr = VecZeroEntries(Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(user->dm, user->Yloc, ADD_VALUES, Y); CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(user->dm, user->Yloc, ADD_VALUES, Y); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the non-linear residual
static PetscErrorCode FormResidual_Ceed(SNES snes, Vec X, Vec Y, void *ptr) {
  PetscErrorCode ierr;
  UserMult user = (UserMult)ptr;

  PetscFunctionBeginUser;
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);
  ierr = DMPlexInsertBoundaryValues(user->dm, PETSC_TRUE, user->Xloc, 0, NULL,
                                    NULL, NULL); CHKERRQ(ierr);
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};

// This function uses libCEED to compute the action of the Jacobian
static PetscErrorCode ApplyJacobian_Ceed(Mat A, Vec X, Vec Y) {
  PetscErrorCode ierr;
  UserMult user;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);
  ierr = VecZeroEntries(user->Xloc); CHKERRQ(ierr);
  ierr = ApplyLocalCeedOp(X, Y, user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};
/*
// This function uses libCEED to compute the action of the prolongation operator
static PetscErrorCode MatMult_Prolong(Mat A, Vec X, Vec Y) {
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
*/
static PetscErrorCode CreateMatrixFreeCtx(MPI_Comm comm, DM dm, Vec vloc,
    CeedData ceeddata, Ceed ceed, UserMult residualCtx, UserMult jacobianCtx) {

  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  residualCtx->comm = comm;
  residualCtx->dm = dm;
  residualCtx->Xloc = vloc;
  ierr = VecDuplicate(vloc, &residualCtx->Yloc); CHKERRQ(ierr);
  residualCtx->xceed = ceeddata->xceed;
  residualCtx->yceed = ceeddata->yceed;
  residualCtx->op = ceeddata->op_apply;
  residualCtx->ceed = ceed;
  ierr = PetscMemcpy(jacobianCtx, residualCtx, sizeof(*residualCtx));
  CHKERRQ(ierr);
  jacobianCtx->op = ceeddata->op_jacob;
  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Jacobian setup
// -----------------------------------------------------------------------------
static PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat Jpre,void *ctx) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  // Nothing to do for J, which is of type MATSHELL and uses information from Newton.
  // Jpre might be AIJ (e.g., when using coloring), so we need to assemble it.
  ierr = MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
};

// -----------------------------------------------------------------------------
// Boundary Functions
// -----------------------------------------------------------------------------
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

  u[0] = exp(2*x)*sin(3*y)*cos(4*z)/1e8;
  u[1] = exp(3*y)*sin(4*z)*cos(2*x)/1e8;
  u[2] = exp(4*z)*sin(2*x)*cos(3*y)/1e8;

  PetscFunctionReturn(0);
};

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
PetscErrorCode BCBend2_ss(PetscInt dim, PetscReal time,
                          const PetscReal coords[],
                          PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscInt *faceID = (PetscInt *)ctx;
  PetscFunctionBeginUser;

  switch (*faceID) {
  case 999: // left side of the cyl-hol
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
    break;
  case 998: //right side of the cyl-hol
    u[0] = 0;
    u[1] = -1; //bend in the -y direction
    u[2] = 0;
    break;
  }
  PetscFunctionReturn(0);
};

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
PetscErrorCode BCBend1_ss(PetscInt dim, PetscReal time,
                          const PetscReal coords[],
                          PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  PetscFunctionReturn(0);
};

PetscErrorCode BCCube(PetscInt dim, PetscReal time,
                          const PetscReal coords[],
                          PetscInt ncompu, PetscScalar *u, void *ctx) {

  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  PetscFunctionReturn(0);
};
#endif //setup_h
