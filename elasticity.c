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
  AppCtx      appCtx;
  Physics     phys;
  DM          dm;
  PetscInt    ncompu = 3; // 3 dofs in 3D



  ierr = PetscInitialize(&argc, &argv, NULL, help);
  if(ierr)
      return ierr;

  comm = PETSC_COMM_WORLD;
  //set mesh-file, polynomial degree, problem type
  ierr = processCommandLineOptions(comm, &appCtx);CHKERRQ(ierr);
  //set Poison's ratio, Young's Modulus
  ierr = processPhysics(comm, &phys);CHKERRQ(ierr);
  //create distributed DM from mesh file (interpolate if polynomial degree > 1)
  ierr = createDistributedDM(comm, &appCtx, &dm);CHKERRQ(ierr);
  //setup DM by polinomial degree
  ierr = SetupDMByDegree(dm, &appCtx, ncompu);


  //Free objects
  ierr = DMDestroy(&dm); CHKERRQ(ierr);




  return PetscFinalize();
}
