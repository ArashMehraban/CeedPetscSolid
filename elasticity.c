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
  AppCtx      appCtx; //contains degree, problem choice & mesh filename
  Physics     phys;   //contains nu and E
  DM          dm;
  PetscInt    ncompu = 3;  // 3 dofs in 3D
  Vec         U,Uloc, R,Rloc; //loc: Local R:Residual
  PetscInt    Usz,Ulsz;    // sz: size



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
  // Create Unknonw vector U and Residual Vector R
  ierr = DMCreateGlobalVector(dm, &U); CHKERRQ(ierr);
  ierr = VecGetLocalSize(U, &Ulsz); CHKERRQ(ierr);
  ierr = VecGetSize(U, &Usz); CHKERRQ(ierr);
  ierr = DMCreateLocalVector(dm, &Uloc); CHKERRQ(ierr);
  ierr = VecGetSize(Uloc, &Ulsz); CHKERRQ(ierr);
  ierr = VecDuplicate(U, &R); CHKERRQ(ierr);
  ierr = VecDuplicate(Uloc, &Rloc); CHKERRQ(ierr);
  ierr = VecZeroEntries(Rloc); CHKERRQ(ierr);


  //Free objects
  ierr = DMDestroy(&dm); CHKERRQ(ierr);
  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = VecDestroy(&Uloc); CHKERRQ(ierr);
  ierr = VecDestroy(&R); CHKERRQ(ierr);
  ierr = VecDestroy(&Rloc); CHKERRQ(ierr);




  return PetscFinalize();
}
