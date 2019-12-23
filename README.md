## libCEED + PETSc Examples
This code is able to solve boundary value linear elasticity, hyperelasticity at small strain and 
hyperelasticity at finite strain (large deformation) using libCEED and PETSc. The hyperlesticity 
at finite strain formulation is Total Lagrangian.

### CEED/PETSc Linear Elascticity problem 

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E]`

Example: `./elasticity ./meshes/beamss16.exo  -degree 2 -nu .3 -E 10e6`

### CEED/PETSc Hyperelasticity at small strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperSS]`

Example: `./elasticity ./meshes/beamss16.exo  -degree 2 -nu .3 -E 10e6 -problem hyperSS`

### CEED/PETSc Hyperelasticity at small strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperFS]`

Example: `./elasticity ./meshes/beamss16.exo  -degree 2 -nu .3 -E 10e6 -problem hyperFS`
