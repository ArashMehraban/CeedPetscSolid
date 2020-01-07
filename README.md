## libCEED + PETSc Examples
This code is able to solve boundary value linear elasticity, hyperelasticity at small strain and
hyperelasticity at finite strain (large deformation) using libCEED and PETSc. The hyperlesticity
at finite strain formulation is Total Lagrangian. This code was tested in valgrind with the
following options: `--track-origins=yes` `--leak-check=full` `--show-leak-kinds=all`

### Valgrind Errors for this branch:
Valgrind verifies: All heap blocks were freed -- no leaks are possible

Valgrind says 6 error, Conditional jump or move depends on uninitialised value(s) in SNES (see err.delme file)

`valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 10e6 -boundary wall -forcing none`


### CEED/PETSc Linear Elasticity problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -boundary [boundary] -forcing [forcing]`

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_4ss_us.exo -degree 2 -nu .3 -E 10e6 -boundary wall -forcing manufactured`

See figure `\meshes\surface999-9.png`:

`4ss` means all 4 sides of `surface999-9.png`. `manufactured` is used with `4ss`

`2ss` means left and right sides of `surface999-9.png`. `none` or `constant` could be used as forcing functions


### CEED/PETSc Hyperelasticity at small strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperSS] -boundary [boundary] -forcing [forcing]`

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 10e6 -problem hyperSS
-boundary wall -forcing none`

See figure `\meshes\surface999-9.png`.

`2ss` means left and right sides of `surface999-9.png`. `none` or `constant` could be used as forcing functions  

### CEED/PETSc Hyperelasticity at finite strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperFS] -boundary [boundary] -forcing [forcing]`

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 10e6 -problem hyperFS
-boundary wall -forcing none`

See figure `\meshes\surface999-9.png`.

`1ss` means left side of `surface999-9.png`. `none` or `constant` could be used as forcing functions  
