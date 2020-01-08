## libCEED + PETSc Examples
This code is able to solve boundary value linear elasticity, hyperelasticity at small strain and
hyperelasticity at finite strain (large deformation) using libCEED and PETSc. The hyperlesticity
at finite strain formulation is Total Lagrangian. This code was tested in valgrind with the
following options: `--track-origins=yes` `--leak-check=full` `--show-leak-kinds=all`

### CEED/PETSc Linear Elasticity problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -boundary [boundary] -forcing [forcing]`

**Note 1:** `mms` stands for Method of Manufactured Solutions. In our case `mms` is:

`u[0] = exp(2x)sin(3y)cos(4z)`\
`u[1] = exp(3y)sin(4z)cos(2x)`\
`u[2] = exp(4z)sin(2x)cos(3y)`

**Note 2:** For the `mms` to work correctly with the meshes in `\meshes\`, mesh files with `_4ss_` must be used. Also, see figure `\meshes\surface999-9.png`. `4ss` refers to the left, right, inner and outer walls of the cylinder with a hole in that figure. (mesh journal file for `mms` is also provided in the same folder.)

Also, `mms` must be used with the `-boundary mms` and `-forcing mms` options.

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_4ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary mms -forcing mms`

`_632e` in the mesh file name means `632` elements and `_us` means `unstructured mesh`

**Note 3:** Two other boundary and forcing functions may be used with this code:

See figure `\meshes\surface999-9.png`.

**1)** left side of the cyl-hol object is attached to a wall:

mesh files with `_1ss` must be used.

`-boundary wall_none` must be used.

forcing function on that could be `none` (no force) or `constant` (constant force in `-y` direction)

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_1ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary wall_none -forcing constant`

**2)** left side of the cyl-hol object is attached to a wall **and** the right side of cyl-hole object has a dead wight on it:

mesh files with `_2ss` must be used.

`-boundary wall_weight` must be used.

forcing function on that could be `none` (no force) or `constant` (constant force in `-y` direction)

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary wall_weight -forcing constant`

### CEED/PETSc Hyperelasticity at small strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperSS] -boundary [boundary] -forcing [forcing]`

Example: `./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 1e6 -problem hyperSS
-boundary wall -forcing none`

See figure `\meshes\surface999-9.png`.

details will be updated when ready.

### CEED/PETSc Hyperelasticity at finite strain problem

To build, run `make`

To run, `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -problem [hyperFS] -boundary [boundary] -forcing [forcing]`

See figure `\meshes\surface999-9.png`.

details will be updated when ready.
