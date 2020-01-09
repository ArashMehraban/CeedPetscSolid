## libCEED + PETSc Examples
This code is able to solve boundary value linear elasticity, hyperelasticity at small strain and
hyperelasticity at finite strain (large deformation) using libCEED and PETSc. The hyperlesticity
at finite strain formulation is Total Lagrangian. This code was tested in valgrind with the
following options: `--track-origins=yes` `--leak-check=full` `--show-leak-kinds=all`

**Boundary:**

Setting boundary is mesh dependent in every FEM problem. As a result, the examples we have provided here depend on the  mesh files in `\mesh\` folder. However, this code is capable of importing any structred or unstructred ExodusII (.exo) mesh file. In such cases, the user is responsible for providing boundary functions in `setup.h`. We have used Trelis/Cubit software to generate mesh. The *journal file*, `.jou` file is provided in the `\meshes\` directory. We have employed the `sideset` feature from Trelis\Cubit software to choose different regions of the geometry. These regions are utilized in the boundary functions in `setup.h` to place *essential* (Dirichlet) boundary values in the sollution vector. `nodeset` must be avoided for the puposes of choosing boundary regions in the mesh as this code runs with high-order polynomials. Everything else about the code is general:

![Image of finger](https://github.com/ArashMehraban/CeedPetscSolid/pictures/finger.png)

**General Notes about mesh naming convention:**

`mms` stands for Method of Manufactured Solutions.\
In cyl-hole_632e_4ss_us.exo file name:\
   `_4ss` refers to the left, right, inner and outer *walls* of the image above.\
   `_2ss` refers to the left and right *walls* of the image above.\
   `_1ss` refers to the left *wall* of the image above.\
   `_632e` in the mesh file name means `632` elements.\
   `_us` means `unstructured mesh`

### CEED/PETSc Linear Elasticity problem

To build, run `make`

To run:\
 `./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -boundary [boundary] -forcing [forcing]`\
 or\
  `mpirun -n [n] ./elasticity -mesh [.exo file]  -degree [degree] -nu [nu] -E [E] -boundary [boundary] -forcing [forcing]`

In our case `mms` is based on the following contrived solution:

`u[0] = exp(2x)sin(3y)cos(4z)`\
`u[1] = exp(3y)sin(4z)cos(2x)`\
`u[2] = exp(4z)sin(2x)cos(3y)`

**Note 1:** For the `mms` to work correctly, you must use: \
            mesh files with `_4ss_` in their name from `\meshes` directory.\
            `-boundary mms` and `-forcing mms` options.

Example:\
 `./elasticity -mesh ./meshes/cyl-hole_632e_4ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary mms -forcing mms`

**Note 2:** Two other boundary and forcing functions may be used with this mesh files provided in `\meshes\`:

**1)** left side of the cyl-hol object is attached to a wall:\
       mesh files with `_1ss` must be used.\
       `-boundary wall_none` must be used.\
       forcing function on that could be `none` (no force) or `constant` (constant force in `-y` direction)

Example:
 `./elasticity -mesh ./meshes/cyl-hole_632e_1ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary wall_none -forcing constant`

**2)** left side of the cyl-hol object is attached to a wall **and** the right side of cyl-hole object has a dead wight hanging off it:\
   mesh files with `_2ss` must be used.\
   `-boundary wall_weight` must be used.\
   forcing function on that could be `none` (no force) or `constant` (constant force in `-y` direction)

Example:
 `./elasticity -mesh ./meshes/cyl-hole_632e_2ss_us.exo -degree 2 -nu .3 -E 1e6 -boundary wall_weight -forcing constant`

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
