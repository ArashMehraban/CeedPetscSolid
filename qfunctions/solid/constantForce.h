#ifndef CONSTANT_H
#define CONSTANT_H

#ifndef __CUDACC__
#  include <math.h>
#endif

#ifndef PHYSICS_STRUCT
#define PHYSICS_STRUCT
typedef struct Physics_private *Physics;
struct Physics_private {
  PetscScalar   nu;      // Poisson's ratio
  PetscScalar   E;       // Young's Modulus
};
#endif

// -----------------------------------------------------------------------------
CEED_QFUNCTION(SetupConstantForce)(void *ctx, const CeedInt Q,
                                   const CeedScalar *const *in,
                                   CeedScalar *const *out) {
  // Inputs
  const CeedScalar *qdata = in[1];

  // Outputs
  CeedScalar *force = out[0];

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    CeedScalar rho = qdata[i];

    // Forcing function
    // -- Component 1
    force[i+0*Q] = 0;

    // -- Component 2
    force[i+1*Q] = -rho;

    // -- Component 3
    force[i+2*Q] = 0;

  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif // End of CONSTANT_H
