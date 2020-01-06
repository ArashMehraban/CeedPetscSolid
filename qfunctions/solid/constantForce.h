#ifndef CONSTANT_H
#define CONSTANT_H

#ifndef __CUDACC__
#  include <math.h>
#endif

#ifndef PHYSICS_STRUCT
#define PHYSICS_STRUCT
typedef struct Physics_private *Physics;
struct Physics_private{
  PetscScalar   nu;      //Poisson's ratio
  PetscScalar   E;       //Young's Modulus
};
#endif

// -----------------------------------------------------------------------------
CEED_QFUNCTION(SetupConstantForce)(void *ctx, const CeedInt Q,
                                   const CeedScalar *const *in,
                                   CeedScalar *const *out) {
  // Inputs
  const CeedScalar *J = in[0], *w = in[1]; // *coords = in[2] - N/A for constant

  // Outputs
  CeedScalar *force = out[0]; // *true_soln = out[1] - N/A for constant force

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    const CeedScalar det = (J[i+Q*0]*(J[i+Q*4]*J[i+Q*8] - J[i+Q*5]*J[i+Q*7]) -
                            J[i+Q*1]*(J[i+Q*3]*J[i+Q*8] - J[i+Q*5]*J[i+Q*6]) +
                            J[i+Q*2]*(J[i+Q*3]*J[i+Q*7] - J[i+Q*4]*J[i+Q*6]));

    // Forcing function
    // -- Component 1
    force[i+0*Q] = det * w[i];
    // -- Component 2
    force[i+1*Q] = force[i+0*Q];
    // -- Component 3
    force[i+2*Q] = force[i+0*Q];
  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif
