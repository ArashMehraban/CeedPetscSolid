#ifndef MANUFACTURED_TRUE_H
#define MANUFACTURED_TRUE_H

#ifndef __CUDACC__
#  include <math.h>
#endif

// -----------------------------------------------------------------------------
CEED_QFUNCTION(MMSTrueSoln)(void *ctx, const CeedInt Q,
                            const CeedScalar *const *in,
                            CeedScalar *const *out) {
  // Inputs
  const CeedScalar *coords = in[0];

  // Outputs
  CeedScalar *true_soln = out[0];

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Setup
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q], z = coords[i+2*Q];

    // True solution
    // -- Component 1
    true_soln[i+0*Q] = exp(2*x)*sin(3*y)*cos(4*z)/1e8;

    // -- Component 2
    true_soln[i+1*Q] = exp(3*y)*sin(4*z)*cos(2*x)/1e8;

    // -- Component 3
    true_soln[i+2*Q] = exp(4*z)*sin(2*x)*cos(3*y)/1e8;

  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif // End MANUFACTURED_TRUE_H
