#ifndef MANUFACTURED_H
#define MANUFACTURED_H

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
CEED_QFUNCTION(SetupMMSForce)(void *ctx, const CeedInt Q,
                              const CeedScalar *const *in,
                              CeedScalar *const *out) {
  // Inputs
  const CeedScalar *coords = in[0], *qdata = in[1];

  // Outputs
  CeedScalar *force = out[0];

  // Context
  const Physics context = ctx;
  const CeedScalar E  = context->E;
  const CeedScalar nu = context->nu;

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q], z = coords[i+2*Q];
    CeedScalar wJ = qdata[i];

    // Forcing function
    force[i+0*Q] = wJ*(-(E*(cos(x*2.0)*cos(y*3.0)*exp(z*4.0)*4.0-cos(z*4.0)*sin(y*3.0)*exp(x*2.0)*8.0)*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*(cos(z*4.0)*sin(y*3.0)*exp(x*2.0)*(4.5)+sin(x*2.0)*sin(z*4.0)*exp(y*3.0)*3.0)*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*nu*cos(x*2.0)*cos(y*3.0)*exp(z*4.0)*8.0)/((nu*2.0-1.0)*(nu+1.0))-(E*nu*sin(x*2.0)*sin(z*4.0)*exp(y*3.0)*6.0)/((nu*2.0-1.0)*(nu+1.0))-
                        (E*cos(z*4.0)*sin(y*3.0)*exp(x*2.0)*(nu-1.0)*4.0)/((nu*2.0-1.0)*(nu+1.0)));

    force[i+1*Q] = wJ*(-(E*(cos(y*3.0)*cos(z*4.0)*exp(x*2.0)*3.0-cos(x*2.0)*sin(z*4.0)*exp(y*3.0)*2.0)*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*(cos(x*2.0)*sin(z*4.0)*exp(y*3.0)*8.0+sin(x*2.0)*sin(y*3.0)*exp(z*4.0)*6.0)*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*nu*cos(y*3.0)*cos(z*4.0)*exp(x*2.0)*6.0)/((nu*2.0-1.0)*(nu+1.0))-(E*nu*sin(x*2.0)*sin(y*3.0)*exp(z*4.0)*12.0)/((nu*2.0-1.0)*(nu+1.0))-
                        (E*cos(x*2.0)*sin(z*4.0)*exp(y*3.0)*(nu-1.0)*9.0)/((nu*2.0-1.0)*(nu+1.0)));

    force[i+2*Q] = wJ*(-(E*(cos(x*2.0)*cos(z*4.0)*exp(y*3.0)*6.0-cos(y*3.0)*sin(x*2.0)*exp(z*4.0)*(4.5))*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*(cos(y*3.0)*sin(x*2.0)*exp(z*4.0)*2.0+sin(y*3.0)*sin(z*4.0)*exp(x*2.0)*4.0)*(nu-0.5))/((nu*2.0-1.0)*(nu+1.0))+
                        (E*nu*cos(x*2.0)*cos(z*4.0)*exp(y*3.0)*12.0)/((nu*2.0-1.0)*(nu+1.0))-(E*nu*sin(y*3.0)*sin(z*4.0)*exp(x*2.0)*8.0)/((nu*2.0-1.0)*(nu+1.0))-
                        (E*cos(y*3.0)*sin(x*2.0)*exp(z*4.0)*(nu-1.0)*16.0)/((nu*2.0-1.0)*(nu+1.0)));
  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif
