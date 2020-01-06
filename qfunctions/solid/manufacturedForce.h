#ifndef MANUFACTURED_H
#define MANUFACTURED_H

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
CEED_QFUNCTION(SetupMMSForce)(void *ctx, const CeedInt Q,
                              const CeedScalar *const *in,
                              CeedScalar *const *out) {
  // Inputs
  const CeedScalar *J = in[0], *w = in[1], *coords = in[2];

  // Outputs
  CeedScalar *force = out[0], *true_soln = out[1];

  // Context
  const Physics context = ctx;
  const CeedScalar E  = context->E;
  const CeedScalar nu = context->nu;

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    CeedScalar x = coords[i+0*Q], y = coords[i+1*Q], z = coords[i+2*Q];
    const CeedScalar det = (J[i+Q*0]*(J[i+Q*4]*J[i+Q*8] - J[i+Q*5]*J[i+Q*7]) -
                            J[i+Q*1]*(J[i+Q*3]*J[i+Q*8] - J[i+Q*5]*J[i+Q*6]) +
                            J[i+Q*2]*(J[i+Q*3]*J[i+Q*7] - J[i+Q*4]*J[i+Q*6]));

    // Forcing function
    CeedScalar rho = det * w[i];
    // -- Component 1
    force[i+0*Q] = rho*(E*(cos(x*2)*cos(y*3)*exp(z*4)*4-cos(z*4)*sin(y*3)*exp(x*2)*8)*(nu-1/2))/((nu*2-1)*(nu+1))+(E*(cos(z*4)*sin(y*3)*exp(x*2)*(9/2)+sin(x*2)*sin(z*4)*exp(y*3)*3)*(nu-1/2))/((nu*2-1)*(nu+1))+(E*nu*cos(x*2)*cos(y*3)*exp(z*4)*8)/((nu*2-1)*(nu+1))-(E*nu*sin(x*2)*sin(z*4)*exp(y*3)*6)/((nu*2-1)*(nu+1))-(E*cos(z*4)*sin(y*3)*exp(x*2)*(nu-1)*4)/((nu*2-1)*(nu+1));
    // -- Component 2
    force[i+1*Q] =rho*(E*(cos(y*3)*cos(z*4)*exp(x*2)*3-cos(x*2)*sin(z*4)*exp(y*3)*2)*(nu-1/2))/((nu*2-1)*(nu+1))+(E*(cos(x*2)*sin(z*4)*exp(y*3)*8+sin(x*2)*sin(y*3)*exp(z*4)*6)*(nu-1/2))/((nu*2-1)*(nu+1))+(E*nu*cos(y*3)*cos(z*4)*exp(x*2)*6)/((nu*2-1)*(nu+1))-(E*nu*sin(x*2)*sin(y*3)*exp(z*4)*1.2e1)/((nu*2-1)*(nu+1))-(E*cos(x*2)*sin(z*4)*exp(y*3)*(nu-1)*9)/((nu*2-1)*(nu+1));
    // -- Component 3
    force[i+2*Q] = rho*(E*(cos(x*2)*cos(z*4)*exp(y*3)*6-cos(y*3)*sin(x*2)*exp(z*4)*(9/2))*(nu-1/2))/((nu*2-1)*(nu+1))+(E*(cos(y*3)*sin(x*2)*exp(z*4)*2+sin(y*3)*sin(z*4)*exp(x*2)*4)*(nu-1/2))/((nu*2-1)*(nu+1))+(E*nu*cos(x*2)*cos(z*4)*exp(y*3)*1.2e1)/((nu*2-1)*(nu+1))-(E*nu*sin(y*3)*sin(z*4)*exp(x*2)*8)/((nu*2-1)*(nu+1))-(E*cos(y*3)*sin(x*2)*exp(z*4)*(nu-1)*1.6e1)/((nu*2-1)*(nu+1));

    // True solution
    // -- Component 1
    true_soln[i+0*Q] = exp(2*x)*sin(3*y)*cos(4*z);
    // -- Component 2
    true_soln[i+1*Q] = exp(3*y)*sin(4*z)*cos(2*x);
    // -- Component 3
    true_soln[i+2*Q] = exp(4*z)*sin(2*x)*cos(3*y);
  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif
