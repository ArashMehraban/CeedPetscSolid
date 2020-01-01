// *****************************************************************************
    // This QFunction sets up the geometric factors required for integration and
    //   coordinate transformations
    //
    // Reference (parent) coordinates: X
    // Physical (current) coordinates: x
    // Change of coordinate matrix: dxdX_{i,j} = x_{i,j} (indicial notation)
    // Inverse of change of coordinate matrix: dXdx_{i,j} = (detJ^-1) * X_{i,j}
    //
    // All quadrature data is stored in 10 field vector of quadrature data.
    //
    // We require the transpose of the inverse of the Jacobian to properly compute
    //   integrals of the form: int( gradv u )
    //
    // Inverse of Jacobian:
    //   dXdx_i,j = Aij / detJ
    //
    // Stored: Aij / detJ
    //   in qdata[1:9] as
    //   (detJ^-1) * [A11 A12 A13]
    //               [A21 A22 A23]
    //               [A31 A32 A33]
    //
    // *****************************************************************************
#ifndef __Physics__
#define __Physics__

    typedef struct Physics_private *Physics;
    struct Physics_private{
      PetscScalar   nu;      //Poisson's ratio
      PetscScalar   E;       //Young's Modulus
    };
#endif



    CEED_QFUNCTION(SetupLinElasGeo)(void *ctx, CeedInt Q,
                          const CeedScalar *const *in, CeedScalar *const *out) {
      // *INDENT-OFF*
      // Inputs
      const CeedScalar (*J)[3][Q] = (CeedScalar(*)[3][Q])in[0],
                       (*w) = in[1];

      // Outputs
      CeedScalar (*qdata)[Q] = (CeedScalar(*)[Q])out[0];
      // *INDENT-ON*

 PetscPrintf(PETSC_COMM_WORLD, "Test: size of Q in SetupLinElasGeo (linElas.h file): %D\n\n",Q);
      CeedPragmaSIMD
      // Quadrature Point Loop
      for (CeedInt i=0; i<Q; i++) {
        // Setup
        const CeedScalar J11 = J[0][0][i];
        const CeedScalar J21 = J[0][1][i];
        const CeedScalar J31 = J[0][2][i];
        const CeedScalar J12 = J[1][0][i];
        const CeedScalar J22 = J[1][1][i];
        const CeedScalar J32 = J[1][2][i];
        const CeedScalar J13 = J[2][0][i];
        const CeedScalar J23 = J[2][1][i];
        const CeedScalar J33 = J[2][2][i];
        const CeedScalar A11 = J22*J33 - J23*J32;
        const CeedScalar A12 = J13*J32 - J12*J33;
        const CeedScalar A13 = J12*J23 - J13*J22;
        const CeedScalar A21 = J23*J31 - J21*J33;
        const CeedScalar A22 = J11*J33 - J13*J31;
        const CeedScalar A23 = J13*J21 - J11*J23;
        const CeedScalar A31 = J21*J32 - J22*J31;
        const CeedScalar A32 = J12*J31 - J11*J32;
        const CeedScalar A33 = J11*J22 - J12*J21;
        const CeedScalar detJ = J11*A11 + J21*A12 + J31*A13;

        // Qdata
        // -- Interp-to-Interp qdata
        qdata[0][i] = w[i] * detJ;
        // -- Interp-to-Grad qdata
        // Inverse of change of coordinate matrix: X_i,j
        qdata[1][i] = A11 / detJ;
        qdata[2][i] = A12 / detJ;
        qdata[3][i] = A13 / detJ;
        qdata[4][i] = A21 / detJ;
        qdata[5][i] = A22 / detJ;
        qdata[6][i] = A23 / detJ;
        qdata[7][i] = A31 / detJ;
        qdata[8][i] = A32 / detJ;
        qdata[9][i] = A33 / detJ;

      } // End of Quadrature Point Loop

      return 0;
    }
// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinElas)(void *ctx, CeedInt Q,
                         const CeedScalar *const *in,
                         CeedScalar *const *out) {
   // Inputs
   const CeedScalar *ug = in[0], *qdata = in[1];

   // Outputs
   CeedScalar *vg = out[0];

   // Context
    const Physics context = ctx;
   const CeedScalar E  = context->E;
   const CeedScalar nu = context->nu;

   // Quadrature Point Loop
     CeedPragmaSIMD
     for (CeedInt i=0; i<Q; i++) {
       // Read spatial derivatives of u
       const CeedScalar du[3][3]   = {{ug[i+(0+0*3)*Q],
                                       ug[i+(0+1*3)*Q],
                                       ug[i+(0+2*3)*Q]},
                                      {ug[i+(1+0*3)*Q],
                                       ug[i+(1+1*3)*Q],
                                       ug[i+(1+2*3)*Q]},
                                      {ug[i+(2+0*3)*Q],
                                       ug[i+(2+1*3)*Q],
                                       ug[i+(2+2*3)*Q]}
                                     };
       // -- Qdata
       const CeedScalar wJ         =    qdata[0*Q+i];
       // *INDENT-OFF*
       const CeedScalar dXdx[3][3] =  {{qdata[1*Q+i],
                                        qdata[2*Q+i],
                                        qdata[3*Q+i]},
                                       {qdata[4*Q+i],
                                        qdata[5*Q+i],
                                        qdata[6*Q+i]},
                                       {qdata[7*Q+i],
                                        qdata[8*Q+i],
                                        qdata[9*Q+i]}
                                      };


     // Apply dXdx^-1
     CeedScalar gradu[3][3];
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         gradu[j][k] = 0;
         for (int m=0; m<3; m++)
           gradu[j][k] += dXdx[j][m]*gradu[k][m];
       }

     // Compute e = 1/2 (grad u + (grad u)^T)
     const CeedScalar e[3][3]     =  {{(gradu[0][0] + gradu[0][0])/2,
                                       (gradu[0][1] + gradu[1][0])/2,
                                       (gradu[0][2] + gradu[2][0])/2},
                                      {(gradu[1][0] + gradu[0][1])/2,
                                       (gradu[1][1] + gradu[1][1])/2,
                                       (gradu[1][2] + gradu[2][1])/2},
                                      {(gradu[2][0] + gradu[0][2])/2,
                                       (gradu[2][1] + gradu[1][2])/2,
                                       (gradu[2][2] + gradu[2][2])/2}
                                     };

     // Sigma = S e
     const CeedScalar ss          =  E/((1+nu)*(1-2*nu));
     const CeedScalar sigma[3][3] =
      { {ss*((1-nu)*e[1][1] + nu*e[2][2] +nu*e[3][3]),
          ss*(1-2*nu)*e[1][2]/2, ss*(1-2*nu)*e[1][3]/2},
         {ss*(1-2*nu)*e[2][1]/2, ss*(nu*e[1][1] + (1-nu)*e[2][2] +nu*e[3][3]),
          ss*(1-2*nu)*e[2][3]/2},
          {ss*(1-2*nu)*e[3][1]/2, ss*(1-2*nu)*e[3][2]/2,
            ss*(nu*e[1][1] + nu*e[2][2] +(1-nu)*e[3][3])}
};

     // *INDENT-ON*

     // Apply dXdx^-T
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         vg[i+(j*3+k)*Q] = 0;
         for (int m=0; m<3; m++)
           vg[i+(j*3+k)*Q] += dXdx[m][j] * sigma[k][m] * wJ;

       }

   } // End of Quadrature Point Loop

   // Return
   return 0;
}
//
