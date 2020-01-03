#ifndef __LIN_ELAS__H
#define __LIN_ELAS__H

    typedef struct Physics_private *Physics;
    struct Physics_private{
      PetscScalar   nu;      //Poisson's ratio
      PetscScalar   E;       //Young's Modulus
    };


// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinElas)(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {
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

     //Compute gradu
     // Apply dXdx^-1 to du = gradu
     CeedScalar gradu[3][3];
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         gradu[j][k] = 0;
         for (int m=0; m<3; m++)
           gradu[j][k] += dXdx[j][m]*du[k][m];
       }

     // Compute Strain : e (epsilon)
     // e = 1/2 (grad u + (grad u)^T)
     const CeedScalar e[3][3]     =  {{(gradu[0][0] + gradu[0][0])*0.5,
                                       (gradu[0][1] + gradu[1][0])*0.5,
                                       (gradu[0][2] + gradu[2][0])*0.5},
                                      {(gradu[1][0] + gradu[0][1])*0.5,
                                       (gradu[1][1] + gradu[1][1])*0.5,
                                       (gradu[1][2] + gradu[2][1])*0.5},
                                      {(gradu[2][0] + gradu[0][2])*0.5,
                                       (gradu[2][1] + gradu[1][2])*0.5,
                                       (gradu[2][2] + gradu[2][2])*0.5}
                                     };
    //strain (epsilon)
    //    and
    //stress (sigma) in Voigt notation:
    //           [e11]              [sigma11]
    //           [e22]              [sigma22]
    // epsilon = [e33]  ,   sigma = [sigma33]
    //           [e23]              [sigma23]
    //           [e13]              [sigma13]
    //           [e12]              [sigma12]
    //
    // Sigma = S * epsilon
    //                         [1-nu   nu    nu                                    ]
    //                         [ nu   1-nu   nu                                    ]
    //                         [ nu    nu   1-nu                                   ]
    // S = E/((1+nu)*(1-2*nu)) [                  (1-2*nu)/2                       ]
    //                         [                             (1-2*nu)/2            ]
    //                         [                                        (1-2*nu)/2 ]

    //Above Voigt Notation is placed in a 3x3 matrix:
     const CeedScalar ss          =  E/((1+nu)*(1-2*nu));
     const CeedScalar sigma[3][3] =
      { {ss*((1-nu)*e[1][1] + nu*e[2][2] +nu*e[3][3]), ss*(1-2*nu)*e[1][2]/2, ss*(1-2*nu)*e[1][3]/2},
        {ss*(1-2*nu)*e[2][1]/2, ss*(nu*e[1][1] + (1-nu)*e[2][2] +nu*e[3][3]), ss*(1-2*nu)*e[2][3]/2},
        {ss*(1-2*nu)*e[3][1]/2, ss*(1-2*nu)*e[3][2]/2, ss*(nu*e[1][1] + nu*e[2][2] +(1-nu)*e[3][3])}
      };

     // *INDENT-ON*

     // Apply dXdx^-T and weight
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         vg[i+(j*3+k)*Q] = 0;
         for (int m=0; m<3; m++)
           vg[i+(j*3+k)*Q] += dXdx[m][j] * sigma[k][m] * wJ;
       }
     } // End of Quadrature Point Loop

   return 0;
}

#endif //End of __LIN_ELAS__H
