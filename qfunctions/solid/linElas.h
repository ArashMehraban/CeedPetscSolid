#ifndef __LIN_ELAS__H
#define __LIN_ELAS__H

#ifndef PHYSICS_STRUCT
#define PHYSICS_STRUCT
    typedef struct Physics_private *Physics;
    struct Physics_private{
      PetscScalar   nu;      //Poisson's ratio
      PetscScalar   E;       //Young's Modulus
    };
#endif

// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinElasF)(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {
   // *INDENT-OFF*
   // Inputs
   const CeedScalar (*ug)[3][Q] = (CeedScalar(*)[3][Q])in[0],
                    (*qdata)[Q] = (CeedScalar(*)[Q])in[1];

   // Outputs
   CeedScalar (*vg)[3][Q] = (CeedScalar(*)[3][Q])out[0];
              // gradu not used for linear elasticity
              // (*gradu)[3][Q] = (CeedScalar(*)[3][Q])out[1];
   // *INDENT-ON*

   // Context
   const Physics context = ctx;
   const CeedScalar E  = context->E;
   const CeedScalar nu = context->nu;


//PetscPrintf(PETSC_COMM_WORLD, "LinElasF in F\n");

   // Quadrature Point Loop
     CeedPragmaSIMD
     for (CeedInt i=0; i<Q; i++) {
       // Read spatial derivatives of u
       // *INDENT-OFF*
       const CeedScalar du[3][3]   = {{ug[0][0][i],
                                       ug[0][1][i],
                                       ug[0][2][i]},
                                      {ug[1][0][i],
                                       ug[1][1][i],
                                       ug[1][2][i]},
                                      {ug[2][0][i],
                                       ug[2][1][i],
                                       ug[2][2][i]}
                                     };
       // -- Qdata
       const CeedScalar wJ         =    qdata[0][i];
       const CeedScalar dXdx[3][3] =  {{qdata[1][i],
                                        qdata[2][i],
                                        qdata[3][i]},
                                       {qdata[4][i],
                                        qdata[5][i],
                                        qdata[6][i]},
                                       {qdata[7][i],
                                        qdata[8][i],
                                        qdata[9][i]}
                                      };
     // *INDENT-ON*

     //Compute gradu
     // Apply dXdx^-T to du = gradu
     CeedScalar gradu[3][3];
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         gradu[j][k] = 0;
         for (int m=0; m<3; m++)
           gradu[j][k] += dXdx[j][m]*du[m][k];
       }

     // Compute Strain : e (epsilon)
     // e = 1/2 (grad u + (grad u)^T)
     // *INDENT-OFF*
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

    // *INDENT-ON*
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
     const CeedScalar ss      =  E/((1+nu)*(1-2*nu));
     const CeedScalar sigma00 = ss*((1-nu)*e[0][0] + nu*e[1][1] +nu*e[2][2]),
                      sigma11 = ss*(nu*e[0][0] + (1-nu)*e[1][1] +nu*e[2][2]),
                      sigma22 = ss*(nu*e[0][0] + nu*e[1][1] +(1-nu)*e[2][2]),
                      sigma12 = ss*(1-2*nu)*e[1][2]/2,
                      sigma02 = ss*(1-2*nu)*e[0][2]/2,
                      sigma01 = ss*(1-2*nu)*e[0][1]/2;
     const CeedScalar sigma[3][3] =
      { {sigma00, sigma01, sigma02},
        {sigma01, sigma11, sigma12},
        {sigma02, sigma12, sigma22}
      };

     // Apply dXdx^-1 and weight
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         vg[j][k][i] = 0;
         for (int m=0; m<3; m++)
           vg[j][k][i] += dXdx[m][j] * sigma[m][k] * wJ;
       }
     } // End of Quadrature Point Loop

   return 0;
}

// -----------------------------------------------------------------------------
CEED_QFUNCTION(LinElasdF)(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {
   // *INDENT-OFF*
   // Inputs
   const CeedScalar (*deltaug)[3][Q] = (CeedScalar(*)[3][Q])in[0],
                    (*qdata)[Q] = (CeedScalar(*)[Q])in[1];
                    // gradu not used for linear elasticity
                    // (*gradu)[3][Q] = (CeedScalar(*)[3][Q])in[2];


//PetscPrintf(PETSC_COMM_WORLD, "LinElasdF in dF\n");

   // Outputs
   CeedScalar (*deltavg)[3][Q] = (CeedScalar(*)[3][Q])out[0];
   // *INDENT-ON*

   // Context
   const Physics context = ctx;
   const CeedScalar E  = context->E;
   const CeedScalar nu = context->nu;

   // Quadrature Point Loop
     CeedPragmaSIMD
     for (CeedInt i=0; i<Q; i++) {
       // Read spatial derivatives of u
       // *INDENT-OFF*
       const CeedScalar deltadu[3][3] = {{deltaug[0][0][i],
                                          deltaug[0][1][i],
                                          deltaug[0][2][i]},
                                         {deltaug[1][0][i],
                                          deltaug[1][1][i],
                                          deltaug[1][2][i]},
                                         {deltaug[2][0][i],
                                          deltaug[2][1][i],
                                          deltaug[2][2][i]}
                                        };
       // -- Qdata
       const CeedScalar wJ         =    qdata[0][i];
       const CeedScalar dXdx[3][3] =  {{qdata[1][i],
                                        qdata[2][i],
                                        qdata[3][i]},
                                       {qdata[4][i],
                                        qdata[5][i],
                                        qdata[6][i]},
                                       {qdata[7][i],
                                        qdata[8][i],
                                        qdata[9][i]}
                                      };
     // *INDENT-ON*

     //Compute graddeltau
     // Apply dXdx^-1 to deltadu = graddeltau
     CeedScalar graddeltau[3][3];
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         graddeltau[j][k] = 0;
         for (int m=0; m<3; m++)
           graddeltau[j][k] += dXdx[j][m]*deltadu[m][k];
       }

     // Compute Strain : e (epsilon)
     // e = 1/2 (grad u + (grad u)^T)
     // *INDENT-OFF*
     const CeedScalar e[3][3]     =  {{(graddeltau[0][0] + graddeltau[0][0])*0.5,
                                       (graddeltau[0][1] + graddeltau[1][0])*0.5,
                                       (graddeltau[0][2] + graddeltau[2][0])*0.5},
                                      {(graddeltau[1][0] + graddeltau[0][1])*0.5,
                                       (graddeltau[1][1] + graddeltau[1][1])*0.5,
                                       (graddeltau[1][2] + graddeltau[2][1])*0.5},
                                      {(graddeltau[2][0] + graddeltau[0][2])*0.5,
                                       (graddeltau[2][1] + graddeltau[1][2])*0.5,
                                       (graddeltau[2][2] + graddeltau[2][2])*0.5}
                                     };

    // *INDENT-ON*
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
    const CeedScalar ss      =  E/((1+nu)*(1-2*nu));
    const CeedScalar sigma00 = ss*((1-nu)*e[0][0] + nu*e[1][1] +nu*e[2][2]),
                     sigma11 = ss*(nu*e[0][0] + (1-nu)*e[1][1] +nu*e[2][2]),
                     sigma22 = ss*(nu*e[0][0] + nu*e[1][1] +(1-nu)*e[2][2]),
                     sigma12 = ss*(1-2*nu)*e[1][2]/2,
                     sigma02 = ss*(1-2*nu)*e[0][2]/2,
                     sigma01 = ss*(1-2*nu)*e[0][1]/2;
    const CeedScalar sigma[3][3] =
     { {sigma00, sigma01, sigma02},
       {sigma01, sigma11, sigma12},
       {sigma02, sigma12, sigma22}
     };

     // Apply dXdx^-T and weight
     for (int j=0; j<3; j++)
       for (int k=0; k<3; k++) {
         deltavg[j][k][i] = 0;
         for (int m=0; m<3; m++)
           deltavg[j][k][i] += dXdx[m][j] * sigma[m][k] * wJ;
       }
     } // End of Quadrature Point Loop

   return 0;
}

#endif //End of __LIN_ELAS__H
