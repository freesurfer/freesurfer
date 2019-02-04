/*

  Performing the continuous max-flow algorithm to solve the
  3D continuous-cut problem with multiple labels (Potts Model)

  Usage: u = CMF3D_ML(bound, iNbIters, fError, cc, steps);

  Inputs:

    - bound: point to the capacities of sink flows pt(x, i=1...nlab);
    - iNbIters: the maximum iteration number
    - fError: the error bound for convergence
    - cc: cc for the step-size of augmented Lagrangian method
    - steps: the step-size for the graident-projection step to the
             total-variation function. Its optimal range is [0.06, 0.12].

  Outputs:

    - u: the computed continuous labeling function u(x,i=1...nlab) in [0,1].
         As the following paper [2], the final label function is
         given by the maximum of u(x,i=1...nlab) at each x.

  The original continuous max-flow algorithm was proposed in the following papers:

  [1] Yuan, J.; Bae, E.;  Tai, X.-C.
     A Study on Continuous Max-Flow and Min-Cut Approaches
     CVPR, 2010

  [2] Yuan, J.; Bae, E.; Tai, X.-C.; Boycov, Y.
     A Continuous Max-Flow Approach to Potts Model
     ECCV, 2010

  The mimetic finite-difference discretization method was proposed for
  the total-variation function in the paper:

  [1] Yuan, J.; Schn{\"o}rr, C.; Steidl, G.
     Simultaneous Optical Flow Estimation and Decomposition
     SIAM J.~Scientific Computing, 2007, vol. 29, page 2283-2304, number 6

  This software can be used only for research purposes, you should cite ALL of
  the aforementioned papers in any resulting publication.

  Please email Jing Yuan (cn.yuanjing@gmail.com) for any questions, suggestions and bug reports

  The Software is provided "as is", without warranty of any kind.

  Version 1.0
  https://sites.google.com/site/wwwjingyuan/

  Copyright 2011 Jing Yuan (cn.yuanjing@gmail.com)

  Ported to python by Andrew

*/

#include <exception>
#include <stdlib.h>
#include <math.h>

#include "labelfusion.h"


#define YES 0
#define NO 1
#define PI 3.1415926
#define MAX(a,b) std::max(a, b)
#define MIN(a,b) std::min(a, b)
#define ABS(x)   std::abs(x)
#define SQRT(x)  std::sqrt(x)
#define SQR(x)   (x)*(x)
#define SIGN(x)  ( x >= 0.0 ? 1.0 : -1.0 )


pyarrayf CMF3D_ML(const pyarrayf &bound, int iNbIters, float fError, float cc, float steps)
{
    float *pfu, *pfbx, *pfby, *pfbz, *pfps, *pfpt, *pfgk, *pfft, *pfdiv;
    float fpt, fps, pfcvg;
    unsigned int iNy, iNx, iNz, iLab, ix, iy, iz, id;
    int ik, iNI;
    int szImg, S3D, S2D;
    int index, index1, indz, indd, indy;

    if (bound.ndim() != 4) throw std::runtime_error("input must have 4 dimensions");
    
    iNy  = bound.shape(0);
    iNx  = bound.shape(1);
    iNz  = bound.shape(2);
    iLab = bound.shape(3);

    S2D = iNy * iNx;
    S3D = S2D * iNz;
    szImg = S3D * iLab;

    const float *pfCt = bound.data(0);

    float * const u = new float[bound.size()];
    pfu = u;

    // allocate the memory for px1
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz*iLab), sizeof(float) );
    if (!pfbx) throw std::runtime_error("memory allocation failure");

    // allocate the memory for py1
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz*iLab), sizeof(float) );
    if (!pfby) throw std::runtime_error("memory allocation failure");

    // allocate the memory for pz1
    pfbz = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)*iLab), sizeof(float) );
    if (!pfbz) throw std::runtime_error("memory allocation failure");

    // allocate the memory for ps
    pfdiv = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfdiv) throw std::runtime_error("memory allocation failure");

    // allocate the memory for ps
    pfps = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps) throw std::runtime_error("memory allocation failure");

    // allocate the memory for pt
    pfpt = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfpt) throw std::runtime_error("memory allocation failure");

    // allocate the memory for gk
    pfgk = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfgk) throw std::runtime_error("memory allocation failure");

    // allocate the memory for ft
    pfft = (float *) calloc( (unsigned)(iLab), sizeof(float) );
    if (!pfft) throw std::runtime_error("memory allocation failure");

    // Preprocessing initial values
    for (iz=0; iz < iNz; iz++){
        indz = iz * S2D;
        for (ix=0; ix < iNx; ix++){
            indy = ix*iNy + indz;
            for (iy=0; iy < iNy; iy++){
                index = indy + iy;

                fpt = pfCt[index];
                ik = 0;

                for (id = 1; id < iLab; id++){
                    index1 = index + id*S3D;
                    if (fpt >= pfCt[index1]){
                        fpt = pfCt[index1];
                        ik = id;
                    }
                }

                pfps[index] = fpt;
                pfu[index+ik*S3D] = 1;

                for (id = 0; id < iLab; id++){
                    pfpt[index+id*S3D] = fpt;
                }
            }
        }
    }

    //  Main iterations

    iNI = 0;

    while (iNI < iNbIters)
    {

        // update the flow fields p(x,i=1...nlab) and pt(x,i=1...nlab) at each layer

        for (id = 0; id < iLab; id++){

            indd = id*S3D;

            // update the spatial flow field p(x,id) = (bx(x,id),by(x,id),bz(x,id))
            // bx, by, bz define the three components of the vector field p(x,id)

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = indy + iy;
                        index1 = indd + index;
                        pfgk[index] = pfdiv[index1] - (pfps[index] - pfpt[index1] + pfu[index1]/cc);
                    }
                }
            }

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx-1; ix++){
                    indy = (ix+1)*iNy+indz;
                    for (iy=0; iy< iNy; iy++){
                        index = indy + iy;
                        index1 = index + indd;
                        pfbx[index1] = steps*(pfgk[index] - pfgk[index-iNy]) + pfbx[index1];
                    }
                }
            }

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz;
                    for(iy = 0; iy < iNy-1; iy ++){
                        index = (iy+1) + indy;
                        index1 = index + indd;
                        pfby[index1] = steps*(pfgk[index] - pfgk[index-1]) + pfby[index1];
                    }
                }
           }

            for (iz=0; iz< iNz-1; iz++){
                indz = (iz+1)*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz;
                    for(iy = 0; iy < iNy; iy ++){
                        index = indy + iy;
                        index1 = index + indd;
                        pfbz[index1] = steps*(pfgk[index] - pfgk[index-S2D]) + pfbz[index1];
                    }
                }
           }

            // projection step

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;

                        fpt = SQRT((SQR(pfbx[index1+iNy]) + SQR(pfbx[index1]) + SQR(pfby[index1+1]) +
                              SQR(pfby[index1]) + SQR(pfbz[index1+S2D]) + SQR(pfbz[index1]))*0.5);

                        if (fpt > 1.0) {
                            pfgk[index] = 1.0 / fpt;
                        } else {
                            pfgk[index] = 1;
                        }
                    }
                }
            }

            // update the component bx(x,id)

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx-1; ix++){
                    indy = (ix+1)*iNy + indz;
                    for (iy=0; iy< iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;
                        pfbx[index1] = (pfgk[index] + pfgk[index-iNy])*0.5*pfbx[index1];
                    }
                }
            }

            // update the component by(x,id)

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz;
                    for (iy=0; iy<iNy-1; iy++){
                       index = (iy+1) + indy;
                       index1 = index + indd;
                       pfby[index1] = 0.5*(pfgk[index-1] + pfgk[index])*pfby[index1];
                    }
                }
            }

            // update the component bz(x,id)

            for (iz=0; iz< iNz-1; iz++){
                indz = (iz+1)*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz;
                    for (iy=0; iy<iNy; iy++){
                       index = iy + indy;
                       index1 = index + indd;
                       pfbz[index1] = 0.5*(pfgk[index-S2D] + pfgk[index])*pfbz[index1];
                    }
                }
            }

            // update the sink flow field pt(x,id)

            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;
                        // update the divergence field div(x,id)
                        pfdiv[index1] = pfbx[index1+iNy] - pfbx[index1] + pfby[index1+1] -
                                pfby[index1] + pfbz[index1+S2D] - pfbz[index1];
                        fpt = pfps[index] + pfu[index1]/cc - pfdiv[index1];
                        pfpt[index1] = MIN(fpt, pfCt[index1]);
                    }
                }
            }
        }

        // update the source flow field ps(x)

        fps = 0;

        for (iz=0; iz< iNz; iz++){
            indz = iz*S2D;
            for (ix=0; ix< iNx; ix++){
                indy = ix*iNy + indz;
                for (iy = 0; iy < iNy; iy++){
                    index = iy + indy;

                    fpt = 0;
                    for (id = 0; id < iLab; id++){
                        index1 = index + id*S3D;
                        pfft[id] = pfdiv[index1] + pfpt[index1];
                        fpt += (pfft[id] - pfu[index1]/cc);
                    }

                    pfps[index] = fpt/iLab + 1/(cc*iLab);

                    // update the multipliers

                    for (id = 0; id < iLab; id++){
                        fpt = cc*(pfft[id] - pfps[index]);
                        pfu[index+id*S3D] -= fpt;
                        fps += ABS(fpt);
                    }
                }
            }
        }

        pfcvg = fps / szImg;
        if (pfcvg <= fError) break;
        iNI ++;
    }

    // Free memory
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfbz );
    free( (float *) pfdiv );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfgk );
    free( (float *) pfft );

    return makeArray({iNy, iNx, iNz, iLab}, MemoryOrder::Fortran, u);
}
