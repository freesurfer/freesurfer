/**
 * @brief test routines
 *
 */
/*
 * Original Author: Y. Tosa
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iostream>
#include <iomanip>


extern "C"
{
#include "error.h"
#include "mri.h"
#include "matrix.h"
const char *Progname = "surftest";
}

#include "transform.h"

using namespace std;

#define V4_LOAD(v, x, y, z, r)  (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z, VECTOR_ELT(v,4)=r) ;

void MPrint(MATRIX *m)
{
  cout.setf(ios::fixed ,ios::floatfield);
  cout.precision(10);
  cout << endl;
  cout << setw(15) << *MATRIX_RELT(m,1,1)
  << setw(15) << *MATRIX_RELT(m,1,2)
  << setw(15) << *MATRIX_RELT(m,1,3)
  << setw(15) << *MATRIX_RELT(m,1,4) << '\n'
  << setw(15) << *MATRIX_RELT(m,2,1)
  << setw(15) << *MATRIX_RELT(m,2,2)
  << setw(15) << *MATRIX_RELT(m,2,3)
  << setw(15) << *MATRIX_RELT(m,2,4) << '\n'
  << setw(15) << *MATRIX_RELT(m,3,1)
  << setw(15) << *MATRIX_RELT(m,3,2)
  << setw(15) << *MATRIX_RELT(m,3,3)
  << setw(15) << *MATRIX_RELT(m,3,4) << '\n'
  << setw(15) << *MATRIX_RELT(m,4,1)
  << setw(15) << *MATRIX_RELT(m,4,2)
  << setw(15) << *MATRIX_RELT(m,4,3)
  << setw(15) << *MATRIX_RELT(m,4,4) << endl;
}

#if 0
# now defined in utils lib
int MRIRASToSurfaceRAS(MRI *mri, Real xr, Real yr, Real zr, Real *sxr, Real *syr, Real *szr)
{
  double m11, m12, m13, m14;
  double m21, m22, m23, m24;
  double m31, m32, m33, m34;
  // MATRIX *conVoxelFromRAS;
  //  MATRIX *surfaceRASFromConVoxel;
  MATRIX *surfaceRASFromRAS;
  VECTOR *vr, *sr;

  // inv = MatrixInverse(conRASFromVoxel, NULL);
  // explicitly calculated
  // ended up with the following
  /////////////////////////////////////////////////////////////
  surfaceRASFromRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  m11= 1.;
  m12= 0.;
  m13= 0.;
  m14= -mri->c_r;
  m21= 0.;
  m22= 1.;
  m23= 0.;
  m24= -mri->c_a;
  m31= 0.;
  m32= 0.;
  m33= 1.;
  m34= -mri->c_s;
  stuff_four_by_four(surfaceRASFromRAS,
                     m11, m12, m13, m14,
                     m21, m22, m23, m24,
                     m31, m32, m33, m34,
                     0.0, 0.0, 0.0, 1.0);

  vr = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(vr, xr, yr, zr, 1.);

  sr = MatrixMultiply(surfaceRASFromRAS, vr, NULL);
  *sxr = V3_X(sr);
  *syr = V3_Y(sr);
  *szr = V3_Z(sr);

  MatrixFree(&surfaceRASFromRAS);
  VectorFree(&vr);
  VectorFree(&sr);

  return (NO_ERROR);
}
#endif

using namespace std;

bool zero(const double &a, const double &epsilon)
{
  return ( a < epsilon ) ? true : false;
}

int main(int argc, char *argv[])
{
//      orig   ---- >   RAS
//        |              |id
//        V              |
//  conformed  ---- >   RAS
//        |              | MRIRASToSurfaceRAS
//        V              V
//  conformed  ---- >   SurfaceRAS (c_(r,a,s) = 0)
//        |id            |id
//        |              |
//     surface ---- >   SurfaceRAS
//
  if (argc < 2)
  {
    cout << "Usage: surftest (volume) [(accuracy)]" << endl;
    exit(-1);
  }

  cout << "Reading " << argv[1] << "...." << endl;
  MRI *mri = MRIread(argv[1]);
  if (mri == 0)
  {
    cerr << "could not read " << argv[1] << endl;
    exit(-1);
  }
  double ep = 0.00001;
  if (argc == 3)
    ep = atof(argv[2]);
  cout << "using accurarcy error to be " << ep << endl;

  int x, y, z;
  double xr1, yr1, zr1;
  double XR1, YR1, ZR1;
  double XR2, YR2, ZR2;

  for (z=0; z < mri->depth; ++z)
  {
    for (y=0; y < mri->height; ++y)
    {
      for (x=0; x < mri->width; ++x)
      {
        // forward path
        MRIvoxelToWorld(mri, x, y, z, &xr1, &yr1, &zr1);
        MRIRASToSurfaceRAS(mri, xr1, yr1, zr1, &XR1, &YR1, &ZR1);
        // vs.
        MRIvoxelToSurfaceRAS(mri, x, y, z, &XR2, &YR2, &ZR2);
        if (!zero(XR1-XR2, ep) || !zero(YR1-YR2, ep) || !zero(ZR1-ZR2, ep))
        {
          cerr << "error : (XR1, YR1, ZR1)= (" << XR1 <<", " << YR1 << ", " << ZR1 << ")" << endl;
          cerr << "error : (XR2, YR2, ZR2)= (" << XR2 <<", " << YR2 << ", " << ZR2 << ")" << endl;
          cerr << "diff  :  XR1 - XR2 = " << XR1 - XR2
          << "  YR1 - YR2 = " << YR1 - YR2
          << "  ZR1 - ZR2 = " << ZR1 - ZR2 << endl;
        }
        // backward path
        double dx, dy, dz;
        MRIsurfaceRASToVoxel(mri, XR2, YR2, ZR2, &dx, &dy, &dz);
        if (!zero(x - dx, ep) || !zero(y - dy, ep) || !zero(z - dz, ep))
        {
          cout << "going back to get back the same vaue" << endl;
          cerr << "error : (  x,   y,   z)= (" << x <<", " << y << ", " << z << ")" << endl;
          cerr << "error : ( dx,  dy,  dz)= (" << dx <<", " << dy << ", " << dz << ")" << endl;
          cerr << "diff  : (" << (x-dx) << ", " << (y-dy) <<", " << (z-dz) << ")" << endl;
        }
      }
    }
    cout << ".";
    cout.flush();
  }
  cout << endl;
}
