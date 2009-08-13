//
// MyMatrix is a class for several matrix operations
//    as used for registration (rigid, affine maps)
//    conversion, halfway spaces,...
//
// written by Martin Reuter
// Aug. 12th ,2009
//

#ifndef MyMatrix_H
#define MyMatrix_H

#ifdef __cplusplus
extern "C"
{
#endif
#include "matrix.h"
#include "mri.h"
#ifdef __cplusplus
}
#endif

#include <utility>
#include <string>
#include <vector>

class MyMatrix
{
public:
  // distances
  static double RigidTransDistSq(MATRIX *a, MATRIX *b = NULL);
  static double AffineTransDistSq(MATRIX *a, MATRIX *b = NULL, double r=100);
  static double getFrobeniusDiff(MATRIX *m1, MATRIX *m2);

  // opeartions
  static MATRIX * MatrixSqrt(MATRIX * m, MATRIX * sqrtm=NULL);

  // conversions
  static MATRIX* getMatrix(std::vector < double > d, int r, int c=-1, MATRIX* m=NULL);
  static MATRIX * aff2mat(MATRIX * aff, MATRIX *outM);
  static MATRIX * getHalfRT (MATRIX * m, MATRIX * mhalf=NULL);
  static double RotMatrixLogNorm(MATRIX * m);
  static double RotMatrixGeoDist(MATRIX * a, MATRIX *b = NULL);

};


#endif
