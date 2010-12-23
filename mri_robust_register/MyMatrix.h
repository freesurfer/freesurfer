/**
 * @file MyMatrix.cpp
 * @brief A static class with Matrix operations
 *
 *    as used for registration (rigid, affine maps)
 *    conversion, halfway spaces,...
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/12/23 00:08:40 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2008-2009
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

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
#include "transform.h"
#ifdef __cplusplus
}
#endif

#include <utility>
#include <string>
#include <vector>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

class MyMatrix
{
public:

  // conversion
  static MATRIX* convertVNL2MATRIX(const vnl_vector <double > & v, MATRIX* outM = NULL);
  static MATRIX* convertVNL2MATRIX(const vnl_matrix <double > & v, MATRIX* outM = NULL);
  static vnl_vector <double > convertVECTOR2VNL(VECTOR* m);
  static vnl_matrix <double > convertMATRIX2VNL(MATRIX* m);

//======== VNL STUFF ===========================================================================

  // operations
  static vnl_matrix < double >  MatrixSqrtIter(const vnl_matrix < double >& m);
  static vnl_matrix < double >  MatrixSqrt(const vnl_matrix < double >& m);
  static std::pair < vnl_matrix < double > , vnl_matrix < double > >
	   MatrixSqrtAndInv(const vnl_matrix < double >& m);

  // distances
  static double RigidTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b  = vnl_matrix<double>());
  static double AffineTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b = vnl_matrix<double>(), double r=100);
  static double getFrobeniusDiff(const vnl_matrix < double >&m1, const vnl_matrix < double >&m2);

  // conversions
  static vnl_matrix < double > getVNLMatrix(std::vector < double > d, int r);
  static double RotMatrixLogNorm(const vnl_matrix_fixed < double, 4, 4 > &m);
  static LTA* VOXmatrix2LTA(const vnl_matrix_fixed < double, 4, 4 >&m, MRI* src, MRI* dst);
  static LTA* RASmatrix2LTA(const vnl_matrix_fixed < double, 4, 4 >&m, MRI* src, MRI* dst);
  static void getRTfromM(const vnl_matrix_fixed < double , 4 , 4 > &m,
	                             vnl_matrix_fixed < double , 3 , 3 > &r, 
															 vnl_vector_fixed < double, 3 >      &t);
  static vnl_matrix_fixed < double , 4 , 4 >  getMfromRT(
	                             const vnl_matrix_fixed < double , 3 , 3 > &r, 
															 const vnl_vector_fixed < double, 3 >      &t);

//========= MATRIX STUFF ========================================================================
	
  // distances
  static double RigidTransDistSq(MATRIX *a, MATRIX *b = NULL);
  static double AffineTransDistSq(MATRIX *a, MATRIX *b = NULL, double r=100);
  static double getFrobeniusDiff(MATRIX *m1, MATRIX *m2);

  // operations
  static MATRIX * MatrixSqrt(MATRIX * m, MATRIX * sqrtm=NULL);

  // conversions
  static MATRIX* getMatrix(std::vector < double > d, int r, int c=-1, MATRIX* m=NULL);
  static MATRIX * aff2mat(MATRIX * aff, MATRIX *outM);
  static MATRIX * getHalfRT (MATRIX * m, MATRIX * mhalf=NULL);
  static double RotMatrixLogNorm(MATRIX * m);
  static double RotMatrixGeoDist(MATRIX * a, MATRIX *b = NULL);
  static std::pair < MATRIX *, VECTOR * > getRTfromM(MATRIX * M, MATRIX * R, VECTOR * T);
  static MATRIX * getMfromRT(MATRIX * R, VECTOR * T, MATRIX * M);
  static LTA* VOXmatrix2LTA(MATRIX * m, MRI* src, MRI* dst);
  static LTA* RASmatrix2LTA(MATRIX * m, MRI* src, MRI* dst);

};


#endif
