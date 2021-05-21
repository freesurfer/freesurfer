/**
 * @brief A static class with Matrix operations
 *
 *    as used for registration (rigid, affine maps)
 *    conversion, halfway spaces,...
 */

/*
 * Original Author: Martin Reuter
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

// written by Martin Reuter
// Aug. 12th ,2009
//
#ifndef MyMatrix_H
#define MyMatrix_H

#include "matrix.h"
#include "mri.h"
#include "transform.h"

#include <utility>
#include <string>
#include <vector>
#include <iostream>

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#undef export

#include <vcl_complex.h>

/** \class MyMatrix
 * \brief A static class with Matrix operations (conversion, sqrt, decomposition, distances etc.)
 */
class MyMatrix
{
public:

  //! Convert VNL vector to FS MATRIX
  static MATRIX* convertVNL2MATRIX(const vnl_vector<double> & v, MATRIX* outM =
      NULL);
      
  //! Convert VNL matrix to FS MATRIX
  static MATRIX* convertVNL2MATRIX(const vnl_matrix<double> & v, MATRIX* outM =
      NULL);
      
  //! Convert FS VECTOR to VNL vector
  static vnl_vector<double> convertVECTOR2VNL(VECTOR* m);
  
  //! Convert FS MATRIX to VNL matrix
  static vnl_matrix<double> convertMATRIX2VNL(MATRIX* m);


//======== VNL MATRIX STUFF ===========================================================================

  //! Matrix Square Root (Denman and Beavers)
  static vnl_matrix<double> MatrixSqrtIter(const vnl_matrix<double>& m);
  
  //! Matrix Square Root and Inverse (Denman and Beavers)
  static std::pair<vnl_matrix<double>, vnl_matrix<double> >
  MatrixSqrtAndInvIter(const vnl_matrix<double>& m);
  
  //! Matrix Square Root for Diagonalizable M (using SVD)
  static vnl_matrix<double> MatrixSqrtEigs(const vnl_matrix<double>& m);
  
  //! Matrix Square Root (using Complex Schur)
  static vnl_matrix<double> MatrixSqrt(const vnl_matrix<double>& m);
  
  //! Matrix Square Root (splitting off the translation)
  static vnl_matrix<double> MatrixSqrtAffine(const vnl_matrix<double>& m);
  

  //! Returns rotation matrix (multiples of 90 degrees)
  static vnl_matrix_fixed<double, 4, 4> getRot(int i, int j, int k);
  
  //! Returns rotation matrix (multiples of 90 degrees)
  static vnl_matrix_fixed<double, 4, 4> getRot(int i);
  
  //! Mean Rotation
  static vnl_matrix_fixed < double, 3, 3> RotationMean(const std::vector < vnl_matrix_fixed < double, 3, 3 > > &vm);

  //! Geometric mean
  static vnl_matrix < double > GeometricMean(const std::vector < vnl_matrix < double > > &vm, int n=-1);
  static vnl_matrix < double > GeometricMean(const std::vector < vnl_matrix_fixed < double, 3, 3 > > &vm, int n=-1);
  
  //! Polar Decomposition: A = R * S  (R orthogonal, S pos. semi def, symmetric)
  static void PolarDecomposition(const vnl_matrix<double> &A,
      vnl_matrix<double> &R, vnl_matrix<double> &S);
      
  //! Advanced Polar Decomposition A = Rot * Shear * Scale  (where scale is diag)
  static void Polar2Decomposition(const vnl_matrix<double> &A,
      vnl_matrix<double> &R, vnl_matrix<double> &S, vnl_diag_matrix<double> &D);

  //! Complex Schur Decomposition
  static void SchurComplex(const vnl_matrix<double> & A,
      vnl_matrix<vcl_complex< double > > & U,
      vnl_matrix < vcl_complex < double > > & T);
      
  //! Test if diagonal (eps)
  static bool isDiag( const vnl_matrix < double > & A, double eps = 0.00000000000001 );
  
  //! Matrix Exponent (see matlab)
  static vnl_matrix < double > MatrixExp(const vnl_matrix < double >& A);
  
  //! Matrix Logarithm (see matlab)
  static vnl_matrix < double > MatrixLog(const vnl_matrix < double >& A, int maxlogiter = 100);
  
  //! Matrix Power
  static vnl_matrix < double > MatrixPow(const vnl_matrix < double >& A, double d);

  //! Squared distance based on rigid transformation
  static double RigidTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b = vnl_matrix<double>());

  //! Squared distance based on affine transformation
  static double AffineTransDistSq(const vnl_matrix < double >&a, const vnl_matrix < double >&b = vnl_matrix<double>(), double r=100);

  //! Distance based on frobenius norm
  static double getFrobeniusDiff(const vnl_matrix < double >&m1, const vnl_matrix < double >&m2);

  //! Quantify tri-linear smoothing induce by this LTA transform
  static double getResampSmoothing(const LTA*);

  //! Construct vnl matrix from double array using r rows
  static vnl_matrix < double > getVNLMatrix(std::vector < double > d, int r);
  
  //! Computes Frobenius norm of log of rotation matrix
  static double RotMatrixLogNorm(const vnl_matrix_fixed < double, 4, 4 > &m);
  
  //! Construct LTA from vox2vox matrix and src, dst geometries
  static LTA* VOXmatrix2LTA(const vnl_matrix_fixed < double, 4, 4 >&m, MRI* src, MRI* dst);
  
  //! Construct LTA from ras2ras matrix and src, dst geometries
  static LTA* RASmatrix2LTA(const vnl_matrix_fixed < double, 4, 4 >&m, MRI* src, MRI* dst);
  
  //! Get vox2vox matrix from LTA
  static vnl_matrix < double > LTA2VOXmatrix (LTA * lta);
  
  //! Get ras2ras matrix from LTA
  static vnl_matrix < double > LTA2RASmatrix (LTA * lta);
  
  //! Compute rotation and translation from rigid 4x4 matrix 
  static void getRTfromM(const vnl_matrix_fixed < double , 4 , 4 > &m,
              vnl_matrix_fixed < double , 3 , 3 > &r,
              vnl_vector_fixed < double, 3 > &t);
              
  //! Construct 4x4 rigid matrix from rotation and translation
  static vnl_matrix_fixed < double , 4 , 4 > getMfromRT(
              const vnl_matrix_fixed < double , 3 , 3 > &r,
              const vnl_vector_fixed < double, 3 > &t);


//========= FS MATRIX STUFF ========================================================================

  //! (FS MATRIX) Squared distance based on rigid transformation
  static double RigidTransDistSq(MATRIX *a, MATRIX *b = NULL);

  //! (FS MATRIX) Squared distance based on affine transformation
  static double AffineTransDistSq(MATRIX *a, MATRIX *b = NULL, double r=100);

  //! (FS MATRIX) Distance based on frobenius norm
  static double getFrobeniusDiff(MATRIX *m1, MATRIX *m2);

  //! (FS MATRIX) Matrix Square Root (Denman and Beavers)
  static MATRIX * MatrixSqrtIter(MATRIX * m, MATRIX * sqrtm=NULL);

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

private:

  // helpers for exp
  static std::vector < double > getPadeCoefficients(unsigned int m);
  static vnl_matrix < double > PadeApproximantOfDegree(const vnl_matrix < double > & A , unsigned int m );
  static void expmchk(const vnl_matrix < double >& A, std::vector < unsigned int > &m_vals, std::vector < double > &theta);

  // helpers for log
  static void gauss_legendre(int n, vnl_vector < double >& x, vnl_vector < double > & w );
  static vnl_matrix < vcl_complex <double > > MatrixLog_pf(const vnl_matrix < vcl_complex <double > >& A , unsigned int m);
  static vnl_matrix < vcl_complex <double > > MatrixLog_isst(const vnl_matrix < vcl_complex < double > >& A, int maxlogiter=100);
  static std::vector < int > blocking(const vnl_matrix < vcl_complex < double > > &A, double delta = 0.1);
  static void swapping(const std::vector < int > &m, std::vector < int > &mm, std::vector < std::vector < int > > &ind);
  static std::vector < int > cumsum0(const std::vector < int > &v, const std::vector < std::pair < double , int > > & w);
  static vnl_matrix < vcl_complex < double > > sylv_tri(
              const vnl_matrix < vcl_complex < double > > & T,
              const vnl_matrix < vcl_complex < double > > & U,
              const vnl_matrix < vcl_complex < double > > & B);
  static vnl_matrix < vcl_complex < double > > getSubMatrix(const vnl_matrix < vcl_complex < double > > & A,
              const std::vector < int > & rows,
              const std::vector < int > & cols);
  static void setSubMatrix(vnl_matrix < vcl_complex < double > > & A,
              const std::vector < int > & rows,
              const std::vector < int > & cols,
              const vnl_matrix < vcl_complex < double > > & B);
  static void OrdSchurComplexLogical (const vnl_matrix < vcl_complex < double > > & U,
              const vnl_matrix < vcl_complex < double > > & T,
              const std::vector < int >& select,
              vnl_matrix < vcl_complex < double > > & US,
              vnl_matrix < vcl_complex < double > > & TS);
  static void OrdSchurComplex (const vnl_matrix < vcl_complex < double > > & U,
              const vnl_matrix < vcl_complex < double > > & T,
              const std::vector < int >& clusters,
              vnl_matrix < vcl_complex < double > > & US,
              vnl_matrix < vcl_complex < double > > & TS);

};

// // helper code to print vectors for debugging
// template <class T> void Print(std::vector<T> & Vec)
// {
//    typename std::vector<T>::iterator p;
// 
//    std::cout << "[ ";
//    for (p = Vec.begin(); p != Vec.end(); p++)
//       std::cout << *p << " ";
//    std::cout << "]" << std::endl << std::endl;
// }
// 
// template <class T, class U> void Print(std::vector<std::pair < T, U > > & Vec)
// {
//    typename std::vector< std::pair < T,U> >::iterator p;
// 
//    std::cout << "[ ";
//    for (p = Vec.begin(); p != Vec.end(); p++)
//       std::cout << "(" << (*p).first << ":" << (*p).second << ") ";
//    std::cout << "]" << std::endl << std::endl;
// }

#endif
