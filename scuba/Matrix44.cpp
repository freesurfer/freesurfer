/**
 * @file  Matrix44.cpp
 * @brief A 4x4 real matrix
 *
 * A wrapper around the libutils MATRIX class for a 4x4 real number
 * matrix. Uses Matrix* functions internally when possible but also
 * includes some geometry functions.Note that all 16 element matrices
 * used in this class are in openGL style format:
 *
 *  [ 0   4   8  12 ]
 *  [ 1   5   9  13 ]
 *  [ 2   6  10  14 ]
 *  [ 3   7  11  15 ]
 *
 * Column/row arguments are in the following format (i,j):
 *
 *  [ 0,0  1,0  2,0  3,0 ]
 *  [ 0,1  1,1  2,1  3,1 ]
 *  [ 0,2  1,2  2,2  3,2 ]
 *  [ 0,3  1,3  2,3  3,3 ]
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.21 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "string_fixed.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <assert.h>

extern "C" {
#include "macros.h"
}

#include "Matrix44.h"
#include "VectorOps.h"

using namespace std;

Matrix44::Matrix44() {

  // Allocate our MATRIX object.
  m = MatrixIdentity( 4, NULL );
}

Matrix44::Matrix44 ( float i0j0, float i1j0, float i2j0, float i3j0,
                     float i0j1, float i1j1, float i2j1, float i3j1,
                     float i0j2, float i1j2, float i2j2, float i3j2,
                     float i0j3, float i1j3, float i2j3, float i3j3 ) {

  // Allocate our MATRIX object.
  m = MatrixIdentity( 4, NULL );

  // Set all our elements.
  SetCR(0,0,i0j0);
  SetCR(1,0,i1j0);
  SetCR(2,0,i2j0);
  SetCR(3,0,i3j0);
  SetCR(0,1,i0j1);
  SetCR(1,1,i1j1);
  SetCR(2,1,i2j1);
  SetCR(3,1,i3j1);
  SetCR(0,2,i0j2);
  SetCR(1,2,i1j2);
  SetCR(2,2,i2j2);
  SetCR(3,2,i3j2);
  SetCR(0,3,i0j3);
  SetCR(1,3,i1j3);
  SetCR(2,3,i2j3);
  SetCR(3,3,i3j3);
}

Matrix44::Matrix44 ( MATRIX const* iMatrix ) {

  // Allocate our MATRIX object.
  m = MatrixIdentity( 4, NULL );

  // Copy the matrix.
  MatrixCopy( const_cast<MATRIX*>(iMatrix), m );
}

Matrix44::Matrix44 ( Matrix44 const& iMatrix ) {

  // Allocate our MATRIX object.
  m = MatrixIdentity( 4, NULL );

  // Copy the matrix.
  SetMatrix( iMatrix );
}

Matrix44::~Matrix44() {

  // Free the matrix.
  assert( m );
  MatrixFree( &m );
}

void
Matrix44::SetMatrix ( float i0j0, float i1j0, float i2j0, float i3j0,
                      float i0j1, float i1j1, float i2j1, float i3j1,
                      float i0j2, float i1j2, float i2j2, float i3j2,
                      float i0j3, float i1j3, float i2j3, float i3j3 ) {

  // Set the matrix elements.
  SetCR(0,0,i0j0);
  SetCR(1,0,i1j0);
  SetCR(2,0,i2j0);
  SetCR(3,0,i3j0);
  SetCR(0,1,i0j1);
  SetCR(1,1,i1j1);
  SetCR(2,1,i2j1);
  SetCR(3,1,i3j1);
  SetCR(0,2,i0j2);
  SetCR(1,2,i1j2);
  SetCR(2,2,i2j2);
  SetCR(3,2,i3j2);
  SetCR(0,3,i0j3);
  SetCR(1,3,i1j3);
  SetCR(2,3,i2j3);
  SetCR(3,3,i3j3);
}

void
Matrix44::SetMatrix ( MATRIX const* iMatrix ) {

  // Copy the matrix.
  MatrixCopy( const_cast<MATRIX*>(iMatrix), m );
}

void
Matrix44::SetMatrix ( Matrix44 const& iMatrix ) {

  // Copy the matrix.
  MatrixCopy( const_cast<MATRIX*>(iMatrix.GetMatrix()), m );
}

Matrix44&
Matrix44::operator= ( Matrix44 const& iMatrix ) {

  // Check for self-assigment.
  if( this == &iMatrix )
    return *this;
  
  // Copy the matrix.
  SetMatrix( iMatrix  );

  // Return ourselves.
  return *this;
}

void
Matrix44::MakeIdentity () {

  // Make our matrix the identity.
  MatrixIdentity( 4, m );
}

void
Matrix44::MakeRotation ( float const iCenterPoint[3],
                         float const iRotationVector[3],
                         float const iRadians ) {

  Point3<float> p( iCenterPoint );
  Point3<float> v( iRotationVector );

  // Rotate around y so that it lies on the x/y plane.
  double radsAroundY = atan2 ( v[2], v[0] );
  Matrix44 yRotation;
  yRotation.MakeYRotation( radsAroundY );

  // Rotate around z so that it lies on the x axis.
  Point3<float> vectorOnXY = yRotation * v;
  double radsAroundZ = -atan2 ( vectorOnXY[1], vectorOnXY[0] );
  Matrix44 zRotation;
  zRotation.MakeZRotation( radsAroundZ );

  // Translation matrix with the center point.
  Matrix44 trans;
  trans.SetMatrix ( 1, 0, 0, p[0],
                    0, 1, 0, p[1],
                    0, 0, 1, p[2],
                    0, 0, 0, 1 );

  // Rotation.
  Matrix44 rotation;
  rotation.MakeXRotation( iRadians );

  // Make all the inverse transformations.
  Matrix44 zRotationInv;
  zRotationInv.MakeInverseZRotation( radsAroundZ );

  Matrix44 yRotationInv;
  yRotationInv.MakeInverseYRotation( radsAroundY );

  Matrix44 transInv;
  transInv.SetMatrix ( 1, 0, 0, -p[0],
                       0, 1, 0, -p[1],
                       0, 0, 1, -p[2],
                       0, 0, 0, 1 );

  // Composition.
  Matrix44 composed = trans * yRotation * zRotation * rotation *
    zRotationInv * yRotationInv * transInv;

  // Set the result.
  SetMatrix( composed );
}

void
Matrix44::MakeXRotation ( float iRadians ) {

  SetMatrix( 1, 0, 0, 0,
             0, cos(iRadians), -sin(iRadians), 0,
             0, sin(iRadians), cos(iRadians), 0,
             0, 0, 0, 1 );
}

void
Matrix44::MakeYRotation ( float iRadians ) {

  SetMatrix( cos(iRadians), 0, sin(iRadians), 0,
             0, 1, 0, 0,
             -sin(iRadians), 0, cos(iRadians), 0,
             0, 0, 0, 1 );
}

void
Matrix44::MakeZRotation ( float iRadians ) {

  SetMatrix( cos(iRadians), -sin(iRadians), 0, 0,
             sin(iRadians), cos(iRadians), 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1 );
}

void
Matrix44::MakeInverseXRotation ( float iRadians ) {

  MakeXRotation( -iRadians );
}

void
Matrix44::MakeInverseYRotation ( float iRadians ) {

  MakeYRotation( -iRadians );
}

void
Matrix44::MakeInverseZRotation ( float iRadians ) {

  MakeZRotation( -iRadians );
}

Matrix44
Matrix44::ExtractRotation () const {

  // Create a matrix identital to this one but without the translation.
  Matrix44 thisWithoutTranslate;
  thisWithoutTranslate.SetMatrix( *this );
  thisWithoutTranslate.SetCR(3,0, 0);
  thisWithoutTranslate.SetCR(3,1, 0);
  thisWithoutTranslate.SetCR(3,2, 0);

  // We need to cancel the scale portion. Extract the scale from this
  // new matrix and divide the scale factors. If it was 0 to begin
  // with, set it to 1.
  Matrix44 scale = thisWithoutTranslate.ExtractScale();
  Matrix44 thisWithoutTranslateOrScale;
  for ( int n = 0; n <= 3; n++ ) {
    if ( FEQUAL( scale(n,n), 0 ) ) {
      thisWithoutTranslateOrScale.SetCR( n, n, 1.0 );
    } else {
      float unscaled = 1.0 / scale(n,n);
      thisWithoutTranslateOrScale.SetCR( n, n, unscaled );
    }
  }

  // Our rotation matrix is thisWithoutTranslateOrScale composed with
  // thisWithoutTranslate.
  return thisWithoutTranslate * thisWithoutTranslateOrScale;
}

Matrix44
Matrix44::ExtractScale () const {

  // Create a matrix identital to this one but without the translation.
  Matrix44 thisWithoutTranslate;
  thisWithoutTranslate.SetMatrix( *this );
  thisWithoutTranslate.SetCR(3,0, 0);
  thisWithoutTranslate.SetCR(3,1, 0);
  thisWithoutTranslate.SetCR(3,2, 0);

  // Create a 1,0,0 vector and multiply it with the new matrix. Get
  // the length of the result. This is the x scale factor.
  Point3<float> v( 1, 0, 0 );
  Point3<float> w;
  thisWithoutTranslate.MultiplyVector3( v.xyz(), w.xyz() );
  float xFactor = VectorOps::Length( w );

  // Do the same for y and z.
  v.Set( 0, 1, 0 );
  thisWithoutTranslate.MultiplyVector3( v.xyz(), w.xyz() );
  float yFactor = VectorOps::Length( w );

  v.Set( 0, 0, 1 );
  thisWithoutTranslate.MultiplyVector3( v.xyz(), w.xyz() );
  float zFactor = VectorOps::Length( w );

  // Now build the result.
  return Matrix44( xFactor, 0, 0, 0,
                   0, yFactor, 0, 0,
                   0, 0, zFactor, 0,
                   0, 0, 0, 1 );
}

Matrix44
Matrix44::ExtractTranslation () const {

  return Matrix44( 1, 0, 0, GetCR(3,0),
		   0, 1, 0, GetCR(3,1),
		   0, 0, 1, GetCR(3,2),
		   0, 0, 0, 1 );
}

void
Matrix44::ApplyTransformMatrix ( Matrix44 const& iTransform ) {

  // Get our transformation components.
  Matrix44 scale     = ExtractScale();
  Matrix44 translate = ExtractTranslation();
  Matrix44 rotate    = ExtractRotation();

  // Find the inverse of the translation and rotation.
  Matrix44 translateInv = translate.Inverse();
  Matrix44 rotateInv    = rotate.Inverse();

  // Extract the transforms to apply.
  Matrix44 scaleApply     = iTransform.ExtractScale();
  Matrix44 translateApply = iTransform.ExtractTranslation();
  Matrix44 rotateApply    = iTransform.ExtractRotation();

  // The new translation is composed with the rotation and inverse rotation.
  Matrix44 tmp1 = rotate * translateApply;
  Matrix44 translateNew = tmp1 * rotateInv;

  // Same with the rotation.
  Matrix44 tmp2 = rotate * rotateApply;
  Matrix44 rotateNew = tmp2 * rotateInv;

  // Now compose everything together.
  Matrix44 composed = translateNew * translate * 
    scaleApply * scale * rotateNew * rotate;

  // Set us.
  SetMatrix( composed );
}

void
Matrix44::MultiplyVector3 ( float const iVector[3], float oVector[3] ) const {

  // This is explicitly written out in an attempt to speed stuff
  // up. Otherwise there's nothing really tricky about it.
  float iX = iVector[0];
  float iY = iVector[1];
  float iZ = iVector[2];

  float m11 = GetCR(0,0);
  float m12 = GetCR(1,0);
  float m13 = GetCR(2,0);
  float m14 = GetCR(3,0);

  float m21 = GetCR(0,1);
  float m22 = GetCR(1,1);
  float m23 = GetCR(2,1);
  float m24 = GetCR(3,1);

  float m31 = GetCR(0,2);
  float m32 = GetCR(1,2);
  float m33 = GetCR(2,2);
  float m34 = GetCR(3,2);

  float a = m11 * iX;
  float b = m12 * iY;
  float c = m13 * iZ;
  float sum0 = a + b + c + m14;

  float d = m21 * iX;
  float e = m22 * iY;
  float f = m23 * iZ;
  float sum1 = d + e + f + m24;

  float g = m31 * iX;
  float h = m32 * iY;
  float i = m33 * iZ;
  float sum2 = g + h + i + m34;

  oVector[0] = sum0;
  oVector[1] = sum1;
  oVector[2] = sum2;
}

void
Matrix44::MultiplyVector3 ( int const iVector[3], float oVector[3] ) const {

  // This is explicitly written out in an attempt to speed stuff
  // up. Otherwise there's nothing really tricky about it.
  float iVectorF[3];
  iVectorF[0] = iVector[0];
  iVectorF[1] = iVector[1];
  iVectorF[2] = iVector[2];

  oVector[0] =
    GetCR(0,0) * iVectorF[0] +
    GetCR(1,0) * iVectorF[1] +
    GetCR(2,0) * iVectorF[2] +
    GetCR(3,0);
  oVector[1] =
    GetCR(0,1) * iVectorF[0] +
    GetCR(1,1) * iVectorF[1] +
    GetCR(2,1) * iVectorF[2] +
    GetCR(3,1);
  oVector[2] =
    GetCR(0,2) * iVectorF[0] +
    GetCR(1,2) * iVectorF[1] +
    GetCR(2,2) * iVectorF[2] +
    GetCR(3,2);
}

void
Matrix44::MultiplyVector3 ( float const iVector[3], int oVector[3] ) const {

  float vectorF[3];
  MultiplyVector3( iVector, vectorF );

  // This rounding is consistent with some other libutils stuff. Just
  // accept it.
  oVector[0] = (int) floor( vectorF[0] + 0.5 );
  oVector[1] = (int) floor( vectorF[1] + 0.5 );
  oVector[2] = (int) floor( vectorF[2] + 0.5 );
}

Matrix44
Matrix44::Inverse() const {

  // Allocate a temp copy. The temp is because MatrixInverse
  // might be unsafe and this function is const.
  MATRIX* temp = MatrixIdentity( 4, NULL );

  // Copy in the main matrix.
  MatrixCopy( m, temp );

  // Perform the inverse and check the result. This allocates mInv.
  MATRIX* inverse = MatrixInverse( temp, NULL );
  if ( NULL == inverse ) {
    
    // Free our temps.
    MatrixFree( &temp );
    MatrixFree( &inverse );

    // Wasn't invertible.
    cerr << "Couldn't invert matrix: " << endl << *this << endl;
    throw runtime_error("Couldn't invert matrix");
  }

  // Construct a return object.
  Matrix44 inverseM( inverse );

  // Free our temps.
  MatrixFree( &temp );
  MatrixFree( &inverse );

  // Return the result.
  return inverseM;
}

MATRIX const*
Matrix44::GetMatrix () const {

  return m;
}

inline Matrix44 operator*( Matrix44 const& m2,
                           Matrix44 const& m1 ) {


  // This is explicitly written out in an attempt to speed stuff
  // up. Otherwise there's nothing really tricky about it.
  float m00 =
    m1(0,0) * m2(0,0) +
    m1(0,1) * m2(1,0) +
    m1(0,2) * m2(2,0) +
    m1(0,3) * m2(3,0);
  float m10 =
    m1(1,0) * m2(0,0) +
    m1(1,1) * m2(1,0) +
    m1(1,2) * m2(2,0) +
    m1(1,3) * m2(3,0);
  float m20 =
    m1(2,0) * m2(0,0) +
    m1(2,1) * m2(1,0) +
    m1(2,2) * m2(2,0) +
    m1(2,3) * m2(3,0);
  float m30 =
    m1(3,0) * m2(0,0) +
    m1(3,1) * m2(1,0) +
    m1(3,2) * m2(2,0) +
    m1(3,3) * m2(3,0);

  float m01 =
    m1(0,0) * m2(0,1) +
    m1(0,1) * m2(1,1) +
    m1(0,2) * m2(2,1) +
    m1(0,3) * m2(3,1);
  float m11 =
    m1(1,0) * m2(0,1) +
    m1(1,1) * m2(1,1) +
    m1(1,2) * m2(2,1) +
    m1(1,3) * m2(3,1);
  float m21 =
    m1(2,0) * m2(0,1) +
    m1(2,1) * m2(1,1) +
    m1(2,2) * m2(2,1) +
    m1(2,3) * m2(3,1);
  float m31 =
    m1(3,0) * m2(0,1) +
    m1(3,1) * m2(1,1) +
    m1(3,2) * m2(2,1) +
    m1(3,3) * m2(3,1);

  float m02 =
    m1(0,0) * m2(0,2) +
    m1(0,1) * m2(1,2) +
    m1(0,2) * m2(2,2) +
    m1(0,3) * m2(3,2);
  float m12 =
    m1(1,0) * m2(0,2) +
    m1(1,1) * m2(1,2) +
    m1(1,2) * m2(2,2) +
    m1(1,3) * m2(3,2);
  float m22 =
    m1(2,0) * m2(0,2) +
    m1(2,1) * m2(1,2) +
    m1(2,2) * m2(2,2) +
    m1(2,3) * m2(3,2);
  float m32 =
    m1(3,0) * m2(0,2) +
    m1(3,1) * m2(1,2) +
    m1(3,2) * m2(2,2) +
    m1(3,3) * m2(3,2);

  float m03 =
    m1(0,0) * m2(0,3) +
    m1(0,1) * m2(1,3) +
    m1(0,2) * m2(2,3) +
    m1(0,3) * m2(3,3);
  float m13 =
    m1(1,0) * m2(0,3) +
    m1(1,1) * m2(1,3) +
    m1(1,2) * m2(2,3) +
    m1(1,3) * m2(3,3);
  float m23 =
    m1(2,0) * m2(0,3) +
    m1(2,1) * m2(1,3) +
    m1(2,2) * m2(2,3) +
    m1(2,3) * m2(3,3);
  float m33 =
    m1(3,0) * m2(0,3) +
    m1(3,1) * m2(1,3) +
    m1(3,2) * m2(2,3) +
    m1(3,3) * m2(3,3);

  return Matrix44( m00, m10, m20, m30,
                   m01, m11, m21, m31,
                   m02, m12, m22, m32,
                   m03, m13, m23, m33 );
};

inline Point3<float> operator* ( Matrix44 const& m,
				 Point3<float> const& p ) {

  return Point3<float> ( m(0,0)*p[0] + m(0,1)*p[1] + m(0,2)*p[2] + m(0,3),
                         m(1,0)*p[0] + m(1,1)*p[1] + m(1,2)*p[2] + m(1,3),
                         m(2,0)*p[0] + m(2,1)*p[1] + m(2,2)*p[2] + m(2,3) );
};

ostream&
operator <<  ( ostream& os, Matrix44 const& iMatrix ) {
  os << "Matrix44:" << endl;
  os << setw(6) << iMatrix(0,0) << " " << setw(6) << iMatrix(1,0) << " "
  << setw(6) << iMatrix(2,0) << " " << setw(6) << iMatrix(3,0) << endl;
  os << setw(6) << iMatrix(0,1) << " " << setw(6) << iMatrix(1,1) << " "
  << setw(6) << iMatrix(2,1) << " " << setw(6) << iMatrix(3,1) << endl;
  os << setw(6) << iMatrix(0,2) << " " << setw(6) << iMatrix(1,2) << " "
  << setw(6) << iMatrix(2,2) << " " << setw(6) << iMatrix(3,2) << endl;
  os << setw(6) << iMatrix(0,3) << " " << setw(6) << iMatrix(1,3) << " "
  << setw(6) << iMatrix(2,3) << " " << setw(6) << iMatrix(3,3) << endl;
  return os;
}
