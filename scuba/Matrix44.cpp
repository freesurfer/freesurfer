#include "string_fixed.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "Matrix44.h"
#include "VectorOps.h"
extern "C" {
  #include "macros.h"
}

using namespace std;

Matrix44::Matrix44() {
  m = MatrixIdentity( 4, NULL );
}

Matrix44::~Matrix44() {
  MatrixFree( &m );
}

void
Matrix44::SetMatrix ( float i0j0, float i1j0, float i2j0, float i3j0,
			  float i0j1, float i1j1, float i2j1, float i3j1,
			  float i0j2, float i1j2, float i2j2, float i3j2,
			  float i0j3, float i1j3, float i2j3, float i3j3 ) {

  SetCR(0,0,i0j0);  SetCR(1,0,i1j0);  SetCR(2,0,i2j0);  SetCR(3,0,i3j0);
  SetCR(0,1,i0j1);  SetCR(1,1,i1j1);  SetCR(2,1,i2j1);  SetCR(3,1,i3j1);
  SetCR(0,2,i0j2);  SetCR(1,2,i1j2);  SetCR(2,2,i2j2);  SetCR(3,2,i3j2);
  SetCR(0,3,i0j3);  SetCR(1,3,i1j3);  SetCR(2,3,i2j3);  SetCR(3,3,i3j3);
}

void
Matrix44::SetMatrix ( MATRIX* iMatrix ) {
  
  MatrixCopy( iMatrix, m );
}

void
Matrix44::SetMatrix ( Matrix44& iMatrix ) {

  MatrixCopy( iMatrix.GetMatrix(), m );
}

void 
Matrix44::MakeIdentity () {

  MatrixIdentity( 4, m );
}

void
Matrix44::MakeRotation ( float iCenterPoint[3], 
			       float iRotationVector[3],
			       float iRadians ) {

  Point3<float> p( iCenterPoint );
  Point3<float> v( iRotationVector );

  double radsAroundY = atan2 ( v[2], v[0] );
  Matrix44 yRotation;
  yRotation.MakeYRotation( radsAroundY );
  Point3<float> vectorOnXY = yRotation * v;

  double radsAroundZ = -atan2 ( vectorOnXY[1], vectorOnXY[0] );

  Matrix44 trans;
  trans.SetMatrix ( 1, 0, 0, p[0],
		       0, 1, 0, p[1],
		       0, 0, 1, p[2],
		       0, 0, 0, 1 );

  // rotate around y so that it lies on the x/y plane.
  yRotation.MakeYRotation( radsAroundY );
  
  // rotate around z so that it lies on the x axis.
  Matrix44 zRotation;
  zRotation.MakeZRotation( radsAroundZ );
  
  Matrix44 rotation;
  rotation.MakeXRotation( iRadians );

  Matrix44 zRotationInv;
  zRotationInv.MakeInverseZRotation( radsAroundZ );
  
  Matrix44 yRotationInv;
  yRotationInv.MakeInverseYRotation( radsAroundY );
  
  Matrix44 transInv;
  transInv.SetMatrix ( 1, 0, 0, -p[0],
		       0, 1, 0, -p[1],
			  0, 0, 1, -p[2],
		       0, 0, 0, 1 );

  /*  Matrix44 final( transInv * 
		  yRotationInv * zRotationInv * 
		  rotation * 
		  zRotation * yRotation *
		  trans ); */
  Matrix44 final( trans * 
		  yRotation * zRotation * 
		  rotation * 
		  zRotationInv * yRotationInv *
		  transInv );

  SetMatrix( final.GetMatrix() );
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

Matrix44& 
Matrix44::ExtractTranslation () {

  Matrix44* translation = new Matrix44();
  translation->SetMatrix( 1, 0, 0, GetCR(3,0),
			  0, 1, 0, GetCR(3,1),
			  0, 0, 1, GetCR(3,2),
			  0, 0, 0, 1 );
  return *translation;
}

Matrix44&
Matrix44::ExtractScale () {

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
  float xFactor = Length( w );

  // Do the same for y and z.
  v.Set( 0, 1, 0 );
  thisWithoutTranslate.MultiplyVector3( v.xyz(), w.xyz() );
  float yFactor = Length( w );
  
  v.Set( 0, 0, 1 );
  thisWithoutTranslate.MultiplyVector3( v.xyz(), w.xyz() );
  float zFactor = Length( w );
  
  // Now build the result.
  Matrix44* scale = new Matrix44();
  scale->SetMatrix( xFactor, 0, 0, 0,
		    0, yFactor, 0, 0,
		    0, 0, zFactor, 0,
		    0, 0, 0, 1 );
  return *scale;
}


Matrix44&
Matrix44::ExtractRotation () {

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
  for( int n = 0; n <= 3; n++ ) {
    if( FEQUAL( scale(n,n), 0 ) ) {
      thisWithoutTranslateOrScale.SetCR( n, n, 1.0 );
    } else {
      float unscaled = 1.0 / scale(n,n);
      thisWithoutTranslateOrScale.SetCR( n, n, unscaled );
    }
  }

  // Our rotation matrix is thisWithoutTranslateOrScale composed with
  // thisWithoutTranslate.
  Matrix44& rotation = thisWithoutTranslate * thisWithoutTranslateOrScale;

  return rotation;
}


void 
Matrix44::MultiplyVector3 ( float const iVector[3], float oVector[3] ) {

#if 0
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
  int sumI0 = (int) sum0;
  
  float d = m21 * iX;
  float e = m22 * iY;
  float f = m23 * iZ;
  float sum1 = d + e + f + m24;
  int sumI1 = (int) sum1;

  float g = m31 * iX;
  float h = m32 * iY;
  float i = m33 * iZ;
  float sum2 = g + h + i + m34;
  int sumI2 = (int) sum2;

  oVector[0] = sumI0;
  oVector[1] = sumI1;
  oVector[2] = sumI2;

#else
  
  oVector[0] =
    GetCR(0,0) * iVector[0] +
    GetCR(1,0) * iVector[1] +
    GetCR(2,0) * iVector[2] +
    GetCR(3,0);
  oVector[1] =
    GetCR(0,1) * iVector[0] +
    GetCR(1,1) * iVector[1] +
    GetCR(2,1) * iVector[2] +
    GetCR(3,1);
  oVector[2] =
    GetCR(0,2) * iVector[0] +
    GetCR(1,2) * iVector[1] +
    GetCR(2,2) * iVector[2] +
    GetCR(3,2);
#endif
}

void
Matrix44::ApplyTransformMatrix ( Matrix44& iTransform ) {

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

  //  cerr << "----------------------------------------" << endl;

  // The new translation is composed with the rotation and inverse rotation.
  Matrix44 translateNew = rotate * translateApply * rotateInv;

  // Same with the rotation.
  Matrix44 rotateNew = rotate * rotateApply * rotateInv;

#if 0
  cerr << "this " << *this << endl;
  cerr << "applying " << iTransform << endl;

  cerr << "scale " << scale << endl;
  cerr << "translate " << translate << endl;
  cerr << "rotate " << rotate << endl;
  cerr << "scaleApply " << scaleApply << endl;
  cerr << "translateApply " << translateApply << endl;
  cerr << "rotateApply " << rotateApply << endl;
  cerr << "translateNew " << translateNew << endl;
  cerr << "rotateNew " << rotateNew << endl;
#endif

  // Now compose everything together.
  Matrix44 t = translateNew * translate * scaleApply * scale * rotateNew * rotate;

#if 0
  cerr << "rotate";
  cerr << rotate << endl;

  cerr << "rotateNew * rotate";
  cerr << rotateNew * rotate << endl;

  cerr << "scale * rotateNew * rotate";
  cerr << scale * rotateNew * rotate << endl;

  cerr << "scaleApply * scale * rotateNew * rotate";
  cerr << scaleApply * scale * rotateNew * rotate << endl;

  cerr << "translate * scaleApply * scale * rotateNew * rotate";
  cerr << translate * scaleApply * scale * rotateNew * rotate << endl;

  cerr << "translateNew * translate * scaleApply * scale * rotateNew * rotate";
  cerr << translateNew * translate * scaleApply * scale * rotateNew * rotate << endl;
#endif

  // Set us.
  SetMatrix( t );
}

void 
Matrix44::MultiplyVector3 ( int const iVector[3], float oVector[3] ) {

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
Matrix44::MultiplyVector3 ( float const iVector[3], int oVector[3] ) {

#if 1
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
  float sum00 = a + b;
  float sum01 = c + m14;
#ifdef RINT
  int sumI0 = (int) rint(sum00 + sum01);
#else
  int sumI0 = (int) (sum00 + sum01);
#endif

  float d = m21 * iX;
  float e = m22 * iY;
  float f = m23 * iZ;
  float sum10 = d + e;
  float sum11 = f + m24;
#ifdef RINT
  int sumI1 = (int) rint(sum10 + sum11);
#else
  int sumI1 = (int) (sum10 + sum11);
#endif

  float g = m31 * iX;
  float h = m32 * iY;
  float i = m33 * iZ;
  float sum20 = g + h;
  float sum21 = i + m34;
#ifdef RINT
  int sumI2 = (int) rint(sum20 + sum21);
#else
  int sumI2 = (int) (sum20 + sum21);
#endif

  oVector[0] = sumI0;
  oVector[1] = sumI1;
  oVector[2] = sumI2;

#else
  oVector[0] = (int) ( GetCR(0,0) * iVector[0] +
		       GetCR(1,0) * iVector[1] +
		       GetCR(2,0) * iVector[2] +
		       GetCR(3,0) );
  oVector[1] = (int) ( GetCR(0,1) * iVector[0] +
		       GetCR(1,1) * iVector[1] +
		       GetCR(2,1) * iVector[2] +
		       GetCR(3,1) );
  oVector[2] = (int) ( GetCR(0,2) * iVector[0] +
		       GetCR(1,2) * iVector[1] +
		       GetCR(2,2) * iVector[2] +
		       GetCR(3,2) );
#endif
}

inline Matrix44& operator*( Matrix44& m2, 
			    Matrix44& m1 ) {


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

  //  MATRIX* mult = MatrixMultiply( m1.GetMatrix(), m2.GetMatrix(), NULL );
  Matrix44* result = new Matrix44();
  //  result->SetMatrix( mult );
  //  MatrixFree( &mult );
  result->SetMatrix( m00, m10, m20, m30,
		     m01, m11, m21, m31,
		     m02, m12, m22, m32,
		     m03, m13, m23, m33 );
  return *result;
};

Matrix44&
Matrix44::Inverse() {

  Matrix44* inverse = new Matrix44();
  MATRIX* inv = MatrixInverse( m, NULL );
  if( NULL == inv ) { 
    cerr << "Couldn't invert matrix: " << endl << *this << endl;
    throw runtime_error("Couldn't invert matrix"); 
  }
  inverse->SetMatrix( inv );
  MatrixFree( &inv );
  return *inverse;
}

// This works for Point3<float>s as vectors or points.
inline Point3<float>& operator*(Matrix44& m,
				Point3<float>& p) {

  Point3<float>* result = 
    new Point3<float> ( m(0,0)*p[0] + m(0,1)*p[1] + m(0,2)*p[2] + m(0,3),
			m(1,0)*p[0] + m(1,1)*p[1] + m(1,2)*p[2] + m(1,3),
			m(2,0)*p[0] + m(2,1)*p[1] + m(2,2)*p[2] + m(2,3) );
  return *result;
};


ostream& 
operator <<  ( ostream& os, Matrix44& iMatrix ) { 
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
