#ifndef Transform44_h
#define Transform44_h

#include "string_fixed.h"
#include "Matrix44.h"
#include "Point3.h"

// Note that all 16 element matrices used in this class are in openGL
// style format:
// [ 0   4   8  12 ]
// [ 1   5   9  13 ]
// [ 2   6  10  14 ]
// [ 3   7  11  15 ]

class Transform44 : public DebugReporter {

  friend class Transform44Tester;

 public:

  Transform44();
  virtual ~Transform44();

  void SetMainTransform ( float i0j0, float i1j0, float i2j0, float i3j0,
			  float i0j1, float i1j1, float i2j1, float i3j1,
			  float i0j2, float i1j2, float i2j2, float i3j2,
			  float i0j3, float i1j3, float i2j3, float i3j3 );

  void SetMainTransform ( MATRIX* iMatrix );
  void SetMainTransform ( Matrix44& iMatrix );
  void SetMainTransform ( Transform44&  iTransform );

  void MakeIdentity ();

  void MakeRotation ( float iCenterPoint[3], 
		      float iRotationVector[3],
		      float iRadians );

  void LoadFromLTAFile ( std::string ifnLTA );

  void ApplyTransform ( Transform44& iTransform );

  void MultiplyVector3 ( float const iVector[3], float oVector[3] );
  void MultiplyVector3 ( int   const iVector[3], float oVector[3] );
  void MultiplyVector3 ( float const iVector[3], int   oVector[3] );
  void InvMultiplyVector3 ( float const iVector[3], float oVector[3] );
  void InvMultiplyVector3 ( int   const iVector[3], float oVector[3] );
  void InvMultiplyVector3 ( float const iVector[3], int   oVector[3] );

  float operator()( int iCol, int iRow ) {
    return m(iCol, iRow);
  }

  Transform44& operator=(Transform44& iTransform) {
    m.SetMatrix( iTransform.GetMainMatrix() );
    ValuesChanged();
    return *this;
  }

  Matrix44& GetMainMatrix () { return m; }

  Transform44& Inverse ();
  
 protected:

  virtual void ValuesChanged ();

  void CalculateInverse ();

  Matrix44 m;
  Matrix44 mInv;
};

// C = A * B
// Transform44 a;
// Transform44 b;
// Transform44 c = a * b;
Transform44& operator*(Transform44& m1, Transform44& m2);

std::ostream& operator << ( std::ostream&, Transform44& iTransform  );

#endif

