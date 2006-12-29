/**
 * @file  Matrix44.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
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


#ifndef Matrix44_h
#define Matrix44_h

#include "string_fixed.h"
extern "C" {
#include "matrix.h"
}
#include "DebugReporter.h"
#include "Point3.h"

// Note that all 16 element matrices used in this class are in openGL
// style format:
// [ 0   4   8  12 ]
// [ 1   5   9  13 ]
// [ 2   6  10  14 ]
// [ 3   7  11  15 ]

class Matrix44 : public DebugReporter {

  friend class Matrix44Tester;

public:

  Matrix44();
  Matrix44( float i0j0, float i1j0, float i2j0, float i3j0,
            float i0j1, float i1j1, float i2j1, float i3j1,
            float i0j2, float i1j2, float i2j2, float i3j2,
            float i0j3, float i1j3, float i2j3, float i3j3 );
  Matrix44 ( MATRIX* iMatrix );
  virtual ~Matrix44();

  void SetMatrix ( float i0j0, float i1j0, float i2j0, float i3j0,
                   float i0j1, float i1j1, float i2j1, float i3j1,
                   float i0j2, float i1j2, float i2j2, float i3j2,
                   float i0j3, float i1j3, float i2j3, float i3j3 );

  void SetMatrix ( MATRIX* iMatrix );

  void SetMatrix ( Matrix44& iMatrix );

  void MakeIdentity ();

  void MakeRotation ( float iCenterPoint[3],
                      float iRotationVector[3],
                      float iRadians );

  void MakeXRotation ( float iRadians );
  void MakeYRotation ( float iRadians );
  void MakeZRotation ( float iRadians );
  void MakeInverseXRotation ( float iRadians );
  void MakeInverseYRotation ( float iRadians );
  void MakeInverseZRotation ( float iRadians );


  Matrix44 ExtractRotation ();
  Matrix44 ExtractScale ();
  Matrix44 ExtractTranslation ();


  void ApplyTransformMatrix ( Matrix44& iMatrix );

  void MultiplyVector3 ( float const iVector[3], float oVector[3] );
  void MultiplyVector3 ( int   const iVector[3], float oVector[3] );
  void MultiplyVector3 ( float const iVector[3], int   oVector[3] );

  MATRIX* GetMatrix () {
    return m;
  }

  float operator()( int iCol, int iRow ) {
    return GetCR(iCol, iRow);
  }

  Matrix44& operator=(Matrix44& m) {
    SetMatrix( m );
    return *this;
  }

  Matrix44 Inverse ();

  inline void SetCR ( int iCol, int iRow, float iValue ) {
    *MATRIX_RELT(m,(iRow+1),(iCol+1)) = iValue;
  }

  inline float GetCR ( int iCol, int iRow ) {
    return *MATRIX_RELT(m,(iRow+1),(iCol+1));
  }
protected:

  MATRIX* m;
  MATRIX* mTmp;
};

// C = A * B
// Matrix44 a;
// Matrix44 b;
// Matrix44 c = a * b;
Matrix44 operator*(Matrix44& m1, Matrix44& m2);
inline Point3<float> operator*(Matrix44& m, Point3<float>& p);

std::ostream& operator << ( std::ostream&, Matrix44& iMatrix   );

#endif

