/**
 * @file  Transform44.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.5 $
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
  Transform44 ( float i0j0, float i1j0, float i2j0, float i3j0,
                float i0j1, float i1j1, float i2j1, float i3j1,
                float i0j2, float i1j2, float i2j2, float i3j2,
                float i0j3, float i1j3, float i2j3, float i3j3 );
  Transform44 ( MATRIX* iMatrix );
  Transform44 ( Matrix44& iMatrix );
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

  Matrix44& GetMainMatrix () {
    return m;
  }

  Transform44 Inverse ();

protected:

  virtual void ValuesChanged ();

  void CalculateInverse ();

  Matrix44 m;
  Matrix44 mInv;

  MATRIX* mTmp;
};

// C = A * B
// Transform44 a;
// Transform44 b;
// Transform44 c = a * b;
Transform44 operator*(Transform44& m1, Transform44& m2);

std::ostream& operator << ( std::ostream&, Transform44& iTransform  );

#endif

