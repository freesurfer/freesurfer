/**
 * @file  Transform44.h
 * @brief A 4x4 real linear transform object
 *
 * A linear transform object that provides a simplified interface to
 * Matrix44, mainly allowing vector mulitplication and automatically
 * caching inverse transforms.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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


#ifndef Transform44_h
#define Transform44_h

#include "string_fixed.h"
extern "C" {
#include "matrix.h"
}

#include "DebugReporter.h"
#include "Matrix44.h"

class Transform44 : public DebugReporter {

  friend class Transform44Tester;

public:
  
  // Initializes the transform to the identity.
  Transform44();

  // Explicity set the main transform matrix.
  Transform44 ( float i0j0, float i1j0, float i2j0, float i3j0,
                float i0j1, float i1j1, float i2j1, float i3j1,
                float i0j2, float i1j2, float i2j2, float i3j2,
                float i0j3, float i1j3, float i2j3, float i3j3 );

  // Initialize the main transform matrix contents from an existing
  // MATRIX, Matrix44, or Transform44.
  Transform44 ( MATRIX const* iMatrix );
  Transform44 ( Matrix44 const& iMatrix );
  Transform44 ( Transform44 const& iTransform );

  virtual ~Transform44();

  // Accessor.
  inline float const& operator() ( int iCol, int iRow ) const {
    return m(iCol,iRow);
  }

  // Settors.
  // Explicitly set element of main transform matirx.
  void SetMainTransform ( float i0j0, float i1j0, float i2j0, float i3j0,
                          float i0j1, float i1j1, float i2j1, float i3j1,
                          float i0j2, float i1j2, float i2j2, float i3j2,
                          float i0j3, float i1j3, float i2j3, float i3j3 );

  // Copy contents from an existing MATRIX, Matrix44, or Transform44.
  void SetMainTransform ( MATRIX const* iMatrix );
  void SetMainTransform ( Matrix44 const& iMatrix );
  void SetMainTransform ( Transform44 const&  iTransform );
  Transform44& operator= (Transform44 const& iTransform );

  // Set individual elements. Does not check bounds!
  inline float& operator() ( int iCol, int iRow ) {
    return m(iCol,iRow);
  }

  // Set to the identity.
  void MakeIdentity ();

  // Make this the rotation matrix for a given center point, rotation
  // vector, and number of radians. 
  void MakeRotation ( float const iCenterPoint[3],
                      float const iRotationVector[3],
                      float const iRadians );

  // Load contents from an LTA file.
  void LoadFromLTAFile ( std::string const& ifnLTA );

  // Applies a transform to this transform.
  void ApplyTransform ( Transform44 const& iTransform );

  // Transform a 3 component vector or point and put the results in
  // the output argument. Can also do the inverse transform.
  void MultiplyVector3 ( float const iVector[3], float oVector[3] ) const;
  void MultiplyVector3 ( int   const iVector[3], float oVector[3] ) const;
  void MultiplyVector3 ( float const iVector[3], int   oVector[3] ) const;
  void InvMultiplyVector3 ( float const iVector[3], float oVector[3] ) const;
  void InvMultiplyVector3 ( int   const iVector[3], float oVector[3] ) const;
  void InvMultiplyVector3 ( float const iVector[3], int   oVector[3] ) const;

  // Return the inverse of the transform. Doesn't modify this transform.
  Transform44 Inverse () const;

  // Direct access to the main matrix of the transform.
  Matrix44 const& GetMainMatrix () const;

protected:

  // Called when values have changed.
  virtual void ValuesChanged ();

  // Calculates the inverse from the main matrix.
  void CalculateInverse ();

  Matrix44 m;     // main matrix
  Matrix44 mInv;  // inverse
};

// Multiplication operator for matrices. You can use it this way:
// C = A * B
// Transform44 a;
// Transform44 b;
// Transform44 c = a * b;
Transform44 operator* ( Transform44 const& m1, Transform44 const& m2 );

std::ostream& operator << ( std::ostream&, Transform44 const& iTransform  );

#endif

