/**
 * @file  Matrix44.h
 * @brief A 4x4 real matrix
 *
 * A wrapper around the libutils MATRIX class for a 4x4 real number
 * matrix. Uses Matrix* functions internally when possible but also
 * includes some geometry functions. Note that all 16 element matrices
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
 *    $Revision: 1.6 $
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


#ifndef Matrix44_h
#define Matrix44_h

#include "string_fixed.h"
extern "C" {
#include "matrix.h"
}

#include "DebugReporter.h"

template <typename T> class Point3;

class Matrix44 : public DebugReporter {

  friend class Matrix44Tester;

public:

  // Initializes the matrix to the identity.
  Matrix44();

  // Explicity set all elements.
  Matrix44( float i0j0, float i1j0, float i2j0, float i3j0,
            float i0j1, float i1j1, float i2j1, float i3j1,
            float i0j2, float i1j2, float i2j2, float i3j2,
            float i0j3, float i1j3, float i2j3, float i3j3 );

  // Initialize with the contents from an existing MATRIX or Matrix44.
  Matrix44 ( MATRIX const* iMatrix );
  Matrix44 ( Matrix44 const& );

  virtual ~Matrix44();

  // Accessors.
  // Get individual elements. Does not check bounds!
  inline float const& GetCR ( int iCol, int iRow ) const {
    return *MATRIX_RELT(m,(iRow+1),(iCol+1));
  }
  float const& operator() ( int iCol, int iRow ) const {
    return *MATRIX_RELT(m,(iRow+1),(iCol+1));
  }

  // Settors.
  // Explicitly set all elements.
  void SetMatrix ( float i0j0, float i1j0, float i2j0, float i3j0,
                   float i0j1, float i1j1, float i2j1, float i3j1,
                   float i0j2, float i1j2, float i2j2, float i3j2,
                   float i0j3, float i1j3, float i2j3, float i3j3 );

  // Copy contents from an existing MATRIX or Matrix44.
  void SetMatrix ( MATRIX const* iMatrix );
  void SetMatrix ( Matrix44 const& iMatrix );
  Matrix44& operator= ( Matrix44 const& iMatrix );
  
  // Set individual elements. Does not check bounds!
  inline void SetCR ( int iCol, int iRow, float iValue ) {
    *MATRIX_RELT(m,(iRow+1),(iCol+1)) = iValue;
  }
  float& operator()( int iCol, int iRow ) {
    return *MATRIX_RELT(m,(iRow+1),(iCol+1));
  }

  // Set to the identity.
  void MakeIdentity ();

  // Make this the rotation matrix for a given center point, rotation
  // vector, and number of radians. 
  void MakeRotation ( float const iCenterPoint[3],
                      float const iRotationVector[3],
                      float const iRadians );

  // Make the rotation matrix around specific axes.
  void MakeXRotation ( float iRadians );
  void MakeYRotation ( float iRadians );
  void MakeZRotation ( float iRadians );
  void MakeInverseXRotation ( float iRadians );
  void MakeInverseYRotation ( float iRadians );
  void MakeInverseZRotation ( float iRadians );

  // Operations.
  // Extract various components from matrices and return a new matrix.
  Matrix44 ExtractRotation () const;
  Matrix44 ExtractScale () const;
  Matrix44 ExtractTranslation () const;

  // Applies a transform matrix to this matrix.
  void ApplyTransformMatrix ( Matrix44 const& iMatrix );

  // Multiply a 3 component vector with this matrix and put the
  // results in the output argument.
  void MultiplyVector3 ( float const iVector[3], float oVector[3] ) const;
  void MultiplyVector3 ( int   const iVector[3], float oVector[3] ) const;
  void MultiplyVector3 ( float const iVector[3], int   oVector[3] ) const;

  // Return the inverse of the matrix. Doesn't modify this matrix.
  Matrix44 Inverse () const;

  // Should only be used when we need the MATRIX for libutils
  // functions. 
  MATRIX const* GetMatrix () const;

protected:

  MATRIX* m; 
};

// Multiplication operator for matrices. You can use it this way:
// C = A * B
// Matrix44 a;
// Matrix44 b;
// Matrix44 c = a * b;
Matrix44 operator* ( Matrix44 const& m1, Matrix44 const& m2 );

// This works for Point3<float>s as vectors or points.
inline Point3<float> operator* (Matrix44 const& m, Point3<float> const& p );

std::ostream& operator << ( std::ostream&, Matrix44 const& iMatrix   );

#endif

