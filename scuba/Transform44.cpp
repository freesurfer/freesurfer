/**
 * @file  Transform44.cpp
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
 *    $Revision: 1.12 $
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
#include <stdexcept>

#include "Transform44.h"

extern "C" {
#include "registerio.h"
#include "transform.h"
}


using namespace std;


Transform44::Transform44() {
}

Transform44::Transform44 ( float i0j0,float i1j0,float i2j0,float i3j0,
                           float i0j1,float i1j1,float i2j1,float i3j1,
                           float i0j2,float i1j2,float i2j2,float i3j2,
                           float i0j3,float i1j3,float i2j3,float i3j3 ) :
  m( i0j0, i1j0, i2j0, i3j0,
     i0j1, i1j1, i2j1, i3j1,
     i0j2, i1j2, i2j2, i3j2,
     i0j3, i1j3, i2j3, i3j3 ) {

  // We set the main matrix so calculate our inverse.
  CalculateInverse();
}

Transform44::Transform44 ( MATRIX const* iMatrix ) :
  m( iMatrix ) {

  // We set the main matrix so calculate our inverse.
  CalculateInverse();
}

Transform44::Transform44 ( Matrix44 const& iMatrix ) :
  m( iMatrix ) {

  // We set the main matrix so calculate our inverse.
  CalculateInverse();
}

Transform44::~Transform44() {
}

void
Transform44::SetMainTransform ( float i0j0,float i1j0,float i2j0,float i3j0,
                                float i0j1,float i1j1,float i2j1,float i3j1,
                                float i0j2,float i1j2,float i2j2,float i3j2,
                                float i0j3,float i1j3,float i2j3,float i3j3) {

  // Set the main matrix.
  m.SetMatrix( i0j0, i1j0, i2j0, i3j0,
               i0j1, i1j1, i2j1, i3j1,
               i0j2, i1j2, i2j2, i3j2,
               i0j3, i1j3, i2j3, i3j3 );

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::SetMainTransform ( MATRIX const* iMatrix ) {

  // Set the main matrix.
  m.SetMatrix( iMatrix );

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::SetMainTransform ( Matrix44 const& iMatrix ) {

  // Set the main matrix.
  m.SetMatrix( iMatrix );

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::SetMainTransform ( Transform44 const& iTransform ) {

  // Set the main matrix.
  m.SetMatrix( iTransform.GetMainMatrix() );

  // Our values have changed.
  ValuesChanged();
}

Transform44& 
Transform44::operator= ( Transform44 const& iTransform ) {

  // Check for self assignment.
  if( this == &iTransform )
    return *this;

  // Set our main transform from this transform.
  SetMainTransform( iTransform );

  // Return a reference to ourselves.
  return *this;
}

void
Transform44::MakeIdentity () {

  // Make the main matrix the identity.
  m.MakeIdentity();

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::MakeRotation ( float const iCenterPoint[3],
                            float const iRotationVector[3],
                            float const iRadians ) {

  // Make the main matrix the identity.
  m.MakeRotation( iCenterPoint, iRotationVector, iRadians );

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::LoadFromLTAFile ( string const& ifnLTA ) {

  // Hack. Since LTAreadEx converts register.dat files to a
  // LINEAR_VOX_TO_VOX in the wrong coordinate space, we'll read
  // register.dat files in manually, and use LTAreadEx for the rest.

  // If we found "register.dat" in the file name...
  string::size_type rFind = ifnLTA.find( "register.dat", 0 );
  if ( rFind != string::npos ) {
					   
    // Attempt to read in the info with regio_read_register. We're
    // ignoring all those other values; we just want the matrix.
    char* sSubject;
    float inPlaneResolution;
    float betweenPlaneResolution;
    float intensity;
    MATRIX* registrationMatrix = NULL;
    int intConversionMethod;

    char* fnLTA = strdup( ifnLTA.c_str() );
    regio_read_register( fnLTA,
                         &sSubject, &inPlaneResolution,
                         &betweenPlaneResolution, &intensity,
                         &registrationMatrix, &intConversionMethod );
    free( fnLTA );
    if ( NULL == registrationMatrix ) {
      throw runtime_error( "Couldn't load registration." );
    }

    // Set our main matrix from the registration matrix we got.
    SetMainTransform( registrationMatrix );

    // Free up the stuff we got from regio_read_register.
    MatrixFree( &registrationMatrix );
    free( sSubject );

  } else {

    // Read in the LTA object.
    LTA* lta = LTAreadEx( ifnLTA.c_str() );
    if ( NULL == lta ) {
      throw runtime_error( "Couldn't load LTA." );
    }

    // Switch on the type and output it.
    switch ( lta->type ) {
    case LINEAR_VOX_TO_VOX:
      cout << "LTA type is LINEAR_VOX_TO_VOX" << endl;
      break;
    case LINEAR_RAS_TO_RAS:
      cout << "LTA type is LINEAR_RAS_TO_RAS" << endl;
      break;
    default:
      cout << "LTA type is unkown" << endl;
      break;
    }
    
    // Get the first transform and matrix (always default in a linear
    // transform).
    LT* transform = &lta->xforms[0];
    MATRIX* matrix = transform->m_L;

    // Set our main matrix from it.
    SetMainTransform( matrix );

    // Free the LTA.
    LTAfree( &lta );
  }

}

void
Transform44::ApplyTransform ( Transform44 const& iTransform ) {

  // Apply this transform to our main matrix.
  m.ApplyTransformMatrix( iTransform.GetMainMatrix() );

  // Our values have changed.
  ValuesChanged();
}

void
Transform44::MultiplyVector3 ( float const iVector[3], 
			       float oVector[3] ) const {

  // Transform by our main matrix.
  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::MultiplyVector3 ( int const iVector[3],
			       float oVector[3] ) const {

  // Transform by our main matrix.
  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::MultiplyVector3 ( float const iVector[3], 
			       int oVector[3] ) const {

  // Transform by our main matrix.
  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( float const iVector[3],
				  float oVector[3] ) const {

  // Transform by our inverse matrix.
  mInv.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( int const iVector[3], 
				  float oVector[3] ) const {

  // Transform by our inverse matrix.
  mInv.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( float const iVector[3], 
				  int oVector[3] ) const {

  // Transform by our inverse matrix.
  mInv.MultiplyVector3( iVector, oVector );
}

Transform44
Transform44::Inverse() const {

  // Get the inverse of our main matrix and make a new transform using
  // that as its main transform.
  return Transform44( m.Inverse() );
}

Matrix44 const&
Transform44::GetMainMatrix () const {

  return m;
}

void
Transform44::ValuesChanged () {

  // If our values changed, recalc the inverse matrix.
  CalculateInverse();
}

void
Transform44::CalculateInverse () {

  // Set our inverse matrix to the inverse of the main matrix.
  // it.
  mInv.SetMatrix( m.Inverse() );
}

Transform44 operator* ( Transform44 const& t1,
			Transform44 const& t2 ) {

  // Use the Matrix44 operator* to multiply them and return a new
  // Transform44 using the result.
  return Transform44( t1.GetMainMatrix() * t2.GetMainMatrix() );
}

ostream&
operator <<  ( ostream& os, Transform44 const& iTransform ) {
  os << "Transform44:" << endl;
  os << setw(6) << iTransform(0,0) << " " << setw(6) << iTransform(1,0) << " "
  << setw(6) << iTransform(2,0) << " " << setw(6) << iTransform(3,0) << endl;
  os << setw(6) << iTransform(0,1) << " " << setw(6) << iTransform(1,1) << " "
  << setw(6) << iTransform(2,1) << " " << setw(6) << iTransform(3,1) << endl;
  os << setw(6) << iTransform(0,2) << " " << setw(6) << iTransform(1,2) << " "
  << setw(6) << iTransform(2,2) << " " << setw(6) << iTransform(3,2) << endl;
  os << setw(6) << iTransform(0,3) << " " << setw(6) << iTransform(1,3) << " "
  << setw(6) << iTransform(2,3) << " " << setw(6) << iTransform(3,3) << endl;
  return os;
}
