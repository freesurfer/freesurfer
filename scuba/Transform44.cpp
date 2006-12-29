/**
 * @file  Transform44.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.10 $
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


#include "string_fixed.h"
#include <iomanip>
#include <stdexcept>
#include "Transform44.h"
extern "C" {
#include "transform.h"
#include "registerio.h"
}

using namespace std;


Transform44::Transform44() {
  mTmp = MatrixIdentity( 4, NULL );
}

Transform44::Transform44 ( float i0j0,float i1j0,float i2j0,float i3j0,
                           float i0j1,float i1j1,float i2j1,float i3j1,
                           float i0j2,float i1j2,float i2j2,float i3j2,
                           float i0j3,float i1j3,float i2j3,float i3j3) {
  mTmp = MatrixIdentity( 4, NULL );
  m.SetMatrix( i0j0, i1j0, i2j0, i3j0,
               i0j1, i1j1, i2j1, i3j1,
               i0j2, i1j2, i2j2, i3j2,
               i0j3, i1j3, i2j3, i3j3 );
  CalculateInverse();
}

Transform44::Transform44 ( MATRIX* iMatrix ) {
  mTmp = MatrixIdentity( 4, NULL );
  m.SetMatrix( iMatrix );
  CalculateInverse();
}

Transform44::Transform44 ( Matrix44& iMatrix ) {
  mTmp = MatrixIdentity( 4, NULL );
  m.SetMatrix( iMatrix );
  CalculateInverse();
}

Transform44::~Transform44() {
  MatrixFree( &mTmp );
}

void
Transform44::SetMainTransform ( float i0j0,float i1j0,float i2j0,float i3j0,
                                float i0j1,float i1j1,float i2j1,float i3j1,
                                float i0j2,float i1j2,float i2j2,float i3j2,
                                float i0j3,float i1j3,float i2j3,float i3j3) {

  m.SetMatrix( i0j0, i1j0, i2j0, i3j0,
               i0j1, i1j1, i2j1, i3j1,
               i0j2, i1j2, i2j2, i3j2,
               i0j3, i1j3, i2j3, i3j3 );

  ValuesChanged();
}

void
Transform44::SetMainTransform ( MATRIX* iMatrix ) {

  m.SetMatrix( iMatrix );
  ValuesChanged();
}

void
Transform44::SetMainTransform ( Matrix44& iMatrix ) {

  m.SetMatrix( iMatrix );
  ValuesChanged();
}

void
Transform44::SetMainTransform ( Transform44&  iTransform ) {

  m.SetMatrix( iTransform.GetMainMatrix() );
  ValuesChanged();
}

void
Transform44::MakeIdentity () {

  m.MakeIdentity();
  ValuesChanged();
}

void
Transform44::MakeRotation ( float iCenterPoint[3],
                            float iRotationVector[3],
                            float iRadians ) {

  m.MakeRotation( iCenterPoint, iRotationVector, iRadians );
  ValuesChanged();
}

void
Transform44::LoadFromLTAFile ( string ifnLTA ) {

  /* Hack. Since LTAreadEx converts register.dat files to a
     LINEAR_VOX_TO_VOX in the wrong coordinate space, we'll read
     register.dat files in manually, and use LTAreadEx for the
     rest. */
  string::size_type rFind = ifnLTA.find( "register.dat", 0 );
  if ( rFind != string::npos ) {

    cerr << "found register.dat in " << ifnLTA << endl;

    char fnLTA[1000];
    char* sSubject;
    float inPlaneResolution;
    float betweenPlaneResolution;
    float intensity;
    MATRIX* registrationMatrix = NULL;
    int intConversionMethod;

    strcpy( fnLTA, ifnLTA.c_str() );
    regio_read_register( fnLTA,
                         &sSubject, &inPlaneResolution,
                         &betweenPlaneResolution, &intensity,
                         &registrationMatrix, &intConversionMethod );
    if ( NULL == registrationMatrix ) {
      throw runtime_error( "Couldn't load registration." );
    }

    SetMainTransform( registrationMatrix );

    MatrixFree( &registrationMatrix );
    free( sSubject );

  } else {

    LTA* lta = LTAreadEx( ifnLTA.c_str() );
    if ( NULL == lta ) {
      throw runtime_error( "Couldn't load LTA." );
    }
    switch ( lta->type ) {
    case LINEAR_VOX_TO_VOX:
      cerr << "LTA type is LINEAR_VOX_TO_VOX" << endl;
      break;
    case LINEAR_RAS_TO_RAS:
      cerr << "LTA type is LINEAR_RAS_TO_RAS" << endl;
      break;
    default:
      cerr << "LTA type is unkown" << endl;
      break;
    }
    LT* transform = &lta->xforms[0];

    MATRIX* matrix = transform->m_L;

    SetMainTransform( matrix );

    LTAfree( &lta );
  }

}

void
Transform44::ApplyTransform ( Transform44& iTransform ) {

  m.ApplyTransformMatrix( iTransform.GetMainMatrix() );
  ValuesChanged();
}

void
Transform44::MultiplyVector3 ( float const iVector[3], float oVector[3] ) {

  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::MultiplyVector3 ( int const iVector[3], float oVector[3] ) {

  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::MultiplyVector3 ( float const iVector[3], int oVector[3] ) {

  m.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( float const iVector[3], float oVector[3] ) {

  mInv.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( int const iVector[3], float oVector[3] ) {

  mInv.MultiplyVector3( iVector, oVector );
}

void
Transform44::InvMultiplyVector3 ( float const iVector[3], int oVector[3] ) {

  mInv.MultiplyVector3( iVector, oVector );
}

void
Transform44::ValuesChanged () {

  CalculateInverse();
}

void
Transform44::CalculateInverse () {
  Matrix44 inverse = m.Inverse();
  mInv.SetMatrix( inverse );
}

Transform44
Transform44::Inverse() {

  // Make a new transform and set its main transform to _our_ inverse
  // transform.
  mTmp = MatrixInverse( m.GetMatrix(), mTmp );
  return Transform44( mTmp );
}

Transform44 operator*( Transform44& t1,
                       Transform44& t2 ) {

  Matrix44& m1 = t1.GetMainMatrix();
  Matrix44& m2 = t2.GetMainMatrix();
  Matrix44 mult = m1 * m2;
  return Transform44( mult );
}

ostream&
operator <<  ( ostream& os, Transform44& iTransform ) {
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
