#include <iomanip>
#include <stdexcept>
#include "Transform44.h"
extern "C" {
#include "transform.h"
#include "registerio.h"
}

using namespace std;

Transform44::Transform44() {
}

Transform44::~Transform44() {

}

void
Transform44::SetMainTransform ( float i0j0,float i1j0,float i2j0,float i3j0,
				float i0j1,float i1j1,float i2j1,float i3j1,
				float i0j2,float i1j2,float i2j2,float i3j2,
				float i0j3,float i1j3,float i2j3,float i3j3){

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
  if( rFind != string::npos ) {

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
    if( NULL == registrationMatrix ) {
      throw runtime_error( "Couldn't load registration." );
    }

    SetMainTransform( registrationMatrix );
    
    MatrixFree( &registrationMatrix );
    free( sSubject );

  } else {

    LTA* lta = LTAreadEx( ifnLTA.c_str() );
    if( NULL == lta ) {
      throw runtime_error( "Couldn't load LTA." );
    }
    switch( lta->type ) {
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

  mInv = m.Inverse();
}

Transform44&
Transform44::Inverse() {

  // Make a new transform and set its main transform to _our_ inverse
  // transform.
  MATRIX* inverseM = MatrixInverse( m.GetMatrix(), NULL );
  Transform44* inverse = new Transform44();
  inverse->SetMainTransform( inverseM );
  MatrixFree( &inverseM );
  return *inverse;
}

Transform44& operator*( Transform44& m1, Transform44& m2 ) {

  Matrix44 mult = m1.GetMainMatrix() * m2.GetMainMatrix();

  Transform44* result = new Transform44();
  result->SetMainTransform( mult );
  return *result;
};


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
