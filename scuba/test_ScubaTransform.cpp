#include "ScubaTransform.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

char* Progname = "test_ScubaTransform";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
	     << "Command: " << sCommand << endl \
	     << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


#define VFEQUAL(v,a,b,c) \
   (FEQUAL(((v)[0]),a) && FEQUAL(((v)[1]),b) && FEQUAL(((v)[2]),c))

class ScubaTransformTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
ScubaTransformTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    ScubaTransform transform;
    
    // Make sure starts out as identity.
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	if( r == c ) { 
	  Assert((transform.GetCR(c,r) == 1),
		 "Didn't start out as identity");
	} else {
	  Assert((transform.GetCR(c,r) == 0),
		 "Didn't start out as identity");
	}
      }
    }

    // Set a bunch of values, make sure they got set.
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	transform.SetCR( c, r, (r * 4) + c);
      }
    }
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	float value = transform.GetCR( c, r );
	ssError << "Check failed for " << c << ", " << r;
	Assert((value - (r*4) == c &&
		(value - c) / 4 == r),
	       ssError.str() );
      }
    }

    // Set to identity, make sure multing a vector returns same vector.
    transform.MakeIdentity();
    float in[3], out[3];
    in[0] = 5;  in[1] = 35.67; in[2] = 1000;
    transform.MultiplyVector3( in, out );
    Assert((in[0] == out[0] && in[1] == out[1] && in[2] == out[2]),
	   "Identity mult check failed");


    // Set to a scale matrix, make sure multing a vector returns
    // correct response.
    transform.SetTransform( 5, 0, 0, 0,
			    0, 5, 0, 0,
			    0, 0, 5, 0,
			    0, 0, 0, 1 );

    in[0] = 5;  in[1] = 6; in[2] = 7;
    transform.MultiplyVector3( in, out );
    Assert((in[0]*5 == out[0] && in[1]*5 == out[1] && in[2]*5 == out[2]),
	   "Scale mult check failed");

    // Check the inverse.
    transform.InvMultiplyVector3( out, in );
    Assert((FEQUAL( in[0], 5.0 ) && 
	    FEQUAL( in[1], 6.0 ) && FEQUAL( in[2], 7.0) ),
	   "Inv scale mult check failed");


    Point3<float> p( 1, 0, 0 );
    Point3<float> q;

    transform.MakeXRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,1,0,0), "X Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,-1,0), "X Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,-1), "X Rotation failed");

    transform.MakeYRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,-1,0,0), "Y Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,1,0), "Y Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,-1), "Y Rotation failed");

    transform.MakeZRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,-1,0,0), "Z Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,-1,0), "Z Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,1), "Z Rotation failed");

    ScubaTransform m;
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	m.SetCR( c, r, (r * 4) + c);
      }
    }
    ScubaTransform id;
    id.MakeIdentity();
    ScubaTransform n = m * id;
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	Assert( (FEQUAL(n(c,r),m(c,r))), "Operator* failed" );
      }
    }

    p.Set( 0, 1, 0 );
    Point3<float> v( 1, 0, 0 );
    m.MakeRotation( p.xyz(), v.xyz(), M_PI );
    p.Set( 1, 0, 0 );
    m.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,1,2,0), "Rotation failed");

    p.Set( 0, 0, 1 );
    v.Set( 0, 1, 0 );
    m.MakeRotation( p.xyz(), v.xyz(), M_PI );
    p.Set( 0, 0, 0 );
    m.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,2), "Rotation failed");


    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "SetTransformValues %d {0 1 2 3 4 5 6 7 8 9 "
	     "10 11 12 13 14 15}", transform.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	float value = transform.GetCR( c, r );
	ssError << "TCL set check failed for " << c << ", " << r;
	Assert((value - (r*4) == c &&
		(value - c) / 4 == r),
	       ssError.str() );
      }
    }

    sprintf( sCommand, "GetTransformValues %d", transform.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    char* sTclResult = Tcl_GetStringResult( iInterp );
    stringstream ssResult( sTclResult );
    for( int r = 0; r < 4; r++ ) {
      for( int c = 0; c < 4; c++ ) {
	float value;
	ssResult >> value;
	ssError << "TCL set check failed for " << c << ", " << r;
	Assert((value - (r*4) == c &&
		(value - c) / 4 == r),
	       ssError.str() );
      }
    }
    


  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );
  
    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );
    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );


    ScubaTransformTester tester0;
    tester0.Test( interp );

 
  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
  
  cerr << "Success" << endl;

  exit( 0 );
}

