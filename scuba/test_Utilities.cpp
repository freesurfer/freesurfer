#include <stdexcept>
#include <sstream>

#include "Utilities.h"
#include "Volume3.h"
#include "Point2.h"
#include "Point3.h"
#include "Scuba-impl.h"

extern "C" {
#include "macros.h"
}


char* Progname = "test_Utilities";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

class UtilitiesTester {
public:
  void Test ();
};


void 
UtilitiesTester::Test () {

  stringstream ssError;

  try {

    Volume3<bool> bVol( 10, 15, 20, false );
    try {
      for( int nZ = 0; nZ < 20; nZ++ ) {
	for( int nY = 0; nY < 15; nY++ ) {
	  for( int nX = 0; nX < 10; nX++ ) {

	    ssError << "get failed for " << nX << ", " << nY << ", " << nZ;
	    bool bVal = bVol.Get( nX, nY, nZ );
	    ssError << "initial value was not false  for " << nX 
		    << ", " << nY << ", " << nZ;
	    Assert( (false == bVal), ssError.str() );

	    ssError << "set failed for " << nX << ", " << nY << ", " << nZ;
	    bVol.Set( nX, nY, nZ, true );

	    bVal = bVol.Get( nX, nY, nZ );
	    ssError << "set value was not true for " << nX 
		    << ", " << nY << ", " << nZ;
	    Assert( (true == bVal), ssError.str() );
	  }
	}
      }
    }
    catch(...) {
      throw runtime_error( ssError.str() );
    }

    try {
      bVol.Get( 11, 0, 0 );
      throw runtime_error( "didn't throw for illegal x access" );
    }
    catch(...) {}
    try {
      bVol.Get( 0, 16, 0 );
      throw runtime_error( "didn't throw for illegal x access" );
    }
    catch(...) {}
    try {
      bVol.Get( 0, 0, 21 );
      throw runtime_error( "didn't throw for illegal x access" );
    }
    catch(...) {}


    
    Point2<int> pt2( 1, 2 );
    cerr << "point is " << pt2 << endl;

    Point3<int> pt3( 1, 2, 3 );
    cerr << "point is " << pt3 << endl;


    int begin[2] = {0, 0};
    int end[2] = {1, 1};
    list<Point2<int> > points;
    Utilities::FindPointsOnLine2d( begin, end, 1, points );
    
      
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

    UtilitiesTester tester0;
    tester0.Test();

 
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

