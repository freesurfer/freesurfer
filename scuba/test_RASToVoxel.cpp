#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "Transform44.h"
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

using namespace std;

char* Progname = "test_RASToVoxel";



void TestCoord ( Transform44& rasToVoxel, Point3<float>& ras ) {


    Point3<int> idxi;
    Point3<float> idxf;
    Point3<float> rasFromIdxICheck;
    Point3<float> rasFromIdxFCheck;
    Point3<int> idxiFromRASFromIdxICheck;
    Point3<float> idxfFromRASFromIdxICheck;

    rasToVoxel.MultiplyVector3( ras.xyz(), idxi.xyz() );
    rasToVoxel.MultiplyVector3( ras.xyz(), idxf.xyz() );
    
    if( idxi.x() != (int) idxf.x() ||
	idxi.y() != (int) idxf.y() ||
	idxi.z() != (int) idxf.z() ) {
      
      stringstream ssError;
      ssError << "RAS -> index float didn't match int: "
	      << ras << " -> " << idxf << ", " << idxi;
      throw runtime_error( ssError.str() );
    }
    
    rasToVoxel.InvMultiplyVector3( idxf.xyz(), rasFromIdxFCheck.xyz() );
    
    if( fabs (rasFromIdxFCheck.x() - ras.x()) > 0.00001 ||
	fabs (rasFromIdxFCheck.y() - ras.y()) > 0.00001 ||
	fabs (rasFromIdxFCheck.z() - ras.z()) > 0.00001 ) {
      
      stringstream ssError;
      ssError << "RAS -> index float -> RAS didn't match orig RAS: "
	      << ras << " -> " << idxf << " -> " << rasFromIdxFCheck;
      throw runtime_error( ssError.str() );
    }
    
    rasToVoxel.InvMultiplyVector3( idxi.xyz(), rasFromIdxICheck.xyz() );
    
    rasToVoxel.MultiplyVector3( rasFromIdxICheck.xyz(), 
				idxiFromRASFromIdxICheck.xyz() );
    rasToVoxel.MultiplyVector3( rasFromIdxICheck.xyz(), 
				idxfFromRASFromIdxICheck.xyz() );
    
    if( idxiFromRASFromIdxICheck.x() != 
	(int) idxfFromRASFromIdxICheck.x() ||
	idxiFromRASFromIdxICheck.y() != 
	(int) idxfFromRASFromIdxICheck.y() ||
	idxiFromRASFromIdxICheck.z() != 
	(int) idxfFromRASFromIdxICheck.z() ) {
      
      stringstream ssError;
      ssError << "RAS -> index int -> RAS -> index float "
	      << "didn't match int: "
	      << ras << " -> " << idxf << " -> " << rasFromIdxICheck
	      << " -> " << idxiFromRASFromIdxICheck << ", " 
	      << idxfFromRASFromIdxICheck;
      throw runtime_error( ssError.str() );
    }
}

int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    Transform44 rasToVoxel;
    rasToVoxel.SetMainTransform( -9.78514, 0.0825722, 1.09133, 189.96,
				 -3.04974e-08, -9.81809, 0.742855, 286.028,
				 1.11155, 0.749787, 9.90971, 113.646,
				 0,      0,      0,      1 );

    Point3<float> ras;

    for( float z = -10.0; z < 10.0; z += 0.0490005 ) {
      for( float y = -10.0; y < 10.0; y += 0.0490003 ) {
	for( float x = -10.0; x < 10.0; x += 0.0490001 ) {
	  ras.Set( x, y, z );
	  TestCoord( rasToVoxel, ras );
	}
      }
    }

    for( int nTrial = 0; nTrial < 1000; nTrial++ ) {

      do {
	ras.Set( (float)random()/(float)random(),
		 (float)random()/(float)random(),
		 (float)random()/(float)random() );
      } while ( ras.x() < -50.0 || ras.x() > 50.0 ||
		ras.y() < -50.0 || ras.y() > 50.0 ||
		ras.z() < -50.0 || ras.z() > 50.0 );

      TestCoord( rasToVoxel, ras );
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

  cerr << "Success" << endl;

  exit( 0 );
}
