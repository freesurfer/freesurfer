#include "ScubaROIVolume.h"

char* Progname = "test_ScubaROIVolume";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


class ScubaROIVolumeTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
ScubaROIVolumeTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    ScubaROIVolume roi;

    int bounds[3];
    try {
      bounds[0] = -1; bounds[1] = 1; bounds[2] = 2;
      roi.SetROIBounds( bounds );
      throw runtime_error( "Illegal SetROIBounds didn't throw" );
    }
    catch(...) {}

    bounds[0] = 10; bounds[1] = 11; bounds[2] = 12;
    roi.SetROIBounds( bounds );
    Assert( (10 == roi.mBounds[0] && 11 == roi.mBounds[1] &&
	     12 == roi.mBounds[2]), "SetROIBounds didn't work" );
    Assert( (NULL != roi.mVoxels), "SetROIBounds didn't create mVoxels" );

    int voxel[3];
    for( int nZ = 0; nZ < bounds[2]; nZ++ ) {
      for( int nY = 0; nY < bounds[1]; nY++ ) {
	for( int nX = 0; nX < bounds[0]; nX++ ) {
	  voxel[0] = nX; voxel[1] = nY; voxel[2] = nZ;
	  if( nZ % 2 && nY % 2 && nX % 2 ) {
	    roi.SelectVoxel( voxel );
	  } else {
	    roi.UnselectVoxel( voxel );
	  }
	}
      }
    }

    for( int nZ = 0; nZ < bounds[2]; nZ++ ) {
      for( int nY = 0; nY < bounds[1]; nY++ ) {
	for( int nX = 0; nX < bounds[0]; nX++ ) {
	  voxel[0] = nX; voxel[1] = nY; voxel[2] = nZ;
	  if( nZ % 2 && nY % 2 && nX % 2 ) {
	    Assert( (roi.IsVoxelSelected( voxel )), "voxel not selected" );
	  } else {
	    Assert( (!roi.IsVoxelSelected( voxel )), "voxel selected" );
	  }
	}
      }
    }

    try {
      voxel[0] = -1;
      roi.SelectVoxel( voxel );
      throw runtime_error( "SelectVoxel with x=-1 didn't throw" );
    }
    catch(...) {}
    
    try {
      voxel[0] = bounds[0];
      roi.SelectVoxel( voxel );
      throw runtime_error( "SelectVoxel with x=bounds[0] didn't throw" );
    }
    catch(...) {}
    
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

    ScubaROIVolumeTester tester0;
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

