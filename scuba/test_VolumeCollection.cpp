#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "VolumeCollection.h"
#include "DataManager.h"
extern "C" {
#include "mri.h"
}
#include "Scuba-impl.h"

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


using namespace std;

char* Progname = "test_VolumeCollection";



class VolumeCollectionTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
VolumeCollectionTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;
  
  try {
  
    string fnMRI = "/Users/kteich/work/subjects/bert/mri/T1";
    
    char* sSubjectsDir = getenv("SUBJECTS_DIR");
    
    if( NULL != sSubjectsDir ) {
      fnMRI = string(sSubjectsDir) + "/bert/mri/T1";
    }

    VolumeCollection vol;
    vol.SetFileName( fnMRI );
    MRI* mri = vol.GetMRI();
    

    Assert( (vol.GetTypeDescription() == "Volume"),
	     "GetTypeDescription didn't return Volume" );

    DataManager dataMgr = DataManager::GetManager();
    MRILoader mriLoader = dataMgr.GetMRILoader();
    Assert( 1 == mriLoader.CountLoaded(), 
	    "CountLoaded didn't return 1" );
    Assert( 1 == mriLoader.CountReferences(mri),
	    "CountReferences didn't return 1" );

    char* fnMRIC = strdup( fnMRI.c_str() );
    MRI* mriComp = MRIread( fnMRIC );
    
    Assert( (MRImatch( mriComp, mri )), "MRImatch failed for load" );
    
    MRIfree( &mriComp );


    // Save it in /tmp, load it, and match it again.
    string fnSave( "/tmp/test.mgz" );
    vol.Save( fnSave );
    
    VolumeCollection testVol;
    testVol.SetFileName( fnSave );
    MRI* testMri = testVol.GetMRI();
    Assert( (MRImatch( testMri, mri )), "MRImatch failed for load after save");



    // Make an ROI and make sure it's a volume ROI.
    try {
      int roiID = vol.NewROI();
      ScubaROIVolume* roi = 
	dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( roiID ));
    }
    catch(...) {
      throw( runtime_error("typecast failed for NewROI") );
    }


    // Try our conversions.
    Point3<float> world;
    Point3<float> data;
    Point3<int> index;
    world.Set( -50, 0, -80 );
    vol.RASToMRIIndex( world.xyz(), index.xyz() );
    if( index.x() != 178 || index.y() != 208 || index.z() != 128 ) {
      cerr << "RASToMRIIndex failed. world " 
	   << world << " index " << index << endl;
      throw( runtime_error( "failed" ) );
    } 

    // Set a transform that scales the volume up by 2x in the world.
    ScubaTransform dataTransform;
    dataTransform.SetMainTransform( 2, 0, 0, 0,
				    0, 2, 0, 0,
				    0, 0, 2, 0,
				    0, 0, 0, 1 );
    vol.SetDataToWorldTransform( dataTransform.GetID() );

    world.Set( -50, 0, -80 );
    vol.RASToDataRAS( world.xyz(), data.xyz() );
    if( !(FEQUAL(data.x(),-25)) || 
	!(FEQUAL(data.y(),0)) || 
	!(FEQUAL(data.z(),-40)) ) {
      cerr << "RASToDataRAS failed. world " 
	   << world << " data " << data << endl;
      throw( runtime_error( "failed" ) );
    } 
    
    vol.RASToMRIIndex( world.xyz(), index.xyz() );

    if( index.x() != 153 || index.y() != 168 || index.z() != 128 ) {
      cerr << "RASToMRIIndex with data transform failed. world " 
	   << world << " index " << index << endl;
      throw( runtime_error( "failed" ) );
    }


    world.Set( -50, 0, -80 );
    VolumeLocation& loc = 
      (VolumeLocation&) vol.MakeLocationFromRAS( world.xyz() );
    if( !vol.IsInBounds( loc ) ) {
      cerr << "IsInBounds failed. world " << world << endl;
      throw( runtime_error( "failed" ) );
    }
    world.Set( -257, 0, 0 );
    VolumeLocation& loc2 =
      (VolumeLocation&) vol.MakeLocationFromRAS( world.xyz() );
    if( vol.IsInBounds( loc2 ) ) {
      cerr << "IsInBounds failed. world " << world << endl;
      throw( runtime_error( "failed" ) );
    }


    dataTransform.SetMainTransform( 2, 0, 0, 0,
				    0, 2, 0, 0,
				    0, 0, 2, 0,
				    0, 0, 0, 1 );
    vol.SetDataToWorldTransform( dataTransform.GetID() );
    world.Set( 0, -258, -254 );
    VolumeLocation& loc3 =
      (VolumeLocation&) vol.MakeLocationFromRAS( world.xyz() );
    if( vol.IsInBounds( loc3 ) ) {
      vol.RASToMRIIndex( world.xyz(), index.xyz() );
      cerr << "IsRASInMRIBounds failed. world " << world 
	   << " index " << index << endl;
      throw( runtime_error( "failed" ) );
    }


    // Create a new one from template.
    VolumeCollection vol2;
    vol2.MakeUsingTemplate( vol.GetID() );
    Assert( (vol.mVoxelSize[0] == vol2.mVoxelSize[0] &&
	     vol.mVoxelSize[1] == vol2.mVoxelSize[1] &&
	     vol.mVoxelSize[2] == vol2.mVoxelSize[2]),
	    "NewUsingTemplate failed, vol2 didn't match vol's voxelsize" );



    dataTransform.SetMainTransform( 1, 0, 0, 0,
				    0, 1, 0, 0,
				    0, 0, 1, 0,
				    0, 0, 0, 1 );
    vol.SetDataToWorldTransform( dataTransform.GetID() );

    Point3<float> p[4];
    Volume3<bool> v( 5, 5, 5, false );
    p[0].Set( 1, 1, 1 );   p[1].Set( 4, 1, 1 );
    p[2].Set( 4, 4, 1 );   p[3].Set( 1, 4, 1 );
    list<Point3<float> > points;
    vol.FindRASPointsInSquare( p[0].xyz(), p[1].xyz(), p[2].xyz(), p[3].xyz(),
			       0, points );
    list<Point3<float> >::iterator tPoint;
    for( tPoint = points.begin(); tPoint != points.end(); ++tPoint ) {
      Point3<float> p = *tPoint;
      v.Set( (int)p[0], (int)p[1], (int)p[2], true );
    }
    for( int nZ = 0; nZ < 5; nZ++ ) {
      for( int nY = 0; nY < 5; nY++ ) {
	for( int nX = 0; nX < 5; nX++ ) {
	  bool filled = v.Get( nX, nY, nZ );
	  if( nX >= 1 && nX <= 4 &&
	      nY >= 1 && nY <= 4 &&
	      nZ == 1 ) {
	    stringstream ssMsg;
	    ssMsg << "Failed: " << Point3<int>(nX,nY,nZ) << "was not filled";
	    Assert( filled == true, ssMsg.str() );
	  } else {
	    stringstream ssMsg;
	    ssMsg << "Failed: " << Point3<int>(nX,nY,nZ) << "was filled";
	    Assert( filled == false, ssMsg.str() );
	  }
	}
      }
    }


    // Check the tcl commands.
    char sCommand[1024];
    int rTcl;

    int id = vol.GetID();
    string fnTest = "test-name";
    sprintf( sCommand, "SetVolumeCollectionFileName %d test-name", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    
    Assert( (vol.mfnMRI == fnTest), 
	    "Setting file name via tcl didn't work" );

  }
  catch( logic_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
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


    VolumeCollectionTester tester0;
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
