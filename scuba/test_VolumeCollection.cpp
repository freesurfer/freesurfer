#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "VolumeCollection.h"
#include "DataManager.h"
extern "C" {
#include "mri.h"
}

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
    
    Assert( (MRImatch( mriComp, mri )), "MRImatch failed" );
    
    MRIfree( &mriComp );


    // Make an ROI and make sure it's a volume ROI.
    try {
      int roiID = vol.NewROI();
      ScubaROIVolume* roi = 
	dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( roiID ));
    }
    catch(...) {
      throw( runtime_error("typecast failed for NewROI") );
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
