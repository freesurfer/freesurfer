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

#define Assert(x,e)   if(!(x)) { throw e; }

using namespace std;

char* Progname = "test_VolumeCollection";

int main ( int argc, char** argv ) {

  try {

    string fnMRI = "/Users/kteich/work/subjects/bert/mri/T1";
    
    char* sSubjectsDir = getenv("SUBJECTS_DIR");
    
    if( NULL != sSubjectsDir ) {
      fnMRI = string(sSubjectsDir) + "/bert/mri/T1";
    }

    VolumeCollection vol( fnMRI );
    MRI* mri = vol.GetMRI();
    
    DataManager dataMgr = DataManager::GetManager();
    MRILoader mriLoader = dataMgr.GetMRILoader();
    Assert( 1 == mriLoader.CountLoaded(), 
	    logic_error( "CountLoaded didn't return 1" ) );
    Assert( 1 == mriLoader.CountReferences(mri),
	    logic_error( "CountReferences didn't return 1" ) );

    char* fnMRIC = strdup( fnMRI.c_str() );
    MRI* mriComp = MRIread( fnMRIC );
    
    Assert( (MRImatch( mriComp, mri )), logic_error( "MRImatch failed" ) );
    
    MRIfree( &mriComp );

  }
  catch( logic_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
