#include <fstream>
#include "VolumeCollection.h"
extern "C" {
#include "macros.h"
}
#include "Timer.h"
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



int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    char* sSubjectsDir = getenv("SUBJECTS_DIR");
    char* sTestDataDir = getenv("FSDEV_TEST_DATA");
    
    string fnVol;
    
    if( NULL != sTestDataDir ) {
      fnVol = string(sTestDataDir) + "anatomical/bert/mri/T1.mgz";
      ifstream fVol( fnVol.c_str(), ios::in );
      if( !fVol ) {
	if( NULL != sSubjectsDir ) {
	  fnVol = string(sSubjectsDir) + "/bert/mri/T1.mgz";
	  ifstream fVol( fnVol.c_str(), ios::in );
	  if( !fVol ) {
	    throw runtime_error( "Couldn't find necessary test data" );
	  }
	}
      }
      fVol.close();
    }

    VolumeCollection col;
    col.SetFileName( fnVol );
    col.GetMRI();

    float RAS[3];
    float index[3];

    Timer timer;
    float fps[10];

#if 0
    cerr << "Ready?" << endl;
    char blah;
    cin >> blah;
#endif

    for( int i = 0; i < 10; i++ ) {
      
      timer.Start();
      
      for( int d = 0; d < 512; d++ ) {
	for( int w = 0; w < 512; w++ ) {
	  col.RASToMRIIndex( RAS, index );
#if 0
	  index[0] = (int) rint( 1.0 * RAS[0] + 2.0 * RAS[1] +
				 3.0 + RAS[2] + 4.0 );
	  index[1] = (int) rint( 1.0 * RAS[0] + 2.0 * RAS[1] +
				 3.0 + RAS[2] + 4.0 );
	  index[2] = (int) rint( 1.0 * RAS[0] + 2.0 * RAS[1] +
				 3.0 + RAS[2] + 4.0 );
#endif
	}
      }
      
      int msec = timer.TimeNow();
      fps[i] = 1.0 / ((float)msec/1000.0);
    
      cerr << "fps: " << fps[i] << endl;
    }

    float totalfps = 0;
    for( int i = 0; i < 10; i++ ) {
      totalfps += fps[i];
    }

    cerr << "avg fps: " << totalfps/10.0 << endl;

  }
  catch( runtime_error& e ) {
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

