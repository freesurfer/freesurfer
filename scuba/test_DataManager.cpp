#include <stdlib.h>
#include <string>
#include <iostream>
#include <tcl.h>
#include "DataManager.h"

#define Assert(x,s)   if(!(x)) { throw (char const*) s; }

using namespace std;

char* Progname = "test_DataManager";


#if 0
void TestMRI ( char const* ifnMRI ) {

  DataManager& dataMgr = DataManager::GetManager();

  // Load a sample data. Check the count of loaded data. Make sure the
  // MRI we have is the one we wanted.
  cerr << "GetMRI( " << ifnMRI << " )" << endl;
  MRI* mri = dataMgr.GetMRI( ifnMRI );

  Assert( 1 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 1" );

  Assert( 0 == strcmp( mri->fname, ifnMRI ),
	  "Filename and MRI filename don't match" );
  
  // Release the MRI and check the count.
  cerr << "Releasing" << endl;
  dataMgr.ReleaseMRI( &mri );

  Assert( 0 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 0" );
  
  // Load the MRI with multiple references. Make sure we still only
  // loaded it once. Make sure all the MRIs are ones we want.
  cerr << "Loading multiple times" << endl;
  MRI* mri1 = dataMgr.GetMRI( ifnMRI );
  MRI* mri2 = dataMgr.GetMRI( ifnMRI );
  MRI* mri3 = dataMgr.GetMRI( ifnMRI );

  Assert( 1 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 1" );

  Assert( mri1 == mri2,
	  "MRIs don't match" );
  Assert( mri2 == mri3,
	  "MRIs don't match" );
  Assert( 0 == strcmp( mri1->fname, ifnMRI ),
	  "Filename and MRI filename don't match" );
  Assert( 0 == strcmp( mri2->fname, ifnMRI ),
	  "Filename and MRI filename don't match" );
  Assert( 0 == strcmp( mri3->fname, ifnMRI ),
	  "Filename and MRI filename don't match" );
  
  // Release some of the references and make sure the MRI is loaded.
  cerr << "Releasing all but one" << endl;
  dataMgr.ReleaseMRI( &mri1 );
  Assert( 1 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 1" );
  dataMgr.ReleaseMRI( &mri2 );
  Assert( 1 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 1" );
  
  // Final release. Check the count.
  cerr << "Releasing final instance" << endl;
  dataMgr.ReleaseMRI( &mri3 );
  Assert( 0 == dataMgr.CountLoadedMRIs(), 
	  "CountLoadedMRIs didn't return 0" );
}
#endif

void TestMRILoader ( string const& ifnMRI ) {

  DataManager& dataMgr = DataManager::GetManager();
  MRILoader& loader = dataMgr.GetMRILoader();
  
  // Load a sample data. Check the count of loaded data. Make sure the
  // MRI we have is the one we wanted.
  cerr << "GetMRI( " << ifnMRI << " )" << endl;
  MRI* mri = loader.GetData( ifnMRI );

  Assert( 1 == loader.CountLoaded(), 
	  "CountLoadedMRIs didn't return 1" );

  Assert( 0 == strcmp( mri->fname, ifnMRI.c_str() ),
	  "Filename and MRI filename don't match" );
  
  // Release the MRI and check the count.
  cerr << "Releasing" << endl;
  loader.ReleaseData( &mri );

  Assert( 0 == loader.CountLoaded(), 
	  "CountLoaded didn't return 0" );
  
  // Load the MRI with multiple references. Make sure we still only
  // loaded it once. Make sure all the MRIs are ones we want.
  cerr << "Loading multiple times" << endl;
  MRI* mri1 = loader.GetData( ifnMRI );
  MRI* mri2 = loader.GetData( ifnMRI );
  MRI* mri3 = loader.GetData( ifnMRI );

  Assert( 1 == loader.CountLoaded(), 
	  "CountLoaded didn't return 1" );

  Assert( mri1 == mri2,
	  "MRIs don't match" );
  Assert( mri2 == mri3,
	  "MRIs don't match" );
  Assert( 0 == strcmp( mri1->fname, ifnMRI.c_str() ),
	  "Filename and MRI filename don't match" );
  Assert( 0 == strcmp( mri2->fname, ifnMRI.c_str() ),
	  "Filename and MRI filename don't match" );
  Assert( 0 == strcmp( mri3->fname, ifnMRI.c_str() ),
	  "Filename and MRI filename don't match" );
  
  // Release some of the references and make sure the MRI is loaded.
  cerr << "Releasing all but one" << endl;
  loader.ReleaseData( &mri1 );
  Assert( 1 == loader.CountLoaded(), 
	  "CountLoaded didn't return 1" );
  loader.ReleaseData( &mri2 );
  Assert( 1 == loader.CountLoaded(), 
	  "CountLoaded didn't return 1" );
  
  // Final release. Check the count.
  cerr << "Releasing final instance" << endl;
  loader.ReleaseData( &mri3 );
  Assert( 0 == loader.CountLoaded(), 
	  "CountLoaded didn't return 0" );
}


int main ( int argc, char** argv ) {

  string fnMRI = "/home/kteich/test_data/anatomical/bert.mgh";

  if( argc > 1 ) {
    fnMRI = argv[1];
  }

  try { 
 
    DataManager& dataMgr = DataManager::GetManager();

#if 0    
    cerr << "Testing MRI" << endl;
    TestMRI( fnMRI.c_str() );
#endif

    cerr << "Testing MRILoader" << endl;
    TestMRILoader( fnMRI );
  }
  catch( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
