#include <stdlib.h>
#include <string>
#include <iostream>
#include <tcl.h>
#include "DataManager.h"

#define Assert(x,s)   if(!(x)) { throw (char const*) s; }

using namespace std;

char* Progname = "test_DataManager";



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


int main ( int argc, char** argv ) {

  char* fnMRI;

  try { 
 
    DataManager& dataMgr = DataManager::GetManager();
    
    TestMRI( "/home/kteich/test_data/anatomical/bert.mgh" );
  }
  catch( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
