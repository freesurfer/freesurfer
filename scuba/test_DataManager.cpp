#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tcl.h>
#include "DataManager.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream sError; \
  sError << "Line " << __LINE__ << ": " << s; \
  cerr << sError.str().c_str() << endl; \
  throw logic_error( sError.str() ); \
  }

using namespace std;

char* Progname = "test_DataManager";

template <typename LoaderType, typename DataType> 
void TestLoaderStackLoad ( string const& ifnData,
			   DataLoader<DataType>& loader,
			   DataType& iData ) {

  iData = loader.GetData( ifnData );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 1 == loader.CountReferences(iData), "CountReferences didn't return 1" );

}

template <typename LoaderType, typename DataType> 
void TestLoaderStackRelease ( string const& ifnData,
			      DataLoader<DataType>& loader,
			      DataType& iData ) {

  loader.ReleaseData( &iData );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );
  Assert( 0 == loader.CountReferences(iData), "CountReferences didn't return 0" );
}


template <typename LoaderType, typename DataType> 
void TestLoader ( string const& ifnData,
		  DataLoader<DataType>& loader ) {

  loader.SetOutputStreamToCerr();

  cerr << "GetData( " << ifnData << " )" << endl;
  DataType data = loader.GetData( ifnData );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 1 == loader.CountReferences(data), "CountReferences didn't return 1" );

  // Release the Data and check the count.
  cerr << "Releasing" << endl;
  loader.ReleaseData( &data );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );
  Assert( 0 == loader.CountReferences(data), "CountReferences didn't return 0" );
  
  // Load the data with multiple references. Make sure we still only
  // loaded it once. Make sure all the Datas are ones we want.
  cerr << "Loading multiple times" << endl;
  DataType data1 = loader.GetData( ifnData );
  Assert( 1 == loader.CountReferences(data1), "CountReferences didn't return 1" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  DataType data2 = loader.GetData( ifnData );
  Assert( 2 == loader.CountReferences(data2), "CountReferences didn't return 2" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  DataType data3 = loader.GetData( ifnData );
  Assert( 3 == loader.CountReferences(data3), "CountReferences didn't return 3" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  Assert( data1 == data2, "Datas don't match" );
  Assert( data2 == data3, "Datas don't match" );
  
  // Release some of the references and make sure the Data is loaded.
  cerr << "Releasing all but one" << endl;
  loader.ReleaseData( &data1 );
  Assert( 2 == loader.CountReferences(data2), "CountReferences didn't return 2" );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  loader.ReleaseData( &data2 );
  Assert( 1 == loader.CountReferences(data3), "CountReferences didn't return 1" );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  
  // Final release. Check the count.
  cerr << "Releasing final instance" << endl;
  loader.ReleaseData( &data3 );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );


  // Load the data in a function, then check the counts, release the
  // data in a function, and check the counts again.
  cerr << "Loading up one level" << endl;
  TestLoaderStackLoad<LoaderType,DataType>( ifnData, loader, data );
  cerr << "Checking down one level" << endl;
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 1 == loader.CountReferences(data), "CountReferences didn't return 1" );
  cerr << "Releasing up one level" << endl;
  TestLoaderStackRelease<LoaderType,DataType>( ifnData, loader, data );
  cerr << "Checking down one level" << endl;
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 0 == loader.CountReferences(data), "CountReferences didn't return 1" );
}


int main ( int argc, char** argv ) {

  string fnMRI = "/Users/kteich/work/subjects/bert/mri/T1";
  string fnMRIS = "/Users/kteich/work/subjects/bert/surf/lh.white";

  char* sSubjectsDir = getenv("SUBJECTS_DIR");

  if( NULL != sSubjectsDir ) {
    fnMRI = string(sSubjectsDir) + "/bert/mri/T1";
    fnMRIS = string(sSubjectsDir) + "/bert/surf/lh.white";
  }


  try { 
 
    DataManager& dataMgr = DataManager::GetManager();
    dataMgr.SetOutputStreamToCerr();

    cerr << "Testing MRILoader" << endl;
    TestLoader<MRILoader,MRI*>( fnMRI, dataMgr.GetMRILoader() );

    cerr << "Testing MRISLoader" << endl;
    TestLoader<MRISLoader,MRIS*>( fnMRIS, dataMgr.GetMRISLoader() );

  }
  catch( exception e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
