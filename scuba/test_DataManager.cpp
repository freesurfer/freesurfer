#include <stdlib.h>
#include <string>
#include <iostream>
#include <tcl.h>
#include "DataManager.h"

#define Assert(x,s)   if(!(x)) { throw (char const*) s; }

using namespace std;

char* Progname = "test_DataManager";

template <typename LoaderType, typename DataType> 
void TestLoader ( string const& ifnData,
		  DataLoader<DataType>& loader ) {

  cerr << "GetData( " << ifnData << " )" << endl;
  DataType data = loader.GetData( ifnData );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );

  // Release the Data and check the count.
  cerr << "Releasing" << endl;
  loader.ReleaseData( &data );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );
  
  // Load the data with multiple references. Make sure we still only
  // loaded it once. Make sure all the Datas are ones we want.
  cerr << "Loading multiple times" << endl;
  DataType data1 = loader.GetData( ifnData );
  DataType data2 = loader.GetData( ifnData );
  DataType data3 = loader.GetData( ifnData );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  Assert( data1 == data2, "Datas don't match" );
  Assert( data2 == data3, "Datas don't match" );
  
  // Release some of the references and make sure the Data is loaded.
  cerr << "Releasing all but one" << endl;
  loader.ReleaseData( &data1 );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  loader.ReleaseData( &data2 );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  
  // Final release. Check the count.
  cerr << "Releasing final instance" << endl;
  loader.ReleaseData( &data3 );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );
}


int main ( int argc, char** argv ) {

  string fnMRI = "/Users/kteich/work/subjects/bert/mri/T1";
  string fnMRIS = "/Users/kteich/work/subjects/bert/surf/lh.white";

  try { 
 
    DataManager& dataMgr = DataManager::GetManager();

    cerr << "Testing MRILoader" << endl;
    TestLoader<MRILoader,MRI*>( fnMRI, dataMgr.GetMRILoader() );

    cerr << "Testing MRISLoader" << endl;
    TestLoader<MRISLoader,MRIS*>( fnMRIS, dataMgr.GetMRISLoader() );

  }
  catch( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
