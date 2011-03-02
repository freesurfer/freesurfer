/**
 * @file  test_DataManager.cpp
 * @brief test DataManager class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.16 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <stdlib.h>
#include <fstream>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
extern "C" {
#define USE_NON_CONST
#include <tcl.h>
#undef USE_NON_CONST
}
#include "DataManager.h"
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream sError; \
  sError << "Line " << __LINE__ << ": " << s; \
  cerr << sError.str().c_str() << endl; \
  throw logic_error( sError.str() ); \
  }

using namespace std;

const char* Progname = "test_DataManager";

template <typename LoaderType, typename DataType>
void TestLoaderStackLoad ( string const& ifnData,
                           DataLoader<DataType>& loader,
                           DataType** ioData ) {

  *ioData = loader.GetData( ifnData );
}

template <typename LoaderType, typename DataType>
void TestLoaderStackRelease ( string const& ifnData,
                              DataLoader<DataType>& loader,
                              DataType** ioData ) {

  loader.ReleaseData( ioData );
}


template <typename LoaderType, typename DataType>
void TestLoader ( string const& ifnData,
                  DataLoader<DataType>& loader ) {

  // Tell the DebugReporter to generate output to cerr.
  loader.SetOutputStreamToCerr();

  // Load in some data and check the counts.
  DataType* data = loader.GetData( ifnData );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 1 == loader.CountReferences(data), "CountReferences didn't return 1" );

  // Release the Data and check the count.
  loader.ReleaseData( &data );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );


  // Load the data with multiple references. Make sure we still only
  // loaded it once. Make sure all the Datas are ones we want.
  DataType* data1 = loader.GetData( ifnData );
  Assert( 1 == loader.CountReferences(data1), "CountReferences didn't return 1" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  DataType* data2 = loader.GetData( ifnData );
  Assert( 2 == loader.CountReferences(data2), "CountReferences didn't return 2" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  DataType* data3 = loader.GetData( ifnData );
  Assert( 3 == loader.CountReferences(data3), "CountReferences didn't return 3" );
  Assert( 1 == loader.CountLoaded(),  "CountLoaded didn't return 1" );
  Assert( data1 == data2, "Datas don't match" );
  Assert( data2 == data3, "Datas don't match" );

  // Release some of the references and make sure the Data is loaded.
  loader.ReleaseData( &data1 );
  Assert( 2 == loader.CountReferences(data2), "CountReferences didn't return 2" );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  loader.ReleaseData( &data2 );
  Assert( 1 == loader.CountReferences(data3), "CountReferences didn't return 1" );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );

  // Final release. Check the count.
  loader.ReleaseData( &data3 );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );


  // Load the data in a function, then check the counts, release the
  // data in a function, and check the counts again.
  DataType* data4 = NULL;
  TestLoaderStackLoad<LoaderType,DataType>( ifnData, loader, &data4 );
  Assert( 1 == loader.CountLoaded(), "CountLoaded didn't return 1" );
  Assert( 1 == loader.CountReferences(data4), "CountReferences didn't return 1" );

  TestLoaderStackRelease<LoaderType,DataType>( ifnData, loader, &data4 );
  Assert( 0 == loader.CountLoaded(), "CountLoaded didn't return 0" );
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    DataManager& dataMgr = DataManager::GetManager();
    dataMgr.SetOutputStreamToCerr();

    string fnMRI = "test_data/bertT1.mgz";
    TestLoader<MRILoader,MRI>( fnMRI, dataMgr.GetMRILoader() );

    string fnMRIS = "test_data/lh.white";
    TestLoader<MRISLoader,MRIS>( fnMRIS, dataMgr.GetMRISLoader() );

  } catch ( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch ( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }
  catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
