#include <pwd.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>


#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "DataCollection.h"
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

#define AssertTclNotOK(x) \
    if( TCL_OK == (x) ) { \
      ssError << "Tcl_Eval returned TCL_OK when it shouldn't've: " << endl  \
	     << "Command: " << sCommand << endl \
	     << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \

using namespace std;

char* Progname = "test_DataCollection";

class TestROI : public ScubaROI {
public:
  TestROI () {}
  ~TestROI () {}
};

class TestCollection : public DataCollection {

public:
  TestCollection( string isLabel ) :
    DataCollection() { 
    SetLabel( isLabel );
  }
  
  virtual void GetInfoAtRAS( float const iX, float const iY, float const iZ,
			     std::map<std::string,std::string>& iLabelValues ) {
    return;
  }

  virtual ScubaROI* DoNewROI () {
    return new TestROI();
  }
};

class DataCollectionTester {
public:
  void Test ( Tcl_Interp* iInterp );
};

void 
DataCollectionTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    TestCollection col1( "col1" );
    Assert( (col1.GetLabel() == "col1"), "col1 label incorrect" );
    TestCollection col2( "col2" );
    Assert( (col2.GetLabel() == "col2"), "col2 label incorrect" );
    TestCollection col3( "col3" );
    Assert( (col3.GetLabel() == "col3"), "col3 label incorrect" );
    Assert( (col1.GetID() != col2.GetID() != col3.GetID()),
	    "not unique IDs" );


    TestCollection* col4 = new TestCollection( "col4" );
    int col4ID = col4->GetID();
    DataCollection& col4comp = DataCollection::FindByID( col4ID );
    Assert( (col4ID == col4comp.GetID()), 
	    "Didn't get correct collection" );
    Assert( (col4->GetLabel() == col4comp.GetLabel()), 
	    "Didn't get correct label" );

    delete col4;
    try {
      DataCollection& col4comp = DataCollection::FindByID( col4ID );
      throw( logic_error("Didn't throw on deleted collection") );
    }
    catch( ... ) {}
    
    int roiID = col1.NewROI ();
    try {
      ScubaROI* roi = col1.mROIMap[roiID];
      Assert( (roi != NULL), "New ROI was null");
    } 
    catch(...) {
      throw( runtime_error( "Couldn't find the ROI that NewROI() should've created") );
    }

    try {
      col1.SelectROI( roiID );
    }
    catch(...) {
      throw( runtime_error( "SelectROI didn't find the new ROI" ));
    }

    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "NewCollectionROI %d", col1.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    char* sTclResult = Tcl_GetStringResult( iInterp );
    roiID = strtol( sTclResult, (char**)NULL, 10);
    try {
      ScubaROI* roi = col1.mROIMap[roiID];
      Assert( (roi != NULL), "TCl NewCollectionROI id was null");
    } 
    catch(...) {
      throw( runtime_error( "Couldn't find the ROI that NewCollectionROI "
			    "should've created") );
    }

    sprintf( sCommand, "SelectCollectionROI %d %d", col1.GetID(), roiID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (roiID == col1.mSelectedROIID), "Tcl SelectCollectionROI failed" );

    sprintf( sCommand, "SelectCollectionROI %d -1", col1.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclNotOK( rTcl );

    // Make a bunch of ROIs. Make a few in a different collection to
    // mix things up.
    col1.NewROI();
    col1.NewROI();
    col2.NewROI();
    col1.NewROI();
    col2.NewROI();
    col1.NewROI();
    col1.NewROI();

    // Get the ID list and for each ID, make sure it exists in the col.
    sprintf( sCommand, "GetROIIDListForCollection %d", col1.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    stringstream ssROIID( sTclResult );
    ssROIID >> roiID;
    while( !ssROIID.eof() ) {
      
      try {
	ScubaROI* roi = col1.mROIMap[roiID];
	Assert( (roi != NULL), "TCL GetROIIDListForCollection returned an "
		"roiID that wasn't in collection.");
      } 
      catch(...) {
      throw( runtime_error( "Couldn't find the ROI from an ID that "
			    "Tcl GetROIIDListForCollection returned.") );
      }

      ssROIID >> roiID;
    }
    
    // Now do the other way; get the list and for every ID in the col
    // roiID list, make sure the ID shows up in the list.
    sprintf( sCommand, "GetROIIDListForCollection %d", col1.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sTclResult = Tcl_GetStringResult( iInterp );
    
    map<int,ScubaROI*>::iterator tIDROI;
    for( tIDROI = col1.mROIMap.begin();
	 tIDROI != col1.mROIMap.end(); ++tIDROI ) {
      int roiID = (*tIDROI).first;
      stringstream ssID;
      ssID << roiID;
      Assert( (NULL != strstr( sTclResult, ssID.str().c_str() )),
	      "Didn't find an ROI in the collection that should've been "
	      "in list returned by Tcl GetROIIDListForCollection" );
    }

    
  }
  catch( exception e ) {
    cerr << "failed: " << e.what() << endl;
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


    DataCollectionTester tester0;
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


