#include "ToglManager.h"
#include "ScubaDataCollectionFactory.h"
#include "VolumeCollection.h"

char* Progname = "test_ScubaDataCollections";

using namespace std;

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


// Non-togl tester --------------------------------------------------------

class ScubaDataCollectionFactoryTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
ScubaDataCollectionFactoryTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    // Make the factory.
    ScubaDataCollectionFactory& colFactory =
      ScubaDataCollectionFactory::GetFactory();
    colFactory.SetOutputStreamToCerr();

    // Try making each kind of collection.
    string sType = "Volume";
    DataCollection& col = colFactory.MakeDataCollection( sType );
    string sCheckType = col.GetTypeDescription();
    Assert( ( sCheckType == "Volume" ),
	    "Data collection factory didn't return correct layer type" );


    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "MakeDataCollection Volume" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    char* sTclResult = Tcl_GetStringResult( iInterp );
    int colID = strtol( sTclResult, (char**)NULL, 10 );
    Assert( (ERANGE != errno), "Error converting return ID from MakeDataCollection" );
    try {
      DataCollection& col2 = DataCollection::FindByID( colID );
      Assert( (col2.GetTypeDescription() == "Volume" ),
	    "Data collection factory didn't return correct collection type via Tcl" );
    }
    catch(...) {
      ssError <<  "Couldn't find collection with return ID from "
	      << "MakeDataCollection, id = " << colID;
      throw new runtime_error( ssError.str() );
    }
    
  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
};  



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


    ScubaDataCollectionFactoryTester tester0;
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
