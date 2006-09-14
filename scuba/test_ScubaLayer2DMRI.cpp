#include "ToglManager.h"
#include "ScubaFrame.h"
#include "ScubaView.h"
#include "ScubaLayer2DMRI.h"
#include "Scuba-impl.h"

char* Progname = "test_ScubaLayers";

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

class ScubaLayer2DMRITester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
ScubaLayer2DMRITester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    string fnMRI = "test_data/bertT1.mgz";
    VolumeCollection vol;
    vol.SetFileName( fnMRI );
    vol.GetMRI();
    vol.GetID();

    ScubaLayer2DMRI layer;
    layer.SetVolumeCollection( vol );
    Assert( (&vol == layer.mVolume), "Didn't set volume collection properly" );
 
    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;
    
    sprintf( sCommand, "Set2DMRILayerVolumeCollection 99 99" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

    int layerID = layer.GetID();
    int volID = vol.GetID();
    sprintf( sCommand, "Set2DMRILayerVolumeCollection %d %d", layerID, volID );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

    
  }
  catch( runtime_error& e ) {
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

    ScubaLayer2DMRITester tester;
    tester.Test( interp );
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


