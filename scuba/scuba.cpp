#include "tk.h"
#include "tix.h"
#include "ScubaGlobalPreferences.h"
#include "PreferencesManager.h"
#include "ToglManager.h"
#include "ScubaFrame.h"
#include "ScubaView.h"
#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"
#include "ScubaDataCollectionFactory.h"
#include "Scuba-impl.h"

char* Progname = "scuba";

// Togl tester ---------------------------------------------------------

extern "C" {
int Scuba_Init ( Tcl_Interp* iInterp ) {
    
    try {
    ToglManager& toglMgr = ToglManager::GetManager();
    toglMgr.InitializeTogl( iInterp );
    toglMgr.SetFrameFactory( new ScubaFrameFactory );

    ScubaFrame::SetViewFactory( new ScubaViewFactory );

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( iInterp );

    ScubaLayerFactory& layerFactory = 
      ScubaLayerFactory::GetFactory();
    layerFactory.SetOutputStreamToCerr();

    ScubaDataCollectionFactory& dataFactory = 
      ScubaDataCollectionFactory::GetFactory();
    dataFactory.SetOutputStreamToCerr();

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    //    prefsMgr.SetOutputStreamToCerr();
    prefsMgr.UseFile( ".scuba" );

    ScubaGlobalPreferences preferences =
      ScubaGlobalPreferences::GetPreferences();
  }
  catch( ... ) {
    return TCL_ERROR;
  }

  return TCL_OK;
}
}





int main ( int argc, char** argv ) {

  try {
    Tcl_Interp* interp = Tcl_CreateInterp();
    if( NULL == interp ) {
      throw runtime_error( "Tcl_CreateInterp returned null" );
    }
    
    int rTcl = Tcl_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tcl_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
    rTcl = Tk_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tk_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
    rTcl = Tix_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tix_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );
    commandMgr.SetCommandLineParameters( argc, argv );

    ToglManager& toglMgr = ToglManager::GetManager();
    toglMgr.InitializeTogl( interp );
    toglMgr.SetFrameFactory( new ScubaFrameFactory );

    ScubaFrame::SetViewFactory( new ScubaViewFactory );

    ScubaLayerFactory& layerFactory = 
      ScubaLayerFactory::GetFactory();
    layerFactory.SetOutputStreamToCerr();

    ScubaDataCollectionFactory& dataFactory = 
      ScubaDataCollectionFactory::GetFactory();
    dataFactory.SetOutputStreamToCerr();

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    //    prefsMgr.SetOutputStreamToCerr();
    prefsMgr.UseFile( ".scuba" );

    ScubaGlobalPreferences preferences =
      ScubaGlobalPreferences::GetPreferences();


    rTcl = Tcl_EvalFile( interp, "scuba.tcl" );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Reading scuba.tcl returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
 
    while( 1 ) {
      Tcl_DoOneEvent( TCL_ALL_EVENTS | TCL_DONT_WAIT );
    }

  }
  catch( runtime_error e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch( exception e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }


  exit( 0 );

}
