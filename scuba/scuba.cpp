extern "C" {
#include "tk.h"
#include "tix.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
}
#include "ScubaGlobalPreferences.h"
#include "PreferencesManager.h"
#include "ToglManager.h"
#include "PathManager.h"
#include "TclProgressDisplayManager.h"
#include "ScubaFrame.h"
#include "ScubaView.h"
#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"
#include "ScubaDataCollectionFactory.h"
#include "Scuba-impl.h"
#include "SegmentationVolumeReport.h"
#include "TclScubaKeyCombo.h"
#include "ScubaVolumeROIIntensityChart.h"
#include "TclChartWindow.h"
#include "ScubaMultiFrameVolumeChart.h"

char* Progname = "scuba";

// Togl tester ---------------------------------------------------------

// extern "C" {
// int Scuba_Init ( Tcl_Interp* iInterp ) {
    
//     try {
//     ToglManager& toglMgr = ToglManager::GetManager();
//     toglMgr.InitializeTogl( iInterp );
//     toglMgr.SetFrameFactory( new ScubaFrameFactory );

//     ScubaFrame::SetViewFactory( new ScubaViewFactory );

//     TclCommandManager& commandMgr = TclCommandManager::GetManager();
//     commandMgr.Start( iInterp );

//     ScubaLayerFactory::GetFactory();
    
//     ScubaDataCollectionFactory::GetFactory();

//     PathManager::GetManager();

//     PreferencesManager& prefsMgr = PreferencesManager::GetManager();
//     prefsMgr.UseFile( ".scuba" );

//     ScubaGlobalPreferences preferences =
//       ScubaGlobalPreferences::GetPreferences();

//   }
//   catch( ... ) {
//     return TCL_ERROR;
//   }

//   return TCL_OK;
// }
//}





int main ( int argc, char** argv ) {

  try {
    Tcl_Interp* interp = Tcl_CreateInterp();
    if( NULL == interp ) {
      throw runtime_error( "Tcl_CreateInterp returned null" );
    }
    
    int rTcl = Tcl_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      const char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tcl_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
    rTcl = Tk_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      const char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tk_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }
    
    rTcl = Tix_Init( interp );
    if( TCL_OK != rTcl ) {
      stringstream ssError;
      const char* sResult = Tcl_GetStringResult( interp );
      ssError <<  "Tix_Init returned not TCL_OK: " << sResult;
      throw runtime_error( ssError.str() );
    }

    
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.Start( interp );
    commandMgr.SetCommandLineParameters( argc, argv );

    ToglManager& toglMgr = ToglManager::GetManager();
    toglMgr.InitializeTogl( interp );
    toglMgr.SetFrameFactory( new ScubaFrameFactory );

    ScubaKeyCombo::SetFactory( new TclScubaKeyComboFactory() );

    ScubaFrame::SetViewFactory( new ScubaViewFactory );

    ScubaLayerFactory::GetFactory();
    
    ScubaDataCollectionFactory::GetFactory();

    PathManager::GetManager();

    ProgressDisplayManager::SetManager( new TclProgressDisplayManager );

    ScubaVolumeROIIntensityChartFactory::GetFactory();

    ScubaMultiFrameVolumeChartFactory::GetFactory();

    ChartWindow::SetFactory( new TclChartWindowFactory );

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    prefsMgr.UseFile( ".scuba" );

    ScubaGlobalPreferences preferences =
      ScubaGlobalPreferences::GetPreferences();

    SegmentationVolumeReport::GetReport();

    TclScubaKeyComboStaticTclListener::GetListener();

    // Look for the script, first in the local dir, then in
    // ../scripts, then in $FREESURFER_HOME/lib/tcl.
   struct stat info;
   string fnScuba( "./scuba.tcl" );
   int rStat = stat( fnScuba.c_str(), &info );
   if( 0 != rStat ||
       !S_ISREG(info.st_mode) ) {
     fnScuba = "../scripts/scuba.tcl";
     rStat = stat( fnScuba.c_str(), &info );
     if( 0 != rStat ||
	 !S_ISREG(info.st_mode) ) {
       char* sFressurferHome = getenv( "FREESURFER_HOME" );
       if( NULL != sFressurferHome ) {
	 fnScuba = sFressurferHome + string("/lib/tcl/scuba.tcl");
	 rStat = stat( fnScuba.c_str(), &info );
       }
     }
   }

   // If we haven't found one by now, bail.
   if( 0 != rStat ||
       !S_ISREG(info.st_mode) ) {
      stringstream ssError;
      ssError <<  "Couldn't find scuba.tcl file.";
      throw runtime_error( ssError.str() );
   }

   char* fnScubaC = strdup( fnScuba.c_str() );
   rTcl = Tcl_EvalFile( interp, fnScubaC );
   if( TCL_OK != rTcl ) {
     stringstream ssError;
     const char* sResult = Tcl_GetStringResult( interp );
     ssError <<  "Reading " << fnScuba << " returned not TCL_OK: " << sResult;
     throw runtime_error( ssError.str() );
   }
   free( fnScubaC );

   cout << "Using " << fnScuba << endl;

    
   while( 1 ) {
     Tcl_DoOneEvent( TCL_ALL_EVENTS );
   }
   
  }
  catch( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }


  exit( 0 );

}
