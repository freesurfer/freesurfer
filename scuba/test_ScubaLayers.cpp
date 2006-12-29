/**
 * @file  test_ScubaLayers.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "ToglManager.h"
#include "ScubaFrame.h"
#include "ScubaView.h"
#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"
#include "ScubaDataCollectionFactory.h"
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


// Togl tester ---------------------------------------------------------

#if BUILD_TCL_TEST
extern "C" {
  int Test_scubalayers_Init ( Tcl_Interp* iInterp ) {

    ToglManager& toglMgr = ToglManager::GetManager();

    try {
      toglMgr.InitializeTogl( iInterp );
      toglMgr.SetFrameFactory( new ScubaFrameFactory );
      ScubaFrame::SetViewFactory( new ScubaViewFactory );

      TclCommandManager& commandMgr = TclCommandManager::GetManager();
      commandMgr.SetOutputStreamToCerr();
      commandMgr.Start( iInterp );

      ScubaLayerFactory& layerFactory =
        ScubaLayerFactory::GetFactory();

      ScubaDataCollectionFactory& dataFactory =
        ScubaDataCollectionFactory::GetFactory();
    } catch ( ... ) {
      return TCL_ERROR;
    }

    return TCL_OK;
  }
}
#endif


// Non-togl tester --------------------------------------------------------

class ScubaLayerFactoryTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
ScubaLayerFactoryTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    // Make the factory.
    ScubaLayerFactory& layerFactory = ScubaLayerFactory::GetFactory();
    layerFactory.SetOutputStreamToCerr();

    // Try making each kind of layer.
    string sType = "2DMRI";
    Layer& layer = layerFactory.MakeLayer( sType );
    Assert( (layer.GetTypeDescription() == "2DMRI" ),
            "Layer factory didn't return correct layer type" );


    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "MakeLayer 2DMRI" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    int layerID = strtol( sTclResult, (char**)NULL, 10 );
    Assert( (ERANGE != errno), "Error converting return ID from MakeLayer" );
    try {
      Layer& layer2 = Layer::FindByID( layerID );
      Assert( (layer2.GetTypeDescription() == "2DMRI" ),
              "Layer factory didn't return correct layer type via Tcl" );
    } catch (...) {
      ssError <<  "Couldn't find layer with return ID from MakeLayer, id ="
      << layerID;
      throw runtime_error( ssError.str() );
    }

  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
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


    ScubaLayerFactoryTester tester0;
    tester0.Test( interp );


  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
