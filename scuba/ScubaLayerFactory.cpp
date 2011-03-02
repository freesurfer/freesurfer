/**
 * @file  ScubaLayerFactory.cpp
 * @brief A factory for Scuba specific Layers
 *
 * Makes layers  for specific types.  This is a little  different that
 * our  other factories  because it  creates  layers based  on a  type
 * (2DMRI and 2DMRIS). Also handles Tcl commands for doing the same.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.9 $
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


#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"
#include "ScubaLayer2DMRIS.h"

using namespace std;

bool ScubaLayerFactory::mbAddedTclCommands = false;

ScubaLayerFactory&
ScubaLayerFactory::GetFactory() {

  static ScubaLayerFactory sFactory;

  if ( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeLayer", 1, "layerType",
                           "Makes a new layer of the given type and returns "
                           "the layerID.");
    commandMgr.AddCommand( sFactory, "DeleteLayer", 1, "layerID",
                           "Deletes a layer. Does not delete any data "
                           "collections associated with it." );

  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaLayerFactory::DoListenToTclCommand( char* isCommand,
    int, char** iasArgv ) {

  // MakeLayer <layerType>   -- returns layer ID
  if ( 0 == strcmp( isCommand, "MakeLayer" ) ) {
    string sType = iasArgv[1];

    try {
      Layer& layer = MakeLayer( sType );
      int id = layer.GetID();
      stringstream ssResult;
      ssResult << id;
      sReturnFormat = "i";
      sReturnValues = ssResult.str();
    } catch ( runtime_error& e ) {
      DebugOutput( << "Bad layer type name" );
      sResult = "bad layer type";
      return error;
    }
  }

  if ( 0 == strcmp( isCommand, "DeleteLayer" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    try {
      Layer& layer = Layer::FindByID( layerID );
      delete &layer;
    } catch ( runtime_error& e ) {
      DebugOutput( << "Couldn't find layer" << layerID );
      sResult = "Couldn't find layer.";
      return error;
    }
  }

  return ok;
}

Layer&
ScubaLayerFactory::MakeLayer ( string isType ) {

  Layer* layer = NULL;

  if ( isType == "2DMRI" ) {

    layer = new ScubaLayer2DMRI();
    DebugOutputStatic( << "Made a ScubaLayer2DMRI" );

  } else if ( isType == "2DMRIS" ) {

    layer = new ScubaLayer2DMRIS();
    DebugOutputStatic( << "Made a ScubaLayer2DMRIS" );

  } else {

    DebugOutputStatic( << "Unknown layer type " << isType );
    throw runtime_error( "Unknown layer type" );
  }

  return *layer;
}
