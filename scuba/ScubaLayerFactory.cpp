#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"
#include "ScubaLayer2DMRIS.h"

using namespace std;

bool ScubaLayerFactory::mbAddedTclCommands = false;

ScubaLayerFactory& 
ScubaLayerFactory::GetFactory() {

  static ScubaLayerFactory sFactory;

  if( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeLayer", 1, "layerType",
			   "Makes a new layer of the given type and returns "
			   "the layerID.");

  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaLayerFactory::DoListenToTclCommand( char* isCommand, 
					 int, char** iasArgv ) {
  
  // MakeLayer <layerType>   -- returns layer ID
  if( 0 == strcmp( isCommand, "MakeLayer" ) ) {
    string sType = iasArgv[1];
    
    try {
      Layer& layer = MakeLayer( sType );
      int id = layer.GetID();
      stringstream ssResult;
      ssResult << id;
      sReturnFormat = "i";
      sReturnValues = ssResult.str();
    }
    catch( runtime_error& e ) {
      DebugOutput( << "Bad layer type name" );
      sResult = "bad layer type";
      return error;
    }
  }

  return ok;
}

Layer&
ScubaLayerFactory::MakeLayer ( string isType ) {

  Layer* layer = NULL;

  if( isType == "2DMRI" ) {

    layer = new ScubaLayer2DMRI();
    DebugOutputStatic( << "Made a ScubaLayer2DMRI" );

  } else if( isType == "2DMRIS" ) {

    layer = new ScubaLayer2DMRIS();
    DebugOutputStatic( << "Made a ScubaLayer2DMRIS" );

  } else {

    DebugOutputStatic( << "Unknown layer type " << isType );
    throw runtime_error( "Unknown layer type" );
  }

  return *layer;
}
