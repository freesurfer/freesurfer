#include "ScubaLayerFactory.h"
#include "ScubaLayer2DMRI.h"

using namespace std;

bool ScubaLayerFactory::mbAddedTclCommands = false;

ScubaLayerFactory& 
ScubaLayerFactory::GetFactory() {

  static ScubaLayerFactory sFactory;

  if( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeLayer" );

  }

  return sFactory;
}

void
ScubaLayerFactory::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // MakeLayer <layerType>   -- returns layer ID
  if( 0 == strcmp( isCommand, "MakeLayer" ) ) {
    if( 2 == iArgc ) {

      string sType = iasArgv[1];

      try {
	Layer& layer = MakeLayer( sType );
	int id = layer.GetID();
	stringstream ssResult;
	ssResult << id;
	sReturnFormat = "i";
	sReturnValues = ssResult.str();
      }
      catch( runtime_error e ) {
	DebugOutput( << "Bad layer type name" );
	sResult = "bad layer type";
	return;
      }

    } else {
      sResult = "wrong # args: should be \"MakeLayer layerType\"";
      DebugOutput( << sResult );
      return;
    }
  }
}

Layer&
ScubaLayerFactory::MakeLayer ( string isType ) {

  Layer* layer = NULL;

  if( isType == "2DMRI" ) {

    layer = new ScubaLayer2DMRI();
    DebugOutputStatic( << "Made a ScubaLayer2DMRI" );

  } else {

    DebugOutputStatic( << "Unknown layer type " << isType );
    throw runtime_error( "Unknown layer type" );
  }

  return *layer;
}
