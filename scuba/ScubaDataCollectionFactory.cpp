#include "ScubaDataCollectionFactory.h"
#include "VolumeCollection.h"
#include "SurfaceCollection.h"

using namespace std;

bool ScubaDataCollectionFactory::mbAddedTclCommands = false;

ScubaDataCollectionFactory& 
ScubaDataCollectionFactory::GetFactory() {

  static ScubaDataCollectionFactory sFactory;

  if( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "MakeDataCollection", 1, "collectionType",
			   "Make a new data collection of the given type "
			   "and return the collectionID." );

  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaDataCollectionFactory::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // MakeDataCollection <collectionType>   -- returns collection ID
  if( 0 == strcmp( isCommand, "MakeDataCollection" ) ) {

    string sType = iasArgv[1];
    
    try {
      DataCollection& collection = MakeDataCollection( sType );
      int id = collection.GetID();
      stringstream ssResult;
      ssResult << id;
      sReturnFormat = "i";
      sReturnValues = ssResult.str();
    }
    catch( runtime_error e ) {
      DebugOutput( << "Bad collection type name" );
      sResult = "bad collection type";
      return error;
    }
    
  }

  return ok;
}

DataCollection&
ScubaDataCollectionFactory::MakeDataCollection ( string isType ) {

  DataCollection* col = NULL;

  if( isType == "Volume" ) {

    col = new VolumeCollection();
    DebugOutputStatic( << "Made a VolumeCollection" );

  } else if( isType == "Surface" ) {

    col = new SurfaceCollection();
    DebugOutputStatic( << "Made a SurfaceCollection" );

  } else {

    DebugOutputStatic( << "Unknown collection type " << isType );
    throw runtime_error( "Unknown collection type" );
  }

  return *col;
}
