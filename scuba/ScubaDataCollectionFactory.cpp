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
    commandMgr.AddCommand( sFactory, "MakeDataCollection" );

  }

  return sFactory;
}

void
ScubaDataCollectionFactory::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // MakeDataCollection <collectionType>   -- returns collection ID
  if( 0 == strcmp( isCommand, "MakeDataCollection" ) ) {
    if( 2 == iArgc ) {

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
	return;
      }

    } else {
      sResult = "wrong # args: should be \"MakeDataCollection "
	"collectionType\"";
      DebugOutput( << sResult );
      return;
    }
  }
}

DataCollection&
ScubaDataCollectionFactory::MakeDataCollection ( string isType ) {

  DataCollection* col = NULL;

  if( isType == "Volume" ) {

    col = new VolumeCollection();
    DebugOutputStatic( << "Made a VolumeCollection" );

  } else {

    DebugOutputStatic( << "Unknown collection type " << isType );
    throw runtime_error( "Unknown collection type" );
  }

  return *col;
}
