#include <stdexcept>
#include "DataCollection.h"

using namespace std;

DeclareIDTracker(DataCollection);


DataCollection::DataCollection() {
  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetCollectionLabel", 2, "collectionID label",
			 "Set the label for a collection." );
  commandMgr.AddCommand( *this, "GetCollectionLabel", 1, "collectionID",
			 "Return the label for this collection." );
}

DataCollection::~DataCollection() {
}

void
DataCollection::GetInfoAtRAS( float const iX, float const iY, float const iZ,
			      std::map<std::string,std::string>& iLabelValues ) {

  return;
}

TclCommandListener::TclCommandResult
DataCollection::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetCollectionLabel <collectionID> <label>
  if( 0 == strcmp( isCommand, "SetCollectionLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetCollectionLabel <collectionID>
  if( 0 == strcmp( isCommand, "SetCollectionLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      sReturnFormat = "s";
      sReturnValues = GetLabel();
    }
  }

  return ok;
}

 
