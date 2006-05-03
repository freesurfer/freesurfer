#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include "DataCollection.h"

using namespace std;

DeclareIDTracker(DataCollection);


DataCollection::DataCollection() :
  Listener( "DataCollection" ),
  Broadcaster( "DataCollection" ),
  msLabel(""),
  mSelectedROIID(-1),
  mbSuspendDataChangedMessage(false),
  mDataToWorldTransform(NULL) {

  // Try setting our initial transform to the default transform with
  // id 0. If it's not there, create it.
  try { 
    mDataToWorldTransform = &(ScubaTransform::FindByID( 0 ));
    mDataToWorldTransform->AddListener( this );
  }
  catch(...) {

    ScubaTransform* transform = new ScubaTransform();
    transform->SetLabel( "Identity" );

    try {
      mDataToWorldTransform = &(ScubaTransform::FindByID( 0 ));
      mDataToWorldTransform->AddListener( this );
    }
    catch(...) {
      DebugOutput( << "Couldn't make default transform!" );
    }
  }

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetCollectionLabel", 2, "collectionID label",
			 "Set the label for a collection." );
  commandMgr.AddCommand( *this, "GetCollectionLabel", 1, "collectionID",
			 "Return the label for this collection." );
  commandMgr.AddCommand( *this, "GetCollectionType", 1, "colID",
			 "Return the type for this layer." );
  commandMgr.AddCommand( *this, "NewCollectionROI", 1, "colID",
			 "Makes a new ROI for this collection and returns "
			 "the ID." );
  commandMgr.AddCommand( *this, "SelectCollectionROI", 2, "colID roiID",
			 "Selects an ROI for this collection." );
  commandMgr.AddCommand( *this, "DeleteCollectionROI", 2, "colID roiID",
			 "Deletes an ROI for this collection." );
  commandMgr.AddCommand( *this, "GetROIIDListForCollection", 1, "colID",
			 "Returns a lit of roiIDs belonging to this "
			 "collection." );
  commandMgr.AddCommand( *this, "SetDataTransform", 2, "colID transformID",
			 "Set the data to world transform for a data "
			 "collection." );
  commandMgr.AddCommand( *this, "GetDataTransform", 1, "colID",
			 "Returns the transformID of a data collection's "
			 "data to world transform." );
}

DataCollection::~DataCollection() {

  // Delete our ROIs when we get deleted.
  map<int,ScubaROI*>::iterator tIDROI;
  for( tIDROI = mROIMap.begin();
       tIDROI != mROIMap.end(); ++tIDROI ) {
    ScubaROI* roi = (*tIDROI).second;
    if( NULL != roi ) {
      delete roi;
    }
  }

  // Stop listening to whoever is still around.
  if( mDataToWorldTransform )
    mDataToWorldTransform->RemoveListener( this );
}

DataLocation&
DataCollection::MakeLocationFromRAS ( float const iRAS[3] ) {
  DataLocation* loc = new DataLocation( iRAS );
  return *loc;
}

void
DataCollection::GetInfo( DataLocation&,
			 std::map<std::string,std::string>& ) {

  return;
}

void 
DataCollection::SetLabel( string isLabel ) {
  msLabel = isLabel; 
  DataChanged();
}

void
DataCollection::GetDataRASBounds ( float oRASBounds[6] ) {

  oRASBounds[0] = oRASBounds[1] = oRASBounds[2] = 
    oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
}

TclCommandListener::TclCommandResult
DataCollection::DoListenToTclCommand( char* isCommand, int, char** iasArgv ) {

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
  if( 0 == strcmp( isCommand, "GetCollectionLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetCollectionType <colID>
  if( 0 == strcmp( isCommand, "GetCollectionType" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      sReturnFormat = "s";
      sReturnValues = GetTypeDescription();
    }
  }

  // NewCollectionROI <colID>
  if( 0 == strcmp( isCommand, "NewCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {

      int roiID = NewROI();

      stringstream ssReturnValues;
      ssReturnValues << roiID;
      sReturnFormat = "i";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SelectCollectionROI <colID> <roiID
  if( 0 == strcmp( isCommand, "SelectCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {

      int roiID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad roi ID";
	return error;
      }
    
      try {
	SelectROI( roiID );
      }
      catch(...) {
	sResult = "That ROI doesn't belong to this collection";
	return error;
      }
    }
  }

  // DeleteCollectionROI <colID> <roiID
  if( 0 == strcmp( isCommand, "DeleteCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {

      int roiID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad roi ID";
	return error;
      }
    
      try {
	DeleteROI( roiID );
      }
      catch(...) {
	sResult = "That ROI doesn't belong to this collection";
	return error;
      }
    }
  }

  // GetROIIDListForCollection <colID>
  if( 0 == strcmp( isCommand, "GetROIIDListForCollection" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      stringstream ssFormat;
      stringstream ssResult;
      ssFormat << "L";

      map<int,ScubaROI*>::iterator tIDROI;
      for( tIDROI = mROIMap.begin();
	   tIDROI != mROIMap.end(); ++tIDROI ) {
	int roiID = (*tIDROI).first;
	ScubaROI* roi = (*tIDROI).second;
	if( NULL != roi ) {
	  ssFormat << "i";
	  ssResult << roiID << " ";
	}
      }
      ssFormat << "l";
      
      sReturnFormat = ssFormat.str();
      sReturnValues = ssResult.str();
    }
  }

  // SetDataTransform <colID> <transformID>
  if( 0 == strcmp( isCommand, "SetDataTransform" ) ) {
    int colID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == colID ) {
      
      int transformID = strtol( iasArgv[2], (char**)NULL, 10 );
      if( ERANGE == errno ) {
	sResult = "bad transformID";
	return error;
      }
      
      SetDataToWorldTransform( transformID );
    }
  }

  // GetDataTransform <viewID>
  if( 0 == strcmp( isCommand, "GetDataTransform" ) ) {
    int colID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == colID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetDataToWorldTransform();
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}

void
DataCollection::DoListenToMessage ( string, void* ) {
  
}

vector<int>
DataCollection::GetROIList () {
  
  std::vector<int> lROIs;
  map<int,ScubaROI*>::iterator tIDROI;
  for( tIDROI = mROIMap.begin();
       tIDROI != mROIMap.end(); ++tIDROI ) {
    int roiID = (*tIDROI).first;
    lROIs.push_back( roiID );
  }

  return lROIs;
}

bool
DataCollection::IsROIInThisCollection ( int iROIID ) {

  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );

  return ( tIDROI != mROIMap.end() );
}

int
DataCollection::NewROI () {

  // Let the subclass make the appropriate ROI.
  ScubaROI* roi = this->DoNewROI();

  // Add it to the roi map.
  mROIMap[roi->GetID()] = roi;

  // Return the ID.
  return roi->GetID();
}

void
DataCollection::SelectROI ( int iROIID ) {

  // Look for this ID in the roi map. If found, select this ROI by
  // saving its ID. Otherwise throw an error.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if( tIDROI != mROIMap.end() ) {
    mSelectedROIID = iROIID;
  } else {
    throw runtime_error( "ROI doesn't belong to this collection" );
  }
}

void
DataCollection::DeleteROI ( int iROIID ) {

  // Look for this ID in the roi map. If found, delete the
  // ROI. Otherwise throw an error.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if( tIDROI != mROIMap.end() ) {

    ScubaROI* roi = (*tIDROI).second;
    delete roi;

    mROIMap.erase( iROIID );

    if( mSelectedROIID == iROIID ) {
      mSelectedROIID = -1;
    }
  } else {
    throw runtime_error( "ROI doesn't belong to this collection" );
  }
}

ScubaROI*
DataCollection::DoNewROI () {

  return NULL;
}
 
void
DataCollection::SetDataToWorldTransform ( int iTransformID ) {

  try {
    mDataToWorldTransform->RemoveListener( this );
    mDataToWorldTransform = &(ScubaTransform::FindByID( iTransformID ));
    mDataToWorldTransform->AddListener( this );
  }
  catch(...) {
    DebugOutput( << "Couldn't find transform " << iTransformID );
  }
}

int
DataCollection::GetDataToWorldTransform () {

  return mDataToWorldTransform->GetID();
}

float
DataCollection::GetPreferredValueIncrement () {
  
  return 0;
}

void
DataCollection::DataChanged () {
  
  if( !mbSuspendDataChangedMessage ) {
    int id = GetID();
    SendBroadcast( "dataChanged", (void*)&id );
  }
}


void 
DataCollection::BeginBatchChanges () {
  mbSuspendDataChangedMessage = true;
}

void 
DataCollection::EndBatchChanges () {
  mbSuspendDataChangedMessage = false;
  DataChanged();
}

void
DataCollection::BeginBatchROIChanges () {
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( mSelectedROIID );
  if( tIDROI != mROIMap.end() ) {
    (*tIDROI).second->BeginBatchChanges();
  }
}

void
DataCollection::EndBatchROIChanges () {
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( mSelectedROIID );
  if( tIDROI != mROIMap.end() ) {
    (*tIDROI).second->EndBatchChanges();
  }
}
