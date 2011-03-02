/**
 * @file  DataCollection.cpp
 * @brief Base abstract data collection object
 *
 * This class is the base class for all data objects in Scuba. A 'collection'
 * is usually a primary data object, such as a volume or surface, and its
 * associated data, such as ROIs. This file also defines the DataLocation
 * class, which is encapsulates different coordinate spaces such as RAS and
 * data indicies.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.33 $
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


#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include "DataCollection.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

DataCollection::DataCollection() :
    Listener( "DataCollection" ),
    Broadcaster( "DataCollection" ),
    msLabel(""),
    mSelectedROIID(-1),
    mbSuspendDataChangedMessage(false),
    mDataToWorldTransform( NULL ) {

  // Try setting our initial transform to the default transform with
  // id 0. If it's not there, create it. This ensures we always have a
  // default identity transform.
  try {
    mDataToWorldTransform = &(ScubaTransform::FindByID( 0 ));
    mDataToWorldTransform->AddListener( *this );
  } catch (...) {

    // This will be created with ID 0.
    ScubaTransform* transform = new ScubaTransform();
    transform->SetLabel( "Identity" );

    try {
      mDataToWorldTransform = &(ScubaTransform::FindByID( 0 ));
      mDataToWorldTransform->AddListener( *this );
      
    } catch (...) {
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
  commandMgr.AddCommand( *this, "GetCollectionRASBounds", 1, "colID",
			 "Returns the RAS bounds of the data collection. Six "
			 "numbers are returned as a list of floats: xmin xmax "
			 "ymin ymax zmin zmax." );
}

DataCollection::~DataCollection() {

  // Delete our ROIs when we get deleted.
  map<int,ScubaROI*>::iterator tIDROI;
  for ( tIDROI = mROIMap.begin();
        tIDROI != mROIMap.end(); ++tIDROI ) {
    ScubaROI* roi = (*tIDROI).second;
    delete roi;
  }

  // Stop listening to whoever is still around.
  if( mDataToWorldTransform )
    mDataToWorldTransform->RemoveListener( *this );
}

DataLocation
DataCollection::MakeLocationFromRAS ( float const iRAS[3] ) const {
  return DataLocation( iRAS );
}

string const&
DataCollection::GetLabel() const {
  
  return msLabel;
}

void
DataCollection::SetLabel( string const& isLabel ) {
  msLabel = isLabel;
  DataChanged();
}

void
DataCollection::GetDataRASBounds ( float oRASBounds[6] ) const {
  
  oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
    oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
}

TclCommandListener::TclCommandResult
DataCollection::DoListenToTclCommand( char* isCommand, int, char** iasArgv ) {

  // SetCollectionLabel <collectionID> <label>
  if ( 0 == strcmp( isCommand, "SetCollectionLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {
      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetCollectionLabel <collectionID>
  if ( 0 == strcmp( isCommand, "GetCollectionLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetCollectionType <colID>
  if ( 0 == strcmp( isCommand, "GetCollectionType" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {
      sReturnFormat = "s";
      sReturnValues = GetTypeDescription();
    }
  }

  // NewCollectionROI <colID>
  if ( 0 == strcmp( isCommand, "NewCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      int roiID = NewROI();

      stringstream ssReturnValues;
      ssReturnValues << roiID;
      sReturnFormat = "i";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SelectCollectionROI <colID> <roiID
  if ( 0 == strcmp( isCommand, "SelectCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      int roiID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad roi ID";
        return error;
      }

      try {
        SelectROI( roiID );
      } catch (...) {
        sResult = "That ROI doesn't belong to this collection";
        return error;
      }
    }
  }

  // DeleteCollectionROI <colID> <roiID
  if ( 0 == strcmp( isCommand, "DeleteCollectionROI" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      int roiID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad roi ID";
        return error;
      }

      try {
        DeleteROI( roiID );
      } catch (...) {
        sResult = "That ROI doesn't belong to this collection";
        return error;
      }
    }
  }

  // GetROIIDListForCollection <colID>
  if ( 0 == strcmp( isCommand, "GetROIIDListForCollection" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssFormat;
      stringstream ssResult;
      ssFormat << "L";

      map<int,ScubaROI*>::iterator tIDROI;
      for ( tIDROI = mROIMap.begin();
            tIDROI != mROIMap.end(); ++tIDROI ) {
        int roiID = (*tIDROI).first;
        ScubaROI* roi = (*tIDROI).second;
        if ( NULL != roi ) {
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
  if ( 0 == strcmp( isCommand, "SetDataTransform" ) ) {
    int colID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == colID ) {

      int transformID = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad transformID";
        return error;
      }

      SetDataToWorldTransform( transformID );
    }
  }

  // GetDataTransform <viewID>
  if ( 0 == strcmp( isCommand, "GetDataTransform" ) ) {
    int colID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == colID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetDataToWorldTransform();
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetCollectionRASBounds <colID>
  if ( 0 == strcmp( isCommand, "GetCollectionRASBounds" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssResult;

      float bounds[6];
      this->GetDataRASBounds( bounds );

      ssResult << bounds[0] << " " << bounds[1] << " "
	       << bounds[2] << " " << bounds[3] << " "
	       << bounds[4] << " " << bounds[5];

      sReturnFormat = "Lffffffl";
      sReturnValues = ssResult.str();
    }
  }

  return ok;
}

void
DataCollection::DoListenToMessage ( string, void* ) {}

int
DataCollection::NewROI () {

  // Let the subclass make the appropriate ROI.
  ScubaROI* roi = this->DoNewROI();
  if( NULL == roi ) 
    throw runtime_error( "DoNewROI returned a NULL ROI" );

  // Add it to the roi map.
  mROIMap[roi->GetID()] = roi;

  // Return the ID.
  return roi->GetID();
}

void
DataCollection::DeleteROI ( int iROIID ) {

  // Look for this ID in the roi map. If found, delete the
  // ROI. Otherwise throw an error.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if ( tIDROI != mROIMap.end() ) {

    // Delete the ROI.
    ScubaROI* roi = (*tIDROI).second;
    delete roi;

    // Delete the entry from the map.
    mROIMap.erase( iROIID );

    // If this was the selected ROI, unselect it.
    if ( mSelectedROIID == iROIID )
      mSelectedROIID = -1;

  } else {
    throw runtime_error( "ROI doesn't belong to this collection" );
  }
}

vector<int>
DataCollection::GetROIList () const {

  // Put all our ROI IDs into a vector and return it.
  std::vector<int> lROIs;
  map<int,ScubaROI*>::const_iterator tIDROI;
  for ( tIDROI = mROIMap.begin();
        tIDROI != mROIMap.end(); ++tIDROI ) {
    int roiID = (*tIDROI).first;
    lROIs.push_back( roiID );
  }

  return lROIs;
}

int
DataCollection::GetNumberOfROIs () const {

  return mROIMap.size();
}

bool
DataCollection::IsROIInThisCollection ( int iROIID ) const {

  // Try to find the ROI ID in our map and return if it was found.
  map<int,ScubaROI*>::const_iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );

  return ( tIDROI != mROIMap.end() );
}

void
DataCollection::SelectROI ( int iROIID ) {

  // Look for this ID in the roi map. If found, select this ROI by
  // saving its ID. Otherwise throw an error.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if ( tIDROI != mROIMap.end() ) {
    mSelectedROIID = iROIID;
  } else {
    throw runtime_error( "ROI doesn't belong to this collection" );
  }
}

int
DataCollection::GetSelectedROI () const {
  return mSelectedROIID;
}

ScubaROI*
DataCollection::DoNewROI () {

  return NULL;
}

void
DataCollection::SetDataToWorldTransform ( int iTransformID ) {

  // Don't set if we're already using this one.
  if ( iTransformID == mDataToWorldTransform->GetID() )
    return;

  try {

    // Try to find this transform.
    ScubaTransform& transform = ScubaTransform::FindByID( iTransformID );

    // If we found it, stop listening to the current one and start
    // listening to this new one.
    mDataToWorldTransform->RemoveListener( *this );
    mDataToWorldTransform = &transform;
    mDataToWorldTransform->AddListener( *this );

  } catch (...) {
    DebugOutput( << "Couldn't find transform " << iTransformID );
  }
}

int
DataCollection::GetDataToWorldTransform () const {

  return mDataToWorldTransform->GetID();
}

float
DataCollection::GetPreferredValueIncrement () const {

  return 0;
}

void
DataCollection::DataChanged () {

  // When our data changes, broadcast a message with our ID.
  if ( !mbSuspendDataChangedMessage ) {
    int id = GetID();
    SendBroadcast( "dataChanged", (void*)&id );
  }
}


void
DataCollection::BeginBatchChanges () {

  // Don't call our changed function until we're done.
  mbSuspendDataChangedMessage = true;
}

void
DataCollection::EndBatchChanges () {

  // We're done, call our changed function.
  mbSuspendDataChangedMessage = false;
  DataChanged();
}

void
DataCollection::BeginBatchROIChanges () {

  // Same as the other batch changes, but we have to notify each ROI.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( mSelectedROIID );
  if ( tIDROI != mROIMap.end() ) {
    (*tIDROI).second->BeginBatchChanges();
  }
}

void
DataCollection::EndBatchROIChanges () {

  // Same as the other batch changes, but we have to notify each ROI.
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( mSelectedROIID );
  if ( tIDROI != mROIMap.end() ) {
    (*tIDROI).second->EndBatchChanges();
  }
}

DataLocation::DataLocation () {
  mRAS[0] = 0;
  mRAS[1] = 0;
  mRAS[2] = 0;
}

DataLocation::DataLocation ( float const iRAS[3] ) {
  mRAS[0] = iRAS[0];
  mRAS[1] = iRAS[1];
  mRAS[2] = iRAS[2];
}

DataLocation::DataLocation ( DataLocation const& iLoc ) {
  mRAS[0] = iLoc.RAS(0);
  mRAS[1] = iLoc.RAS(1);
  mRAS[2] = iLoc.RAS(2);
}
