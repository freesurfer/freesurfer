/**
 * @file  ScubaDataCollectionFactory.cpp
 * @brief Factory for creating DataCollection subclasse based on type
 *
 * Creates Scuba specfic DataCollections based on type strings (MRI
 * and MRIS). Also handles Tcl commands to do the same.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.10 $
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


#include "ScubaDataCollectionFactory.h"
#include "VolumeCollection.h"
#include "SurfaceCollection.h"

using namespace std;

bool ScubaDataCollectionFactory::mbAddedTclCommands = false;

ScubaDataCollectionFactory&
ScubaDataCollectionFactory::GetFactory() {

  static ScubaDataCollectionFactory sFactory;

  if ( !mbAddedTclCommands ) {

    mbAddedTclCommands = true;
    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( sFactory, "GetDataCollectionIDList", 0, "",
                           "Return a list of all collectionIDs." );
    commandMgr.AddCommand( sFactory, "MakeDataCollection", 1, "collectionType",
                           "Make a new data collection of the given type "
                           "and return the collectionID." );
    commandMgr.AddCommand( sFactory, "DeleteDataCollection", 1, "id",
                           "Delete a data collection." );

  }

  return sFactory;
}

TclCommandListener::TclCommandResult
ScubaDataCollectionFactory::DoListenToTclCommand( char* isCommand,
    int, char** iasArgv ) {

  // GetDataCollectionIDList
  if ( 0 == strcmp( isCommand, "GetDataCollectionIDList" ) ) {
    list<int> idList;
    DataCollection::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::iterator tID;
    for ( tID = idList.begin(); tID != idList.end(); ++tID ) {
      int id = *tID;
      ssFormat << "i";
      ssResult << id << " ";
    }
    ssFormat << "l";

    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  // MakeDataCollection <collectionType>   -- returns collection ID
  if ( 0 == strcmp( isCommand, "MakeDataCollection" ) ) {

    string sType = iasArgv[1];

    try {
      DataCollection& collection = MakeDataCollection( sType );
      int id = collection.GetID();
      stringstream ssResult;
      ssResult << id;
      sReturnFormat = "i";
      sReturnValues = ssResult.str();
    } catch ( runtime_error& e ) {
      DebugOutput( << "Bad collection type name" );
      sResult = "bad collection type";
      return error;
    }
  }

  // DeleteDataCollection <id>
  if ( 0 == strcmp( isCommand, "DeleteDataCollection" ) ) {

    int colID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    try {
      DataCollection& collection = DataCollection::FindByID( colID );
      delete &collection;

    } catch ( runtime_error& e ) {
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

  if ( isType == "Volume" ) {

    col = new VolumeCollection();
    DebugOutputStatic( << "Made a VolumeCollection" );

  } else if ( isType == "Surface" ) {

    col = new SurfaceCollection();
    DebugOutputStatic( << "Made a SurfaceCollection" );

  } else {

    DebugOutputStatic( << "Unknown collection type " << isType );
    throw runtime_error( "Unknown collection type" );
  }

  return *col;
}
