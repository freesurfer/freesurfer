/**
 * @file  PathManager.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.10 $
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


#include <fstream>
#include "PathManager.h"

using namespace std;

PathManager&
PathManager::GetManager() {

  static PathManager* sManager = NULL;
  if ( NULL == sManager ) {
    sManager = new PathManager();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sManager, "WritePathFile", 1, "fileName",
                           "Write paths to a file." );
    commandMgr.AddCommand( *sManager, "ReadPathFile", 1, "fileName",
                           "Read paths from a file." );
  }

  return *sManager;
}

PathManager::PathManager () :
    Broadcaster( "PathManager" ),
    Listener( "PathManager" ) {
  mbSendUpdates = true;
}

void
PathManager::ManagePath ( Path<float>& iPath ) {

  iPath.AddListener( this );
  mPaths.push_back( &iPath );

  // Notify listeners of changee.
  int pathID = iPath.GetID();
  SendBroadcast( "pathChanged", (void*)&pathID );
}

void
PathManager::UnmanagePath ( Path<float>& iPath ) {

  list<Path<float>*>::iterator tPath;
  for ( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    if ( path->GetID() == iPath.GetID() ) {
      int pathID = path->GetID();
      mPaths.erase( tPath );

      // Notify listeners of changee.
      SendBroadcast( "pathChanged", (void*)&pathID );
      break;
    }
  }
}

list<Path<float>*>&
PathManager::GetPathList () {
  return mPaths;
}

TclCommandManager::TclCommandResult
PathManager::DoListenToTclCommand ( char* isCommand, int,
                                    char** iasArgv ) {

  // WritePathFile
  if ( 0 == strcmp( isCommand, "WritePathFile" ) ) {

    string fnPaths = iasArgv[1];
    WritePathFile( fnPaths );
  }

  // ReadPathFile
  if ( 0 == strcmp( isCommand, "ReadPathFile" ) ) {

    string fnPaths = iasArgv[1];
    ReadPathFile( fnPaths );
  }


  return ok;
}

void
PathManager::EnableUpdates () {
  mbSendUpdates = true;
}

void
PathManager::DisableUpdates () {
  mbSendUpdates = false;
}

void
PathManager::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "pathChanged" ||
       isMessage == "pathVertexAdded" ) {
    if ( mbSendUpdates ) {
      SendBroadcast( isMessage, iData );
    }
  }
}

void
PathManager::ReadPathFile ( string ifnPaths ) {

  ifstream fPaths( ifnPaths.c_str(), ios::in );
  if ( !fPaths || fPaths.bad() ) {
    throw runtime_error( "Couldn't open paths file." );
  }

  int cPaths;
  fPaths >> cPaths;

  for ( int nPath = 0; nPath < cPaths; nPath++ ) {
    Path<float>* path = new Path<float>;
    path->ReadFromStream( fPaths );
    ManagePath( *path );

    // Notify listeners of changee.
    int pathID = path->GetID();
    SendBroadcast( "pathChanged", (void*)&pathID );
  }

  fPaths.close();
}

void
PathManager::WritePathFile ( string ifnPaths ) {

  ofstream fPaths( ifnPaths.c_str(), ios::out );
  if ( !fPaths || fPaths.bad() ) {
    throw runtime_error( "Couldn't write paths file." );
  }

  fPaths << mPaths.size() << endl;

  list<Path<float>*>::iterator tPath;
  for ( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    path->WriteToStream( fPaths );
  }

  fPaths.close();
}
