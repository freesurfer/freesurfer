#include "PathManager.h"

using namespace std;

PathManager& 
PathManager::GetManager() {

  static PathManager* sManager = NULL;
  if( NULL == sManager ) {
    sManager = new PathManager();
  }

  return *sManager;
}

PathManager::PathManager () {
}

void
PathManager::ManagePath ( Path<float>& iPath ) {

  iPath.AddListener( this );
  mPaths.push_back( &iPath );
}

void
PathManager::UnmanagePath ( Path<float>& iPath ) {

  list<Path<float>*>::iterator tPath;
  for( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    if( path->GetID() == iPath.GetID() ) {
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
PathManager::DoListenToTclCommand ( char* isCommand, int iArgc, 
				    char** iasArgv ) {

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

  if( isMessage == "pathChanged" ||
      isMessage == "pathVertexAdded" ) {
    if( mbSendUpdates ) {
      SendBroadcast( isMessage, iData );
    }
  }
}

