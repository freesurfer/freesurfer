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

Path<float>*
PathManager::NewPath () {

  Path<float>* path = new Path<float>;
  path->AddListener( this );
  mPaths.push_back( path );
  return path;
}

void
PathManager::DeletePath ( Path<float>* iPath ) {

  list<Path<float>*>::iterator tPath;
  for( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    if( path == iPath ) {
      mPaths.erase( tPath );
      delete path;
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
