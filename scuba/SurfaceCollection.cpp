#include <stdexcept>
#include "SurfaceCollection.h"
#include "DataManager.h"

using namespace std;


SurfaceCollection::SurfaceCollection () :
  DataCollection() {

  mMRIS = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetSurfaceCollectionFileName", 2, 
			 "collectionID fileName", 
			 "Sets the file name for a given surface collection.");
}

SurfaceCollection::~SurfaceCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRISLoader mrisLoader = dataMgr.GetMRISLoader();
  try { 
    mrisLoader.ReleaseData( &mMRIS );
  } 
  catch(...) {
    cerr << "Couldn't release data"  << endl;
  }
}


void
SurfaceCollection::SetSurfaceFileName ( string& ifnMRIS ) {

  mfnMRIS = ifnMRIS;
}

MRIS*
SurfaceCollection::GetMRIS () { 

  if( NULL == mMRIS ) {
    
    DataManager dataMgr = DataManager::GetManager();
    MRISLoader mrisLoader = dataMgr.GetMRISLoader();
    
    mMRIS = NULL;
    try { 
      mMRIS = mrisLoader.GetData( mfnMRIS );
    }
    catch( exception e ) {
      throw logic_error( "Couldn't load MRIS" );
    }
  }

  if( msLabel == "" ) {
    SetLabel( mfnMRIS );
  }

  return mMRIS;
}

TclCommandListener::TclCommandResult 
SurfaceCollection::DoListenToTclCommand ( char* isCommand,
					 int iArgc, char** iasArgv ) {

  // SetSurfaceCollectionFileName <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "SetSurfaceCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      string fnSurface = iasArgv[2];
      SetSurfaceFileName( fnSurface );
    }
  }
  
  return DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

