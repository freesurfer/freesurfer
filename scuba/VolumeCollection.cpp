#include <stdexcept>
#include "VolumeCollection.h"
#include "DataManager.h"

using namespace std;


VolumeCollection::VolumeCollection () :
  DataCollection() {
  mMRI = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollectionFileName" );
}

VolumeCollection::~VolumeCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();
  try { 
    mriLoader.ReleaseData( &mMRI );
  } 
  catch(...) {
    cerr << "Couldn't release data"  << endl;
  }
}

void
VolumeCollection::SetFileName ( string& ifnMRI ) {

  mfnMRI = ifnMRI;
}

MRI*
VolumeCollection::GetMRI() { 

  if( NULL == mMRI ) {
    
    DataManager dataMgr = DataManager::GetManager();
    MRILoader mriLoader = dataMgr.GetMRILoader();

    try { 
      mMRI = mriLoader.GetData( mfnMRI );
    }
    catch( exception e ) {
      throw logic_error( "Couldn't load MRI" );
    }

    if( msLabel == "" ) {
      SetLabel( mfnMRI );
    }
  }

  return mMRI; 
}

void 
VolumeCollection::DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv ) {

  // SetVolumeCollectionFileName <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "SetVolumeCollectionFileName" ) ) {
    if( 3 == iArgc ) {
      int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad collection ID";
	return;
      }

      if( mID == collectionID ) {

	string fnVolume = iasArgv[2];
	SetFileName( fnVolume );
      }
    } else {
      sResult = "wrong # args: should be \"SetVolumeCollectionFileName "
	"collectionID fileName\"";
      DebugOutput( << sResult );
      return;
    }
  }

  DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

