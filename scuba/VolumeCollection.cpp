#include <stdexcept>
#include "VolumeCollection.h"
#include "DataManager.h"

using namespace std;


VolumeCollection::VolumeCollection ( std::string& fnMRI ) :
  DataCollection( "" ) {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();

  mMRI = NULL;
  try { 
    mMRI = mriLoader.GetData( fnMRI );
  }
  catch( exception e ) {
    throw logic_error( "Couldn't load MRI" );
  }

  SetLabel( fnMRI );
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
