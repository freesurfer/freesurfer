#include <stdexcept>
#include "SurfaceCollection.h"
#include "DataManager.h"

using namespace std;


SurfaceCollection::SurfaceCollection () :
  DataCollection() {

  mMRIS = NULL;
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
