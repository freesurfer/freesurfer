#include <stdexcept>
#include "SurfaceCollection.h"
#include "DataManager.h"

using namespace std;


SurfaceCollection::SurfaceCollection ( std::string& fnMRIS ) :
  DataCollection( "" ) {
  
  DataManager dataMgr = DataManager::GetManager();
  MRISLoader mrisLoader = dataMgr.GetMRISLoader();

  mMRIS = NULL;
  try { 
    mMRIS = mrisLoader.GetData( fnMRIS );
  }
  catch( exception e ) {
    throw logic_error( "Couldn't load MRIS" );
  }

  SetLabel( fnMRIS );
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
