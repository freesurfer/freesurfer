[#include <stdexcept>
#include "DataCollection.h"

using namespace std;

DataCollection::ID DataCollection::mNextID = 0;
DataCollection::CollectionIDMap DataCollection::mCollectionIDs;

DataCollection::DataCollection( string isLabel ) {

  mID = GetNextID();
  mCollectionIDs[mID] = this;

  msLabel = isLabel;
}

DataCollection::~DataCollection() {

  mCollectionIDs[mID] = NULL;
}

void
DataCollection::GetInfoAtRAS( float const iX, float const iY, float const iZ,
			      std::list<std::string> olLabels,
			      std::list<std::string> olValues ) const {

  return;
}

DataCollection& 
DataCollection::GetDataCollection( ID const iID ) {

  DataCollection* dataCln = mCollectionIDs[iID];
  if( NULL == dataCln ) {
    DebugOutputStatic( << "GetDataCollection ID=" << iID << " was null" );
    throw domain_error( "ID no longer valid" );
  }
  
  return *dataCln;
}
