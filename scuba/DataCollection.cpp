#include <stdexcept>
#include "DataCollection.h"

using namespace std;

DeclareIDTracker(DataCollection);


DataCollection::DataCollection() {
}

DataCollection::~DataCollection() {
}

void
DataCollection::GetInfoAtRAS( float const iX, float const iY, float const iZ,
			      std::map<std::string,std::string>& iLabelValues ) {

  return;
}

void
DataCollection::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

}

 
