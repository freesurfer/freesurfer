#include <iostream>
#include "Layer.h"

using namespace std;

template IDTracker<Layer>;
int IDTracker<Layer>::mNextID = 0;
map<int,Layer*> IDTracker<Layer>::mIDMap;

Layer::Layer() {
  mOpacity = 1.0;
  msLabel = "";
}

Layer::~Layer() {
}

void 
Layer::DrawIntoBuffer( GLbyte* iBuffer, ViewState& iViewState ) {
  cerr << "Layer " << msLabel << " is drawing into buffer" << endl;
}

void 
Layer::GetInfoAtRAS ( float inX, float inY, float inZ,
		      std::map<std::string,std::string>& iLabelValues ) {

  string sLabel;
  if( msLabel != "" ) {
    sLabel = msLabel;
  } else {
    stringstream ssLabel;
    ssLabel << mID;
    sLabel = ssLabel.str();
  }

  iLabelValues[sLabel] = "Hello world";
}
 
