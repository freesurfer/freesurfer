#include <iostream>
#include "Layer.h"

using namespace std;

DeclareIDTracker(Layer);

bool Layer::mbRegisteredStaticListener = false;
LayerStaticTclListener Layer::mStaticListener;

Layer::Layer() {
  mOpacity = 1.0;
  msLabel = "";

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetLayerLabel", 2, "layerID label",
			 "Set the label for a layer." );
  commandMgr.AddCommand( *this, "GetLayerLabel", 1, "layerID",
			 "Return the label for this layer." );
  commandMgr.AddCommand( *this, "GetLayerType", 1, "layerID",
			 "Return the type for this layer." );
  commandMgr.AddCommand( *this, "GetLayerOpacity", 1, "layerID",
			 "Return the opacity for this layer." );
  commandMgr.AddCommand( *this, "SetLayerOpacity", 2, "layerID opacity",
			 "Set the opacity for this layer. opacity should be "
			 "a float from 0 to 1." );
  if( !mbRegisteredStaticListener ) {
    commandMgr.AddCommand( mStaticListener, "GetLayerIDList", 0, "", 
			   "Return a list of all layerIDs." );
    mbRegisteredStaticListener = true;
  }
}

Layer::~Layer() {
}

void 
Layer::DrawIntoBuffer( GLubyte* iBuffer, int iWidth, int iHeight,
		       ViewState& iViewState,
		       ScubaWindowToRASTranslator& iTranslator ) {
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
 
TclCommandListener::TclCommandResult
Layer::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetLayerLabel <layerID> <label>
  if( 0 == strcmp( isCommand, "SetLayerLabel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetLayerLabel <layerID>
  if( 0 == strcmp( isCommand, "GetLayerLabel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "s";
      sReturnValues = GetLabel();
    }
  }

  // GetLayerType <layerID>
  if( 0 == strcmp( isCommand, "GetLayerType" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "s";
      sReturnValues = GetTypeDescription();
    }
  }

  // GetLayerOpacity <layerID>
  if( 0 == strcmp( isCommand, "GetLayerOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetOpacity();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetLayerOpacity <layerID> <opacity>
  if( 0 == strcmp( isCommand, "SetLayerOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float opacity = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad opacity";
	return error;
      }
      
      SetOpacity( opacity );
    }
  }
  
  return ok;
}

TclCommandListener::TclCommandResult
LayerStaticTclListener::DoListenToTclCommand ( char* isCommand, 
					       int iArgc, char** iasArgv ) {

  // GetLayerIDList
  if( 0 == strcmp( isCommand, "GetLayerIDList" ) ) {
    list<int> idList;
    Layer::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::iterator tID;
    for( tID = idList.begin(); tID != idList.end(); ++tID ) {
      int id = *tID;
      ssFormat << "i";
      ssResult << id << " ";
    }
    ssFormat << "l";
    
    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  return ok;
}

LayerStaticTclListener::~LayerStaticTclListener () {
}
