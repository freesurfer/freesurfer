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
  commandMgr.AddCommand( *this, "SetLayerLabel" );
  commandMgr.AddCommand( *this, "GetLayerLabel" );
  commandMgr.AddCommand( *this, "GetLayerType" );
  commandMgr.AddCommand( *this, "GetLayerOpacity" );
  commandMgr.AddCommand( *this, "SetLayerOpacity" );
  if( !mbRegisteredStaticListener ) {
    commandMgr.AddCommand( mStaticListener, "GetLayerIDList" );
    mbRegisteredStaticListener = true;
  }
}

Layer::~Layer() {
}

void 
Layer::DrawIntoBuffer( GLubyte* iBuffer, ViewState& iViewState,
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
 
void
Layer::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetLayerLabel <layerID> <label>
  if( 0 == strcmp( isCommand, "SetLayerLabel" ) ) {
    if( 3 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {

	string sLabel = iasArgv[2];
	SetLabel( sLabel );
      }
    } else {
      sResult = "wrong # args: should be \"SetLayerLabel layerID label\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // GetLayerLabel <layerID>
  if( 0 == strcmp( isCommand, "GetLayerLabel" ) ) {
    if( 2 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {
	sReturnFormat = "s";
	sReturnValues = GetLabel();
      }
    } else {
      sResult = "wrong # args: should be \"GetLayerLabel layerID\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // GetLayerType <layerID>
  if( 0 == strcmp( isCommand, "GetLayerType" ) ) {
    if( 2 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {
	sReturnFormat = "s";
	sReturnValues = GetTypeDescription();
      }
    } else {
      sResult = "wrong # args: should be \"GetLayerType layerID\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // GetLayerOpacity <layerID>
  if( 0 == strcmp( isCommand, "GetLayerOpacity" ) ) {
    if( 2 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {
	sReturnFormat = "f";
	stringstream ssReturnValues;
	ssReturnValues << GetOpacity();
	sReturnValues = ssReturnValues.str();
      }
    } else {
      sResult = "wrong # args: should be \"GetLayerOpacity layerID\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // SetLayerOpacity <layerID> <opacity>
  if( 0 == strcmp( isCommand, "SetLayerOpacity" ) ) {
    if( 3 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {

	float opacity = strtof( iasArgv[2], (char**)NULL );
	if( ERANGE == errno ) {
	  sResult = "bad opacity";
	  return;
	}
	  
	SetOpacity( opacity );
      }
    } else {
      sResult = "wrong # args: should be \"SetLayerOpacity layerID opacity\"";
      DebugOutput( << sResult );
      return;
    }
  }
}

void 
LayerStaticTclListener::DoListenToTclCommand ( char* isCommand, 
					       int iArgc, char** iasArgv ) {

  // SetViewInPlane <frameID> <inPlane>
  if( 0 == strcmp( isCommand, "GetLayerIDList" ) ) {
    if( 1 == iArgc ) {

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
      
    } else {
      sResult = "wrong # args: should be \"GetLayerIDList\"";
      DebugOutputStatic( << sResult );
      return;
    }
  }
}

LayerStaticTclListener::~LayerStaticTclListener () {
}
