#include "ScubaROI.h"
#include "ScubaColorLUT.h"

using namespace std;

DeclareIDTracker(ScubaROI);

ScubaROIStaticTclListener ScubaROI::mStaticListener;


ScubaROI::ScubaROI () :
  Broadcaster( "ScubaROI" ) ,
  msLabel( "" ),
  mType( Free ),
  mLUTID( 0 ),
  mStructure( 0 ),
  mbSuspendROIChangedMessage( false )
{
  mColor[0] = mColor[1] = mColor[2] = 0;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetROILabel", 2, "roiID label",
			 "Set the label for a ROI." );
  commandMgr.AddCommand( *this, "GetROILabel", 1, "roiID",
			 "Returns the label for a ROI." );
  commandMgr.AddCommand( *this, "SetROIType", 2, "roiID type",
			 "Sets ROI type. type should be structure or free." );
  commandMgr.AddCommand( *this, "GetROIType", 1, "roiID",
			 "Returns ROI type of structure or free." );
  commandMgr.AddCommand( *this, "SetROILUTID", 2, "roiID lutID",
			 "Sets the LUT ID for an ROI." );
  commandMgr.AddCommand( *this, "GetROILUTID", 1, "roiID",
			 "Returns LUT ID for an ROI." );
  commandMgr.AddCommand( *this, "SetROIStructure", 2, "roiID structure",
			 "Sets ROI structure index." );
  commandMgr.AddCommand( *this, "GetROIStructure", 1, "roiID",
			 "Returns ROI structure index." );
  commandMgr.AddCommand( *this, "SetROIColor", 4, "roiID red green blue",
			 "Sets ROI color. red, green, and blue should be"
			 "0-255 integers." );
  commandMgr.AddCommand( *this, "GetROIColor", 1, "roiID",
			 "Returns ROI color as a list of red, green, and "
			 "blue integers from 0-255." );

}

ScubaROI::~ScubaROI () {

}


TclCommandManager::TclCommandResult
ScubaROI::DoListenToTclCommand ( char* isCommand,
				 int, char** iasArgv ) {

  // SetROILabel <transformID> <label>
  if( 0 == strcmp( isCommand, "SetROILabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetROILabel <transformID>
  if( 0 == strcmp( isCommand, "GetROILabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }
  
  // SetROIType <roiID> <type>
  if( 0 == strcmp( isCommand, "SetROIType" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      Type type;
      if( 0 == strcmp( iasArgv[2], "structure" ) ) {
	type = Structure;
      } else if( 0 == strcmp( iasArgv[2], "free" ) ) {
	type = Free;
      } else {
	sResult = "bad type \"" + string(iasArgv[2]) +
	  "\", should be structure or free";
	return error;
      }
      
      SetType( type );
    }
  }

  // GetROIType <roiID>
  if( 0 == strcmp( isCommand, "GetROIType" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      switch( mType ) {
      case Structure:
	sReturnValues = "structure";
	break;
      case Free:
	sReturnValues = "free";
	break;
      }
      sReturnFormat = "s";
    }
  }

  // SetROILUTID <roiID> <lutID>
  if( 0 == strcmp( isCommand, "SetROILUTID" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {

      int lutID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad lut ID";
	return error;
      }
    
      try {
	ScubaColorLUT& lut = ScubaColorLUT::FindByID( lutID );
	if( lut.GetID() != lutID ) {
	  sResult = "Got wrong lut";
	  return error;
	}
      }
      catch(...) {
	sResult = "bad lut ID, doesn't exist";
	return error;
      }
      
      SetColorLUT( lutID );
    }
  }

  // GetROILUTID <roiID>
  if( 0 == strcmp( isCommand, "GetROILUTID" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << GetColorLUT();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // SetROIStructure <roiID> <structure>
  if( 0 == strcmp( isCommand, "SetROIStructure" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      int structure = strtol( iasArgv[2], (char**)NULL, 10 );
      if( ERANGE == errno ) {
	sResult = "bad structure";
	return error;
      }

      SetStructure( structure );
    }
  }

  // GetROIStructure <roiID>
  if( 0 == strcmp( isCommand, "GetROIStructure" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << mStructure;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // SetROIColor <roiID> <red> <green> <blue>
  if( 0 == strcmp( isCommand, "SetROIColor" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      int red = strtol( iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad red";
	return error;
      }
      
      int green = strtol( iasArgv[3], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad green";
	return error;
      }
      
      int blue = strtol( iasArgv[4], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad blue";
	return error;
      }
      
      int color[3];
      color[0] = red;
      color[1] = green;
      color[2] = blue;
      SetColor( color );
    }
  }

  // GetROIColor <roiID>
  if( 0 == strcmp( isCommand, "GetROIColor" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }
    
    if( mID == roiID ) {
      
      sReturnFormat = "Liiil";
      stringstream ssReturnValues;
      ssReturnValues << mColor[0] << " " << mColor[1] << " "
		     << mColor[2];
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}

void 
ScubaROI::SetColor( int iColor[3] ) {

  mColor[0] = iColor[0];
  mColor[1] = iColor[1];
  mColor[2] = iColor[2];
}

void 
ScubaROI::GetColor( int oColor[3] ) {

  oColor[0] = mColor[0];
  oColor[1] = mColor[1];
  oColor[2] = mColor[2];
}

void 
ScubaROI::GetDrawColor( int oColor[3] ) {

  switch( mType ) {
  case Free:
    oColor[0] = mColor[0];
    oColor[1] = mColor[1];
    oColor[2] = mColor[2];
    break;
  case Structure: {
    try {
      ScubaColorLUT& lut = ScubaColorLUT::FindByID( mLUTID );
      lut.GetColorAtIndex( mStructure, oColor );
    }
    catch(...) {
      // Couldn't find LUT, just use red.
      oColor[0] = 255;
      oColor[1] = 0;
      oColor[2] = 0;
    }
  }
    break;
  }
}

void
ScubaROI::ROIChanged() {

  if( !mbSuspendROIChangedMessage ) {
    int id = GetID();
    SendBroadcast( "roiChanged", (void*)&id );
  }
}

void 
ScubaROI::BeginBatchChanges () {
  mbSuspendROIChangedMessage = true;
}

void 
ScubaROI::EndBatchChanges () {
  mbSuspendROIChangedMessage = false;
  ROIChanged();
}

ScubaROIStaticTclListener::ScubaROIStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetROIIDList", 0, "", 
			 "Return a list of all ROIs." );
}

ScubaROIStaticTclListener::~ScubaROIStaticTclListener () {

}

TclCommandManager::TclCommandResult
ScubaROIStaticTclListener::DoListenToTclCommand ( char* isCommand, 
						  int, char** ) {

  // GetROIIDList
  if( 0 == strcmp( isCommand, "GetROIIDList" ) ) {
    list<int> idList;
    ScubaROI::GetIDList( idList );
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
