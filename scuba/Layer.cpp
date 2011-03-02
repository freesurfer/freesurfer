/**
 * @file  Layer.cpp
 * @brief Basic Layer superclass
 *
 * This is a superclass for a visual layer that can be displayed by a
 * View. Contains the basic functionality for drawing, returning info
 * at an RAS point, responding to Tcl and Broadcaster messages, and
 * tool UI stuff.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.46 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <iostream>
#include "Layer.h"
#include "macros.h"
#include "Utilities.h"
#include "Point2.h"

using namespace std;

bool Layer::kbDefaultReportInfo = true;

LayerStaticTclListener Layer::mStaticListener;

Layer::Layer() :
    Listener( "Layer" ),
    Broadcaster( "Layer" ),
    mWidth(0),
    mHeight(0),
    msLabel(""),
    mOpacity(1.0),
    mbVisible(true),
    mbPostRedisplay(false),
    mbReportInfoAtRAS(kbDefaultReportInfo),
    mBytesPerPixel(4) {

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
  commandMgr.AddCommand( *this, "GetLayerPreferredBrushRadiusIncrement", 1,
                         "layerID", "Return a preferrerd brush radius "
                         "increment based on the data." );
  commandMgr.AddCommand( *this, "GetLayerPreferredValueIncrement", 1,
                         "layerID", "Return a preferrerd value "
                         "increment based on the data." );
  commandMgr.AddCommand( *this, "ProcessLayerOptionList", 2,
                         "layerID layerOptionList", "Process a string of "
                         "options in the format "
                         "option[=value][:option[=value]]..." );
  commandMgr.AddCommand( *this, "GetLayerReportInfo", 1, "layerID",
                         "Return whether or not a layer is reporting info." );
  commandMgr.AddCommand( *this, "SetLayerReportInfo", 2, "layerID report",
                         "Set whether or not a layer should report info." );
  commandMgr.AddCommand( *this, "GetLayerReportItemInfo", 2, 
			 "layerID itemLabel",
			 "Return whether or not a layer is reporting info "
			 "for the item with the given label" );
  commandMgr.AddCommand( *this, "SetLayerReportItemInfo", 3, 
			 "layerID itemLabel report",
			 "Set whether or not a layer should report info "
			 "for the item with the given label" );
  commandMgr.AddCommand( *this, "GetLayerReportableInfo", 1, "layerID",
			 "Returns a list of reportable info labels as a list "
			 "of strings." );
  commandMgr.AddCommand( *this, "GetLayerMainDataCollection", 1, "layerID",
                         "Returns the collection ID of the main collection "
                         "for this layer." );
}

Layer::~Layer() {
  int id = GetID();
  SendBroadcast( "layerDeleted", (void*)&id );
}

void
Layer::DrawIntoBuffer( GLubyte*, int, int,
                       ViewState&,
                       ScubaWindowToRASTranslator& ) {
  cerr << "Layer " << msLabel << " is drawing into buffer" << endl;
}

void
Layer::DrawIntoGL ( ViewState& iViewState,
                    ScubaWindowToRASTranslator& iTranslator ) {
  cerr << "Layer " << msLabel << " is drawing into GL" << endl;
}

void
Layer::GetInfoAtRAS ( float[3], list<InfoAtRAS>& ioInfo ) {

  if ( !mbReportInfoAtRAS )
    return;

  string sLabel;
  if ( msLabel != "" ) {
    sLabel = msLabel;
  } else {
    stringstream ssLabel;
    ssLabel << mID;
    sLabel = ssLabel.str();
  }

  InfoAtRAS info;
  info.SetLabel( sLabel );
  info.SetValue( "Hello world" );
  ioInfo.push_back( info );
}

void
Layer::SetReportInfo ( bool ibReport ) {

  mbReportInfoAtRAS = ibReport;
  SendBroadcast( "layerInfoSettingsChanged", NULL );
}

void
Layer::SetReportInfo ( string isInfoLabel, bool ibReport ) {

  if( maReportableInfo.find( isInfoLabel ) == maReportableInfo.end() )
    throw runtime_error( string("Info item with label ") + 
			 isInfoLabel + " not found." );

  maReportableInfo[isInfoLabel] = ibReport;
  SendBroadcast( "layerInfoSettingsChanged", NULL );
}

bool
Layer::GetReportInfo ( string isInfoLabel ) {

  if( maReportableInfo.find( isInfoLabel ) == maReportableInfo.end() )
    throw runtime_error( string("Info item with label ") + 
			 isInfoLabel + " not found." );

  return maReportableInfo[isInfoLabel];
}

void
Layer::AddReportableInfo ( string isInfoLabel, bool ibReport ) {

  if( maReportableInfo.find( isInfoLabel ) != maReportableInfo.end() )
    throw runtime_error( string("Info item with label ") + 
			 isInfoLabel + " already in map." );

  // Add it to the map.
  maReportableInfo[isInfoLabel] = ibReport;
}

void
Layer::GetReportableInfo ( vector<string>& iolItems ) {

  // Go through the map and push all items into the list.
  for( map<string,bool>::iterator tReportableInfo = maReportableInfo.begin();
       tReportableInfo != maReportableInfo.end();
       ++tReportableInfo ) {
    
    iolItems.push_back( tReportableInfo->first );
  }
}

TclCommandListener::TclCommandResult
Layer::DoListenToTclCommand( char* isCommand, int, char** iasArgv ) {

  // SetLayerLabel <layerID> <label>
  if ( 0 == strcmp( isCommand, "SetLayerLabel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetLayerLabel <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerLabel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetLayerType <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerType" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "s";
      sReturnValues = GetTypeDescription();
    }
  }

  // GetLayerOpacity <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetOpacity();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetLayerOpacity <layerID> <opacity>
  if ( 0 == strcmp( isCommand, "SetLayerOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float opacity = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad opacity";
        return error;
      }

      SetOpacity( opacity );
    }
  }

  // GetLayerPreferredBrushRadiusIncrement <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerPreferredBrushRadiusIncrement" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetPreferredBrushRadiusIncrement();
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetLayerPreferredValueIncrement <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerPreferredValueIncrement" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetPreferredValueIncrement();
      sReturnValues = ssReturnValues.str();
    }
  }

  // ProcessLayerOptionList <layerID> <layerOptions>
  if ( 0 == strcmp( isCommand, "ProcessLayerOptionList" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        ProcessOptionsList( string(iasArgv[2]) );
      } catch ( runtime_error& e ) {
        sResult = "bad options \"" + string(iasArgv[2]) + "\", \n" + e.what();
        return error;
      }
    }
  }

  // GetLayerReportInfo <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerReportInfo" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( GetReportInfo() );
      sReturnFormat = "i";
    }
  }

  // SetLayerReportInfo <layerID> <report>
  if ( 0 == strcmp( isCommand, "SetLayerReportInfo" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bReport =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetReportInfo( bReport );
      } catch ( runtime_error& e ) {
        sResult = "bad report \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // GetLayerReportItemInfo <layerID> <itemLabel>
  if ( 0 == strcmp( isCommand, "GetLayerReportItemInfo" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      bool bReport = GetReportInfo( iasArgv[2] );
      sReturnValues =
	TclCommandManager::ConvertBooleanToReturnValue( bReport );
      sReturnFormat = "i";
    }
  }

  // SetLayerReportItemInfo <layerID> <itemLabel> <report>
  if ( 0 == strcmp( isCommand, "SetLayerReportItemInfo" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bReport =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
        SetReportInfo( iasArgv[2], bReport );
      } catch ( runtime_error& e ) {
        sResult = "bad report \"" + string(iasArgv[3]) + "\"," + e.what();
        return error;
      }
    }
  }

  // GetLayerReportableInfo <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerReportableInfo" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssFormat;
      stringstream ssResult;
      ssFormat << "L";
     
      for( map<string,bool>::iterator tReportableInfo = 
	     maReportableInfo.begin();
	   tReportableInfo != maReportableInfo.end();
	   ++tReportableInfo ) {
	
	ssFormat << "s";
	ssResult << "\"" << tReportableInfo->first << "\" ";
      }

      ssFormat << "l"; 

      sReturnFormat = ssFormat.str();
      sReturnValues = ssResult.str();
    }
  }

  // GetLayerMainDataCollection <layerID>
  if ( 0 == strcmp( isCommand, "GetLayerMainDataCollection" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      DataCollection* collection = GetMainDataCollection();
      if ( collection ) {
        stringstream ssReturnValues;
        ssReturnValues << collection->GetID();
        sReturnValues = ssReturnValues.str();
        sReturnFormat = "i";
      } else {
        stringstream ssReturnValues;
        ssReturnValues << -1;
        sReturnValues = ssReturnValues.str();
        sReturnFormat = "i";
      }
    }
  }

  return ok;
}

void
Layer::ProcessOptionsList ( string isOptionList ) {

  vector<string> lsOptions;
  int cDelims =
    Utilities::SplitString( isOptionList, ":", lsOptions );

  // If cDelims is zero, there might still be an optionvalue here,
  // it's just that the string we got was "option=value" so there's no
  // delimiter. So just push the whole thing.
  if ( cDelims == 0 ) {
    lsOptions.push_back( isOptionList );
  }

  bool bError = false;
  stringstream ssErrors;
  vector<string>::iterator tOption;
  for ( tOption = lsOptions.begin();
        tOption != lsOptions.end(); ++tOption ) {
    string sOptionValue = *tOption;

    int nEquals = sOptionValue.find( '=' );
    if ( nEquals < 0 ) {

      try {
        this->ProcessOption( sOptionValue, "" );
      } catch ( runtime_error& e) {
        ssErrors << e.what() << ". ";
        bError = true;
      }

    } else {

      try {
        int z = sOptionValue.length();
        this->ProcessOption( sOptionValue.substr( 0, nEquals ),
                             sOptionValue.substr( nEquals+1, z-nEquals+1 ) );
      } catch ( runtime_error& e) {
        ssErrors << e.what() << ". ";
        bError = true;
      }

    }

  }

  if ( bError ) {
    throw runtime_error( ssErrors.str() );
  }
}

void
Layer::ProcessOption ( string isOption, string isValue ) {

  char sValue[1024];
  strcpy( sValue, isValue.c_str() );

  if ( 0 == isOption.compare( "label" ) ) {
    SetLabel( sValue );

  } else if ( 0 == isOption.compare( "opacity" ) ) {
    float opacity = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad opacity value" );
    }
    SetOpacity( opacity );

  } else if ( 0 == isOption.compare( "visible" ) ) {
    int visible = strtol( sValue, (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad visible value" );
    }
    SetVisible( (bool)visible );

  } else {

    stringstream ssError;
    ssError << "Unrecognized option " << isOption;
    throw runtime_error( ssError.str() );
  }
}

void
Layer::DoListenToMessage ( string isMessage, void* ) {

  if ( isMessage == "dataChanged" ) {
    DataChanged ();
    int id = GetID();
    SendBroadcast( "layerChanged", (void*)&id );
  }
}

DataCollection*
Layer::GetMainDataCollection() {
  return NULL;
}


int
Layer::GetSelectedROI () {

  DataCollection* col = this->GetMainDataCollection();
  if ( NULL != col ) {
    return col->GetSelectedROI();
  } else {
    return -1;
  }
}

void
Layer::DataChanged () {
  // Default behavior of data changing is to request a redisplay.
  cerr << "Layer::DataChanged()" << endl;
  RequestRedisplay();
}

void
Layer::SetLabel ( string isLabel ) {
  msLabel = isLabel;
}

string
Layer::GetLabel () {
  return msLabel;
}

void
Layer::SetWidth( int iWidth ) {
  mWidth = iWidth;
}

void
Layer::SetHeight( int iHeight ) {
  mHeight = iHeight;
}

void
Layer::HandleTool ( float[3], ViewState&,
                    ScubaWindowToRASTranslator&,
                    ScubaToolState&, InputState& ) {}

void
Layer::Timer () {
  this->DoTimer ();
}

void
Layer::DrawPixelIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                             int iWindow[2], int iColor[3], float iOpacity ) {

  if ( iWindow[0] >= 0 && iWindow[0] < iWidth &&
       iWindow[1] >= 0 && iWindow[1] < iHeight ) {
    GLubyte* dest = iBuffer + (iWindow[0] * mBytesPerPixel) +
                    (iWindow[1] * iWidth * mBytesPerPixel);
    DrawPixelIntoBuffer( dest, iColor, iOpacity );
  }
}

void
Layer::DrawPixelIntoBuffer ( GLubyte* iAddress,
                             int iColor[3], float iOpacity ) {

  iAddress[0] = (GLubyte) (((float)iAddress[0] * (1.0 - iOpacity)) +
                           ((float)iColor[0] * iOpacity));
  iAddress[1] = (GLubyte) (((float)iAddress[1] * (1.0 - iOpacity)) +
                           ((float)iColor[1] * iOpacity));
  iAddress[2] = (GLubyte) (((float)iAddress[2] * (1.0 - iOpacity)) +
                           ((float)iColor[2] * iOpacity));
  iAddress[3] = 255;
}



void
Layer::DrawAALineIntoBuffer ( GLubyte* iBuffer, int iWidth, int,
                              int iFromWindow[2], int iToWindow[2],
                              int iColor[3], int iThickness, float ) {

  // Set up x1, x2, y1, y2 so that x1 is always less than x2;
  int x1 = MIN( iFromWindow[0], iToWindow[0] );
  int x2 = MAX( iFromWindow[0], iToWindow[0] );
  int y1 = iFromWindow[1];
  int y2 = iToWindow[1];

  // Calculate deltas.
  int dx = x2 - x1;
  int dy = y2 - y1;

  // If the dy is negative, remember it and reverse dy.
  bool bNegativeY = false;
  if ( dy < 0 ) {
    bNegativeY = true;   // This affects
    dy = -dy;
  }

  // If the line is steep, remember it and swap dx and dy so dx >
  // dy. This makes sure that x will be our major axis and y our
  // minor.
  bool bSteep = false;
  if ( dy > dx ) {
    bSteep = true;
    int temp = dx;
    dx = dy;
    dy = temp;
  }

  // Calculate our step increments based on the negative y-ness and
  // steep-ness.
  int adjInc, diagInc, orthInc;
  if ( !bSteep && !bNegativeY ) {
    adjInc  = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    orthInc = ( 0 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
  } else if ( bSteep && !bNegativeY ) {
    adjInc  = ( 0 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    orthInc = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
  } else if ( !bSteep && bNegativeY ) {
    adjInc  = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + (-1 * iWidth * mBytesPerPixel);
    orthInc = ( 0 * mBytesPerPixel) + (-1 * iWidth * mBytesPerPixel);
  } else { //  bSteep && bNegativeY
    adjInc  = ( 0 * mBytesPerPixel) + (-1 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + (-1 * iWidth * mBytesPerPixel);
    orthInc = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
  }

  // Calculate the slope. In the divide by zero case, set it to a big
  // number.
  float slope;
  if ( 0 == dx ) {
    slope = 100000;
  } else {
    slope = dy / dx;
  }
  float Poinc = sqrt( slope );
  if ( 0 == Poinc ) Poinc = 1;
  float Painc = slope * Poinc;
  float Pdinc = Painc - Poinc;
  float Pmid = 0;
  float Pnow = 0;

  float Bainc = dy * 2.0;
  float Bdinc = (dy-dx) * 2.0;
  float Bvar = Bainc - dx;

  GLubyte* midPixel = iBuffer + ( x1 * mBytesPerPixel ) + ( y1 * iWidth * mBytesPerPixel );
  GLubyte* curPixel;

  do {
    // Draw at the middle pixel with our current Pmid.
    DrawPixelIntoBuffer( midPixel, iColor, Pmid );

    // Walk up the minor axis.
    for ( Pnow = Poinc - Pmid,   curPixel = midPixel + (GLubyte)orthInc;
          Pnow < iThickness;
          Pnow += Poinc,         curPixel += (GLubyte)orthInc ) {
      DrawPixelIntoBuffer( curPixel, iColor, Pnow );
    }

    // Walk down the minor axis.
    for ( Pnow = Poinc - Pmid,   curPixel = midPixel - (GLubyte)orthInc;
          Pnow < iThickness;
          Pnow += Poinc,         curPixel -= (GLubyte)orthInc ) {
      DrawPixelIntoBuffer( curPixel, iColor, Pnow );
    }

    if ( Bvar < 0 ) {
      Bvar += Bainc;
      midPixel += adjInc;
      Pmid += Painc;
    } else {
      Bvar += Bdinc;
      midPixel += diagInc;
      Pmid += Pdinc;
    }

    dx--;
  } while ( dx >= 0 );

}

void
Layer::DrawLineIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                            int iFromWindow[2], int iToWindow[2],
                            int iColor[3], int iThickness, float iOpacity ) {


  list<Point2<int> > points;
  Utilities::FindPointsOnLine2d( iFromWindow, iToWindow, points );

  list<Point2<int> >::iterator tPoints;
  for ( tPoints = points.begin(); tPoints != points.end(); ++tPoints ) {

    Point2<int>& point = *tPoints;
    int window[2];
    window[0] = point.x();
    window[1] = point.y();
    DrawPixelIntoBuffer( iBuffer, iWidth, iHeight, window,
                         iColor, iOpacity );
  }
}

void
Layer::GetPreferredThroughPlaneIncrements ( float oIncrements[3] ) {

  oIncrements[0] = 1.0;
  oIncrements[1] = 1.0;
  oIncrements[2] = 1.0;
}

float
Layer::GetPreferredBrushRadiusIncrement () {

  return 1.0;
}

float
Layer::GetPreferredValueIncrement () {

  return 1.0;
}

void
Layer::DoTimer () {}

LayerStaticTclListener::LayerStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetLayerIDList", 0, "",
                         "Return a list of all layerIDs." );
}

LayerStaticTclListener::~LayerStaticTclListener () {}

TclCommandListener::TclCommandResult
LayerStaticTclListener::DoListenToTclCommand ( char* isCommand,
    int, char** ) {

  // GetLayerIDList
  if ( 0 == strcmp( isCommand, "GetLayerIDList" ) ) {
    list<int> idList;
    Layer::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::iterator tID;
    for ( tID = idList.begin(); tID != idList.end(); ++tID ) {
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



Layer::InfoAtRAS::InfoAtRAS()
    : msLabel(""), msValue(""), msTclCallback(""), msInputFilter(""),
    mbShortenHint(false) {}

void
Layer::InfoAtRAS::Clear () {
  msLabel = msValue = msTclCallback = msInputFilter = "";
  mbShortenHint = false;
}
