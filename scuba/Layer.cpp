#include <iostream>
#include "Layer.h"
#include "macros.h"
#include "Utilities.h"
#include "Point2.h"

using namespace std;

DeclareIDTracker(Layer);

LayerStaticTclListener Layer::mStaticListener;

Layer::Layer() {
  mOpacity = 1.0;
  msLabel = "";
  mbPostRedisplay = false;
  mBytesPerPixel = 4;

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
Layer::GetInfoAtRAS ( float iRAS[3],
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
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
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
      
      float opacity = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad opacity";
	return error;
      }
      
      SetOpacity( opacity );
    }
  }
  
  // GetLayerPreferredBrushRadiusIncrement <layerID>
  if( 0 == strcmp( isCommand, "GetLayerPreferredBrushRadiusIncrement" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetPreferredBrushRadiusIncrement();
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}


void
Layer::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "dataChanged" ) {
    DataChanged ();
    SendBroadcast( "layerChanged", NULL );
  }
}

void
Layer::DataChanged () {
  // Default behavior of data changing is to request a redisplay.
  cerr << "Layer::DataChanged()" << endl;
  RequestRedisplay();
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
Layer::HandleTool ( float iRAS[3], ViewState& iViewState,
		    ScubaWindowToRASTranslator& iTranslator,
		    ScubaToolState& iTool, InputState& iInput ) {

}


void 
Layer::DrawPixelIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
			     int iWindow[2], int iColor[3], float iOpacity ){

  if( iWindow[0] >= 0 && iWindow[0] < iWidth && 
      iWindow[1] >= 0 && iWindow[1] < iHeight ) {
    GLubyte* dest = iBuffer + (iWindow[0] * mBytesPerPixel) + 
      (iWindow[1] * iWidth * mBytesPerPixel);
    DrawPixelIntoBuffer( dest, iColor, iOpacity );
  }
}

void 
Layer::DrawPixelIntoBuffer ( GLubyte* iAddress,
			     int iColor[3], float iOpacity ){

  iAddress[0] = (GLubyte) (((float)iAddress[0] * (1.0 - iOpacity)) +
			   ((float)iColor[0] * iOpacity));
  iAddress[1] = (GLubyte) (((float)iAddress[1] * (1.0 - iOpacity)) +
			   ((float)iColor[1] * iOpacity));
  iAddress[2] = (GLubyte) (((float)iAddress[2] * (1.0 - iOpacity)) +
			   ((float)iColor[2] * iOpacity));
  iAddress[3] = 255; 
}



void 
Layer::DrawAALineIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
			      int iFromWindow[2], int iToWindow[2],
			      int iColor[3], int iThickness, float iOpacity ) {

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
  if( dy < 0 ) {
    bNegativeY = true;   // This affects
    dy = -dy;
  }

  // If the line is steep, remember it and swap dx and dy so dx >
  // dy. This makes sure that x will be our major axis and y our
  // minor.
  bool bSteep = false;
  if( dy > dx ) {
    bSteep = true;
    int temp = dx; dx = dy; dy = temp;
  }

  // Calculate our step increments based on the negative y-ness and
  // steep-ness.
  int adjInc, diagInc, orthInc;
  if( !bSteep && !bNegativeY ) {
    adjInc  = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    orthInc = ( 0 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
  } else if( bSteep && !bNegativeY ) {
    adjInc  = ( 0 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    diagInc = ( 1 * mBytesPerPixel) + ( 1 * iWidth * mBytesPerPixel);
    orthInc = ( 1 * mBytesPerPixel) + ( 0 * iWidth * mBytesPerPixel);
  } else if( !bSteep && bNegativeY ) {
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
  if( 0 == dx ) {
    slope = 100000;
  } else {
    slope = dy / dx;
  }
  float Poinc = sqrt( slope );
  if( 0 == Poinc ) Poinc = 1;
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
    for( Pnow = Poinc - Pmid,   curPixel = midPixel + (GLubyte)orthInc;
	 Pnow < iThickness;
	 Pnow += Poinc,         curPixel += (GLubyte)orthInc ) {
      DrawPixelIntoBuffer( curPixel, iColor, Pnow );
    }

    // Walk down the minor axis.
    for( Pnow = Poinc - Pmid,   curPixel = midPixel - (GLubyte)orthInc;
	 Pnow < iThickness;
	 Pnow += Poinc,         curPixel -= (GLubyte)orthInc ) {
      DrawPixelIntoBuffer( curPixel, iColor, Pnow );
    }
    
    if( Bvar < 0 ) {
      Bvar += Bainc;
      midPixel += adjInc;
      Pmid += Painc;
    } else {
      Bvar += Bdinc;
      midPixel += diagInc;
      Pmid += Pdinc;
    }
    
    dx--;
  } while( dx >= 0 );

}

void
Layer::DrawLineIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
			    int iFromWindow[2], int iToWindow[2],
			    int iColor[3], int iThickness, float iOpacity ) {
  

  list<Point2<int> > points;
  Utilities::FindPointsOnLine2d( iFromWindow, iToWindow, iThickness,
				 points );

  list<Point2<int> >::iterator tPoints;
  for( tPoints = points.begin(); tPoints != points.end(); ++tPoints ) {

    Point2<int>& point = *tPoints;
    int window[2];
    window[0] = point.x();
    window[1] = point.y();
    DrawPixelIntoBuffer( iBuffer, iWidth, iHeight, window,
			 iColor, iOpacity );
  }
}

void 
Layer::GetPreferredInPlaneIncrements ( float oIncrements[3] ) {
  
  oIncrements[0] = 1.0;
  oIncrements[1] = 1.0;
  oIncrements[2] = 1.0;
}

float
Layer::GetPreferredBrushRadiusIncrement () {
  
  return 1.0;
}



LayerStaticTclListener::LayerStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetLayerIDList", 0, "", 
			 "Return a list of all layerIDs." );
}

LayerStaticTclListener::~LayerStaticTclListener () {
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



