#include "ScubaFrame.h"
#include "PreferencesManager.h"
extern "C" {
#include "glut.h"
}

using namespace std;

ViewFactory* ScubaFrame::mFactory = NULL;

ScubaFrame::ScubaFrame( ToglFrame::ID iID ) 
  : ToglFrame( iID ) {

  DebugOutput( << "Created ScubaFrame " << iID );
  SetOutputStreamToCerr();

  mnSelectedViewCol = 0;
  mnSelectedViewRow = 0;

  SetViewConfiguration( c1 );

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetFrameViewConfiguration", 2, 
			 "frameID, configuration", 
			 "Sets a frame's view configuration. Supported "
			 "configurations, where each number is the number "
			 "of columns in a row: c1 c22 c44 c13" );
  commandMgr.AddCommand( *this, "GetViewIDFromFrameColRow", 3, 
			 "frameID col row",
			 "Return the viewID from a view at a certain "
			 "location. col and row must be valid for the current "
			 "view configuration." );
  commandMgr.AddCommand( *this, "GetSelectedViewID", 1, "frameID",
			 "Return the viewID of the selected view." );
  commandMgr.AddCommand( *this, "GetNumberOfRowsInFrame", 1, "frameID",
			 "Return the number of rows in a frame." );
  commandMgr.AddCommand( *this, "GetNumberOfColsAtRowInFrame", 2, 
			 "frameID row", "Return the number of columns in a"
			 "row." );
  commandMgr.AddCommand( *this, "GetViewIDAtFrameLocation", 3, 
			 "frameID windowX windowY", 
			 "Return the view ID at a window location." );
  commandMgr.AddCommand( *this, "GetColumnOfViewInFrame", 2, "frameID viewID", 
			 "Return the column of the view ID in a frame." );
  commandMgr.AddCommand( *this, "GetRowOfViewInFrame", 2, "frameID viewID", 
			 "Return the row of the view ID in a frame." );
  commandMgr.AddCommand( *this, "RedrawFrame", 1, "frameID", 
			 "Tells a frame to redraw." );

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.UseFile( ".scuba" );
  PreferencesManager::StringPrefValue cycleKey( "q" );
  prefsMgr.RegisterValue( "key-CycleViewsInFrame", 
			  "Key to cycle view in a frame.", cycleKey );
  msCycleKey = prefsMgr.GetValue( "key-CycleViewsInFrame" );
}

ScubaFrame::~ScubaFrame() {
}

TclCommandListener::TclCommandResult
ScubaFrame::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetFrameViewConfiguration <frameID> <configuration>
  if( 0 == strcmp( isCommand, "SetFrameViewConfiguration" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      ViewConfiguration config;
      if( 0 == strcmp( iasArgv[2], "c1" ) ) {
	config = c1;
      } else if( 0 == strcmp( iasArgv[2], "c22" ) ) {
	config = c22;
      } else if( 0 == strcmp( iasArgv[2], "c44" ) ) {
	config = c44;
      } else if( 0 == strcmp( iasArgv[2], "c13" ) ) {
	config = c13;
      } else {
	sResult = "bad configuration \"" + string(iasArgv[2]) + 
	  "\", should be c1, c22, c44, or c13";
	return error;
      }
      
      SetViewConfiguration( config );
    }
  }

  // GetViewIDFromFrameColRow <frameID> <col> <row>
  if( 0 == strcmp( isCommand, "GetViewIDFromFrameColRow" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      int nCol = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad row";
	return error;
      }
      
      int nRow = strtol(iasArgv[3], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad col";
	return error;
      }
      
      if( nRow >= 0 && nRow < mcRows ) {
	
	if( nCol >= 0 && nCol < mcCols[nRow] ) {
	  
	  try { 
	    View* view = GetViewAtColRow( nCol, nRow );
	    int id = view->GetID();
	    stringstream sID;
	    sID << id;
	    sReturnFormat = "i";
	    sReturnValues = sID.str();
	    return ok;
	  }
	  catch(...) {
	    stringstream sError;
	    sError << "couldn't get view at col " << nCol <<
	      " row " << nRow;
	    sResult = sError.str();
	    DebugOutput( << sResult );
	    return error;
	  }
	  
	} else {
	  stringstream sError;
	  sError << "bad col \"" << string(iasArgv[1]) << 
	    "\", should be between 0 and " << mcCols[nRow];
	  sResult = sError.str();
	  DebugOutput( << sResult );
	  return error;
	}
	
      } else {
	stringstream sError;
	sError << "bad row \"" << string(iasArgv[1]) <<
	  "\", should be between 0 and " << mcRows;
	sResult = sError.str();
	DebugOutput( << sResult );
	return error;
      }
      
    }
  }


  // GetSelectedViewID <frameID>
  if( 0 == strcmp( isCommand, "GetSelectedViewID" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      try { 
	View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
	int id = view->GetID();
	stringstream sID;
	sID << id;
	sReturnFormat = "i";
	sReturnValues = sID.str();
	return ok;
      }
      catch(...) {
	stringstream sError;
	sError << "couldn't get selected view at col " 
	       << mnSelectedViewCol << " row " << mnSelectedViewRow;
	sResult = sError.str();
	return error;
      }
    }
  }
  
  // GetNumberOfRowsInFrame <frameID>
  if( 0 == strcmp( isCommand, "GetNumberOfRowsInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    if( mID == frameID ) {
      stringstream cRows;
      cRows << mcRows;
      sReturnFormat = "i";
      sReturnValues = cRows.str();
      return ok;
    }
  }

  // GetNumberOfColsAtRowInFrame <frameID> <row>
  if( 0 == strcmp( isCommand, "GetNumberOfColsAtRowInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      int row = strtol(iasArgv[2], (char**)NULL, 10);
      if( row >= 0 && row < mcRows ) {
	stringstream cCols;
	cCols << mcCols[row];
	sReturnFormat = "i";
	sReturnValues = cCols.str();
	return ok;
      } else {
	sResult = "bad row";
	DebugOutput( << sResult );
	return error;
      }
      
    }
  }

  // GetColumnOfViewInFrame <frameID> <viewID> 
  if( 0 == strcmp( isCommand, "GetColumnOfViewInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return error;
      }
      
      for( int nRow = 0; nRow < mcRows; nRow++ ) {
	int cCols = mcCols[nRow];
	for( int nCol = 0; nCol < cCols; nCol++ ) {
	  
	  try {
	    View* view = GetViewAtColRow( nCol, nRow );
	    if( view->GetID() == viewID ) {
	      stringstream ssReturn;
	      ssReturn << nCol;
	      sReturnFormat = "i";
	      sReturnValues = ssReturn.str();
	      return ok;
	    }
	  } 
	  catch(...) {
	  }
	}
      }

      sResult = "bad view ID";
      return error;
    }
  }

  // GetRowOfViewInFrame <frameID> <viewID> 
  if( 0 == strcmp( isCommand, "GetRowOfViewInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad view ID";
	return error;
      }
      
      for( int nRow = 0; nRow < mcRows; nRow++ ) {
	int cCols = mcCols[nRow];
	for( int nCol = 0; nCol < cCols; nCol++ ) {
	  
	  try {
	    View* view = GetViewAtColRow( nCol, nRow );
	    if( view->GetID() == viewID ) {
	      stringstream ssReturn;
	      ssReturn << nRow;
	      sReturnFormat = "i";
	      sReturnValues = ssReturn.str();
	      return ok;
	    }
	  } 
	  catch(...) {
	  }
	}
      }

      sResult = "bad view ID";
      return error;
    }
  }

  // GetViewIDAtFrameLocation <frameID> <windowX> <windowY>
  if( 0 == strcmp( isCommand, "GetViewIDAtFrameLocation" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      int windowX = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad x";
	return error;
      }
      
      int windowY = strtol(iasArgv[3], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad y";
	return error;
      }
      
      if( windowX < 0 || windowX >= GetWidth() ||
	  windowY < 0 || windowY >= GetHeight() ) {
	sResult = "location is out of bounds";
	return error;
      }
	
      try { 
	// We need to y flip this since we're getting the coords right
	// from Tcl, just like we y flip them in ToglFrame.
	View* view = FindViewAtWindowLoc( windowX, (mHeight - windowY),
					  NULL, NULL );
	if( NULL != view ) {
	  int id = view->GetID();
	  stringstream sID;
	  sID << id;
	  sReturnFormat = "i";
	  sReturnValues = sID.str();
	  return ok;
	}
      }
      catch(...) {
	stringstream sError;
	sError << "couldn't find view at x " << windowX <<
	  " y " << windowY;
	sResult = sError.str();
	DebugOutput( << sResult );
	return error;
      }
    }
  }

  // RedrawFrame <frameID>
  if( 0 == strcmp( isCommand, "RedrawFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    
    if( mID == frameID ) {
      
      RequestRedisplay();
    }
  }

  return ok;
}

void
ScubaFrame::SizeViewsToConfiguration() {

  for( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for( int nCol = 0; nCol < cCols; nCol++ ) {
      
      View* view;
      try {
	view = GetViewAtColRow( nCol, nRow );
	view->Reshape( mWidth / cCols, mHeight / mcRows );
      } 
      catch(...) {
	DebugOutput( << "Couldn't create new view because factory "
		     << "has not been set" );
      }
      
    }
  }  
}

void
ScubaFrame::DoDraw() {
  
  for( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for( int nCol = 0; nCol < cCols; nCol++ ) {

      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
      //      glOrtho( 0, mWidth, 0, mHeight, -1.0, 1.0 );
      glOrtho( 0, mWidth/cCols, 0, mHeight/mcRows, -1.0, 1.0 );
      glMatrixMode( GL_MODELVIEW );

      // We change the y position so that the 0,0 view is in the top
      // left corner and the cCols-1,mcRows-1 view is in the bottom
      // right, even tho the GL port's 0,0 is in the lower left
      // corner. 
      int x = (mWidth/cCols) * nCol;
      int y = (mHeight/mcRows) * (mcRows - nRow-1);

      // Use glViewport to change the origin of the gl context so our
      // views can start drawing at 0,0. However leave the width and
      // height at the frame's width and height, otherwise we'll mess
      // up the proportions.
      //      glViewport( x, y, mWidth, mHeight );
      glViewport( x, y, mWidth/cCols, mHeight/mcRows );
      glRasterPos2i( 0, 0 );

      // Tell the view to draw.
      try {
	View* view = GetViewAtColRow( nCol, nRow );
	view->Draw();
      }
      catch(...){
      }

      // If this is our selected view, draw a green box around it.
      if( nRow == mnSelectedViewRow && 
	  nCol == mnSelectedViewCol ) {

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( 0, mWidth, 0, mHeight, -1.0, 1.0 );
	glMatrixMode( GL_MODELVIEW );
	glColor3f ( 0.0, 1.0, 0.0 );
	glViewport( 0, 0, mWidth, mHeight );
	glBegin( GL_LINE_STRIP );
	glVertex2d( x, y );
	glVertex2d( x + (mWidth/cCols)-1, y );
	glVertex2d( x + (mWidth/cCols)-1, y + (mHeight/mcRows)-1 );
	glVertex2d( x, y + (mHeight/mcRows)-1 );
	glVertex2d( x, y );
	glEnd ();
      }
    }
  }
}

void
ScubaFrame::DoReshape() {

  SizeViewsToConfiguration();

}

void
ScubaFrame::DoTimer() {

  // In our timer function we scan our views and ask if they want
  // redisplays.
  for( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for( int nCol = 0; nCol < cCols; nCol++ ) {
      
      View* view;
      try {
	view = GetViewAtColRow( nCol, nRow );
	if( view->WantRedisplay() ) {
	  RequestRedisplay();
	  view->RedisplayPosted();
	}
      } 
      catch(...) {
      }
    }
  }  
}

void
ScubaFrame::DoMouseMoved( int inX, int inY, InputState& iState ) {

  try {
    int nRow, nCol;
    View* view = FindViewAtWindowLoc( inX, inY, &nCol, &nRow );
    int width = mWidth / mcCols[nRow];
    int height = mHeight / mcRows;
    view->MouseMoved( inX - (nCol * width), 
		      inY - ((mcRows-1 - nRow) * height), iState );
  }
  catch(...) {
  }
}

void
ScubaFrame::DoMouseUp( int inX, int inY, InputState& iState ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->MouseUp( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoMouseDown( int inX, int inY, InputState& iState ) {

  try {
    int nRow;
    int nCol;
    View* view = FindViewAtWindowLoc( inX, inY, &nCol, &nRow );

    // Select this view and request a redisplay that we can draw our
    // frame around it.
    mnSelectedViewCol = nCol;
    mnSelectedViewRow = nRow;

    RequestRedisplay();

    view->MouseDown( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoKeyDown( int inX, int inY, InputState& iState ) {

  try {

    // Tab key cycles selected view.
    if( iState.Key() == msCycleKey ) {
      if( mnSelectedViewCol < mcCols[mnSelectedViewRow] - 1 ) {
	mnSelectedViewCol++;
      } else if( mnSelectedViewRow < mcRows - 1 ) {
	mnSelectedViewRow++;
	mnSelectedViewCol = 0;
      } else {
	mnSelectedViewRow = 0;
	mnSelectedViewCol = 0;
      }

      // Redraw the frame.
      RequestRedisplay();  
    }

    View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
    //    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->KeyDown( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoKeyUp( int inX, int inY, InputState& iState ) {

  try {
    View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
    //    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->KeyUp( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::SetViewConfiguration( ScubaFrame::ViewConfiguration iConfig ) {

  mViewConfiguration = iConfig;

  switch( mViewConfiguration ) {
    case c1:
      mcRows = 1;
      mcCols[0] = 1;
      break;
    case c22:
      mcRows = 2;
      mcCols[0] = mcCols[1] = 2;
      break;
    case c44:
      mcRows = 4;
      mcCols[0] = mcCols[1] = mcCols[2] = mcCols[3] = 4;
      break;
    case c13:
      mcRows = 2;
      mcCols[0] = 1;
      mcCols[1] = 3;
      break;
  }

  for( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for( int nCol = 0; nCol < cCols; nCol++ ) {
      
      View* view;
      try {
	view = GetViewAtColRow( nCol, nRow );
      } 
      catch(...) {
	if( NULL != mFactory ) {
	  
	  view = mFactory->NewView();
	  SetViewAtColRow( nCol, nRow, view );

	  stringstream sID;
	  sID << nCol << ", " << nRow;
	  view->SetLabel( sID.str() );
	  
	} else {
	  DebugOutput( << "Couldn't create new view because factory "
		       << "has not been set" );
	}
      }
    }
  }  

  SizeViewsToConfiguration();

  // Make sure the selected col/row are in bounds. 
  if( mnSelectedViewRow >= mcRows ) {
    mnSelectedViewRow = mcRows - 1;
  }
  if( mnSelectedViewCol >= mcCols[mnSelectedViewRow] ) {
    mnSelectedViewRow = mcCols[mnSelectedViewRow] - 1;
  }


  RequestRedisplay();
}


View*
ScubaFrame::GetViewAtColRow( int iCol, int iRow ) {
  try { 
    View* view = (mViews[iRow])[iCol];
    if( NULL == view ) {
      throw logic_error( "no view" );
    }
    return view;
  }
  catch(...) {
    stringstream sError;
    sError << "Requested view that doesn't exist at " 
	   << iCol << ", " << iRow;
    throw new out_of_range( sError.str() );;
  }
}

void 
ScubaFrame::SetViewAtColRow( int iCol, int iRow, View* const iView ) {
  (mViews[iRow])[iCol] = iView;
}

View*
ScubaFrame::FindViewAtWindowLoc( int iWindowX, int iWindowY,
				 int* onCol, int* onRow ) {

  try {

    int nRow = (int)floor( (float)(mHeight - iWindowY) / 
			   ((float)mHeight / (float)mcRows) );
    int nCol = (int)floor( (float)iWindowX / 
			   ((float)mWidth / (float)mcCols[nRow]) );
    
    if( NULL != onCol ) 
      *onCol = nCol;
    if( NULL != onRow )
      *onRow = nRow;

    return GetViewAtColRow( nCol, nRow );
  }
  catch(...) {
    return NULL;
  }
}
