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

  SetViewConfiguration( c11 );

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetFrameViewConfiguration" );
  commandMgr.AddCommand( *this, "GetViewIDFromFrameColRow" );
  commandMgr.AddCommand( *this, "GetSelectedViewID" );

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.UseFile( ".scuba" );
  PreferencesManager::StringPrefValue cycleKey( "q" );
  prefsMgr.RegisterValue( "key-CycleViewsInFrame", 
			  "Key to cycle view in a frame.", cycleKey );
  msCycleKey = prefsMgr.GetValue( "key-CycleViewsInFrame" );
}

ScubaFrame::~ScubaFrame() {
}

void
ScubaFrame::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetFrameViewConfiguration <frameID> <configuration>
  if( 0 == strcmp( isCommand, "SetFrameViewConfiguration" ) ) {
    if( 3 == iArgc ) {
      int frameID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad frame ID";
	return;
      }
      
      if( mID == frameID ) {
	
	ViewConfiguration config;
	if( 0 == strcmp( iasArgv[2], "c11" ) ) {
	  config = c11;
	} else if( 0 == strcmp( iasArgv[2], "c22" ) ) {
	  config = c22;
	} else if( 0 == strcmp( iasArgv[2], "c44" ) ) {
	  config = c44;
	} else if( 0 == strcmp( iasArgv[2], "c13" ) ) {
	  config = c13;
	} else {
	  sResult = "bad configuration \"" + string(iasArgv[2]) + 
	    "\", should be c11, c22, c44, or c13";
	  return;
	}
      
	SetViewConfiguration( config );
      }
    } else {
      sResult = "wrong # args: should be \"SetViewConfiguration "
	"frameID configuration\"";
      DebugOutput( << sResult );
      return;
    }
  }

  // GetViewIDFromFrameColRow <frameID> <col> <row>
  if( 0 == strcmp( isCommand, "GetViewIDFromFrameColRow" ) ) {
    if( 4 == iArgc ) {
      int frameID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad frame ID";
	return;
      }
      
      if( mID == frameID ) {
	
	int nCol = strtol(iasArgv[2], (char**)NULL, 10);
	if( ERANGE == errno ) {
	  sResult = "bad row";
	  return;
	}
	
	int nRow = strtol(iasArgv[3], (char**)NULL, 10);
	if( ERANGE == errno ) {
	  sResult = "bad col";
	  return;
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
	      return;
	    }
	    catch(...) {
	      stringstream sError;
	      sError << "couldn't get view at col " << nCol <<
		" row " << nRow;
	      sResult = sError.str();
	      DebugOutput( << sResult );
	      return;
	    }

	  } else {
	    stringstream sError;
	    sError << "bad col \"" << string(iasArgv[1]) << 
	      "\", should be between 0 and " << mcCols[nRow];
	    sResult = sError.str();
	    DebugOutput( << sResult );
	    return;
	  }
	  
	} else {
	  stringstream sError;
	  sError << "bad row \"" << string(iasArgv[1]) <<
	    "\", should be between 0 and " << mcRows;
	  sResult = sError.str();
	  DebugOutput( << sResult );
	  return;
	}
      
      }
    } else {
      sResult = "wrong # args: should be \"GetViewIDFromFrameColRow "
	"frameID row col\"";
      DebugOutput( << sResult << " (iArgc = " << iArgc << ")" );
      return;
    }
  }


  // GetSelectedViewID <frameID>
  if( 0 == strcmp( isCommand, "GetSelectedViewID" ) ) {
    if( 2 == iArgc ) {
      int frameID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad frame ID";
	return;
      }
      
      if( mID == frameID ) {
	
	try { 
	  View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
	  int id = view->GetID();
	  stringstream sID;
	  sID << id;
	  sReturnFormat = "i";
	  sReturnValues = sID.str();
	  return;
	}
	catch(...) {
	  stringstream sError;
	  sError << "couldn't get selected view at col " 
		 << mnSelectedViewCol << " row " << mnSelectedViewRow;
	  sResult = sError.str();
	  return;
	}
      }
    }
  }
}

void
ScubaFrame::DoDraw() {
  
  glClear( GL_COLOR_BUFFER_BIT );

  for( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for( int nCol = 0; nCol < cCols; nCol++ ) {

      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
      glOrtho( 0, mWidth, 0, mHeight, -1.0, 1.0 );
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
      glViewport( x, y, mWidth, mHeight );

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

  ViewRowMap::iterator tRow;
  for( tRow = mViews.begin(); tRow != mViews.end(); ++tRow ) {
    
    ViewColMap::iterator tCol;
    ViewColMap viewCols = (*tRow).second;
    for( tCol = viewCols.begin(); tCol != viewCols.end(); ++tCol ) {

      View* view = (*tCol).second;
      view->Reshape( 0, 0 );
    }
  }
}

void
ScubaFrame::DoTimer() {

}

void
ScubaFrame::DoMouseMoved( int inX, int inY, InputState& iState ) {

  try {
    
    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->MouseMoved( inX, inY, iState );
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

    // Select this view.
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
    }
    
    RequestRedisplay();  

    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->KeyDown( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoKeyUp( int inX, int inY, InputState& iState ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY, NULL, NULL );
    view->KeyUp( inX, inY, iState );
  }
  catch(...) {
  } 
}

void
ScubaFrame::SetViewConfiguration( ScubaFrame::ViewConfiguration iConfig ) {

  mViewConfiguration = iConfig;

  switch( mViewConfiguration ) {
    case c11:
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

      if( NULL != view ) {
	view->SetWidth( mWidth / cCols );
	view->SetHeight( mHeight / mcRows );
      }
      
    }
  }  

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

    int nRow = (int)floor( (float)iWindowY / 
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
