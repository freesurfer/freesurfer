#include "ScubaFrame.h"
extern "C" {
#include "glut.h"
}

using namespace std;

ViewFactory* ScubaFrame::mFactory = NULL;

ScubaFrame::ScubaFrame( ToglFrame::ID iID ) 
  : ToglFrame( iID ) {

  DebugOutput( << "Created ScubaFrame " << iID );
  SetOutputStreamToCerr();

  SetViewConfiguration( c11 );

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetViewConfiguration" );
}

ScubaFrame::~ScubaFrame() {
}

void
ScubaFrame::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  if( 0 == strcmp( isCommand, "SetViewConfiguration" ) ) {
    if( 3 == iArgc ) {
      
      int frameID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad frame ID";
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
	  sResult = "bad configuration \"" + string(iasArgv[1]) + 
	    "\", should be c11, c22, c44, or c13";
	  return;
	}
      
	SetViewConfiguration( config );
      }

    } else {
      sResult = "wrong # args: should be \"SetViewConfiguration "
	"frameID configuration\"";
      return;
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
      glViewport( x, y, mWidth, mHeight );

      try {
	View* view = GetViewAtColRow( nCol, nRow );
	view->Draw();
      }
      catch(...){
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
ScubaFrame::DoMouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY );
    view->MouseMoved( inX, inY, iButton, iModifiers );
  }
  catch(...) {
  }
}

void
ScubaFrame::DoMouseUp( int inX, int inY, int iButton, int iModifiers ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY );
    view->MouseUp( inX, inY, iButton, iModifiers );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoMouseDown( int inX, int inY, int iButton, int iModifiers ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY );
    view->MouseDown( inX, inY, iButton, iModifiers );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoKeyDown( int inX, int inY, string isKey, int iModifiers ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY );
    view->KeyDown( inX, inY, isKey, iModifiers );
  }
  catch(...) {
  } 
}

void
ScubaFrame::DoKeyUp( int inX, int inY, string isKey, int iModifiers ) {

  try {
    View* view = FindViewAtWindowLoc( inX, inY );
    view->KeyUp( inX, inY, isKey, iModifiers );
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
ScubaFrame::FindViewAtWindowLoc( int iWindowX, int iWindowY ) {

  try {

    int nRow = iWindowY / (mHeight / mcRows);
    int nCol = iWindowX / (mWidth / mcCols[nRow]);

    return GetViewAtColRow( nCol, nRow );
  }
  catch(...) {
    return NULL;
  }
}


