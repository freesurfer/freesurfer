/**
 * @file  ScubaFrame.cpp
 * @brief Scuba specific subclass of WindowFrame that manages Views
 *
 * This class adds the capability to display Views, which are panels
 * in the GL context that have their own display state and event
 * handling. The ScubaFrame creates and arranges the views according
 * to ViewConfigurations. It also listens to Tcl commands and
 * Broadcast messages. It also owns the tool state.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.45 $
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


#include "ScubaFrame.h"
#include <math.h>
#include "PreferencesManager.h"
extern "C" {
#include "glut.h"
#include "rgb_image.h"
#include "tiffio.h"
}
#include "ScubaView.h"

using namespace std;

ViewFactory* ScubaFrame::mFactory = NULL;

ScubaFrame::ScubaFrame() :
    WindowFrame() ,
    Listener( "ScubaFrame" ) {

  SetOutputStreamToCerr();

  mcRows = 0;
  mnSelectedViewCol = 0;
  mnSelectedViewRow = 0;

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
  commandMgr.AddCommand( *this, "SetSelectedViewID", 2, "frameID viewID",
                         "Sets the select view in a frame." );
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
  commandMgr.AddCommand( *this, "CopyViewLayersToAllViewsInFrame", 2,
                         "frameID viewID", "Copies the layer settings "
                         "in a view to all other views in a frame." );
  commandMgr.AddCommand( *this, "GetToolIDForFrame", 1, "frameID",
                         "Returns the ID of the tool for this frame." );
  commandMgr.AddCommand( *this, "CycleCurrentViewInFrame", 1, "frameID",
                         "Selects the next view in a frame." );
  commandMgr.AddCommand( *this, "CaptureFrameToFile", 2, "frameID fileName",
                         "Make a screen capture of the frame." );
}

ScubaFrame::~ScubaFrame() {

  // Stop listening to our views.
  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      try {
        View* view;
        view = GetViewAtColRow( nCol, nRow );
        if ( view )
          view->RemoveListener( *this );
      } catch (...) {}
    }
  }
}

TclCommandListener::TclCommandResult
ScubaFrame::DoListenToTclCommand( char* isCommand, int, char** iasArgv ) {

  // SetFrameViewConfiguration <frameID> <configuration>
  if ( 0 == strcmp( isCommand, "SetFrameViewConfiguration" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      ViewConfiguration config;
      if ( 0 == strcmp( iasArgv[2], "c1" ) ) {
        config = c1;
      } else if ( 0 == strcmp( iasArgv[2], "c22" ) ) {
        config = c22;
      } else if ( 0 == strcmp( iasArgv[2], "c44" ) ) {
        config = c44;
      } else if ( 0 == strcmp( iasArgv[2], "c13" ) ) {
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
  if ( 0 == strcmp( isCommand, "GetViewIDFromFrameColRow" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int nCol = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad row";
        return error;
      }

      int nRow = strtol(iasArgv[3], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad col";
        return error;
      }

      if ( nRow >= 0 && nRow < mcRows ) {

        if ( nCol >= 0 && nCol < mcCols[nRow] ) {

          try {
            View* view = GetViewAtColRow( nCol, nRow );
            int id = view->GetID();
            stringstream sID;
            sID << id;
            sReturnFormat = "i";
            sReturnValues = sID.str();
            return ok;
          } catch (...) {
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

  // SetSelectedViewID <frameID> <viewID>
  if ( 0 == strcmp( isCommand, "SetSelectedViewID" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad view ID";
        return error;
      }

      try {
        SetSelectedView( viewID );
      } catch (...) {
        stringstream ssError;
        ssError << "View ID " << viewID << " is not currently on screen.";
        sResult = ssError.str();
        return error;
      }
    }
  }


  // GetSelectedViewID <frameID>
  if ( 0 == strcmp( isCommand, "GetSelectedViewID" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      try {
        View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
        int id = view->GetID();
        stringstream sID;
        sID << id;
        sReturnFormat = "i";
        sReturnValues = sID.str();
        return ok;
      } catch (...) {
        stringstream sError;
        sError << "couldn't get selected view at col "
        << mnSelectedViewCol << " row " << mnSelectedViewRow;
        sResult = sError.str();
        return error;
      }
    }
  }

  // GetNumberOfRowsInFrame <frameID>
  if ( 0 == strcmp( isCommand, "GetNumberOfRowsInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }
    if ( mID == frameID ) {
      stringstream cRows;
      cRows << mcRows;
      sReturnFormat = "i";
      sReturnValues = cRows.str();
      return ok;
    }
  }

  // GetNumberOfColsAtRowInFrame <frameID> <row>
  if ( 0 == strcmp( isCommand, "GetNumberOfColsAtRowInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {
      int row = strtol(iasArgv[2], (char**)NULL, 10);
      if ( row >= 0 && row < mcRows ) {
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
  if ( 0 == strcmp( isCommand, "GetColumnOfViewInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad view ID";
        return error;
      }

      for ( int nRow = 0; nRow < mcRows; nRow++ ) {
        int cCols = mcCols[nRow];
        for ( int nCol = 0; nCol < cCols; nCol++ ) {

          try {
            View* view = GetViewAtColRow( nCol, nRow );
            if ( view->GetID() == viewID ) {
              stringstream ssReturn;
              ssReturn << nCol;
              sReturnFormat = "i";
              sReturnValues = ssReturn.str();
              return ok;
            }
          } catch (...) {
            stringstream sError;
            sError << "couldn't get column of view " << viewID;
            sResult = sError.str();
            return error;
          }
        }
      }

      sResult = "bad view ID";
      return error;
    }
  }

  // GetRowOfViewInFrame <frameID> <viewID>
  if ( 0 == strcmp( isCommand, "GetRowOfViewInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad view ID";
        return error;
      }

      for ( int nRow = 0; nRow < mcRows; nRow++ ) {
        int cCols = mcCols[nRow];
        for ( int nCol = 0; nCol < cCols; nCol++ ) {

          try {
            View* view = GetViewAtColRow( nCol, nRow );
            if ( view->GetID() == viewID ) {
              stringstream ssReturn;
              ssReturn << nRow;
              sReturnFormat = "i";
              sReturnValues = ssReturn.str();
              return ok;
            }
          } catch (...) {
            stringstream sError;
            sError << "couldn't get row of view " << viewID;
            sResult = sError.str();
            return error;
          }
        }
      }

      sResult = "bad view ID";
      return error;
    }
  }

  // GetViewIDAtFrameLocation <frameID> <windowX> <windowY>
  if ( 0 == strcmp( isCommand, "GetViewIDAtFrameLocation" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int windowX = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad x";
        return error;
      }

      int windowY = strtol(iasArgv[3], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad y";
        return error;
      }

      // This actually happens quite a lot but it's ok, so just return
      // nothing.
      if ( windowX < 0 || windowX >= GetWidth() ||
           windowY < 0 || windowY >= GetHeight() ) {
        sReturnFormat = "i";
        sReturnValues = "-1";
        return ok;
      }

      try {
        // We need to y flip this since we're getting the coords right
        // from Tcl, just like we y flip them in WindowFrame.
        int windowCoords[2];
        windowCoords[0] = windowX;
        windowCoords[1] = (mHeight - windowY);
        View* view = FindViewAtWindowLoc( windowCoords, NULL, NULL );
        if ( NULL != view ) {
          int id = view->GetID();
          stringstream sID;
          sID << id;
          sReturnFormat = "i";
          sReturnValues = sID.str();
          return ok;
        }
      } catch (...) {
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
  if ( 0 == strcmp( isCommand, "RedrawFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      RequestRedisplay();
    }
  }

  // CopyViewLayersToAllViewsInFrame <frameID> <viewID>
  if ( 0 == strcmp( isCommand, "CopyViewLayersToAllViewsInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      int viewID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad view ID";
        return error;
      }

      try {
        CopyLayerSettingsToAllViews( viewID );
      } catch (...) {
        sResult = "Couldn't auto-configure.";
        return error;
      }

      RequestRedisplay();
    }
  }

  // GetToolIDForFrame <frameID>
  if ( 0 == strcmp( isCommand, "GetToolIDForFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      stringstream sID;
      sID << mTool.GetID();
      ;
      sReturnFormat = "i";
      sReturnValues = sID.str();
    }
  }

  // CycleCurrentViewInFrame <frameID>
  if ( 0 == strcmp( isCommand, "CycleCurrentViewInFrame" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad frame ID";
      return error;
    }

    if ( mID == frameID ) {

      if ( mnSelectedViewCol < mcCols[mnSelectedViewRow] - 1 ) {
        mnSelectedViewCol++;
      } else if ( mnSelectedViewRow < mcRows - 1 ) {
        mnSelectedViewRow++;
        mnSelectedViewCol = 0;
      } else {
        mnSelectedViewRow = 0;
        mnSelectedViewCol = 0;
      }

    }
  }

  // CaptureFrameToFile <frameID> <fileName>
  if ( 0 == strcmp( isCommand, "CaptureFrameToFile" ) ) {
    int frameID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad view ID";
      return error;
    }

    if ( mID == frameID ) {

      string fn( iasArgv[2] );
      CaptureToFile( fn );
    }
  }



  return ok;
}

void
ScubaFrame::DoListenToMessage ( string isMessage, void* ) {

  if ( isMessage == "viewChanged" ) {
    RequestRedisplay();
  }
}

View*
ScubaFrame::GetSelectedView () {

  try {
    View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );
    return view;
  } catch (...) {
    throw runtime_error( "Couldn't find selected view." );
  }
}

void
ScubaFrame::SetSelectedView ( int iViewID ) {

  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      try {
        View* view = GetViewAtColRow( nCol, nRow );
        if ( view->GetID() == iViewID ) {
          mnSelectedViewRow = nRow;
          mnSelectedViewCol = nCol;
          return;
        }
      } catch (...) {
        stringstream ssError;
        ssError << "couldn't get view at col "
        << nCol << " row " << nRow;
        throw runtime_error( ssError.str() );
      }
    }
  }
}

void
ScubaFrame::TranslateWindowToView ( int iWindow[2], int inCol, int inRow,
                                    int oView[2] ) {

  oView[0] = iWindow[0] - (inCol * (mWidth / mcCols[inRow]));
  oView[1] = iWindow[1] - ((mcRows-1 - inRow) * (mHeight/mcRows));
}

void
ScubaFrame::SizeViewsToConfiguration() {


  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      View* view;
      try {
        view = GetViewAtColRow( nCol, nRow );
        view->Reshape( mWidth / cCols, mHeight / mcRows );
      } catch (...) {
        DebugOutput( << "Couldn't find a view where there was supposed "
                     << "to be one: " << nCol << ", " << nRow );
      }

    }
  }
}

void
ScubaFrame::DoDraw() {

  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
      glOrtho( 0, mWidth/cCols-1, 0, mHeight/mcRows-1, -1.0, 1.0 );
      glMatrixMode( GL_MODELVIEW );

      // We change the y position so that the 0,0 view is in the top
      // left corner and the cCols-1,mcRows-1 view is in the bottom
      // right, even tho the GL port's 0,0 is in the lower left
      // corner.
      int x = (mWidth/cCols) * nCol;
      int y = (mHeight/mcRows) * (mcRows - nRow-1);

      // Use glViewport to change the origin of the gl context so our
      // views can start drawing at 0,0.
      glViewport( x, y, mWidth/cCols-1, mHeight/mcRows-1 );
      glRasterPos2i( 0, 0 );

      // Tell the view to draw.
      try {
        View* view = GetViewAtColRow( nCol, nRow );
        view->Draw();
      } catch (...) {
        cerr << "Error drawing view at " << nCol << ", " << nRow << endl;
      }

      // If this is our selected view, draw a green box around it.
      if ( mcRows > 1 && cCols > 1 &&
           nRow == mnSelectedViewRow &&
           nCol == mnSelectedViewCol ) {

        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        glOrtho( 0, mWidth-1, 0, mHeight-1, -1.0, 1.0 );
        glMatrixMode( GL_MODELVIEW );
        glColor3f ( 0.0, 1.0, 0.0 );
        glViewport( 0, 0, mWidth-1, mHeight-1 );
        glBegin( GL_LINE_STRIP );
        glVertex2d( x, y );
        glVertex2d( x + (mWidth/cCols)-2, y );
        glVertex2d( x + (mWidth/cCols)-2, y + (mHeight/mcRows)-2 );
        glVertex2d( x, y + (mHeight/mcRows)-2 );
        glVertex2d( x, y );
        glEnd ();
      }
    }
  }
}

void
ScubaFrame::DoReshape() {

  glViewport( 0, 0, mWidth, mHeight );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho( 0, mWidth-1, 0, mHeight-1, -1.0, 1.0 );
  glMatrixMode( GL_MODELVIEW );

  SizeViewsToConfiguration();
}

void
ScubaFrame::DoTimer() {

#if 0
  // In our timer function we scan our views and ask if they want
  // redisplays.
  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      View* view;
      try {
        view = GetViewAtColRow( nCol, nRow );
        view->Timer();
        if ( view->WantRedisplay() ) {
          RequestRedisplay();
          view->RedisplayPosted();
        }
      } catch (...) {
        DebugOutput( << "Couldn't find a view where there was supposed "
                     << "to be one: " << nCol << ", " << nRow );
      }
    }
  }
#endif
}

void
ScubaFrame::DoMouseMoved( int iWindow[2], InputState& iInput ) {

  int nRow, nCol;

  View* view = FindViewAtWindowLoc( iWindow, &nCol, &nRow );
  if ( NULL != view ) {
    int viewCoords[2];
    TranslateWindowToView( iWindow, nCol, nRow, viewCoords );
    view->MouseMoved( viewCoords, iInput, mTool );

    if ( view->WantRedisplay() ) {
      RequestRedisplay();
      view->RedisplayPosted();
    }
  }
}

void
ScubaFrame::DoMouseUp( int iWindow[2], InputState& iInput ) {

  int nRow, nCol;
  View* view = FindViewAtWindowLoc( iWindow, &nCol, &nRow );
  if ( NULL != view ) {
    int viewCoords[2];
    TranslateWindowToView( iWindow, nCol, nRow, viewCoords );
    view->MouseUp( viewCoords, iInput, mTool );

    if ( view->WantRedisplay() ) {
      RequestRedisplay();
      view->RedisplayPosted();
    }
  }
}

void
ScubaFrame::DoMouseDown( int iWindow[2], InputState& iInput ) {

  int nRow, nCol;
  View* view = FindViewAtWindowLoc( iWindow, &nCol, &nRow );
  if ( NULL != view ) {

    // Select this view and request a redisplay that we can draw our
    // frame around it.
    mnSelectedViewCol = nCol;
    mnSelectedViewRow = nRow;

    RequestRedisplay();

    int viewCoords[2];
    TranslateWindowToView( iWindow, nCol, nRow, viewCoords );
    view->MouseDown( viewCoords, iInput, mTool );
  }
}

void
ScubaFrame::DoKeyDown( int iWindow[2], InputState& iInput ) {

  View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );

  int viewCoords[2];
  TranslateWindowToView( iWindow, mnSelectedViewCol, mnSelectedViewRow,
                         viewCoords );

  view->KeyDown( viewCoords, iInput, mTool );

  if ( view->WantRedisplay() ) {
    RequestRedisplay();
    view->RedisplayPosted();
  }
}

void
ScubaFrame::DoKeyUp( int iWindow[2], InputState& iInput ) {

  View* view = GetViewAtColRow( mnSelectedViewCol, mnSelectedViewRow );

  int viewCoords[2];
  TranslateWindowToView( iWindow, mnSelectedViewCol, mnSelectedViewRow,
                         viewCoords );

  view->KeyUp( viewCoords, iInput, mTool );

  if ( view->WantRedisplay() ) {
    RequestRedisplay();
    view->RedisplayPosted();
  }
}

void
ScubaFrame::SetViewConfiguration( ScubaFrame::ViewConfiguration iConfig ) {

  mViewConfiguration = iConfig;

  int cNewRows;
  map<int,int> cNewCols;

  switch ( mViewConfiguration ) {
  case c1:
    cNewRows = 1;
    cNewCols[0] = 1;
    break;
  case c22:
    cNewRows = 2;
    cNewCols[0] = cNewCols[1] = 2;
    break;
  case c44:
    cNewRows = 4;
    cNewCols[0] = cNewCols[1] = cNewCols[2] = cNewCols[3] = 4;
    break;
  case c13:
    cNewRows = 2;
    cNewCols[0] = 1;
    cNewCols[1] = 3;
    break;
  default:
    cNewRows = 0;
  }

  // First disable existing views that won't be in the new
  // configuration.
  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    if ( cNewCols[nRow] < cCols ) {
      for ( int nCol = cNewCols[nRow]; nCol < cCols; nCol++ ) {
        View* view;
        try {
          view = GetViewAtColRow( nCol, nRow );
          view->SetVisibleInFrame( false );
        } catch (...) {
          cerr << "didn't get view col " << nCol << " row " << nRow << endl;
        }
      }
    }
  }

  // Now make sure all the new views are made, saving the new
  // configuration in the process.
  mcRows = cNewRows;
  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow] = cNewCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      View* view;
      try {
        view = GetViewAtColRow( nCol, nRow );
        view->SetVisibleInFrame( true );
      } catch (...) {
        if ( NULL != mFactory ) {

          view = mFactory->NewView();
          SetViewAtColRow( nCol, nRow, view );

          stringstream sID;
          sID << nCol << ", " << nRow;
          view->SetLabel( sID.str() );

          view->SetVisibleInFrame( true );

          view->AddListener( *this );

        } else {
          DebugOutput( << "Couldn't create new view because factory "
                       << "has not been set" );
        }
      }
    }
  }

  SizeViewsToConfiguration();

  // Make sure the selected col/row are in bounds.
  if ( mnSelectedViewRow >= mcRows ) {
    mnSelectedViewRow = mcRows - 1;
  }
  if ( mnSelectedViewCol >= mcCols[mnSelectedViewRow] ) {
    mnSelectedViewRow = mcCols[mnSelectedViewRow] - 1;
  }


  RequestRedisplay();
}

string
ScubaFrame::GetViewConfigurationAsString () {

  switch ( mViewConfiguration ) {
  case c1:
    return "c1";
    break;
  case c22:
    return "c22";
    break;
  case c44:
    return "c44";
    break;
  case c13:
    return "c13";
    break;
  }

  return "Unknown";
}

View*
ScubaFrame::GetViewAtColRow ( int iCol, int iRow ) {

  try {
    View* view = (mViews[iRow])[iCol];
    if ( NULL == view ) {
      throw logic_error( "no view" );
    }
    return view;
  } catch (...) {
    stringstream sError;
    sError << "Requested view that doesn't exist at "
    << iCol << ", " << iRow;
    throw out_of_range( sError.str() );
    ;
  }
}

void
ScubaFrame::SetViewAtColRow( int iCol, int iRow, View* const iView ) {
  (mViews[iRow])[iCol] = iView;
}

void
ScubaFrame::CopyLayerSettingsToAllViews ( int iViewID ) {

  View& srcView = View::FindByID( iViewID );
// ScubaView& scubaView = dynamic_cast<ScubaView&>(view);
  ScubaView& srcScubaView = (ScubaView&)srcView;

  // For each view in this frame...
  for ( int nRow = 0; nRow < mcRows; nRow++ ) {
    int cCols = mcCols[nRow];
    for ( int nCol = 0; nCol < cCols; nCol++ ) {

      try {
        View* destView = GetViewAtColRow( nCol, nRow );
        ScubaView& destScubaView = *(ScubaView*)destView;
        if ( destView->GetID() != iViewID ) {
          srcScubaView.CopyLayerSettingsToView( destScubaView );
        }
      } catch (...) {
        DebugOutput( << "Couldn't find a view where there was supposed "
                     << "to be one: " << nCol << ", " << nRow );
      }
    }
  }
}

int
ScubaFrame::GetNumberOfRows () {
  return mcRows;
}

int
ScubaFrame::GetNumberOfColsAtRow ( int inRow ) {
  return mcCols[inRow];
}

View*
ScubaFrame::FindViewAtWindowLoc( int iWindow[2], int* onCol, int* onRow ) {

  try {

    int nRow = (int)floor( (float)(mHeight - iWindow[1]) /
                           ((float)mHeight / (float)mcRows) );
    int nCol = (int)floor( (float)iWindow[0] /
                           ((float)mWidth / (float)mcCols[nRow]) );

    if ( NULL != onCol )
      *onCol = nCol;
    if ( NULL != onRow )
      *onRow = nRow;

    return GetViewAtColRow( nCol, nRow );
  } catch (...) {
    return NULL;
  }
}

void
ScubaFrame::CaptureToFile ( string ifn ) {

  GLint rowlength, skiprows, skippixels, alignment;
  GLboolean swapbytes, lsbfirst;
  GLubyte* pixelData = NULL;

  try {

    // Read from the front buffer.
    glReadBuffer( GL_FRONT );

    // Save our unpack attributes.
    glGetBooleanv(GL_PACK_SWAP_BYTES, &swapbytes);
    glGetBooleanv(GL_PACK_LSB_FIRST, &lsbfirst);
    glGetIntegerv(GL_PACK_ROW_LENGTH, &rowlength);
    glGetIntegerv(GL_PACK_SKIP_ROWS, &skiprows);
    glGetIntegerv(GL_PACK_SKIP_PIXELS, &skippixels);
    glGetIntegerv(GL_PACK_ALIGNMENT, &alignment);

    // Set them.
    glPixelStorei(GL_PACK_SWAP_BYTES, GL_FALSE);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    // Read RGB pixel data.
    pixelData = (GLubyte*) malloc( mWidth * mHeight * 3 );
    glReadPixels(0, 0, mWidth, mHeight, GL_RGB,
                 GL_UNSIGNED_BYTE, (GLvoid*)pixelData);
    GLenum eGL = glGetError ();
    if ( GL_NO_ERROR != eGL ) {
      throw runtime_error( "Error reading pixels" );
    }

    // Restore the attributes.
    glPixelStorei(GL_PACK_SWAP_BYTES, swapbytes);
    glPixelStorei(GL_PACK_LSB_FIRST, lsbfirst);
    glPixelStorei(GL_PACK_ROW_LENGTH, rowlength);
    glPixelStorei(GL_PACK_SKIP_ROWS, skiprows);
    glPixelStorei(GL_PACK_SKIP_PIXELS, skippixels);
    glPixelStorei(GL_PACK_ALIGNMENT, alignment);

    // Open a TIFF.
    char fn[1000];
    strcpy( fn, ifn.c_str() );
    TIFF* fTIFF = TIFFOpen( fn, "w" );
    if ( NULL == fTIFF ) {
      throw runtime_error( "Couldn't create file." );
    }

    // Set the TIFF info.
    TIFFSetField( fTIFF, TIFFTAG_IMAGEWIDTH, mWidth );
    TIFFSetField( fTIFF, TIFFTAG_IMAGELENGTH, mHeight );
    TIFFSetField( fTIFF, TIFFTAG_SAMPLESPERPIXEL, 3 );
    TIFFSetField( fTIFF, TIFFTAG_BITSPERSAMPLE, 8 );
    TIFFSetField( fTIFF, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
    TIFFSetField( fTIFF, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
    TIFFSetField( fTIFF, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB );

    // Calculate some sizes and allocate a line buffer.
    tsize_t cLineBytes = 3 * mWidth;
    unsigned char* lineBuffer = NULL;
    int scanLineSize = TIFFScanlineSize( fTIFF );
    if ( scanLineSize != cLineBytes ) {
      cerr << "scanLineSize = " << scanLineSize << ", cLineBytes = "
      << cLineBytes << endl;
    }

    lineBuffer = (unsigned char*) _TIFFmalloc( scanLineSize  );
    if ( NULL == lineBuffer ) {
      throw runtime_error( "Couldn't allocate line buffer." );
    }

    // Set the strip size to default.
    int stripSize = TIFFDefaultStripSize( fTIFF, mWidth * 3 );
    TIFFSetField( fTIFF, TIFFTAG_ROWSPERSTRIP, stripSize );

    // Write line by line (bottom to top).
    for ( int nRow = 0; nRow < mHeight; nRow++ ) {
      memcpy( lineBuffer,
              &pixelData[(mHeight-nRow-1) * cLineBytes],
              cLineBytes );
      TIFFWriteScanline( fTIFF, lineBuffer, nRow, 0 );
    }

    // Close the tiff file and free the line buffer.
    TIFFClose( fTIFF );
    _TIFFfree( lineBuffer );
  } catch ( runtime_error& e) {
    throw runtime_error( string("Error saving TIFF: ") + e.what() );
  } catch (...) {
    throw runtime_error( "Error saving TIFF: " );
  }

  if ( NULL != pixelData )
    free( pixelData );
}


