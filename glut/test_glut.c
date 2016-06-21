/**
 * @file  test_glut.c
 * @brief simple glut lib test
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/01/04 19:05:40 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2006-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "../include/xGLutWindow.h"
#include "../include/version.h"

const char* Progname = "test_glut";

/* window data */
static xGLutWindowRef  gWindow = NULL;
static int             gzWidth = 0;
static int             gzHeight = 0;
static GLubyte*        gFrameBuffer = NULL;
static int             gBrushRadius = 1;

#define kzInitialWidth  300
#define kzInitialHeight 300

#define kBytesPerPixel 4

#define ThrowIf(x) if( (x) ) goto error;

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

void EventCallback ( void* iWindow, xGWin_tEventRef iEvent );
int BuildWindow ();
int InitalizeFrameBuffer ();
int Brush ( int iX, int iY );
int EditFrameBuffer( int iX, int iY, int iRed, int iGreen, int iBlue );
int Redraw ();

void EventCallback ( void* iWindow, xGWin_tEventRef iEvent ) {

  switch ( iEvent->mType ) {

  case xGWin_tEventType_KeyDown:
    switch ( iEvent->mKey ) {
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      gBrushRadius = (int)(iEvent->mKey - '1') + 1;
      ;
      break;
    }
    break;

  case xGWin_tEventType_MouseMoved:
    /*
      fprintf( stderr, "moved %d %d(%d)\n",
      iEvent->mWhere.mnX, iEvent->mWhere.mnY,
      (gzHeight - iEvent->mWhere.mnY) );
    */
    Brush( iEvent->mWhere.mnX, (gzHeight - iEvent->mWhere.mnY) );
    Redraw();
    break;

  case xGWin_tEventType_Resize:
    gzWidth  = iEvent->mWhere.mnX;
    gzHeight = iEvent->mWhere.mnY;
    glutReshapeWindow( gzWidth, gzHeight );
    glViewport( 0, 0, gzWidth, gzHeight );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho( 0, gzWidth, 0, gzHeight, -1.0, 1.0 );
    glMatrixMode( GL_MODELVIEW );
    InitalizeFrameBuffer();
    Redraw();
    break;

  case xGWin_tEventType_Draw:
    glRasterPos2i( 0, 0 );
    glDrawPixels( gzWidth, gzHeight, GL_RGBA, GL_UNSIGNED_BYTE, gFrameBuffer );
    glutSwapBuffers();
    break;

  default:
    break;
  }

}

/* Returns 1 if fail. */
int BuildWindow () {

  xGWin_tErr rGWin = xGWin_tErr_NoErr;
  int        rCode = 0;

  rGWin = xGWin_New( &gWindow, kzInitialWidth, kzInitialHeight, "test_glut" );
  ThrowIf( xGWin_tErr_NoErr != rGWin );

  rGWin = xGWin_SetEventHandlerFunc( gWindow, EventCallback, gWindow );
  ThrowIf( xGWin_tErr_NoErr != rGWin );

  rGWin = xGWin_ActivateIdleEvents( gWindow );
  ThrowIf( xGWin_tErr_NoErr != rGWin );

  rGWin = xGWin_ActivatePassiveMotionEvents( gWindow );
  ThrowIf( xGWin_tErr_NoErr != rGWin );

  gzWidth = kzInitialWidth;
  gzHeight = kzInitialHeight;

  rCode = InitalizeFrameBuffer();
  ThrowIf( rCode );

  goto cleanup;
error:
  fprintf( stderr, "Error in BuildWindow()\n" );
cleanup:

  return (xGWin_tErr_NoErr != rGWin);
}

int InitalizeFrameBuffer () {

  int rCode = 0;
  int nPixel = 0;

  if ( NULL != gFrameBuffer ) {
    free( gFrameBuffer );
  }

  gFrameBuffer = malloc( gzWidth * gzHeight * kBytesPerPixel );
  ThrowIf( NULL == gFrameBuffer );

  for ( nPixel = 0; nPixel < gzWidth * gzHeight * kBytesPerPixel; nPixel += kBytesPerPixel ) {
    gFrameBuffer[nPixel] = (GLubyte) 0;
    gFrameBuffer[nPixel+1] = (GLubyte) 0;
    gFrameBuffer[nPixel+2] = (GLubyte) 0;
    gFrameBuffer[nPixel+3] = (GLubyte) 1;
  }

  goto cleanup;
error:
  fprintf( stderr, "Error in InitalizeFrameBuffer()\n" );
  rCode = 1;
cleanup:

  return rCode;
}

int Brush ( int iX, int iY ) {

  int x = 0;
  int y = 0;
  int xMax = 0;
  int yMax = 0;
  int xMin = 0;
  int yMin = 0;

  xMin = MAX( iX - (gBrushRadius - 1), 0 );
  yMin = MAX( iY - (gBrushRadius - 1), 0 );
  xMax = MIN( iX + gBrushRadius, gzWidth - 1 );
  yMax = MIN( iY + gBrushRadius, gzHeight - 1 );

  for ( x = xMin; x <= xMax; x++ ) {
    for ( y = yMin; y <= yMax; y++ ) {
      EditFrameBuffer( x, y, random() % 255, random() % 255, random() % 255 );
    }
  }

  return 0;
}

int EditFrameBuffer( int iX, int iY, int iRed, int iGreen, int iBlue ) {

  int rCode  = 0;
  int nPixel = 0;

  if ( NULL == gFrameBuffer )
    goto error;
  if ( iX < 0 || iX >= gzWidth ||
       iY < 0 || iY >= gzHeight )
    goto error;

  nPixel = (iY * gzWidth * kBytesPerPixel) + (iX * kBytesPerPixel);
  gFrameBuffer[nPixel] = (GLubyte) iRed;
  gFrameBuffer[nPixel+1] = (GLubyte) iGreen;
  gFrameBuffer[nPixel+2] = (GLubyte) iBlue;
  gFrameBuffer[nPixel+3] = (GLubyte) 1;

  goto cleanup;
error:
  fprintf( stderr, "Error in EditFrameBuffer()\n" );
  rCode = 1;
cleanup:

  return rCode;
}

int Redraw () {

  glutPostRedisplay();

  return 0;
}

int main ( int argc, char** argv ) {

  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: test_glut.c,v 1.2 2010/01/04 19:05:40 nicks Exp $", "$Name: stable6 $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  BuildWindow();
  glutMainLoop();

  return 0;
}
