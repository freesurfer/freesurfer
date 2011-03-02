/**
 * @file  test_window_env.c
 * @brief test OpenGL stuff
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.10 $
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

#include <stdlib.h>
#include <stdio.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>
#include "version.h"

const char* Progname = "test_window_env";
Display *gDisplay = NULL;
int gScreen = 0;
Window gRoot = 0;

void InitDisplay ();
void FindVisual ();
void PrintAttributeList ( FILE* iOutput, 
                          char* isPrefix, 
                          int iaAttributes[] );

int main ( int argc, char** argv ) {

  int nargs;
  char *env;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: test_window_env.c,v 1.10 2011/03/02 00:04:40 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  env = getenv("DISPLAY");
  if ( NULL == env ) {
    printf("DISPLAY env var not set.  Skipping test.\n");
    exit (77);
  }

#if 0
  GLenum bSuccess;

  bSuccess = tkoInitWindow( "test" );
  printf( "tkoInitWindow() returned %s\n",
          (GL_TRUE == bSuccess ? "GL_TRUE" : "GL_FALSE") );
#endif

  InitDisplay();
  FindVisual();

  return 0;
}


void InitDisplay () {

  int errorBase;
  int eventBase;

  gDisplay = XOpenDisplay( 0 );
  if ( !gDisplay ) {
    printf( "XOpenDisplay() did not return a display.\n" );
    exit( 1 );
  }
  if ( !glXQueryExtension( gDisplay, &errorBase, &eventBase ) ) {
    printf( "glXQueryExtension returned false, no glx extensions.\n" );
    exit( 1 );
  }

  gScreen = DefaultScreen( gDisplay );
  gRoot = RootWindow( gDisplay, gScreen );

}

#define CLEAR_ATTRIBUTES() \
  nAttribute = 0; \
  aAttributes[0] = (int)None;
#define SET_ATTRIBUTE(a) \
  aAttributes[nAttribute++] = a; \
  aAttributes[nAttribute] = (int)None;
#define SET_ATTRIBUTE_VALUE(a,v) \
  aAttributes[nAttribute++] = a; \
  aAttributes[nAttribute++] = v; \
  aAttributes[nAttribute] = (int)None;
#define RUN_TEST() \
  printf( "\nTesting attribute set:\n" ); \
  PrintAttributeList( stdout, "\t", aAttributes ); \
  pVisualInfo = glXChooseVisual( gDisplay, gScreen, aAttributes ); \
  if( pVisualInfo ) { \
    printf( "\t-> Got a visual\n" ); \
    XFree( pVisualInfo ); \
  } else { \
    printf( "\tXX Didn't get a visual\n" ); \
    exit( 1 ); \
  }

void FindVisual () {

  XVisualInfo* pVisualInfo;
  int aAttributes[32];
  int nAttribute = 0;

  CLEAR_ATTRIBUTES();
  SET_ATTRIBUTE_VALUE( GLX_LEVEL, 0 );
  SET_ATTRIBUTE( GLX_RGBA );
  RUN_TEST();

  /* RGB screen. try the test at different phases. */
  SET_ATTRIBUTE_VALUE( GLX_RED_SIZE, 1 );
  SET_ATTRIBUTE_VALUE( GLX_GREEN_SIZE, 1 );
  SET_ATTRIBUTE_VALUE( GLX_BLUE_SIZE, 1 );
  RUN_TEST();

  SET_ATTRIBUTE_VALUE( GLX_ALPHA_SIZE, 1 );
  RUN_TEST();

  SET_ATTRIBUTE_VALUE( GLX_ACCUM_RED_SIZE, 1 );
  SET_ATTRIBUTE_VALUE( GLX_ACCUM_GREEN_SIZE, 1 );
  SET_ATTRIBUTE_VALUE( GLX_ACCUM_BLUE_SIZE, 1 );
  RUN_TEST();

  SET_ATTRIBUTE_VALUE( GLX_ACCUM_ALPHA_SIZE, 1 );
  RUN_TEST();

  /* try with and without double buffer. */
  SET_ATTRIBUTE( GLX_DOUBLEBUFFER );
  RUN_TEST();

  /* try with and without depth buffer. */
  SET_ATTRIBUTE_VALUE( GLX_DEPTH_SIZE, 1 );
  RUN_TEST();

  /* try with and without stencil buffer. */
  SET_ATTRIBUTE_VALUE( GLX_STENCIL_SIZE, 1 );
  RUN_TEST();

#if 0
  /* indexed screen */
  CLEAR_ATTRIBUTES();
  SET_ATTRIBUTE_VALUE( GLX_LEVEL, 0 );
  SET_ATTRIBUTE_VALUE( GLX_BUFFER_SIZE, 1 );
  RUN_TEST();

  /* try with and without double buffer. */
  SET_ATTRIBUTE( GLX_DOUBLEBUFFER );
  RUN_TEST();

  /* try with and without depth buffer. */
  SET_ATTRIBUTE_VALUE( GLX_DEPTH_SIZE, 1 );
  RUN_TEST();

  /* try with and without stencil buffer. */
  SET_ATTRIBUTE_VALUE( GLX_STENCIL_SIZE, 1 );
  RUN_TEST();
#endif
}

#define PRINT_ATTRIBUTE_VALUE_CASE(a) \
    case a: \
      fprintf( iOutput, "%s"#a" = %d\n", \
               isPrefix, iaAttributes[nAttribute++] ); \
      break;
#define PRINT_ATTRIBUTE_CASE(a) \
    case a: \
      fprintf( iOutput, "%s"#a"\n", isPrefix ); \
      break;


void PrintAttributeList ( FILE* iOutput, char* isPrefix, int iaAttributes[] ) {

  int nAttribute = 0;

  while ( (int)None != iaAttributes[nAttribute] ) {
    switch ( iaAttributes[nAttribute++] ) {
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_LEVEL);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_RED_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_GREEN_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_BLUE_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_ALPHA_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_ACCUM_RED_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_ACCUM_GREEN_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_ACCUM_BLUE_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_ACCUM_ALPHA_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_BUFFER_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_DEPTH_SIZE);
      PRINT_ATTRIBUTE_VALUE_CASE(GLX_STENCIL_SIZE);
      PRINT_ATTRIBUTE_CASE(GLX_DOUBLEBUFFER);
      PRINT_ATTRIBUTE_CASE(GLX_RGBA);
    default:
      fprintf( iOutput, "Unrecognized attribute: %d\n",
               iaAttributes[nAttribute] );
    }
  }
}
