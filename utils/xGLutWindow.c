/**
 * @file  xGLutWindow.c
 * @brief general purpose utils
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.10 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifdef HAVE_GLUT_LIBS

#include <stdlib.h>
#include "xGLutWindow.h"
#include "xDebug.h"

static xGLutWindowRef mLookupList [xGWin_knMaxNumWindows];

static char *ksErrorStrings [xGWin_knNumErrorCodes] =
{
  "No error.",
  "Invalid window ptr (probably NULL).",
  "Invalid window signature.",
  "Memory allocation failed.",
  "Invalid input paramter.",
  "Invalid error code."
};

xGWin_tErr xGWin_New ( xGLutWindowRef* oppWindow,
                       int             inWidth,
                       int             inHeight,
                       char*           isTitle )
{

  xGWin_tErr     eResult = xGWin_tErr_NoErr;
  xGLutWindowRef this    = NULL;

  // check our params.
  if ( NULL == isTitle
       || inWidth < 0
       || inHeight < 0 )
  {
    eResult = xGWin_tErr_InvalidParameter;
    goto error;
  }

  // allocate the window.
  this = (xGLutWindowRef) malloc ( sizeof ( xGLutWindow ) );
  if ( NULL == this )
  {
    eResult = xGWin_tErr_AllocationFailed;
    goto error;
  }

  // set signature.
  this->mSignature = xGWin_kSignature;

  // set the window size.
  glutInitWindowSize ( inWidth, inHeight );

  // create an glut window and save its id
  this->mnGLutWindowID = glutCreateWindow ( isTitle );

  // set the glut window's handlers to our handlers.
  glutDisplayFunc  ( xGWin_GLutDrawCallback );
  glutKeyboardFunc ( xGWin_GLutKeyboardCallback );
  glutSpecialFunc  ( xGWin_GLutSpecialCallback );
  glutMouseFunc    ( xGWin_GLutMouseCallback );
  glutMotionFunc   ( xGWin_GLutMotionCallback );
  glutReshapeFunc  ( xGWin_GLutResizeCallback );

  // set func ptrs and data to nil.
  this->mpHandlerFunc     = NULL;
  this->mpHandlerFuncData = NULL;

  // add this window to the window list.
  xGWin_AddWindowIDToLookupList ( this->mnGLutWindowID, this );

  // return it.
  *oppWindow = this;

  goto cleanup;

error:

  // delete window if it was allocated.
  if ( NULL != this )
  {
    free ( this );
  }

  DebugPrint( ("Error %d in xGWin_New: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;
}

xGWin_tErr xGWin_Delete ( xGLutWindowRef* ioppWindow )
{

  xGWin_tErr     eResult = xGWin_tErr_NoErr;
  xGLutWindowRef this    = NULL;

  // grab window.
  this = *ioppWindow;

  // verify it.
  eResult = xGWin_Verify ( this );
  if ( xGWin_tErr_NoErr != eResult )
  {
    goto error;
  }

  // destroy the glut window.
  glutDestroyWindow ( this->mnGLutWindowID );

  // trash signature.
  this->mSignature = 0x1;

  // delete window.
  free ( this );

  // return nill.
  *ioppWindow = NULL;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xGWin_Delete: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;
}

xGWin_tErr xGWin_SetWindowTitle ( xGLutWindowRef this,
                                  char*          isTitle )
{

  xGWin_tErr eResult = xGWin_tErr_NoErr;

  // verify it.
  eResult = xGWin_Verify ( this );
  if ( xGWin_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* set window title */
  glutSetWindowTitle( isTitle );

  goto cleanup;

error:

  DebugPrint( ("Error %d in xGWin_SetWindowTitle: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;

}

xGWin_tErr xGWin_SetEventHandlerFunc ( xGLutWindowRef          this,
                                       xGWin_tEventHandlerFunc ipFunc,
                                       void*                   ipData )
{

  xGWin_tErr eResult = xGWin_tErr_NoErr;

  // verify it.
  eResult = xGWin_Verify ( this );
  if ( xGWin_tErr_NoErr != eResult )
  {
    goto error;
  }

  // set func and data.
  this->mpHandlerFunc     = ipFunc;
  this->mpHandlerFuncData = ipData;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xGWin_SetEventHandlerFunc: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;
}

xGWin_tErr xGWin_ActivateIdleEvents ( xGLutWindowRef  this )
{

  xGWin_tErr eResult       = xGWin_tErr_NoErr;
  int        nSaveWindowID = 0;

  // verify it.
  eResult = xGWin_Verify ( this );
  if ( xGWin_tErr_NoErr != eResult )
  {
    goto error;
  }

  // save the current window and set the window to this window.
  nSaveWindowID = glutGetWindow ();
  glutSetWindow ( this->mnGLutWindowID );

  // register our idle callback.
  glutIdleFunc ( xGWin_GLutIdleCallback );

  // restore the window.
  glutSetWindow ( nSaveWindowID );

  goto cleanup;

error:

  DebugPrint( ("Error %d in xGWin_ActivateIdleEvents: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;
}

xGWin_tErr xGWin_ActivatePassiveMotionEvents ( xGLutWindowRef  this )
{

  xGWin_tErr eResult       = xGWin_tErr_NoErr;
  int        nSaveWindowID = 0;

  // verify it.
  eResult = xGWin_Verify ( this );
  if ( xGWin_tErr_NoErr != eResult )
  {
    goto error;
  }

  // save the current window and set the window to this window.
  nSaveWindowID = glutGetWindow ();
  glutSetWindow ( this->mnGLutWindowID );

  // register our passive motion callback.
  glutPassiveMotionFunc ( xGWin_GLutPassiveMotionCallback );

  // restore the window.
  glutSetWindow ( nSaveWindowID );

  goto cleanup;

error:

  DebugPrint( ("Error %d in xGWin_ActivatePassiveMotionEvents: %s",
               eResult, xGWin_GetErrorString ( eResult ) ) );

cleanup:

  return eResult;
}

xGWin_tErr xGWin_Verify ( xGLutWindowRef this )
{

  xGWin_tErr eResult = xGWin_tErr_NoErr;

  // check ptr.
  if ( NULL == this )
  {
    eResult = xGWin_tErr_InvalidWindowPtr;
    goto error;
  }

  // check signature.
  if ( xGWin_kSignature != this->mSignature )
  {
    eResult = xGWin_tErr_InvalidWindowSignature;
    goto error;
  }

  goto cleanup;

error:
cleanup:

  return eResult;
}

char* xGWin_GetErrorString ( xGWin_tErr ieCode )
{

  if ( ieCode < 0 || ieCode >= xGWin_knNumErrorCodes )
  {
    ieCode = xGWin_tErr_InvalidErrorCode;
  }

  return ksErrorStrings [ieCode];
}

void xGWin_AddWindowIDToLookupList ( int            inWindowID,
                                     xGLutWindowRef ipWindow )
{

  /* this used to be a list of number/ptr nodes, but the glut specification
     says that numbers always start at 1 and increase sequentially, so i
     guess it's okay to use an arrray if you don't mind a fixed number
     of max windows. */

  // check the bounds.
  if ( inWindowID < 0
       || inWindowID >= xGWin_knMaxNumWindows )
  {
    DebugPrint( ("xGWin_AddWindowIDToLookupList: invalid ID: %d\n",
                 inWindowID ) );
    return;
  }

  // save this ptr.
  mLookupList[inWindowID] = ipWindow;
}

void xGWin_GetWindowFromID ( int             inWindowID,
                             xGLutWindowRef* oppWindow )
{

  // check the bounds.
  if ( inWindowID < 0
       || inWindowID >= xGWin_knMaxNumWindows )
  {
    DebugPrint( ("xGWin_AddWindowIDToLookupList: invalid ID: %d\n",
                 inWindowID ) );
    return;
  }

  // return this window ptr.
  *oppWindow = mLookupList[inWindowID];
}

void xGWin_RemoveWindowIDFromLookupList ( int inWindowID )
{

  // check the bounds.
  if ( inWindowID < 0
       || inWindowID >= xGWin_knMaxNumWindows )
  {
    DebugPrint( ("xGWin_AddWindowIDToLookupList: invalid ID: %d\n",
                 inWindowID ) );
    return;
  }

  // set the ptr to null.
  mLookupList[inWindowID] = NULL;
}

void xGWin_PassEventToCurrentWindow ( xGWin_tEventRef ipEvent )
{

  int            nWindowID = 0;
  xGLutWindowRef pWindow   = NULL;

  // get current window id.
  nWindowID = glutGetWindow ();

  // get a window ptr.
  xGWin_GetWindowFromID ( nWindowID, &pWindow );

  if ( NULL != pWindow )
  {

    // call the window's event handler with the window handler data
    // and this event.
    pWindow->mpHandlerFunc ( pWindow->mpHandlerFuncData, ipEvent );

  }
  else
  {

    DebugPrint( ("xGWin_PassEventToCurrentWindow: Couldn't find current"
                 "window, id is %d", nWindowID ) );
  }
}

void xGWin_PassEventToAllWindows ( xGWin_tEventRef ipEvent )
{

  int            nSaveWindowID = 0;
  int            nWindowID     = 0;
  xGLutWindowRef pWindow       = NULL;

  // save the current window.
  nSaveWindowID = glutGetWindow ();

  // for all windows...
  for ( nWindowID = 0; nWindowID < xGWin_knMaxNumWindows; nWindowID++ )
  {

    // try to get a window ptr.
    xGWin_GetWindowFromID ( nWindowID, &pWindow );

    // if we got one...
    if ( NULL != pWindow )
    {

      // set the window
      glutSetWindow ( nWindowID );

      // call the window's event handler with the window handler data
      // and this event.
      pWindow->mpHandlerFunc ( pWindow->mpHandlerFuncData, ipEvent );
    }
  }

  // restore the window.
  glutSetWindow ( nSaveWindowID );
}

void xGWin_GLutKeyboardCallback ( unsigned char icKey,
                                  int           inX,
                                  int           inY )
{

  xGWin_tEventRef pEvent        = NULL;
  int             nState        = 0;
  unsigned char   ucModifiedKey = 0;

  /* in all these funcs, we just create an event, set the relevant fields,
     and pass it to the event sender */
  xGWin_NewEvent ( &pEvent );

  pEvent->mType      = xGWin_tEventType_KeyDown;
  pEvent->mKey       = icKey;
  pEvent->mWhere.mnX = inX;
  pEvent->mWhere.mnY = inY;

  /* finds the state of modifer keys and set booleans appropriatly */
  nState = glutGetModifiers();
  if ( nState & GLUT_ACTIVE_SHIFT )
  {
    pEvent->mbShiftKey = TRUE;
  }
  if ( nState & GLUT_ACTIVE_CTRL )
  {
    pEvent->mbCtrlKey = TRUE;
  }
  if ( nState & GLUT_ACTIVE_ALT )
  {
    pEvent->mbAltKey = TRUE;
  }

  /* if control is down, the char we got was not an actual char, it was
     a special character. so we'll add the value of a minus the value of
     ctrl-a. */
  if ( pEvent->mbCtrlKey )
  {
    ucModifiedKey = (unsigned char)((int)icKey + ((int)'a' - xGWin_knCtrlA));
    if ( (ucModifiedKey >= 'a' && ucModifiedKey <= 'z')
         || (ucModifiedKey >= 'A' && ucModifiedKey <= 'Z') )
    {
      pEvent->mKey = ucModifiedKey;
    }
    else
    {
      /* we also have special codes for the control-number keys. check
      those. */
      switch ( icKey )
      {
      case xGWin_knCtrl0:
        pEvent->mKey = '0';
        break;
      case xGWin_knCtrl1:
        pEvent->mKey = '1';
        break;
      case xGWin_knCtrl2:
        pEvent->mKey = '2';
        break;
      case xGWin_knCtrl3:
        pEvent->mKey = '3';
        break;
      case xGWin_knCtrl4:
        pEvent->mKey = '4';
        break;
      case xGWin_knCtrl5:
        pEvent->mKey = '5';
        break;
      case xGWin_knCtrl6:
        pEvent->mKey = '6';
        break;
      case xGWin_knCtrl7:
        pEvent->mKey = '7';
        break;
      case xGWin_knCtrl8:
        pEvent->mKey = '8';
        break;
      case xGWin_knCtrl9:
        pEvent->mKey = '9';
        break;
      }
    }
  }

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutSpecialCallback ( int inKey,
                                 int inX,
                                 int inY )
{

  xGWin_tEventRef pEvent = NULL;
  int             nState = 0;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType      = xGWin_tEventType_KeyDown;
  pEvent->mWhere.mnX = inX;
  pEvent->mWhere.mnY = inY;

  /* choose the proper key. */
  switch ( inKey )
  {
  case GLUT_KEY_LEFT:
    pEvent->mKey = xGWin_tKey_LeftArrow;
    break;
  case GLUT_KEY_UP:
    pEvent->mKey = xGWin_tKey_UpArrow;
    break;
  case GLUT_KEY_RIGHT:
    pEvent->mKey = xGWin_tKey_RightArrow;
    break;
  case GLUT_KEY_DOWN:
    pEvent->mKey = xGWin_tKey_DownArrow;
    break;
  case GLUT_KEY_PAGE_UP:
    pEvent->mKey = xGWin_tKey_PageUp;
    break;
  case GLUT_KEY_PAGE_DOWN:
    pEvent->mKey = xGWin_tKey_PageDown;
    break;
  case GLUT_KEY_HOME:
    pEvent->mKey = xGWin_tKey_Home;
    break;
  case GLUT_KEY_END:
    pEvent->mKey = xGWin_tKey_End;
    break;
  case GLUT_KEY_INSERT:
    pEvent->mKey = xGWin_tKey_Insert;
    break;
  }

  /* finds the state of modifer keys and set booleans appropriatly */
  nState = glutGetModifiers();
  if ( nState & GLUT_ACTIVE_SHIFT )
  {
    pEvent->mbShiftKey = TRUE;
  }
  if ( nState & GLUT_ACTIVE_CTRL )
  {
    pEvent->mbCtrlKey = TRUE;
  }
  if ( nState & GLUT_ACTIVE_ALT )
  {
    pEvent->mbAltKey = TRUE;
  }

  xGWin_PassEventToCurrentWindow ( pEvent );
}

/* we use these variables to save the button pressed and modifers when
   the mouse is first pressed down. we use these values when generating
   mouse-moved events, because glut doesn't let us glutGetModifiers on
   a mouse-moved callback. */
static xGWin_tButton mButton = 0;
static tBoolean mbShiftKey = FALSE;
static tBoolean mbCtrlKey  = FALSE;
static tBoolean mbAltKey   = FALSE;

void xGWin_GLutMouseCallback ( int inButton,
                               int inState,
                               int inX,
                               int inY )
{

  xGWin_tEventRef pEvent = NULL;
  int             nState = 0;

  xGWin_NewEvent ( &pEvent );

  /* find out if it was up or down. */
  if ( GLUT_DOWN == inState )
  {
    pEvent->mType = xGWin_tEventType_MouseDown;
  }
  if ( GLUT_UP   == inState )
  {
    pEvent->mType = xGWin_tEventType_MouseUp;
  }

  /* find out which button */
  if ( GLUT_LEFT_BUTTON   == inButton )
  {
    pEvent->mButton = 1;
  }
  if ( GLUT_MIDDLE_BUTTON == inButton )
  {
    pEvent->mButton = 2;
  }
  if ( GLUT_RIGHT_BUTTON  == inButton )
  {
    pEvent->mButton  = 3;
  }
  mButton = pEvent->mButton;

  /* set location */
  pEvent->mWhere.mnX = inX;
  pEvent->mWhere.mnY = inY;

  /* finds the state of modifer keys and set booleans appropriatly */
  nState = glutGetModifiers();
  if ( nState & GLUT_ACTIVE_SHIFT )
  {
    pEvent->mbShiftKey = TRUE;
    mbShiftKey         = TRUE;
  }
  else
  {
    mbShiftKey         = FALSE;
  }
  if ( nState & GLUT_ACTIVE_CTRL )
  {
    pEvent->mbCtrlKey = TRUE;
    mbCtrlKey         = TRUE;
  }
  else
  {
    mbCtrlKey         = FALSE;
  }
  if ( nState & GLUT_ACTIVE_ALT )
  {
    pEvent->mbAltKey = TRUE;
    mbAltKey         = TRUE;
  }
  else
  {
    mbAltKey         = FALSE;
  }

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutMotionCallback ( int inX,
                                int inY )
{

  xGWin_tEventRef pEvent = NULL;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType = xGWin_tEventType_MouseMoved;

  /* get button from last mouse down.*/
  pEvent->mButton = mButton;

  /* set location */
  pEvent->mWhere.mnX = inX;
  pEvent->mWhere.mnY = inY;

  /* get modifiers from last mouse down. */
  pEvent->mbShiftKey = mbShiftKey;
  pEvent->mbCtrlKey  = mbCtrlKey;
  pEvent->mbAltKey   = mbAltKey;

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutPassiveMotionCallback ( int inX,
                                       int inY )
{

  xGWin_tEventRef pEvent = NULL;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType = xGWin_tEventType_MouseMoved;

  /* no mouse down.*/
  pEvent->mButton = 0;

  /* set location */
  pEvent->mWhere.mnX = inX;
  pEvent->mWhere.mnY = inY;

  /* no modifiers */
  pEvent->mbShiftKey = FALSE;
  pEvent->mbCtrlKey  = FALSE;
  pEvent->mbAltKey   = FALSE;

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutResizeCallback ( int inWidth,
                                int inHeight )
{

  xGWin_tEventRef pEvent;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType      = xGWin_tEventType_Resize;
  pEvent->mWhere.mnX = inWidth;
  pEvent->mWhere.mnY = inHeight;

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutDrawCallback ()
{

  xGWin_tEventRef pEvent;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType = xGWin_tEventType_Draw;

  xGWin_PassEventToCurrentWindow ( pEvent );
}

void xGWin_GLutIdleCallback ()
{

  xGWin_tEventRef pEvent;

  xGWin_NewEvent ( &pEvent );

  pEvent->mType = xGWin_tEventType_Idle;

  xGWin_PassEventToCurrentWindow ( pEvent );
}


void xGWin_RedrawAllWindows ()
{

  int            nSaveWindowID = 0;
  int            nWindowID     = 0;
  xGLutWindowRef pWindow       = NULL;

  // save the current window.
  nSaveWindowID = glutGetWindow ();

  // for all windows...
  for ( nWindowID = 0; nWindowID < xGWin_knMaxNumWindows; nWindowID++ )
  {

    // try to get a window ptr.
    xGWin_GetWindowFromID ( nWindowID, &pWindow );

    // if we got one...
    if ( NULL != pWindow )
    {

      // set the window
      glutSetWindow ( nWindowID );

      // post a redisplay for it.
      glutPostRedisplay ();
    }
  }

  // restore the window.
  glutSetWindow ( nSaveWindowID );

}

void xGWin_NewEvent ( xGWin_tEventRef* oppEvent )
{

  xGWin_tEventRef this   = NULL;

  // allocate event storage.
  this = (xGWin_tEventRef) malloc ( sizeof ( xGWin_tEvent ) );

  // set all to defaults.
  this->mType       = xGWin_tEventType_NoEvent;
  this->mWhere.mnX  = -1;
  this->mWhere.mnY  = -1;
  this->mButton     = -1;
  this->mKey        = (xGWin_tKey)' ';
  this->mbCtrlKey   = FALSE;
  this->mbAltKey    = FALSE;
  this->mbShiftKey  = FALSE;

  *oppEvent = this;
}

void xGWin_DeleteEvent ( xGWin_tEventRef* ioppEvent )
{

  xGWin_tEventRef this = NULL;

  this = *ioppEvent;

  // just delete the event.
  free ( this );

  *ioppEvent = NULL;
}

void xGWin_DebugPrintEvent ( xGWin_tEventRef this )
{

  DebugPrint( ("XGWin_tEvent:\n" ) );

  switch ( this->mType )
  {
  case xGWin_tEventType_Draw:
    DebugPrint( ("\tType: Draw\n" ) );
    break;
  case xGWin_tEventType_KeyDown:
    DebugPrint( ("\tType: KeyDown\n" ) );
    break;
  case xGWin_tEventType_MouseDown:
    DebugPrint( ("\tType: MouseDown\n" ) );
    break;
  case xGWin_tEventType_MouseUp:
    DebugPrint( ("\tType: MouseUp\n" ) );
    break;
  case xGWin_tEventType_MouseMoved:
    DebugPrint( ("\tType: MouseMoved\n" ) );
    break;
  case xGWin_tEventType_Resize:
    DebugPrint( ("\tType: Resize\n" ) );
    break;
  default:
    break;
  }

  DebugPrint( ("\tWhere: %d %d\n", this->mWhere.mnX, this->mWhere.mnY
              ) );
  DebugPrint( ("\tButton: %d\n", this->mButton ) );
  DebugPrint( ("\tKey: %c (%d)\n", this->mKey, (int)this->mKey ) );
  DebugPrint( ("\tModifiers: ctrl %d alt %d shift %d\n",
               (int)(this->mbCtrlKey), (int)(this->mbAltKey), (int)(this->mbShiftKey)
              ) );
  DebugPrint( ("\n" ) );
}

#endif // #ifdef HAVE_GLUT_LIBS
