/**
 * @file  xGLutWindow.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.11 $
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


#ifndef xGLutWindow_h
#define xGLutWindow_h

#include <glut.h>
#include "xTypes.h"

typedef enum
{

  xGWin_tErr_NoErr = 0,
  xGWin_tErr_InvalidWindowPtr,
  xGWin_tErr_InvalidWindowSignature,
  xGWin_tErr_AllocationFailed,
  xGWin_tErr_InvalidParameter,
  xGWin_tErr_InvalidErrorCode,
  xGWin_knNumErrorCodes

} xGWin_tErr;

/* the event class */
typedef enum
{

  xGWin_tEventType_NoEvent = 0,
  xGWin_tEventType_KeyDown,
  xGWin_tEventType_MouseDown,
  xGWin_tEventType_MouseUp,
  xGWin_tEventType_MouseMoved,
  xGWin_tEventType_Resize,
  xGWin_tEventType_Draw,
  xGWin_tEventType_Idle,
  xGWin_NumEventTypes

} xGWin_tEventType;

typedef int           xGWin_tButton;
typedef unsigned char xGWin_tKey;

/* special keys */
enum
{

  xGWin_tKey_Tab = 9,  /* 9 myst remain hardcoded */
  xGWin_tKey_UpArrow = 10,
  xGWin_tKey_DownArrow,
  xGWin_tKey_LeftArrow,
  xGWin_tKey_RightArrow,
  xGWin_tKey_PageUp,
  xGWin_tKey_PageDown,
  xGWin_tKey_End,
  xGWin_tKey_Insert,
  xGWin_tKey_Home
};

#define xGWin_knCtrlA 1
#define xGWin_knCtrl0 48
#define xGWin_knCtrl1 49
#define xGWin_knCtrl2 0
#define xGWin_knCtrl3 27
#define xGWin_knCtrl4 28
#define xGWin_knCtrl5 29
#define xGWin_knCtrl6 30
#define xGWin_knCtrl7 31
#define xGWin_knCtrl8 127
#define xGWin_knCtrl9 57

typedef struct
{

  xGWin_tEventType mType;      // type of event that occured.
  xPoint2n         mWhere;     // where mouse was (or what window resized to)
  xGWin_tButton    mButton;    // which button (1,2,3)
  xGWin_tKey       mKey;       // what key (a letter or a special enum)
  tBoolean         mbCtrlKey;  // ctrl key was down?
  tBoolean         mbAltKey;   // alt key was down?
  tBoolean         mbShiftKey; // shift key was down?

}
xGWin_tEvent, *xGWin_tEventRef;

/* create and destroy events */
void xGWin_NewEvent    ( xGWin_tEventRef* oppEvent );
void xGWin_DeleteEvent ( xGWin_tEventRef* ioppEvent );

/* print event data */
void xGWin_DebugPrintEvent ( xGWin_tEventRef this );

/* event handle function signature */
typedef void(*xGWin_tEventHandlerFunc) ( void* ipData, xGWin_tEventRef );
typedef void(*xGWin_tIdleFunc)         ( void );

#define xGWin_kSignature 0x1234ABCD

/* window class */
typedef struct
{

  tSignature              mSignature;
  int                     mnGLutWindowID;    // openGL window id
  xGWin_tEventHandlerFunc mpHandlerFunc;     // subclass event handler
  void*                   mpHandlerFuncData; // ptr to subclass data
}
xGLutWindow, *xGLutWindowRef;

/* creates and deletes windows. */
xGWin_tErr xGWin_New    ( xGLutWindowRef* oppWindow,
                          int             inWidth,
                          int             inHeight,
                          char*           isTitle );
xGWin_tErr xGWin_Delete ( xGLutWindowRef* ioppWindow );

/* sets the window title */
xGWin_tErr xGWin_SetWindowTitle ( xGLutWindowRef ipWindow,
                                  char*          isTitle );

/* set the event handler for this window. when an event is received,
   calls ipFunc( ipData ). */
xGWin_tErr xGWin_SetEventHandlerFunc ( xGLutWindowRef          this,
                                       xGWin_tEventHandlerFunc ipFunc,
                                       void*                   ipData );

/* activate idle events for this window */
xGWin_tErr xGWin_ActivateIdleEvents ( xGLutWindowRef  this );

/* activate passive motion events for this window */
xGWin_tErr xGWin_ActivatePassiveMotionEvents ( xGLutWindowRef  this );


xGWin_tErr xGWin_Verify ( xGLutWindowRef this );

char* xGWin_GetErrorString ( xGWin_tErr );


/* since glut doesn't provide us a way to attach a ptr to a glut window,
   we have to make a lookup list of gl id numbers and window ptrs, so when
   we get an event, we figure out the current window, look up its ptr,
   and pass the event to the right window. */

#define xGWin_knMaxNumWindows 200

/* this adds an id and window to the list. it will also init the list if
   it hasn't been done yet. */
void xGWin_AddWindowIDToLookupList ( int            inWindowID,
                                     xGLutWindowRef ipWindow );

/* this gets a window ptr from an id. returns null if it can't find the
   the id. */
void xGWin_GetWindowFromID ( int             inWindowID,
                             xGLutWindowRef* oppWindow );

/* removes a window from the list. does nothing if it can't find the id. */
void xGWin_RemoveWindowIDFromLookupList ( int inWindowID );

/* passes an event to a window. */
void xGWin_PassEventToCurrentWindow ( xGWin_tEventRef ipEvent );
void xGWin_PassEventToAllWindows    ( xGWin_tEventRef ipEvent );

/* these are the callbacks we register with glut. since glut's event system
   is pretty barebones, and each callback has a different signature, we
   stuff each event into a custom structure and pass it to a single event
   handler on the window side. */
void xGWin_GLutKeyboardCallback        ( unsigned char icKey,
    int           inX,
    int           inY );
void xGWin_GLutSpecialCallback         ( int           inKey,
    int           inX,
    int           inY );
void xGWin_GLutMouseCallback           ( int           inButton,
    int           inState,
    int           inX,
    int           inY );
void xGWin_GLutMotionCallback          ( int           inX,
    int           inY );
void xGWin_GLutPassiveMotionCallback   ( int           inX,
    int           inY );
void xGWin_GLutResizeCallback          ( int           inWidth,
    int           inHeight );
void xGWin_GLutDrawCallback            ();
void xGWin_GLutIdleCallback            ();

/* posts a redisplay for all windows. note that the redisplays go into the
   normal glut event queue. */
void xGWin_RedrawAllWindows ();

#endif
