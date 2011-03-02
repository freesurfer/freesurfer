/**
 * @file  ToglManager.cpp
 * @brief Main C++ object started by Tcl
 *
 * Togl is a bunch of code that implements a Tk GL
 * context. ToglManager is a class that interfaces with togl and acts
 * as the frontline for getting events from Tcl via the togl
 * object. When the ToglManager starts up, it registers a callback for
 * the creation of a togl object, and when the Tcl code creates one,
 * this code makes a WindowFrame, our base Window object. ToglManager
 * takes events and passes them to WindowFrames.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.30 $
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


#include "string_fixed.h"
#include <stdexcept>
#include "ToglManager.h"
#include "Timer.h"
#include "TclScubaKeyCombo.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcpy
#ifdef __cplusplus
}
#endif

using namespace std;

map<int,WindowFrame*> ToglManager::mFrames;
WindowFrameFactory* ToglManager::mFactory = NULL;
InputState ToglManager::mState;
int ToglManager::mCurrentWindowCoords[2];


void
ToglManager::DrawCallback ( struct Togl* iTogl ) {

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];
  frame->Draw();
  Togl_SwapBuffers( iTogl );
}

void
ToglManager::CreateCallback ( struct Togl* iTogl ) {

  if ( NULL != mFactory ) {
    int id = atoi( Togl_Ident( iTogl ));
    WindowFrame* frame = (WindowFrame*) mFactory->NewWindowFrame();
    mFrames[id] = frame;
  }
}

void
ToglManager::DestroyCallback ( struct Togl* iTogl ) {

  char* sIdent = Togl_Ident( iTogl ); // sometimes this is null?
  if ( NULL != sIdent ) {
    int id = atoi( Togl_Ident( iTogl ));
    WindowFrame* frame = mFrames[id];
    delete frame;
    mFrames[id] = NULL;
  }
}

void
ToglManager::ReshapeCallback ( struct Togl* iTogl ) {

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];
  int width = Togl_Width( iTogl );
  int height = Togl_Height( iTogl );
  frame->Reshape( width, height );

  // Post a redisplay if the frame wants one.
  if ( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }
}

void
ToglManager::TimerCallback ( struct Togl* iTogl ) {

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];
  frame->Timer();

  // Post a redisplay if the frame wants one.
  if ( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }
}

int
ToglManager::MouseMotionCallback ( struct Togl* iTogl,
                                   int iArgc, char* iArgv[] ) {

  int eTcl = TCL_OK;

  // widget MouseMotionCallback x y button
  if ( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse is dragging.
  if ( mState.Button() != 0 ) {
    mState.SetButtonDragEvent();
  }

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];

  mState.AddButtonDelta( atoi(iArgv[2]) - mCurrentWindowCoords[0],
                         YFlip(frame, atoi(iArgv[3])) - mCurrentWindowCoords[1] );

  mCurrentWindowCoords[0] = atoi(iArgv[2]);
  mCurrentWindowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  try {
    frame->MouseMoved( mCurrentWindowCoords, mState );
  } catch ( runtime_error& e) {
    char sError[1024];
    strcpy( sError, e.what() );
    Tcl_SetResult( Togl_Interp(iTogl), sError, TCL_VOLATILE );
    eTcl = TCL_ERROR;
  } catch (...) {
    cerr << "Uncaught exception in MouseMotionCallback" << endl;
    eTcl = TCL_ERROR;
  }

  // Post a redisplay if the frame wants one.
  if ( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return eTcl;
}

int
ToglManager::MouseDownCallback ( struct Togl* iTogl,
                                 int iArgc, char* iArgv[] ) {

  int eTcl = TCL_OK;

  // widget MouseDownCallback x y button
  if ( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse down event.
  mState.SetButtonDownEvent( atoi(iArgv[4]) );

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];

  mCurrentWindowCoords[0] = atoi(iArgv[2]);
  mCurrentWindowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  try {
    frame->MouseDown( mCurrentWindowCoords, mState );
  } catch ( runtime_error& e) {
    char sError[1024];
    strcpy( sError, e.what() );
    Tcl_SetResult( Togl_Interp(iTogl), sError, TCL_VOLATILE );
    eTcl = TCL_ERROR;
  } catch (...) {
    cerr << "Uncaught exception in MouseDownCallback" << endl;
    eTcl = TCL_ERROR;
  }

  // Post a redisplay if the frame wants one.
  if ( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return eTcl;
}

int
ToglManager::MouseUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  int eTcl = TCL_OK;

  // widget MouseUpCallback x y button
  if ( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse up event.
  mState.SetButtonUpEvent();

  int id = atoi( Togl_Ident( iTogl ));
  WindowFrame* frame = mFrames[id];

  mCurrentWindowCoords[0] = atoi(iArgv[2]);
  mCurrentWindowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  try {
    frame->MouseUp( mCurrentWindowCoords, mState );
  } catch ( runtime_error& e) {
    char sError[1024];
    strcpy( sError, e.what() );
    Tcl_SetResult( Togl_Interp(iTogl), sError, TCL_VOLATILE );
    eTcl = TCL_ERROR;
  } catch (...) {
    cerr << "Uncaught exception in MouseUpCallback" << endl;
    eTcl = TCL_ERROR;
  }

  // Clear the mouse events.
  mState.ClearEvents();

  // Post a redisplay if the frame wants one.
  if ( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return eTcl;
}

int
ToglManager::KeyDownCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  int eTcl = TCL_OK;

  // widget KeyDownCallback x y key
  if ( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Look for modifiers. If it's shift (Shift_L or Shift_R), alt
  // (Alt_L or Alt_R), or ctrl (Control_L or Control_R) then just mark
  // our modifiers but don't pass the key along.
  string sKey = iArgv[4];
  if ( sKey == "Shift_L" || sKey == "Shift_R" ) {
    mState.SetShiftKey( true );

  } else if ( sKey == "Alt_L" || sKey == "Alt_R" ) {
    mState.SetAltKey( true );

  } else if ( sKey == "Control_L" || sKey == "Control_R" ) {
    mState.SetControlKey( true );

  } else {

    // If shift, change key to lowercase.
    if ( mState.IsShiftKeyDown() ) {
      sKey = tolower(sKey[0]);
    }

    // Record the key.
    ScubaKeyCombo* key = ScubaKeyCombo::MakeKeyCombo();
    key->SetFromString( sKey );
    mState.CopyKey( *key );

    int id = atoi( Togl_Ident( iTogl ));
    WindowFrame* frame = mFrames[id];

    mCurrentWindowCoords[0] = atoi(iArgv[2]);
    mCurrentWindowCoords[1] = YFlip(frame, atoi(iArgv[3]));
    try {
      frame->KeyDown( mCurrentWindowCoords, mState );
    } catch ( runtime_error& e) {
      char sError[1024];
      strcpy( sError, e.what() );
      Tcl_SetResult( Togl_Interp(iTogl), sError, TCL_VOLATILE );
      eTcl = TCL_ERROR;
    } catch (...) {
      cerr << "Uncaught exception in KeyDownCallback" << endl;
      eTcl = TCL_ERROR;
    }

    // Post a redisplay if the frame wants one.
    if ( frame->WantRedisplay() ) {
      Togl_PostRedisplay( iTogl );
      frame->RedisplayPosted();
    }
  }

  return eTcl;
}

int
ToglManager::KeyUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  int eTcl = TCL_OK;

  // widget KeyDownCallback x y key
  if ( iArgc != 5 ) {
    return TCL_ERROR;
  }
  // Look for modifiers. If it's shift (Shift_L or Shift_R), alt
  // (Alt_L or Alt_R), or ctrl (Control_L or Control_R) then just mark
  // our modifiers but don't pass the key along.
  string sKey = iArgv[4];
  if ( sKey == "Shift_L" || sKey == "Shift_R" ) {
    mState.SetShiftKey( false );

  } else if ( sKey == "Alt_L" || sKey == "Alt_R" ) {
    mState.SetAltKey( false );

  } else if ( sKey == "Control_L" || sKey == "Control_R" ) {
    mState.SetControlKey( false );

  } else {

    // Record the key.
    ScubaKeyCombo* key = ScubaKeyCombo::MakeKeyCombo();
    key->SetFromString( sKey );
    mState.CopyKey( *key );

    int id = atoi( Togl_Ident( iTogl ));
    WindowFrame* frame = mFrames[id];

    mCurrentWindowCoords[0] = atoi(iArgv[2]);
    mCurrentWindowCoords[1] = YFlip(frame, atoi(iArgv[3]));
    try {
      frame->KeyUp( mCurrentWindowCoords, mState );
    } catch ( runtime_error& e) {
      char sError[1024];
      strcpy( sError, e.what() );
      Tcl_SetResult( Togl_Interp(iTogl), sError, TCL_VOLATILE );
      eTcl = TCL_ERROR;
    } catch (...) {
      cerr << "Uncaught exception in KeyUpCallback" << endl;
      eTcl = TCL_ERROR;
    }

    // Post a redisplay if the frame wants one.
    if ( frame->WantRedisplay() ) {
      Togl_PostRedisplay( iTogl );
      frame->RedisplayPosted();
    }
  }

  return eTcl;
}

int
ToglManager::ExitCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if ( iArgc != 2 ) {
    return TCL_ERROR;
  }

  // Just clear the modifiers.
  mState.SetShiftKey( false );
  mState.SetAltKey( false );
  mState.SetControlKey( false );
  mState.ClearEvents();

  return TCL_OK;
}

int
ToglManager::EnterCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if ( iArgc != 2 ) {
    return TCL_ERROR;
  }

  // Just clear the modifiers.
  mState.SetShiftKey( false );
  mState.SetAltKey( false );
  mState.SetControlKey( false );
  mState.ClearEvents();

  return TCL_OK;
}


ToglManager&
ToglManager::GetManager() {
  static ToglManager sManager;
  return sManager;
}

void
ToglManager::InitializeTogl ( Tcl_Interp* iInterp ) {

  // Initialize the Togl module.
  int rTogl = Togl_Init( iInterp );
  if ( TCL_ERROR == rTogl ) {
    throw runtime_error( "Couldn't initialize Togl" );
  }

  // Register our Togl callbacks.
  Togl_CreateFunc( ToglManager::CreateCallback );
  Togl_DestroyFunc( ToglManager::DestroyCallback );
  Togl_DisplayFunc( ToglManager::DrawCallback );
  Togl_ReshapeFunc( ToglManager::ReshapeCallback );
  Togl_TimerFunc( ToglManager::TimerCallback );

  // Create our Tcl commands that will be bound to the Togl object as
  // callbacks.
  Togl_CreateCommand( (char*)"MouseMotionCallback", ToglManager::MouseMotionCallback);
  Togl_CreateCommand( (char*)"MouseDownCallback", ToglManager::MouseDownCallback );
  Togl_CreateCommand( (char*)"MouseUpCallback", ToglManager::MouseUpCallback );
  Togl_CreateCommand( (char*)"KeyDownCallback", ToglManager::KeyDownCallback );
  Togl_CreateCommand( (char*)"KeyUpCallback", ToglManager::KeyUpCallback );
  Togl_CreateCommand( (char*)"ExitCallback", ToglManager::ExitCallback );
  Togl_CreateCommand( (char*)"EnterCallback", ToglManager::EnterCallback );
}

