#include "string_fixed.h"
#include <stdexcept>
#include "ToglManager.h"
#include "Timer.h"

using namespace std;

map<ToglFrame::ID,ToglFrame*> ToglManager::mFrames;
WindowFrameFactory* ToglManager::mFactory = NULL;
InputState ToglManager::mState;


void
ToglManager::DrawCallback ( struct Togl* iTogl ) {

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->Draw();
  Togl_SwapBuffers( iTogl );
}

void
ToglManager::CreateCallback ( struct Togl* iTogl ) {

  if( NULL != mFactory ) {
    ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
    ToglFrame* frame = (ToglFrame*) mFactory->NewWindowFrame( id );
    mFrames[id] = frame;
  }
}

void
ToglManager::DestroyCallback ( struct Togl* iTogl ) {

  char* sIdent = Togl_Ident( iTogl ); // sometimes this is null?
  if( NULL != sIdent ) {
    ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
    ToglFrame* frame = mFrames[id];
    delete frame;
    mFrames[id] = NULL;
  }
}

void
ToglManager::ReshapeCallback ( struct Togl* iTogl ) {

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  int width = Togl_Width( iTogl );
  int height = Togl_Height( iTogl );
  frame->Reshape( width, height );

  // Post a redisplay if the frame wants one. 
  if( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }
}

void
ToglManager::TimerCallback ( struct Togl* iTogl ) {

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->Timer();

  // Post a redisplay if the frame wants one. 
  if( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }
}

int
ToglManager::MouseMotionCallback ( struct Togl* iTogl, 
				   int iArgc, char* iArgv[] ) {

  // widget MouseMotionCallback x y button
  if( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse is dragging.
  if( mState.Button() != 0 ) {
    mState.SetButtonDragEvent();
  }

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];

  int windowCoords[2];
  windowCoords[0] = atoi(iArgv[2]);
  windowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  frame->MouseMoved( windowCoords, mState );

  // Post a redisplay if the frame wants one. 
  if( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return TCL_OK;
}

int
ToglManager::MouseDownCallback ( struct Togl* iTogl, 
				 int iArgc, char* iArgv[] ) {

  // widget MouseDownCallback x y button
  if( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse down event.
  mState.SetButtonDownEvent( atoi(iArgv[4]) );

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];

  int windowCoords[2];
  windowCoords[0] = atoi(iArgv[2]);
  windowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  frame->MouseDown( windowCoords, mState );

  // Post a redisplay if the frame wants one. 
  if( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return TCL_OK;
}

int
ToglManager::MouseUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  // widget MouseUpCallback x y button
  if( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Mouse up event.
  mState.SetButtonUpEvent();
  
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];

  int windowCoords[2];
  windowCoords[0] = atoi(iArgv[2]);
  windowCoords[1] = YFlip(frame, atoi(iArgv[3]));
  frame->MouseUp( windowCoords, mState );

  // Clear the mouse events.
  mState.ClearEvents();

  // Post a redisplay if the frame wants one. 
  if( frame->WantRedisplay() ) {
    Togl_PostRedisplay( iTogl );
    frame->RedisplayPosted();
  }

  return TCL_OK;
}

int
ToglManager::KeyDownCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  // widget KeyDownCallback x y key
  if( iArgc != 5 ) {
    return TCL_ERROR;
  }

  // Look for modifiers. If it's shift (Shift_L or Shift_R), alt
  // (Alt_L or Alt_R), or ctrl (Control_L or Control_R) then just mark
  // our modifiers but don't pass the key along.
  string sKey = iArgv[4];
  if( sKey == "Shift_L" || sKey == "Shift_R" ) {
    mState.mbShiftKey = true;

  } else if( sKey == "Alt_L" || sKey == "Alt_R" ) {
    mState.mbAltKey = true;

  } else if( sKey == "Control_L" || sKey == "Control_R" ) {
    mState.mbControlKey = true;

  } else {
  
    // Record the key.
    mState.msKey = sKey;
    
    ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
    ToglFrame* frame = mFrames[id];

    int windowCoords[2];
    windowCoords[0] = atoi(iArgv[2]);
    windowCoords[1] = YFlip(frame, atoi(iArgv[3]));
    frame->KeyDown( windowCoords, mState );
    
    // Post a redisplay if the frame wants one. 
    if( frame->WantRedisplay() ) {
      Togl_PostRedisplay( iTogl );
      frame->RedisplayPosted();
    }
  }

  return TCL_OK;
}

int
ToglManager::KeyUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

   // widget KeyDownCallback x y key
 if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  // Look for modifiers. If it's shift (Shift_L or Shift_R), alt
  // (Alt_L or Alt_R), or ctrl (Control_L or Control_R) then just mark
  // our modifiers but don't pass the key along.
  string sKey = iArgv[4];
  if( sKey == "Shift_L" || sKey == "Shift_R" ) {
    mState.mbShiftKey = false;

  } else if( sKey == "Alt_L" || sKey == "Alt_R" ) {
    mState.mbAltKey = false;

  } else if( sKey == "Control_L" || sKey == "Control_R" ) {
    mState.mbControlKey = false;

  } else {
  
  
    // Record the key.
    mState.msKey = sKey;
    
    ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
    ToglFrame* frame = mFrames[id];

    int windowCoords[2];
    windowCoords[0] = atoi(iArgv[2]);
    windowCoords[1] = YFlip(frame, atoi(iArgv[3]));
    frame->KeyUp( windowCoords, mState );
    
    // Clear this in the keyboard state.
    mState.msKey = "";
    
    // Post a redisplay if the frame wants one.
    if( frame->WantRedisplay() ) {
      Togl_PostRedisplay( iTogl );
      frame->RedisplayPosted();
    }
  }

  return TCL_OK;
}

int
ToglManager::ExitCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if( iArgc != 2 ) {
    return TCL_ERROR;
  }

  // Just clear the modifiers.
  mState.mbShiftKey = false;
  mState.mbAltKey = false;
  mState.mbControlKey = false;
  mState.msKey = "";
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
  if( TCL_ERROR == rTogl ) {
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
  Togl_CreateCommand( "MouseMotionCallback", ToglManager::MouseMotionCallback);
  Togl_CreateCommand( "MouseDownCallback", ToglManager::MouseDownCallback );
  Togl_CreateCommand( "MouseUpCallback", ToglManager::MouseUpCallback );
  Togl_CreateCommand( "KeyDownCallback", ToglManager::KeyDownCallback );
  Togl_CreateCommand( "KeyUpCallback", ToglManager::KeyUpCallback );
  Togl_CreateCommand( "ExitCallback", ToglManager::ExitCallback );
}


void 
ToglFrame::Reshape( int iWidth, int iHeight ) {

  mWidth = iWidth;
  mHeight = iHeight;

  glViewport( 0, 0, mWidth, mHeight );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho( 0, mWidth, 0, mHeight, -1.0, 1.0 );
  glMatrixMode( GL_MODELVIEW );

  WindowFrame::Reshape( iWidth, iHeight );
}
