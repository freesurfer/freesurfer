#include <stdexcept>
#include "ToglManager.h"

using namespace std;

map<ToglFrame::ID,ToglFrame*> ToglManager::mFrames;
ToglFrameFactory* ToglManager::mFactory = NULL;
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
    ToglFrame* frame = mFactory->NewToglFrame( id );
    mFrames[id] = frame;
  }
}

void
ToglManager::DestroyCallback ( struct Togl* iTogl ) {

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  delete frame;
  mFrames[id] = NULL;
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
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseMoved( atoi(iArgv[2]), YFlip(frame, atoi(iArgv[3])), mState );

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

  // Record this in the keyboard state.
  mState.mButton = atoi(iArgv[4]);

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseDown( atoi(iArgv[2]),  YFlip(frame, atoi(iArgv[3])), mState );

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

  // Record this in the keyboard state.
  mState.mButton = atoi(iArgv[4]);
  
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseUp( atoi(iArgv[2]),  YFlip(frame, atoi(iArgv[3])), mState );

  // Clear this in the keyboard state.
  mState.mButton = 0;

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
    frame->KeyDown( atoi(iArgv[2]),  YFlip(frame, atoi(iArgv[3])), mState );
    
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
    frame->KeyUp( atoi(iArgv[2]),  YFlip(frame, atoi(iArgv[3])), mState );
    
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
  mState.mButton = 0;
  mState.msKey = "";

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
    throw logic_error( "Couldn't initialize Togl" );
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



ToglFrame::ToglFrame( ID iID ) {

  mID = iID;
}

ToglFrame::~ToglFrame() {

}

void 
ToglFrame::Draw() {

  this->DoDraw();
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

  this->DoReshape();
}

void 
ToglFrame::Timer() {

  this->DoTimer();
}

void
ToglFrame::MouseMoved( int inX, int inY, InputState& iState ) {

  if( inX > 0 && inX < mWidth-1 &&
      inY > 0 && inY < mHeight-1 ) {
    this->DoMouseMoved( inX, inY, iState );
  }
}

void
ToglFrame::MouseUp( int inX, int inY, InputState& iState ) {

  if( inX > 0 && inX < mWidth-1 &&
      inY > 0 && inY < mHeight-1 ) {
    this->DoMouseUp( inX, inY, iState );
  }
}

void
ToglFrame::MouseDown( int inX, int inY, InputState& iState ) {

  if( inX > 0 && inX < mWidth-1 &&
      inY > 0 && inY < mHeight-1 ) {
    this->DoMouseDown( inX, inY, iState );
  }
}

void
ToglFrame::KeyDown( int inX, int inY, InputState& iState ) {

  if( inX > 0 && inX < mWidth-1 &&
      inY > 0 && inY < mHeight-1 ) {
    this->DoKeyDown( inX, inY, iState );
  }
}

void
ToglFrame::KeyUp( int inX, int inY, InputState& iState ) {

  if( inX > 0 && inX < mWidth-1 &&
      inY > 0 && inY < mHeight-1 ) {
    this->DoKeyUp( inX, inY, iState );
  }
}

void
ToglFrame::DoDraw() {

  DebugOutput( << "ToglFrame " << mID << ": DoDraw()" );
}

void
ToglFrame::DoReshape() {

  DebugOutput( << "ToglFrame " << mID << ": DoReshape()" );
}

void
ToglFrame::DoTimer() {

  DebugOutput( << "ToglFrame " << mID << ": DoTimer()" );
}

void
ToglFrame::DoMouseMoved( int inX, int inY, InputState& iState ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseMoved()" );
}

void
ToglFrame::DoMouseUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseUp()" );
}

void
ToglFrame::DoMouseDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseDown()" );
}

void
ToglFrame::DoKeyDown( int inX, int inY, InputState& iState ) {

  DebugOutput( << "ToglFrame " << mID << ": DoKeyDown()" );
}

void
ToglFrame::DoKeyUp( int inX, int inY, InputState& iState ) {

  DebugOutput( << "ToglFrame " << mID << ": DoKeyUp()" );
}
