#include <stdexcept>
#include "ToglManager.h"

using namespace std;

map<ToglFrame::ID,ToglFrame*> ToglManager::mFrames;
ToglFrameFactory* ToglManager::mFactory = NULL;

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
}

void
ToglManager::TimerCallback ( struct Togl* iTogl ) {

  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->Timer();
}

int
ToglManager::MouseMotionCallback ( struct Togl* iTogl, 
				   int iArgc, char* iArgv[] ) {

  if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseMoved( atoi(iArgv[2]), atoi(iArgv[3]), atoi(iArgv[4]), 0 );
  return TCL_OK;
}

int
ToglManager::MouseDownCallback ( struct Togl* iTogl, 
				 int iArgc, char* iArgv[] ) {

  if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseDown( atoi(iArgv[2]), atoi(iArgv[3]), atoi(iArgv[4]), 0 );
  return TCL_OK;
}

int
ToglManager::MouseUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->MouseUp( atoi(iArgv[2]), atoi(iArgv[3]), atoi(iArgv[4]), 0 );
  return TCL_OK;
}

int
ToglManager::KeyDownCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->KeyDown( atoi(iArgv[2]), atoi(iArgv[3]), string(iArgv[4]), 0 );
  return TCL_OK;
}

int
ToglManager::KeyUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] ) {

  if( iArgc != 5 ) {
    return TCL_ERROR;
  }
  ToglFrame::ID id = atoi( Togl_Ident( iTogl ));
  ToglFrame* frame = mFrames[id];
  frame->KeyUp( atoi(iArgv[2]), atoi(iArgv[3]), string(iArgv[4]), 0 );
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
  Togl_DestroyFunc( ToglManager::CreateCallback );
  Togl_DisplayFunc( ToglManager::DrawCallback );
  Togl_ReshapeFunc( ToglManager::ReshapeCallback );
  Togl_TimerFunc( ToglManager::TimerCallback );

  // Create our Tcl commands that will be bound to the Togl object as
  // callbacks.
  Togl_CreateCommand( "MouseMotionCallback", ToglManager::MouseMotionCallback );
  Togl_CreateCommand( "MouseDownCallback", ToglManager::MouseDownCallback );
  Togl_CreateCommand( "MouseUpCallback", ToglManager::MouseUpCallback );
  Togl_CreateCommand( "KeyDownCallback", ToglManager::KeyDownCallback );
  Togl_CreateCommand( "KeyUpCallback", ToglManager::KeyUpCallback );
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
ToglFrame::MouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  this->DoMouseMoved( inX, inY, iButton, iModifiers );
}

void
ToglFrame::MouseUp( int inX, int inY, int iButton, int iModifers ) {

  this->DoMouseUp( inX, inY, iButton, iModifers );
}

void
ToglFrame::MouseDown( int inX, int inY, int iButton, int iModifers ) {

  this->DoMouseDown( inX, inY, iButton, iModifers );
}

void
ToglFrame::KeyDown( int inX, int inY, string isKey, int iModifers ) {

  this->DoKeyDown( inX, inY, isKey, iModifers );
}

void
ToglFrame::KeyUp( int inX, int inY, string isKey, int iModifers ) {

  this->DoKeyUp( inX, inY, isKey, iModifers );
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
ToglFrame::DoMouseMoved( int inX, int inY, int iButton, int iModifiers ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseMoved()" );
}

void
ToglFrame::DoMouseUp( int inX, int inY, int iButton, int iModifers ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseUp()" );
}

void
ToglFrame::DoMouseDown( int inX, int inY, int iButton, int iModifers ) {

  DebugOutput( << "ToglFrame " << mID << ": DoMouseDown()" );
}

void
ToglFrame::DoKeyDown( int inX, int inY, string isKey, int iModifers ) {

  DebugOutput( << "ToglFrame " << mID << ": DoKeyDown()" );
}

void
ToglFrame::DoKeyUp( int inX, int inY, string isKey, int iModifers ) {

  DebugOutput( << "ToglFrame " << mID << ": DoKeyUp()" );
}


