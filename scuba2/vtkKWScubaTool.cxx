/**
 * @file  vtkKWScubaTool.cxx
 * @brief Base tool class
 *
 * Base tool class that tools should subclass to implement editing,
 * navigating, etc functionality in reponse to mouse
 * clicks. Subclasses should implement the Do* functions to respond to
 * input events. They can also subclass the AddControls and
 * RemoveControls functions.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


#include "vtkKWScubaTool.h"
#include "vtkKWScubaWindow.h"
#include "vtkKWScubaView.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindowInteractor.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaTool );
vtkCxxRevisionMacro( vtkKWScubaTool, "$Revision: 1.2 $" );

vtkKWScubaTool::vtkKWScubaTool () :
    mButton(0),
    mbShiftKeyDown(false),
    mbControlKeyDown(false),
    mbAltKeyDown(false) {
  mTotalDelta[0] = mTotalDelta[1] = 0;
  mCurrentDelta[0] = mCurrentDelta[1] = 0;
}

vtkKWScubaTool::~vtkKWScubaTool () {}

void
vtkKWScubaTool::SetLabel ( const char* isLabel ) {
  msLabel = isLabel;
}

const char*
vtkKWScubaTool::GetLabel () {
  return msLabel.c_str();
}

void
vtkKWScubaTool::PopulateControlPage ( vtkKWWidget* iPanel ) {
  // Add the subclass's controls.
  this->AddControls( iPanel );
}

void
vtkKWScubaTool::DepopulateControlPage () {
  // Tell subclassses to remove controls.
  this->RemoveControls();
}

void
vtkKWScubaTool::MouseMoveEvent ( vtkKWScubaWindow* iWindow,
                                 vtkKWScubaView* iView,
                                 vtkKWScubaLayer* iLayer,
                                 float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  if ( this->WhichButtonDown() == 0 ) {
    this->DoMouseMove( iWindow, iView, iLayer, iRAS );
  } else {
    this->DoMouseDrag( iWindow, iView, iLayer, iRAS );
  }
}

void
vtkKWScubaTool::LeftMouseDownEvent ( vtkKWScubaWindow* iWindow,
                                     vtkKWScubaView* iView,
                                     vtkKWScubaLayer* iLayer,
                                     float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 1 );
  this->DoMouseDown( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::LeftMouseUpEvent ( vtkKWScubaWindow* iWindow,
                                   vtkKWScubaView* iView,
                                   vtkKWScubaLayer* iLayer,
                                   float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 1 ); // Set mouse button now for the event
  this->DoMouseUp( iWindow, iView, iLayer, iRAS );
  this->SetWhichButton( 0 ); // Clear mouse button after event handled
}

void
vtkKWScubaTool::MiddleMouseDownEvent ( vtkKWScubaWindow* iWindow,
                                       vtkKWScubaView* iView,
                                       vtkKWScubaLayer* iLayer,
                                       float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 2 );
  this->DoMouseDown( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::MiddleMouseUpEvent ( vtkKWScubaWindow* iWindow,
                                     vtkKWScubaView* iView,
                                     vtkKWScubaLayer* iLayer,
                                     float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 2 );
  this->DoMouseUp( iWindow, iView, iLayer, iRAS );
  this->SetWhichButton( 0 );
}

void
vtkKWScubaTool::RightMouseDownEvent ( vtkKWScubaWindow* iWindow,
                                      vtkKWScubaView* iView,
                                      vtkKWScubaLayer* iLayer,
                                      float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 3 );
  this->DoMouseDown( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::RightMouseUpEvent ( vtkKWScubaWindow* iWindow,
                                    vtkKWScubaView* iView,
                                    vtkKWScubaLayer* iLayer,
                                    float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->SetWhichButton( 3 );
  this->DoMouseUp( iWindow, iView, iLayer, iRAS );
  this->SetWhichButton( 0 );
}

void
vtkKWScubaTool::KeyDownEvent ( vtkKWScubaWindow* iWindow,
                               vtkKWScubaView* iView,
                               vtkKWScubaLayer* iLayer,
                               float iRAS[3] ) {

  this->SetStateFromRWI( iView->GetRenderWindowInteractor() );
  this->DoKeyDown( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::EnterEvent ( vtkKWScubaWindow* iWindow,
                             vtkKWScubaView* iView,
                             vtkKWScubaLayer* iLayer,
                             float iRAS[3] ) {

  this->ClearState();
  this->DoEnter( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::LeaveEvent ( vtkKWScubaWindow* iWindow,
                             vtkKWScubaView* iView,
                             vtkKWScubaLayer* iLayer,
                             float iRAS[3] ) {

  this->ClearState();
  this->DoLeave( iWindow, iView, iLayer, iRAS );
}

void
vtkKWScubaTool::AddControls ( vtkKWWidget* iPanel ) {}

void
vtkKWScubaTool::RemoveControls () {}

void
vtkKWScubaTool::DoMouseMove ( vtkKWScubaWindow*, vtkKWScubaView*,
                              vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoMouseDrag ( vtkKWScubaWindow*, vtkKWScubaView*,
                              vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoMouseDown ( vtkKWScubaWindow*, vtkKWScubaView*,
                              vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoMouseUp ( vtkKWScubaWindow*, vtkKWScubaView*,
                            vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoKeyDown ( vtkKWScubaWindow*, vtkKWScubaView*,
                            vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoEnter ( vtkKWScubaWindow*, vtkKWScubaView*,
                          vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoLeave ( vtkKWScubaWindow*, vtkKWScubaView*,
                          vtkKWScubaLayer*, float[3] ) {}

void
vtkKWScubaTool::DoClearState () {}

bool
vtkKWScubaTool::IsButtonDown () {
  return !(0 == mButton);
}

bool
vtkKWScubaTool::IsButtonUp () {
  return (0 == mButton);
}

int
vtkKWScubaTool::WhichButtonDown () {
  return mButton;
}

int
vtkKWScubaTool::GetTotalButtonDeltaX () {
  return mTotalDelta[0];
}

int
vtkKWScubaTool::GetTotalButtonDeltaY () {
  return mTotalDelta[1];
}

int
vtkKWScubaTool::GetButtonDeltaX () {
  return mCurrentDelta[0];
}

int
vtkKWScubaTool::GetButtonDeltaY () {
  return mCurrentDelta[1];
}


bool
vtkKWScubaTool::IsShiftKeyDown () {
  return mbShiftKeyDown;
}

bool
vtkKWScubaTool::IsControlKeyDown () {
  return mbControlKeyDown;
}

bool
vtkKWScubaTool::IsAltKeyDown () {
  return mbAltKeyDown;
}

void
vtkKWScubaTool::SetStateFromRWI ( vtkRenderWindowInteractor* iRWI ) {

  mCurrentDelta[0] =
    iRWI->GetEventPosition()[0] - iRWI->GetLastEventPosition()[0];
  mCurrentDelta[1] =
    iRWI->GetEventPosition()[1] - iRWI->GetLastEventPosition()[1];

  mTotalDelta[0] += mCurrentDelta[0];
  mTotalDelta[1] += mCurrentDelta[1];

  this->SetShiftKeyDown( iRWI->GetShiftKey() );
  this->SetControlKeyDown( iRWI->GetControlKey() );
  this->SetAltKeyDown( false ); // How to get this?
}

void
vtkKWScubaTool::ClearState () {

  mButton = 0;
  mTotalDelta[0] = mTotalDelta[1] = 0;
  mCurrentDelta[0] = mCurrentDelta[1] = 0;
  mbShiftKeyDown = false;
  mbControlKeyDown = false;
  mbAltKeyDown = false;

  this->DoClearState();
}

void
vtkKWScubaTool::SetWhichButton ( int iButton ) {
  mButton = iButton;
}

void
vtkKWScubaTool::SetShiftKeyDown ( bool ibDown ) {
  mbShiftKeyDown = ibDown;
}

void
vtkKWScubaTool::SetControlKeyDown ( bool ibDown ) {
  mbControlKeyDown = ibDown;
}

void
vtkKWScubaTool::SetAltKeyDown ( bool ibDown ) {
  mbAltKeyDown = ibDown;
}

