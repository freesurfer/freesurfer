/**
 * @file  InputState.cpp
 * @brief Represents current input state
 *
 * Used by the ToglManager to communicate input state information to
 * the frame, and then downwards through to the view and layer.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.13 $
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


#include "InputState.h"

using namespace std;

InputState::InputState () :
    mbShiftKey(false), 
    mbAltKey(false),
    mbControlKey(false),
    mbButtonDownEvent(false),
    mbButtonDragEvent(false), 
    mbButtonUpEvent(false),
    mButton(0),
    mKey( ScubaKeyCombo::MakeKeyCombo() ) {

  mDelta[0] = mDelta[1] = 0;
}

bool InputState::IsShiftKeyDown () const {
  return mbShiftKey;
}

bool InputState::IsAltKeyDown () const {
  return mbAltKey;
}

bool InputState::IsControlKeyDown () const {
  return mbControlKey;
}

ScubaKeyCombo const& InputState::Key () const {
  return *mKey;
}

int InputState::Button () const {
  return mButton;
}

bool
InputState::IsButtonDown () const {
  return (mbButtonDownEvent || mbButtonDragEvent);
}

bool
InputState::IsButtonUp () const {
  return (mbButtonUpEvent ||
          (!(mbButtonDragEvent) && !(mbButtonDownEvent)));
}

bool
InputState::IsButtonDownEvent () const {
  return mbButtonDownEvent;
}

bool
InputState::IsButtonDragEvent () const {
  return mbButtonDragEvent;
}

bool
InputState::IsButtonUpEvent () const {
  return mbButtonUpEvent;
}

int
InputState::GetTotalButtonDeltaX () const {
  return mDelta[0];
}

int
InputState::GetTotalButtonDeltaY () const {
  return mDelta[1];
}

void
InputState::ClearEvents () {
  mbButtonDownEvent = false;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = false;
  mButton = 0;
  mDelta[0] = mDelta[1] = 0;
  //  mbShiftKey   = false;
  //  mbAltKey     = false;
  //  mbControlKey = false;
}

void
InputState::SetShiftKey ( bool ibShiftKey ) {
  mbShiftKey = ibShiftKey;
  mKey->SetShiftKeyDown( ibShiftKey );
}

void
InputState::SetAltKey ( bool ibAltKey ) {
  mbAltKey = ibAltKey;
  mKey->SetAltKeyDown( ibAltKey );
}

void
InputState::SetControlKey ( bool ibControlKey ) {
  mbControlKey = ibControlKey;
  mKey->SetControlKeyDown( ibControlKey );
}

void
InputState::SetKeyFromString ( std::string const& isKey ) {
  mKey->SetFromString( isKey );
}

void
InputState::SetKeyCode ( int iKey ) {
  mKey->SetKeyCode( iKey );
}

void
InputState::CopyKey ( ScubaKeyCombo const& iKey ) {
  mKey->CopyFrom( iKey );
}

void
InputState::AddButtonDelta ( int iX, int iY ) {
  mDelta[0] += iX;
  mDelta[1] += iY;
}

void
InputState::SetButtonDownEvent ( int iButton ) {
  mbButtonDownEvent = true;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = false;
  mButton = iButton;
  mDelta[0] = mDelta[1] = 0;
}

void
InputState::SetButtonDragEvent () {
  mbButtonDownEvent = false;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = true;
}

void
InputState::SetButtonUpEvent () {
  mbButtonDownEvent = false;
  mbButtonUpEvent   = true;
  mbButtonDragEvent = false;
}

std::ostream& operator << ( std::ostream& os, InputState const& iInput ) {
  if ( iInput.IsShiftKeyDown() )
    os << "shift ";
  if ( iInput.IsAltKeyDown() )
    os << "alt ";
  if ( iInput.IsControlKeyDown() )
    os << "control ";
  if ( iInput.Button() != 0 )
    os << "button " << iInput.Button();
  if ( iInput.Key().GetKeyCode() != -1 )
    os << "key " << iInput.Key();
  return os;
}

