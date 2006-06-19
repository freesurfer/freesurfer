#include "InputState.h"

using namespace std;

InputState::InputState () :
  mbShiftKey(false), mbAltKey(false), mbControlKey(false), 
  mbButtonDownEvent(false), mbButtonDragEvent(false), mbButtonUpEvent(false),
  mButton(0), mKey(NULL) { 
  mKey = ScubaKeyCombo::MakeKeyCombo();
  mDelta[0] = mDelta[1] = 0;
}

bool InputState::IsShiftKeyDown () { 
  return mbShiftKey; 
}

bool InputState::IsAltKeyDown () { 
  return mbAltKey;
}

bool InputState::IsControlKeyDown () { 
  return mbControlKey;
}

int InputState::Button () { 
  return mButton;
}

ScubaKeyCombo* InputState::Key () { 
  return mKey;
}

bool 
InputState::IsButtonDown () { 
  return (mbButtonDownEvent || mbButtonDragEvent); 
}

bool 
InputState::IsButtonUp () { 
  return (mbButtonUpEvent || 
	  (!(mbButtonDragEvent) && !(mbButtonDownEvent))); 
}

bool 
InputState::IsButtonDownEvent () { 
  return mbButtonDownEvent; 
}

bool 
InputState::IsButtonDragEvent () { 
  return mbButtonDragEvent; 
}

bool 
InputState::IsButtonUpEvent () { 
  return mbButtonUpEvent; 
}

void
InputState::AddButtonDelta ( int iX, int iY ) {
  mDelta[0] += iX;
  mDelta[1] += iY;
}

int
InputState::GetTotalButtonDeltaX () {
  return mDelta[0];
}

int
InputState::GetTotalButtonDeltaY () {
  return mDelta[1];
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

void
InputState::ClearEvents () {
  mbButtonDownEvent = false;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = false;
  mButton = 0;
  mDelta[0] = mDelta[1] = 0;
  mbShiftKey   = false;
  mbAltKey     = false;
  mbControlKey = false;
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
InputState::SetKeyFromString ( std::string isKey ) { 
  mKey->SetFromString( isKey ); 
}

void 
InputState::SetKeyCode ( int iKey ) { 
  mKey->SetKeyCode( iKey ); 
}

std::ostream& operator << ( std::ostream& os, InputState& iInput ) { 
  if( iInput.IsShiftKeyDown() )
    os << "shift ";
  if( iInput.IsAltKeyDown() )
    os << "alt ";
  if( iInput.IsControlKeyDown() )
    os << "control ";
  if( iInput.Button() != 0 )
    os << "button " << iInput.Button();
  if( iInput.Key()->GetKeyCode() != -1 )
    os << "key " << iInput.Key();
  return os;
}

