#include "InputState.h"

using namespace std;

InputState::InputState () { 
  mbShiftKey = false; 
  mbAltKey = false; 
  mbControlKey = false; 
  mButton = 0;
  msKey = "";
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

string InputState::Key () { 
  return msKey;
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
InputState::SetButtonDownEvent ( int iButton ) {
  mbButtonDownEvent = true;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = false;
  mButton = iButton;
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
  mButton = 0;
}

void
InputState::ClearEvents () {
  mbButtonDownEvent = false;
  mbButtonUpEvent   = false;
  mbButtonDragEvent = false;
  mButton = 0;
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
  if( iInput.Key() != "" )
    os << "key " << iInput.Key();
  return os;
}

