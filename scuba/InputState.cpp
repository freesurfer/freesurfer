#include "InputState.h"

using namespace std;

InputState::InputState() { 
  mbShiftKey = false; 
  mbAltKey = false; 
  mbControlKey = false; 
  mButton = 0;
  msKey = "";
}

bool InputState::IsShiftKeyDown() { 
  return mbShiftKey; 
}

bool InputState::IsAltKeyDown() { 
  return mbAltKey;
}

bool InputState::IsControlKeyDown() { 
  return mbControlKey;
}

int InputState::Button() { 
  return mButton;
}

string InputState::Key() { 
  return msKey;
}

std::ostream& operator << ( std::ostream& os, InputState& iState ) { 
  if( iState.IsShiftKeyDown() )
    os << "shift ";
  if( iState.IsAltKeyDown() )
    os << "alt ";
  if( iState.IsControlKeyDown() )
    os << "control ";
  if( iState.Button() != 0 )
    os << "button " << iState.Button();
  if( iState.Key() != "" )
    os << "key " << iState.Key();
  return os;
}

