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

