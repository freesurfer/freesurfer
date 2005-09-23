#ifndef InputState_h
#define InputState_h

#include <iostream>
#include "string_fixed.h"
#include "ScubaKeyCombo.h"

class InputState {
  friend class InputStateTester;
 public:
  InputState ();
  bool IsShiftKeyDown ();
  bool IsAltKeyDown ();
  bool IsControlKeyDown ();
  int Button ();
  ScubaKeyCombo* Key ();

  // For polling button state.
  bool IsButtonDown ();
  bool IsButtonUp ();

  // For getting the exact mouse event.
  bool IsButtonDownEvent ();
  bool IsButtonDragEvent ();
  bool IsButtonUpEvent ();

  // For setting the event type.
  void SetButtonDownEvent ( int iButton );
  void SetButtonDragEvent ();
  void SetButtonUpEvent ();
  void ClearEvents ();

  void SetShiftKey ( bool ibShiftKey );
  void SetAltKey ( bool ibAltKey );
  void SetControlKey ( bool ibControlKey );
  void SetKeyFromString ( std::string isKey );
  void SetKeyCode ( int iKey );

 protected:
  bool mbShiftKey;
  bool mbAltKey;
  bool mbControlKey;
  bool mbButtonDownEvent;
  bool mbButtonDragEvent;
  bool mbButtonUpEvent;
  int mButton;
  ScubaKeyCombo* mKey;
};

std::ostream& operator << ( std::ostream& os, InputState& iInput );



#endif
