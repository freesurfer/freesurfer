#ifndef InputState_h
#define InputState_h

#include <iostream>
#include <string>

class InputState {
  friend class ToglManager;
  friend class InputStateTester;
 public:
  InputState ();
  bool IsShiftKeyDown ();
  bool IsAltKeyDown ();
  bool IsControlKeyDown ();
  int Button ();
  std::string Key ();

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

 protected:
  bool mbShiftKey;
  bool mbAltKey;
  bool mbControlKey;
  bool mbButtonDownEvent;
  bool mbButtonDragEvent;
  bool mbButtonUpEvent;
  int mButton;
  std::string msKey;
};

std::ostream& operator << ( std::ostream& os, InputState& iInput );



#endif
