/**
 * @file  InputState.h
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
 *    $Revision: 1.12 $
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


#ifndef InputState_h
#define InputState_h

#include <iostream>
#include "string_fixed.h"
#include "ScubaKeyCombo.h"

class InputState {
  friend class InputStateTester;
public:

  InputState ();

  // Acessors.
  // Key modifiers.
  bool IsShiftKeyDown () const;
  bool IsAltKeyDown () const;
  bool IsControlKeyDown () const;
  ScubaKeyCombo const& Key () const;

  // Mouse button state.
  int Button () const;
  bool IsButtonDown () const;
  bool IsButtonUp () const;
  bool IsButtonDownEvent () const;
  bool IsButtonDragEvent () const;
  bool IsButtonUpEvent () const;

  // Returns the button delta during a drag event. This is in window
  // coords.
  int GetTotalButtonDeltaX () const;
  int GetTotalButtonDeltaY () const;


  // Settors.
  void ClearEvents ();

  // Key modifiers and codes.
  void SetShiftKey ( bool ibShiftKey );
  void SetAltKey ( bool ibAltKey );
  void SetControlKey ( bool ibControlKey );

  // For modifying the ScubaKeyCombo.
  void SetKeyFromString ( std::string const& isKey );
  void SetKeyCode ( int iKey );
  void CopyKey ( ScubaKeyCombo const& iKey );

  // Mouse information.
  void AddButtonDelta( int iX, int iY );
  void SetButtonDownEvent ( int iButton );
  void SetButtonDragEvent ();
  void SetButtonUpEvent ();

protected:

  bool mbShiftKey;
  bool mbAltKey;
  bool mbControlKey;
  bool mbButtonDownEvent;
  bool mbButtonDragEvent;
  bool mbButtonUpEvent;
  int mButton;
  ScubaKeyCombo* mKey;
  int mDelta[2];
};

std::ostream& operator << ( std::ostream& os, InputState const& iInput );



#endif
