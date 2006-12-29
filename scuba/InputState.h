/**
 * @file  InputState.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:13 $
 *    $Revision: 1.10 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
  bool IsShiftKeyDown ();
  bool IsAltKeyDown ();
  bool IsControlKeyDown ();
  int Button ();
  ScubaKeyCombo* Key ();

  // For polling button state.
  bool IsButtonDown ();
  bool IsButtonUp ();

  // Returns the button delta during a drag event. This is in window
  // coords.
  void AddButtonDelta( int iX, int iY );
  int GetTotalButtonDeltaX ();
  int GetTotalButtonDeltaY ();

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
  int mDelta[2];
};

std::ostream& operator << ( std::ostream& os, InputState& iInput );



#endif
