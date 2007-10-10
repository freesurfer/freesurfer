/**
 * @file  ToglManager.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/10 17:34:03 $
 *    $Revision: 1.14 $
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


#ifndef ToglManager_h
#define ToglManager_h

#include <map>
extern "C" {
#include "togl.h"
}
#include "string_fixed.h"
#include "WindowFrame.h"
#include "DebugReporter.h"
#include "InputState.h"

class ToglManager {

public:
  // These are the Togl callbacks that we will register to get event
  // notifications from the Togl system. ToglManager determines to
  // which window they apply and calls ToglFrame functions in a
  // Togl-free way.
  static void DrawCallback ( struct Togl* iTogl );
  static void CreateCallback ( struct Togl* iTogl );
  static void DestroyCallback ( struct Togl* iTogl );
  static void ReshapeCallback ( struct Togl* iTogl );
  static void TimerCallback ( struct Togl* iTogl );

  // These are the Togl frame-specific callbacks that are attached to
  // the Togl Tcl/Tk object with the Tk bind command. ToglManager
  // determines to which window they apply, parses the arguments, and
  // calls ToglFrame functions in a Togl-free way.
  static int MouseMotionCallback ( struct Togl* iTogl, int, char* iArgv[] );
  static int MouseDownCallback ( struct Togl* iTogl, int, char* iArgv[] );
  static int MouseUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int KeyDownCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int KeyUpCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int ExitCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );
  static int EnterCallback ( struct Togl* iTogl, int iArgc, char* iArgv[] );

  // The main entry point should call this function to register all
  // the callbacks.
  static void InitializeTogl ( Tcl_Interp *iInterp );

  // Gets the static instance of the manager.
  static ToglManager& GetManager();

  // Sets the factory to use for creating new frames.
  static void SetFrameFactory( WindowFrameFactory* iFactory ) {
    mFactory = iFactory;
  }

protected:

  static inline int YFlip ( WindowFrame* iFrame, int iY ) {
    return (iFrame->GetHeight() - iY);
  }

  // Maps Window IDs to frames.
  static std::map<WindowFrame::ID,WindowFrame*> mFrames;

  static WindowFrameFactory* mFactory;

  // Our modifers for shift, contrl, and alt keys.
  static InputState mState;

  static int mCurrentWindowCoords[2];
};


#endif
