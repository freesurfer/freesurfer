/**
 * @file  ToglManager.h
 * @brief Main C++ object started by Tcl
 *
 * Togl is a bunch of code that implements a Tk GL
 * context. ToglManager is a class that interfaces with togl and acts
 * as the frontline for getting events from Tcl via the togl
 * object. When the ToglManager starts up, it registers a callback for
 * the creation of a togl object, and when the Tcl code creates one,
 * this code makes a WindowFrame, our base Window object. ToglManager
 * takes events and passes them to WindowFrames.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.16 $
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

  // Maps togl ids to frames.
  static std::map<int,WindowFrame*> mFrames;

  static WindowFrameFactory* mFactory;

  // Our modifers for shift, contrl, and alt keys.
  static InputState mState;

  static int mCurrentWindowCoords[2];
};


#endif
