/**
 * @file  vtkKWScubaApp.h
 * @brief Main app for Scuba
 *
 * The main app code. Parse command line options and does some initial
 * set up.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/05/25 18:18:05 $
 *    $Revision: 1.2 $
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


// .NAME vtkKWScubaApp - Scuba application
// .SECTION Description
// Starts up the up, creates a vtkKWScubaWindow, and does some command
// line option processing.

#ifndef __vtkKWScubaApp_h
#define __vtkKWScubaApp_h

#include <string>
#include "vtkKWApplication.h"

class vtkKWScubaWindow;

class vtkKWScubaApp : public vtkKWApplication {

public:

  static vtkKWScubaApp* New ();
  vtkTypeRevisionMacro ( vtkKWScubaApp, vtkKWApplication );

  // Description:
  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );

  // Description:
  // Override to save/load our prefs.
  virtual void RestoreApplicationSettingsFromRegistry ();
  virtual void SaveApplicationSettingsToRegistry ();

  // Description:
  // Tells the window to load something.
  void LoadVolume ( const char* ifnVolume );
  void LoadSurface ( const char* ifnSurface );
  void LoadDTI ( const char* ifnDTI );
  void LoadPath ( const char* ifnPath );

protected:

  vtkKWScubaApp ();
  ~vtkKWScubaApp ();

  vtkKWScubaWindow* mWindow;

  int mMainWindowX;
  int mMainWindowY;
  int mMainWindowWidth;
  int mMainWindowHeight;
};

#endif
