/**
 * @file  vtkKWOrientMRIApp.h
 * @brief Command line parsing and startup
 *
 * Creates our window, does some setup stuff, parses command line
 * args.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/27 21:24:36 $
 *    $Revision: 1.3 $
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


#ifndef __vtkKWOrientMRIApp_h
#define __vtkKWOrientMRIApp_h

#include <string>
#include "vtkKWApplication.h"

class vtkKWOrientMRIWindow;

class vtkKWOrientMRIApp : public vtkKWApplication {

public:

  static vtkKWOrientMRIApp* New ();
  vtkTypeRevisionMacro ( vtkKWOrientMRIApp, vtkKWApplication );

  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );

  // Override to save/load our prefs.
  virtual void RestoreApplicationSettingsFromRegistry ();
  virtual void SaveApplicationSettingsToRegistry ();

  // Registry key constants.
  static const char* sMainWindowWidthRegKey;
  static const char* sMainWindowHeightRegKey;

  void LoadVolume ( const char* ifnVolume );

protected:

  vtkKWOrientMRIApp ();
  ~vtkKWOrientMRIApp ();

  vtkKWOrientMRIWindow* mWindow;

  int mMainWindowWidth;
  int mMainWindowHeight;
};

#endif
