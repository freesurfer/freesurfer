/**
 * @file  vtkKWOrientMRIApp.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:11 $
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
