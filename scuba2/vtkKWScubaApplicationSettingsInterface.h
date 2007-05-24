/**
 * @file  vtkKWScubaApplicationSettingsInterface.h
 * @brief Preferences for Scuba
 *
 * Subclass of vtkKWApplicationSettingsInterface that adds our own
 * custom settings.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/05/24 20:20:59 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2007,
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

#ifndef __vtkKWScubaApplicationSettingsInterface_h
#define __vtkKWScubaApplicationSettingsInterface_h

#include "vtkKWApplicationSettingsInterface.h"

class vtkKWCheckButton;
class vtkKWScubaWindow;

class vtkKWScubaApplicationSettingsInterface : public vtkKWApplicationSettingsInterface
{
public:
  static vtkKWScubaApplicationSettingsInterface* New();
  vtkTypeRevisionMacro(vtkKWScubaApplicationSettingsInterface,vtkKWApplicationSettingsInterface);

  // Description:
  // Create the widget.
  virtual void Create ();

  // Description:
  // Set/Get the window (do not ref count it since the window will ref count
  // this widget).
  vtkGetObjectMacro ( ScubaWindow, vtkKWScubaWindow );
  virtual void SetScubaWindow ( vtkKWScubaWindow* );

  // Description:
  // Refresh the interface given the current value of the Window and its
  // views/composites/widgets.
  virtual void Update ();

  // Description:
  // Update the "enable" state of the object and its internal parts.
  // Depending on different Ivars (this->Enabled, the application's 
  // Limited Edition Mode, etc.), the "enable" state of the object is updated
  // and propagated to its internal parts/subwidgets. This will, for example,
  // enable/disable parts of the widget UI, enable/disable the visibility
  // of 3D widgets, etc.
  virtual void UpdateEnableState ();

  // Description:
  // Callbacks. Internal, do not use.
  virtual void AutoSizeInfoAreaCallback ( int state );

protected:
  vtkKWScubaApplicationSettingsInterface ();
  ~vtkKWScubaApplicationSettingsInterface ();

  // Interface settings
  vtkKWCheckButton* mChkBtnAutoSizeInfoArea;

  // Window.
  vtkKWScubaWindow* ScubaWindow;

private:
  vtkKWScubaApplicationSettingsInterface(const vtkKWScubaApplicationSettingsInterface&); // Not implemented
  void operator=(const vtkKWScubaApplicationSettingsInterface&); // Not Implemented
};

#endif


