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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


