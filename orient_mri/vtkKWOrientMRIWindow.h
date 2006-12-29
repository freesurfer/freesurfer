/**
 * @file  vtkKWOrientMRIWindow.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:11 $
 *    $Revision: 1.4 $
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


#ifndef vtkKWOrientMRIWindow_h
#define vtkKWOrientMRIWindow_h

#include "vtkKWWindow.h"
#include "vtkKWOrientMRIView.h"

class vtkKWMenu;
class vtkKWPushButton;

class vtkKWOrientMRIWindow : public vtkKWWindow {

public:

  static vtkKWOrientMRIWindow* New ();
  vtkTypeRevisionMacro( vtkKWOrientMRIWindow, vtkKWWindow );

  // Override to create our interior.
  virtual void Create ();

  void LoadVolumeFromDlog ();
  void LoadVolume ( const char* ifnVolume );

  void SaveVolumeWithConfirm ();
  void SaveVolumeAsFromDlog ();

  void RevertToSavedTransform ();

  void TransformVolume ();

  void RestoreView ();

  void ZoomBy ( float iFactor );
  void ZoomIn ();
  void ZoomOut ();

protected:

  vtkKWOrientMRIWindow ();
  virtual ~vtkKWOrientMRIWindow ();

  vtkKWOrientMRIView* mView;


  // Enable or disable buttons and menu items based on program state.
  void UpdateCommandStatus ();

  // For keeping track of buttons and menu items that relate to
  // commands, and whether or not they should be enabled.
  //BTX
  enum Command {
    CmdLoadVolume = 0,
    CmdSaveVolume,
    CmdSaveVolumeAs,
    CmdTransformVolume,
    CmdRevertVolume,
    CmdRestoreView,
    CmdZoomOut,
    CmdZoomIn,
    kcCommands
  };
  bool maCommandEnabled[kcCommands];

  // Struct for associating a menu and an entry item.
  typedef struct {
    vtkKWMenu* menu;
    int nItem;
  }
  MenuItem;

  // The menu items associated with each command.
  MenuItem maMenuItems[kcCommands];

  // The toolbar button associated with each command.
  vtkKWPushButton* maPushButtons[kcCommands];

  // The icons to associate with our commands.
  vtkKWIcon* maIcons[kcCommands];

  //ETX
};

#endif
