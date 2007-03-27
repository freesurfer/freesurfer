/**
 * @file  vtkKWOrientMRIWindow.h
 * @brief Loads and works on data, handles UI commands
 *
 * Populates menus and toolbars with commands. Handles dialog boxes
 * for loading and saving data. Owns the data objects. Calculates the
 * new transform based on the camera orientation and writes it to the
 * volume.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/27 21:24:36 $
 *    $Revision: 1.5 $
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
class vtkFSVolumeSource;
class vtkScalarsToColors;
class vtkMatrix4x4;

class vtkKWOrientMRIWindow : public vtkKWWindow {

public:

  static vtkKWOrientMRIWindow* New ();
  vtkTypeRevisionMacro( vtkKWOrientMRIWindow, vtkKWWindow );

  // Description:
  // Override to create our interior.
  virtual void Create ();

  // Description:
  // Load a volume and display it in the view.
  void LoadVolumeFromDlog ();
  void LoadVolume ( const char* ifnVolume );

  // Description:
  // Save the volume with its new transform.
  void SaveVolumeWithConfirm ();
  void SaveVolumeAsFromDlog ();
  void SaveVolume ( const char* ifnVolume );
  void SaveVolume ();

  // Description:
  // Throw away the user modified transform and go back to the
  // original version.
  void RevertToSavedTransform ();
  
  // Description:
  // Set the volume's transform based on the current camera
  // orientation.
  void TransformVolume ();

  // Description:
  // Go back to the original camera position that shows the entirety
  // of the volume with the X plane parallel to the camera plane.
  void RestoreView ();

  // Description:
  // Zoom.
  void ZoomBy ( float iFactor );
  void ZoomIn ();
  void ZoomOut ();

protected:

  vtkKWOrientMRIWindow ();
  virtual ~vtkKWOrientMRIWindow ();

  // Our view object.
  vtkKWOrientMRIView* mView;

  // Dirty if we've done any transforming that could be saved.
  bool mbDirty;

  // Description:
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
  
  // Data.
  vtkFSVolumeSource* mVolume;
  vtkLookupTable* mLUT;

  // Transform objects.
  vtkMatrix4x4* mOriginalVoxelToRASMatrix;
  vtkMatrix4x4* mOriginalView;
  vtkMatrix4x4* mOriginalViewI;

  //ETX
};

#endif
