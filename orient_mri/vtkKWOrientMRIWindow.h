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
 *    $Date: 2007/03/28 20:04:49 $
 *    $Revision: 1.7 $
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
class vtkKWRadioButton;
class vtkFreesurferLookupTable;
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
  // Load a color table for the volume.
  void LoadLUTFromDlog ();
  void LoadLUT ( const char* ifnLUT );

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

  // Description:
  // Change the color table for the volume.
  void UseGrayScaleColors ();
  void UseLUTColors ();
  

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

  //BTX
  enum Command {
    CmdLoadVolume = 0,
    CmdSaveVolume,
    CmdSaveVolumeAs,
    CmdLoadLUT,
    CmdTransformVolume,
    CmdRevertVolume,
    CmdRestoreView,
    CmdZoomOut,
    CmdZoomIn,
    CmdUseGrayScaleColors,
    CmdUseLUTColors,
    CmdRotateXPos,
    CmdRotateXNeg,
    CmdRotateYPos,
    CmdRotateYNeg,
    CmdRotateZPos,
    CmdRotateZNeg,
    kcCommands
  };

  // Struct for associating a menu and an entry item.
  typedef struct {
    vtkKWMenu* menu;
    int nItem;
  }
  MenuItem;

  // The menu items associated with each command.
  MenuItem mMenuLoadVolume;
  MenuItem mMenuSaveVolume;
  MenuItem mMenuSaveVolumeAs;
  MenuItem mMenuLoadLUT;
  MenuItem mMenuTransformVolume;
  MenuItem mMenuRevertVolume;
  MenuItem mMenuRestoreView;
  MenuItem mMenuZoomOut;
  MenuItem mMenuZoomIn;
  MenuItem mMenuUseGrayScaleColors;
  MenuItem mMenuUseLUTColors;

  // The toolbar buttons associated with each command.
  vtkKWPushButton* mBtnLoadVolume;
  vtkKWPushButton* mBtnSaveVolume;
  vtkKWPushButton* mBtnTransformVolume;
  vtkKWPushButton* mBtnRevertVolume;
  vtkKWPushButton* mBtnRestoreView;
  vtkKWPushButton* mBtnZoomOut;
  vtkKWPushButton* mBtnZoomIn;
  vtkKWRadioButton* mRadBtnUseGrayScaleColors;
  vtkKWRadioButton* mRadBtnUseLUTColors;
  vtkKWPushButton* mBtnRotateXPos;
  vtkKWPushButton* mBtnRotateXNeg;
  vtkKWPushButton* mBtnRotateYPos;
  vtkKWPushButton* mBtnRotateYNeg;
  vtkKWPushButton* mBtnRotateZPos;
  vtkKWPushButton* mBtnRotateZNeg;

  // Data.
  vtkFSVolumeSource* mVolume;
  vtkLookupTable* mGrayScaleColors;
  vtkFreesurferLookupTable* mLUTColors;

  // Transform objects.
  vtkMatrix4x4* mOriginalVoxelToRASMatrix;
  vtkMatrix4x4* mOriginalView;
  vtkMatrix4x4* mOriginalViewI;

  //ETX
};

#endif
