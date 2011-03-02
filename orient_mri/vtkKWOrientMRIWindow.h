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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.11 $
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


#ifndef vtkKWOrientMRIWindow_h
#define vtkKWOrientMRIWindow_h

#include "vtkKWWindow.h"
#include "vtkCommand.h"
#include "vtkImageReslice.h"
#include "vtkSmartPointer.h"

class vtkFSVolumeSource;
class vtkFreesurferLookupTable;
class vtkKWMenu;
class vtkKWOrientMRIView2D;
class vtkKWOrientMRIView3D;
class vtkKWPushButton;
class vtkKWRadioButton;
class vtkKWScale;
class vtkLookupTable;
class vtkMatrix4x4;
class vtkScalarsToColors;

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
  // Set the volume's transform based on the current user transform,
  // rebuilds the RAS transform, and resets the user transform.
  void TransformVolumeWithUserTransform ();

  // Description: 
  // Make a transformation matrix based on the orthogonalized
  // reorthogonalization lines, and compose that with the user
  // transform.
  void ReorthogonalizeVolume ();

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

  // Description:
  // Get the user transform from the 3D view and copy it into
  // mUserTransform, updating the 2D views.
  void UpdateUserTransform ();

  // Description:
  // Our callback for UserTransformChanged events.
  static void UserTransformChanged ( vtkObject* iCaller, 
				     unsigned long iEventId,
				     void* iClientData, void* iCallData );
  
  // Description:
  // Show a dialog from which the user can enter a number of degrees
  // to rotate the user transform. The first argument is the axis and
  // will be passed to RotateUserTransform. The second transform
  // should be -1 or 1, which will be multiplied to the value the user
  // enter.
  void DoRotateUserTransformDialog ( int iAxis, int iMultiplier );

protected:

  vtkKWOrientMRIWindow ();
  virtual ~vtkKWOrientMRIWindow ();

  // Our view objects.
  //BTX
  enum { X = 0, Y = 1, Z = 2 };
  vtkSmartPointer<vtkKWOrientMRIView3D> mView3D;
  vtkSmartPointer<vtkKWOrientMRIView2D> mView2D[3];

  // Our scale widgets to control the through plane for the 2D views.
  vtkSmartPointer<vtkKWScale> mScaleThroughPlane[3];
  
  // Dirty if we've done any transforming that could be saved.
  bool mbDirty;

  // Description:
  // Enable or disable buttons and menu items based on program state.
  void UpdateCommandStatus ();

  // Description:
  // Recalculate our mVolumeToRASTransform with the data in the volume
  // source.
  void RebuildVolumeToRASTransform ();
  
  // Description:
  // Get the current user transform.
  vtkMatrix4x4& GetUserTransform ();

  // Inserts a 5 pixel wide spacer into the toolbar.
  enum { kToolbarSpacerWidth = 5 };
  void AddSpacerToToolbar ( vtkKWToolbar& iToolbar, 
			    int iWidth = kToolbarSpacerWidth );

  // Mini class for creating a menu item and enabling/disabling
  // it. Parameters are self explanatory.
  class MenuItem {
  public:
    MenuItem ();

    // Creates a simple menu command.
    void MakeCommand ( vtkKWMenu* iMenu, int inItem, const char* isText,
		       vtkObject* iCommandObject, const char* isCommand, 
		       const char* isAccelerator, const char* isIconKey );

    // Create a check box item.
    void MakeCheckButton ( vtkKWMenu* iMenu, int inItem, const char* isText,
			   vtkObject* iCommandObject, const char* isCommand, 
			   const char* isAccelerator, const char* isIconKey );
    
    // Disable or enable the menu item.
    void SetStateToDisabled ();
    void SetStateToNormal ();

    // Set/get the check value if it's a check button.
    void SetSelectedState ( int ibOn );
    int GetSelectedState ();
    
    // Get the index of the menu item.
    int GetIndex () const;

  protected:
    vtkSmartPointer<vtkKWMenu> mMenu;
    int mnItem ;
  };

  // The menu items associated with each command.
  MenuItem* mMenuLoadVolume;
  MenuItem* mMenuSaveVolume;
  MenuItem* mMenuSaveVolumeAs;
  MenuItem* mMenuLoadLUT;
  MenuItem* mMenuTransformVolume;
  MenuItem* mMenuReorthogonalizeVolume;
  MenuItem* mMenuRevertVolume;
  MenuItem* mMenuRestoreView;
  MenuItem* mMenuZoomOut;
  MenuItem* mMenuZoomIn;
  MenuItem* mMenuUseGrayScaleColors;
  MenuItem* mMenuUseLUTColors;

  // The toolbar buttons associated with each command.
  vtkSmartPointer<vtkKWPushButton> mBtnLoadVolume;
  vtkSmartPointer<vtkKWPushButton> mBtnSaveVolume;
  vtkSmartPointer<vtkKWPushButton> mBtnTransformVolume;
  vtkSmartPointer<vtkKWPushButton> mBtnRevertVolume;
  vtkSmartPointer<vtkKWPushButton> mBtnRestoreView;
  vtkSmartPointer<vtkKWPushButton> mBtnZoomOut;
  vtkSmartPointer<vtkKWPushButton> mBtnZoomIn;
  vtkSmartPointer<vtkKWRadioButton> mRadBtnUseGrayScaleColors;
  vtkSmartPointer<vtkKWRadioButton> mRadBtnUseLUTColors;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateXPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateXNeg;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateYPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateYNeg;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateZPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateZNeg;

  // Data.
  vtkSmartPointer<vtkFSVolumeSource> mVolume;
  vtkSmartPointer<vtkImageReslice> mVolumeToRASTransform;
  vtkSmartPointer<vtkImageReslice> mUserTransform;
  vtkSmartPointer<vtkLookupTable> mGrayScaleColors;
  vtkSmartPointer<vtkFreesurferLookupTable> mLUTColors;

  // Transform objects.
  vtkSmartPointer<vtkMatrix4x4> mOriginalVoxelToRASMatrix;

  //ETX
};

#endif
