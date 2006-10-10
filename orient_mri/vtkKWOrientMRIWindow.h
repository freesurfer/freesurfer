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
  enum Command { CmdLoadVolume = 0,
		 CmdSaveVolume,
		 CmdTransformVolume,
		 CmdRevertVolume, 
		 CmdRestoreView,
		 CmdZoomOut,
		 CmdZoomIn,
		 kcCommands };
  bool maCommandEnabled[kcCommands];

  // Struct for associating a menu and an entry item.
  typedef struct { vtkKWMenu* menu; int nItem; } MenuItem;

  // The menu items associated with each command.
  MenuItem maMenuItems[kcCommands];

  // The toolbar button associated with each command.
  vtkKWPushButton* maPushButtons[kcCommands];
  //ETX
};
  
#endif
