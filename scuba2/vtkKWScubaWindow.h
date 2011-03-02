/**
 * @file  vtkKWScubaWindow.h
 * @brief The main window
 *
 * Manages the 'current' tool, view, and layer. Holds one or more
 * views and lays them out. Top level document loading and creation of
 * layers. Runs the info area and gets updates from mouseovered
 * objects. Runs the toolbars. Manages the UI panels for settings.
 * 
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.7 $
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


// .NAME vtkKWScubaWindow - Scuba window
// .SECTION Description
// Window subclass that loads a vtkKWScubaView in the view area,
// panels and stuff in the main area, a toolbar, and the info area in
// the secondary area.

#ifndef vtkKWScubaWindow_h
#define vtkKWScubaWindow_h

#include "vtkKWWindow.h"
#include "vtkKWScubaView.h"
#include "vtkKWScubaLayerCollection.h"
#include "Listener.h"
#include "vtkSmartPointer.h"

class vtkKWFrame;
class vtkKWMenu;
class vtkKWMenuButton;
class vtkKWMultiColumnList;
class vtkKWPushButton;
class vtkKWRadioButton;
class vtkKWRadioButtonSet;
class vtkKWScaleWithEntry;
class vtkKWScubaApplicationSettingsInterface;
class vtkKWScubaLayer;
class vtkKWScubaLayerCollection;
class vtkKWScubaTool;
class vtkKWToolbar;

class vtkKWScubaWindow : public vtkKWWindow
      //BTX
      , public Listener
      //ETX
{

public:

  static vtkKWScubaWindow* New ();
  vtkTypeRevisionMacro( vtkKWScubaWindow, vtkKWWindow );

  // Description:
  // Creates our UI.
  virtual void CreateWidget ();

  // Description:
  // Save our window geometry before deleting.
  virtual void PrepareForDelete ();

  // Desciption:
  // Override to return our custom interface class. This will end up
  // in the prefs dialog.
  virtual vtkKWApplicationSettingsInterface* GetApplicationSettingsInterface();

  // Description:
  // Bring up dialog boxes and then pass the results to LoadVolume,
  // LoadSurface, SaveVolume, or LoadDTI.
  void LoadVolumeFromDlog ();
  void LoadSurfaceFromDlog ();
  void SaveVolumeWithConfirm ();
  void LoadDTIFromDlog ();
  void LoadPathFromDlog ();
  void LoadODFFromDlog ();

  // Description:
  // Does all the UI stuff around loading a volume or surface. Creates
  // a proper layer collection for the data type loaded with the
  // filename we got, and sets the collection in the first unused
  // slot.
  void LoadVolume ( const char* ifnVolume );
  void LoadSurface ( const char* ifnSurface );
  void LoadDTI ( const char* ifnDTI );
  void LoadPath ( const char* ifnPath );
  void LoadODF ( const char* ifnODF );

  // Description:
  // Zoom the view. Calls the relevant function on the view.
  void ZoomBy ( float iFactor );
  void ZoomIn ();
  void ZoomOut ();

  //BTX
  enum ViewLayout { NoViewLayout, Single, TwoByTwo };
  void SetViewLayout ( ViewLayout iLayout );
  ViewLayout GetViewLayout () const;
  //ETX
  void SetViewLayoutToSingle ();
  void SetViewLayoutToTwoByTwo ();

  // Description:
  // Set the current tool. Puts this tool's controls in the tool
  // page.
  vtkKWScubaTool* GetCurrentTool ();
  void SetCurrentTool ( vtkKWScubaTool& iTool );
  void SetCurrentToolFromMenu ();

  // Description:
  // Set the current view. Puts this views's controls in the tool
  // page.
  void SetCurrentView ( vtkKWScubaView& iView );
  void SetCurrentViewFromMenu ();

  // Description:
  // Set the current layer collection. Puts this layer collection and
  // current layer's controls in the layer page.
  void SetCurrentLayerCollection ( vtkKWScubaLayerCollection& iCol );
  void SetCurrentLayerCollectionFromMenu ();

  // Description:
  // The interactor style will call this function when an event has
  // been handled. The window responds by asking the view if its info
  // has changed, and will call UpdateInfoArea() if so.
  void EventDone ();

  // Description:
  // We listen to the LayerLabelChanged message and update our menu.
  //BTX
  virtual void DoListenToMessage ( std::string const isMessage, 
				   void* const iData );
  //ETX

  // Description:
  // Get and set the flag for resizing the info area.
  int GetAutoSizeInfoArea ();
  void SetAutoSizeInfoArea ( int ibSize );

protected:

  vtkKWScubaWindow ();
  virtual ~vtkKWScubaWindow ();

  // Description:
  // Rebuild our menus with current available items.
  void UpdateToolMenu ();
  void UpdateViewMenu ();
  void UpdateLayerMenu ();

  // Description:
  // Enable or disable buttons and menu items based on program state.
  void UpdateCommandStatus ();

  // Description:
  // Updates the info area, polling views for data to display.
  void UpdateInfoArea ();

  // Description:
  // Set all view settings to a default for the current layout.
  void InitializeViewSettingsForLayout ();

  // Description:
  // Get a view at a slot. Makes and initializes it if it doesn't
  // exist.
  vtkKWScubaView* GetOrMakeNthView ( int inView );
  
  // Description:
  // Adds the LayerCollection to the views and broadcasts the change.
  void AddLayerCollectionToViews ( vtkKWScubaLayerCollection* col );  

  //BTX
  // The views in the view area.
  std::map<int,vtkSmartPointer<vtkKWScubaView> > maView;

   // Inserts a 5 pixel wide spacer into the toolbar.
  static int const kToolbarSpacerWidth;
  void AddSpacerToToolbar ( vtkKWToolbar* iToolbar, 
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

  // Panels.
  vtkSmartPointer<vtkKWFrame> mToolUIFrame;
  vtkSmartPointer<vtkKWFrame> mViewUIFrame;
  vtkSmartPointer<vtkKWFrame> mLayerUIFrame;

  // Menus in panels.
  vtkSmartPointer<vtkKWMenuButton> mMenuTool;
  vtkSmartPointer<vtkKWMenuButton> mMenuView;
  vtkSmartPointer<vtkKWMenuButton> mMenuLayer;

  // Toolbar areas.
  vtkSmartPointer<vtkKWToolbar> mToolbarTools;
  vtkSmartPointer<vtkKWToolbar> mToolbarWindow;
  vtkSmartPointer<vtkKWToolbar> mToolbarView;

  // Radio button sets for toolbars.
  vtkSmartPointer<vtkKWRadioButtonSet> mRadBtnSetViewLayout;
  vtkSmartPointer<vtkKWRadioButtonSet> mRadBtnSetTool;

  // Lookups for menu items to pointers.
  std::map<int,vtkSmartPointer<vtkKWScubaTool> > mToolMenuIndexToPointerMap;
  std::map<int,vtkSmartPointer<vtkKWScubaView> > mViewMenuIndexToPointerMap;
  std::map<int,vtkSmartPointer<vtkKWScubaLayerCollection> > mLayerMenuIndexToPointerMap;

  // Our menu items.
  MenuItem* mMenuLoadVolume;
  MenuItem* mMenuLoadSurface;
  MenuItem* mMenuLoadDTI;
  MenuItem* mMenuLoadPath;
  MenuItem* mMenuLoadODF;
  MenuItem* mMenuSaveVolume;
  MenuItem* mMenuZoomOut;
  MenuItem* mMenuZoomIn;

  // Our control buttons.
  vtkSmartPointer<vtkKWPushButton> mBtnLoadVolume;
  vtkSmartPointer<vtkKWPushButton> mBtnSaveVolume;

  // Our info tables.
  vtkSmartPointer<vtkKWMultiColumnList> mCursorInfoTable;
  vtkSmartPointer<vtkKWMultiColumnList> mMouseOverInfoTable;

  // The tool associated with this window.
  vtkSmartPointer<vtkKWScubaTool> mCurrentTool;
  vtkSmartPointer<vtkKWScubaView> mCurrentView;
  vtkSmartPointer<vtkKWScubaLayerCollection> mCurrentLayerCollection;

  ViewLayout mCurrentViewLayout;

  // Whether or not to automatically size the info area.
  int mbAutoSizeInfoArea;
  
  // default file extension to filter for
  static const std::string sDefaultVolumeFileExtension;

  // Registry keys.
  static const char* sRegistryKey;
  static const char* sAutoSizeInfoAreaKey;

  //ETX

};

#endif
