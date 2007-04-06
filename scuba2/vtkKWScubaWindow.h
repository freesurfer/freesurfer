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
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:06 $
 *    $Revision: 1.1 $
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

class vtkKWMenu;
class vtkKWPushButton;
class vtkKWRadioButtonSet;
class vtkKWRadioButton;
class vtkKWScaleWithEntry;
class vtkKWMultiColumnList;
class vtkKWScubaLayer;
class vtkKWScubaLayerCollection;
class vtkKWMenuButton;
class vtkKWScubaTool;
class vtkKWToolbar;
class vtkKWFrame;

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
  virtual void Create ();

  // Description:
  // Bring up dialog boxes and then pass the results to LoadVolume,
  // LoadSurface, SaveVolume, or LoadDTI.
  void LoadVolumeFromDlog ();
  void LoadSurfaceFromDlog ();
  void SaveVolumeWithConfirm ();
  void LoadDTIFromDlog ();
  void LoadPathFromDlog ();

  // Description:
  // Does all the UI stuff around loading a volume or surface. Creates
  // a proper layer collection for the data type loaded with the
  // filename we got, and sets the collection in the first unused
  // slot.
  void LoadVolume ( const char* ifnVolume );
  void LoadSurface ( const char* ifnSurface );
  void LoadDTI ( const char* ifnDTI );
  void LoadPath ( const char* ifnPath );

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
  vtkKWScubaView* GetNthView ( int inView );
  
  // Description:
  // Adds the LayerCollection to the views and broadcasts the change.
  void AddLayerCollectionToViews ( vtkKWScubaLayerCollection* col );  

  //BTX
  // The views in the view area.
  std::map<int,vtkKWScubaView*> maView;

  // Struct for associating a menu and an entry item.
  class MenuItem {
  public:
    vtkKWMenu* menu;
    int nItem;
    MenuItem () { menu = NULL; nItem = 0; }
  };

  // Panels.
  vtkKWFrame* mToolUIFrame;
  vtkKWFrame* mViewUIFrame;
  vtkKWFrame* mLayerUIFrame;

  // Menus in panels.
  vtkKWMenuButton* mMenuTool;
  vtkKWMenuButton* mMenuView;
  vtkKWMenuButton* mMenuLayer;

  // Toolbar areas.
  vtkKWToolbar* mToolbarTools;
  vtkKWToolbar* mToolbarWindow;
  vtkKWToolbar* mToolbarView;

  // Radio button set for toolbar tools.
  vtkKWRadioButtonSet* mRadBtnSetTool;

  // Lookups for menu items to pointers.
  std::map<int,vtkKWScubaTool*> mToolMenuIndexToPointerMap;
  std::map<int,vtkKWScubaView*> mViewMenuIndexToPointerMap;
  std::map<int,vtkKWScubaLayerCollection*> mLayerMenuIndexToPointerMap;

  // Our menu items.
  MenuItem mMenuLoadVolume;
  MenuItem mMenuLoadSurface;
  MenuItem mMenuLoadDTI;
  MenuItem mMenuLoadPath;
  MenuItem mMenuSaveVolume;
  MenuItem mMenuZoomOut;
  MenuItem mMenuZoomIn;

  // Our control buttons.
  vtkKWPushButton* mBtnLoadVolume;
  vtkKWPushButton* mBtnSaveVolume;

  // Our info tables.
  vtkKWMultiColumnList* mCursorInfoTable;
  vtkKWMultiColumnList* mMouseOverInfoTable;

  // The tool associated with this window.
  vtkKWScubaTool* mCurrentTool;
  vtkKWScubaView* mCurrentView;
  vtkKWScubaLayerCollection* mCurrentLayerCollection;

  ViewLayout mCurrentViewLayout;
  
  // default file extension to filter for
  static const std::string DEFAULT_VOLUME_FILE_EXTENSION;

  //ETX

};

#endif
