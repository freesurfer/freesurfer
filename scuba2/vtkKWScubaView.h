/**
 * @file  vtkKWScubaView.h
 * @brief A VTK view containing ScubaLayerCollections
 *
 * Implementation of a vtkKWRenderWidget that can have multiple
 * ScubaLayerCollections. Also handles passing events down to them and
 * managing their viewing state.
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

#ifndef vtkKWScubaView_h
#define vtkKWScubaView_h

#include <map>
#include <list>
#include "Broadcaster.h"
#include "IDTracker.h"
#include "Listener.h"
#include "ScubaInfoItem.h"
#include "ScubaViewProperties.h"
#include "vtkKWRenderWidget.h"
#include "vtkSmartPointer.h"

class vtkKWScubaWindow;
class vtkKWScubaLayer;
class vtkKWScubaLayerCollection;
class vtkKWScubaTool;
class vtkKWComboBox;
class vtkKWRadioButton;
class vtkKWRadioButtonSet;
class vtkCamera;
class vtkKWScaleWithEntry;
class vtkKWToolbar;
class vtkKWPushButton;

class vtkKWScubaView : public vtkKWRenderWidget
      //BTX
      , public ScubaViewProperties
      , public IDTracker<vtkKWScubaView>
      , public Listener
      , public Broadcaster
      //ETX
{

public:

  static vtkKWScubaView* New ();
  vtkTypeRevisionMacro( vtkKWScubaView, vtkKWRenderWidget );

  // Description:
  // Create the widget.
  void CreateWidget ();

  // Description:
  // This should be the label that will be used in the GUI.
  virtual void SetLabel ( const char* isLabel );
  virtual const char* GetLabel () const;

  // Description: 
  // Populate a toolbar with controls for this view. Depopulate will
  // be called when the view is no longer displaying its controls.
  virtual void PopulateToolbar ( vtkKWToolbar* iToolbar );
  virtual void DepopulateToolbar ();

  // Description: 
  // Populate a UI page with controls for this view. Depopulate will
  // be called when the view is no longer displaying its controls on
  // the panel.
  virtual void PopulateControlPage ( vtkKWWidget* iPanel );
  virtual void DepopulateControlPage ();

  // Description: 
  // Convenience function to check the collection map if a collection
  // exists at a slot. Just returns a pointer to the collection or
  // NULL, no error is raised if the slot is out of bounds.
  vtkKWScubaLayerCollection* GetCollectionAtSlot ( int inSlot ) const;

  // Description:
  // All layer collections must be set at specific slots. 0 is the
  // 'bottom' slot. Higher numbered slots are drawn on top of lower
  // numbered slots.
  int GetFirstUnusedLayerSlot () const;
  int GetHighestFilledLayerSlot () const;
  void SetLayerCollectionAtSlot ( int inSlot, 
				  vtkKWScubaLayerCollection* iCollection );
  void SetLayerCollectionAtSlotByComboBoxValue ( int inSlot,
						 const char* isValue );

  // Description:
  // Returns the slot of this layer or collection, or -1 if the layer
  // is not in this view.
  int GetSlotOfLayer ( const vtkKWScubaLayer* iLayer ) const;
  int GetSlotOfLayerCollection ( const vtkKWScubaLayerCollection* iCollection ) const;

  // Description:

  void MakeMenuForSlot ( int inSlot );

  // Description:
  // Called when the layers available have changed.
  void LayerListChanged ();

  // Description: 
  // Reset all cameras in all display modes. Useful when a new piece
  // of data is loaded.
  void ResetAllCameras ();

  // Description:
  // Adjust the currently displayed 2DRASZ or 3DRAS{X,Y,Z} slider ranges.
  void AdjustViewPlaneSliders ();

  // Description:
  // Set the display mode for this view. This will get the right layer
  // from the layer collection in each slot and display it.
  //BTX
  enum DisplayMode { NoDisplayMode, TwoDee, ThreeDee };
  void SetDisplayMode ( DisplayMode iMode );
  DisplayMode GetDisplayMode () const;
  //ETX
  virtual void DisplayModeChanged ();
  void SetDisplayModeTo2D ();
  void SetDisplayModeTo3D ();

  // Description:
  // The tool will notify the view when a view update operation is
  // beginning so it can adjust layers' view settings accordingly (set
  // them to FastMode).
  void StartFastMode ();
  void StopFastMode ();
  bool GetFastMode () const;

  // Description:
  // These functions change the camera and RASZ for when the user
  // clicks and drags between points on the screen. A reasonable
  // movement vector is made, and a delta is calculated and applied.
  void PanBetweenWindowCoords ( int iDelta[2] );
  void ScrollBetweenWindowCoords ( int iDelta[2] );
  void ZoomBetweenWindowCoords ( int iDelta[2] );
  void RotateBetweenWindowCoords ( int iDelta[2] );

  // Description:
  // Polled to determine if the view contains layers that have data in
  // a modified and can be saved.
  bool IsDataDirty () const;

  // Description:
  // Get and set properties for the 2D view mode. InPlane determines
  // which plane is being shown in cross section. is 0 for x, 1 for y,
  // and 2 for z.
  void Set2DInPlane ( int iPlane );
  int Get2DInPlane () const;

  // Description:
  // Get and set properties for the 2D view mode. RASZ is the current
  // slice in the in plane we're showing.
  void Set2DRASZ ( float iRASZ );
  float Get2DRASZ () const;

  // Description:
  // The zoom level is an arbitrary number. Zoom level 1 is the level
  // at which all layers' content fill the view.
  void Set2DZoomLevel ( float iZoomLevel );
  float Get2DZoomLevel ();
  void ZoomOut ();
  void ZoomIn ();

  // Description: 
  // Rotate commands. These do a cool little animation of rotating the
  // camera arond the object.
  void AnimateCameraElevateNegative ();
  void AnimateCameraElevatePositive ();
  void AnimateCameraAzimuthNegative ();
  void AnimateCameraAzimuthPositive ();
  void AnimateCameraRollNegative ();
  void AnimateCameraRollPositive ();

  // Description:
  // This returns the range of z values with content. This is union of
  // the RAS bounds of all layers for the current in plane.
  void Get2DRASZRange ( float& oMin, float& oMax ) const;

  // Description:
  // A hint as to what Z increment to use, for example, when clicking
  // a "next slice" button, for the current in plane.
  float Get2DRASZIncrementHint () const;
 
  // Description:
  // This returns the range of z values with content for all
  // orientations.
  void Get3DRASRange ( float oMin[3], float oMax[3] ) const;

  // Description:
  // A hint as to what Z increments to use, for example, when clicking
  // a "next slice" button, for all orientations.
  void Get3DRASIncrementHint ( float oInc[3] ) const;
 
  // Description:
  // Get and set properties for the 3D view mode. These are used only
  // by layers that display volumes in slices in 3D. Changes are
  // broadcast with Layer3DInfoChanged.
  void Set3DRASX ( float i3DRASX );
  float Get3DRASX () const;
  void Set3DRASY ( float i3DRASY );
  float Get3DRASY () const;
  void Set3DRASZ ( float i3DRASZ );
  float Get3DRASZ () const;
  void Get3DRAS ( float& o3DRASX, float& o3DRASY, float& o3DRASZ ) const;

  // Description:
  // Event callbacks. iCoords is in window coords and can be
  // transformed by the view into world coords.
  void MouseMoveEvent        ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void LeftButtonDownEvent   ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void LeftButtonUpEvent     ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void MiddleButtonDownEvent ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void MiddleButtonUpEvent   ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void RightButtonDownEvent  ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void RightButtonUpEvent    ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void KeyDownEvent          ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void KeyUpEvent            ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void EnterEvent            ( vtkKWScubaWindow* iWindow, int iCoords[2] );
  void LeaveEvent            ( vtkKWScubaWindow* iWindow, int iCoords[2] );

  // Description:
  // World coords are the one VTK will use, and will have a 0 z
  // component. This changes the z coordinate to be the RASZ we're on.
  //BTX
  void ConvertWorldCoordsToRAS ( double const iWorldCoords[3],
                                 float  oRAS[3] ) const;
  //ETX

  // Description:
  // Set the cursor for this view.
  //BTX
  void SetCursorFromWorld ( double const iWorldCoords[3] );
  void SetCursorFromRAS   ( float const iRASCoords[3] );
  //ETX

  // Description:
  // Get the info that should current be displayed in the window's
  // info area. This will also be called on all layers.
  //BTX
  void GetMouseOverInfoItems ( std::list<ScubaInfoItem>& ilInfo );
  void GetCursorInfoItems ( std::list<ScubaInfoItem>& ilInfo );
  //ETX

  // Description:
  // The window will call this to determine if the info in the window
  // containing the view should be changed. The view will set this
  // base on its own internal info and by polling all the layers to
  // see if their info has changed.
  bool IsInfoChanged ();
  void InfoUpdated ();

  // Description:
  // PipelineChanged: Render the window
  // LayerLabelChanged: Pass message to window
  // PropListChanged: Rebuild the prop listo
  //BTX
  virtual void DoListenToMessage ( std::string const isMessage, 
				   void* const iData );
  //ETX

  // Description:
  // A way to find an instance of a ScubaView from a vtkRenderWindow.
  static vtkKWScubaView* GetViewFromRenderWindow( vtkRenderWindow* );

protected:

  vtkKWScubaView ();
  virtual ~vtkKWScubaView ();

  // Description:
  // Packs the controls for the current display mode.
  void PackDisplayModeControls ();

  // Description:
  // Fills our slot/layer combo boxes with the current layer names.
  void UpdateLayerComboBoxes ();

  //BTX

  // Description:
  // Provides the basic template code for event handling. Gets the
  // picker, picks a position at the given window coords, finds the
  // actor, then the layer, gets the current tool, and calls the
  // passed event function on the tool. It will set the world coords
  // calculated at the pick point in the passed in float array if not
  // NULL.
  void PickPointAndCallToolEvent ( vtkKWScubaWindow* iWindow, int iCoords[2],
                                   void(vtkKWScubaTool::*iToolFunc)(vtkKWScubaWindow*,vtkKWScubaView*,vtkKWScubaLayer*,float[3]),
                                   double* ioWorldCoords );

  // Description:
  // Generate ScubaInfoItems for a particular RAS coord, and poll all
  // layers for the same.
  void GetInfoItemsAtRAS ( float iRAS[3],
                           std::list<ScubaInfoItem>& ilInfo );

  // Description:
  // Call this function if info has changed.
  void InfoChanged ();

  // Description:
  // Populates the view with props. First removes all existing props
  // and then goes through the layers in order and adds their props.
  void AddLayerPropsToView ();

  //BTX
  // Display mode buttons.
  vtkSmartPointer<vtkKWRadioButton> mRadBtnDisplayMode2D;
  vtkSmartPointer<vtkKWRadioButton> mRadBtnDisplayMode3D;
  
  // Frame for display mode controls.
  vtkSmartPointer<vtkKWFrame> mFrameDisplayModeControls;

  // 2D: In plane radio buttons.
  vtkSmartPointer<vtkKWRadioButtonSet> mRadBtnSet2DInPlane;
  vtkSmartPointer<vtkKWRadioButton> mRadBtn2DInPlaneX;
  vtkSmartPointer<vtkKWRadioButton> mRadBtn2DInPlaneY;
  vtkSmartPointer<vtkKWRadioButton> mRadBtn2DInPlaneZ;

  // 2D: Z scale.
  vtkSmartPointer<vtkKWScaleWithEntry> mScale2DRASZ;

  // 3D: scales for RAS slices.
  vtkSmartPointer<vtkKWScaleWithEntry> mScale3DRASX;
  vtkSmartPointer<vtkKWScaleWithEntry> mScale3DRASY;
  vtkSmartPointer<vtkKWScaleWithEntry> mScale3DRASZ;
  
  // Frame for packing slot menus.
  vtkSmartPointer<vtkKWFrame> mFrameSlotMenus;

  // Zoom and rotate buttons for toolbar.
  vtkSmartPointer<vtkKWPushButton> mBtnZoomOut;
  vtkSmartPointer<vtkKWPushButton> mBtnZoomIn;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateXPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateXNeg;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateYPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateYNeg;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateZPos;
  vtkSmartPointer<vtkKWPushButton> mBtnRotateZNeg;

  std::string msLabel;

  bool mbDirty;
  bool mbInfoChanged;

  // The collection currently loaded into each slot.
  typedef std::map<int,vtkSmartPointer<vtkKWScubaLayerCollection> >
    SlotCollectionMapType;
  SlotCollectionMapType maCol;

  DisplayMode mDisplayMode;

  vtkSmartPointer<vtkCamera> m2DCamera;
  vtkSmartPointer<vtkCamera> m3DCamera;

  bool mbFastMode;

  // Our 2D info.
  float m2DRASZ;
  int m2DInPlane;

  // Our 3D info.
  float m3DRASX, m3DRASY, m3DRASZ;

  // The current mouseover and cursor coords.
  float mMouseOverRASCoords[3];
  float mCursorRASCoords[3];

  // Used by GetViewFromRenderWindow to map scuba views to vtk render
  // windows. 
  static std::map<vtkRenderWindow*,vtkSmartPointer<vtkKWScubaView> > 
    mRenderWindowToViewMap;

  // A pointer to a menu for each slot. The menu contains the labels
  // for all the layer collections that are currently loaded, and
  // currently displays the collection loaded at that slot.
  typedef std::map<int,vtkSmartPointer<vtkKWComboBox> > SlotMenuMapType;
  SlotMenuMapType maSlotMenu;

  // A pointer to a collection for a menu index. The menus in
  // maSlotMenus use this map to go from an entry index to an actual
  // collection.
  typedef std::map<int,vtkSmartPointer<vtkKWScubaLayerCollection> >
    MenuIndexCollectionMapType;
  MenuIndexCollectionMapType maMenuIndexCollection;
  //ETX
};

#endif
