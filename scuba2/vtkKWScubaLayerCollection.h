/**
 * @file  vtkKWScubaLayerCollection.h
 * @brief A collection of layers for different vtkKWScubaView display modes
 *
 * A class that collects multiples layers designed to render one kind
 * of data in multiple display modes. The vtkKWScubaView will ask the
 * collection for a layer for a certain display mode, and the
 * collection is responsible for returning the right kind of layer,
 * and instantiating it and setting it up with initial data if
 * necessary. Subclasses should include function necessary to set
 * filenames for data to be loaded, so they can later create instances
 * of layers and load them. The collection can also keep common data 
 * objects that all layers will use.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.3 $
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

#ifndef vtkKWScubaLayerCollection_h
#define vtkKWScubaLayerCollection_h

#include <string>
#include "vtkKWObject.h"
#include "ScubaCollectionProperties.h"
#include "vtkKWScubaView.h"
#include "vtkKWScubaLayer.h"

class vtkKWWidget;
class vtkKWEntry;
class vtkKWScaleWithEntry;
class vtkKWFrame;

class vtkKWScubaLayerCollection : public vtkKWObject
      //BTX
      , public ScubaCollectionProperties,
      public Broadcaster,
      public Listener,
      public IDTracker<vtkKWScubaLayerCollection>
      //ETX
 {

 public:

  static vtkKWScubaLayerCollection* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollection, vtkKWObject );

  // Description:
  // A label for this layer collection, for use in the UI.
  void SetLabel ( const char* isLabel );
  const char* GetLabel () const;

  // Description:
  // Get a layer for the current display mode and view. The layer for
  // this mode will be made the current layer; when the
  // PopulateControlPage, the current layer's PopulateControlPage will
  // be called as well. When the layer is made, this function will set
  // its collection, application, and add this as a listener to it.
  //BTX
  vtkKWScubaLayer* GetCurrentLayer () const;
  void SetDisplayModeAndView ( vtkKWScubaView::DisplayMode iMode,
			       vtkKWScubaView* iView );
  //ETX

  // Description:
  // Populate a UI page with controls for this layer. Subclasses
  // should override AddControls(). This will also call invidiual
  // layers' PopulateControlPage functions.
  void PopulateControlPage ( vtkKWWidget* iPanel );
  void DepopulateControlPage ();

  // Description:
  // We listen to the messages from the layers and rebroadcast them.
  //BTX
  virtual void DoListenToMessage ( std::string const isMessage,
				   void* const iData );
  //ETX

  // Description:
  // Set overall opacity for this layer. Subclasses should draw all
  // their contents at this opacity (they can override
  // OpacityChanged() to respond).
  void SetOpacity ( float iOpacity );
  float GetOpacity () const; // Implements ScubaCollectionProperties
  virtual void OpacityChanged ();

  // Description:
  // Layers should provide fast drawing props for Fast Mode, which
  // will be used when the user is scrolling through the slices or
  // doing something else that requires fast FPS and low
  // precision. Subclasses can override FastModeChanged() to set up
  // their pipelines accordingly.
  void SetFastMode ( bool ibFastMode );
  bool GetFastMode () const; // Implements ScubaCollectionProperties
  virtual void FastModeChanged ();

  // Description:
  // Set properties for the 2D view mode. Subclasses should override
  // Layer2DInfoChanged() to respond.
  void Set2DInPlane ( int iPlane );
  int Get2DInPlane () const; // Implements ScubaCollectionProperties
  void Set2DRASZ ( float iRASZ );
  float Get2DRASZ () const; // Implements ScubaCollectionProperties
  virtual void Layer2DInfoChanged ();

  // Description:
  // Set properties for the 3D view mode. Subclasses should override
  // Layer3DInfoChanged() to respond.
  void Set3DRASX ( float i3DRASX );
  float Get3DRASX () const; // Implements ScubaCollectionProperties
  void Set3DRASY ( float i3DRASY );
  float Get3DRASY () const; // Implements ScubaCollectionProperties
  void Set3DRASZ ( float i3DRASZ );
  float Get3DRASZ () const; // Implements ScubaCollectionProperties
  virtual void Layer3DInfoChanged ();

 protected:

  vtkKWScubaLayerCollection ();
  ~vtkKWScubaLayerCollection ();

  // Description:
  // Subclasses should override these to add and remove their own
  // controls to the layer control panel.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // The subclass should override to create and setup a layer of the
  // proper type. It should set any special properties in the layer if
  // necessary.
  //BTX
  virtual vtkKWScubaLayer* 
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

 private:

  //BTX
  std::string msLabel;

  // Shared state ---------------------------------------------------------
  // Opacity for the layer.
  float mOpacity;

  // If true, draw stuff fast and imprecise, otherwise draw normally.
  bool mbFastMode;

  // The orientation and location of the current 'slice' in 2D mode.
  float m2DRASZ;
  int mInPlane;

  // The 3D in planes
  float m3DRASX;
  float m3DRASY;
  float m3DRASZ;
  // ----------------------------------------------------------------------

  // Control widgets ------------------------------------------------------
  vtkKWEntry* mEntryLabel;
  vtkKWScaleWithEntry* mScaleOpacity;
  vtkKWFrame* mLayerSettingsPanel;
  // ----------------------------------------------------------------------

  // Our display mode and view.
  vtkKWScubaView::DisplayMode mDisplayMode;
  vtkKWScubaView* mCurrentView;

  // This is a map of layers for each display mode. The first time the
  // display mode is set, a layer will be crated for that mode.
  typedef std::map<vtkKWScubaView*,vtkKWScubaLayer*> ViewLayerMapType;
  typedef std::map<vtkKWScubaView::DisplayMode,ViewLayerMapType> 
    ModeViewLayerMapType;
  ModeViewLayerMapType maLayers;
  //ETX
  
};

#endif
