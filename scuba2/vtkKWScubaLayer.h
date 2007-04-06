/**
 * @file  vtkKWScubaLayer.h
 * @brief Base class for layers
 *
 * Provides base class for layers. Subclasses should load data, create
 * vtk objects, use AddProp() to add them to the view, and handle
 * listened to events.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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

#ifndef vtkKWScubaLayer_h
#define vtkKWScubaLayer_h

#include <list>
#include <map>
#include "Broadcaster.h"
#include "IDTracker.h"
#include "Listener.h"
#include "ScubaCollectionProperties.h"
#include "ScubaViewProperties.h"
#include "vtkKWObject.h"

class vtkKWScubaLayerCollection;
class vtkKWWidget;
class vtkProp;
class vtkPropCollection;
//BTX
class ScubaInfoItem;
//ETX

class vtkKWScubaLayer : public vtkKWObject
      //BTX
      , public Broadcaster,
      public Listener,
      public IDTracker<vtkKWScubaLayer>
      //ETX
{

public:

  static vtkKWScubaLayer* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer, vtkKWObject );

  // Description:
  // Set the collection properties. Properties should provide the
  // layer with accessors to data sources, settings, and shared
  // pipeline objects.  
  //BTX
  void SetCollectionProperties ( ScubaCollectionProperties* const iProperties );
  //ETX

  // Description:
  // Set the view properties. Implemented by the view in which this
  // layer is being display, Properties should provide the layer with
  // accessors to current RASZ settings, in plane, etc.
  //BTX
  void SetViewProperties ( ScubaViewProperties* const iProperties );
  //ETX

  // Description:
  // Called after initializtion and after properties are set so that
  // the layer may init itself.
  virtual void Create () {}

  // Description:
  // Implements Listener. This base class doesn't listen for any
  // messages, but this is provided to be overwritten.  
  //BTX
  virtual void DoListenToMessage ( std::string const isMessage,
				   void* const iData ) {}
  //ETX

  // Description:
  // Just calls GetLayer() from the properties.
  const char* GetLabel () const;

  // Description:
  // The prop collection for this layer.
  vtkPropCollection* GetPropCollection () const;

  // Description:
  // Populate a UI page with controls for this layer. Depopulate will
  // be called when the layer is no longer displaying its controls on
  // the panel. Subclasses should override AddControls() and
  // RemoveControls().
  void PopulateControlPage ( vtkKWWidget* iPanel );
  void DepopulateControlPage ();

  // Description:
  // Get the RAS bounds for this layer.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Get a suggested increment for the Z plane movement. This is
  // usually the pixel size in the in plane direction.
  virtual void Get2DRASZIncrementHint ( float ioHint[3]) const;

  // Description:
  // Get the info that should current be displayed in the window's
  // info area for these coords.
  //BTX
  virtual void GetInfoItems ( float iRAS[3],
                              std::list<ScubaInfoItem>& ilInfo ) const {}
  //ETX

  // Description:
  // The view and window will use IsInfoChanged() to see if info has
  // changed and call GetInfoItems to get the new info, and then
  // InfoUpdated() will be called to clear the flag.
  bool IsInfoChanged () const;
  void InfoUpdated ();

  // Description:
  // A static function to get the layer associated with an prop. The
  // view will use this to get the layer associated with a poked
  // prop.
  static vtkKWScubaLayer* GetLayerFromProp ( vtkProp* iProp );

protected:

  // Description:
  // Subclasses can take arguments such as file names and create vtk
  // objects.
  vtkKWScubaLayer ();
  virtual ~vtkKWScubaLayer ();

  // Description:
  // Call to add or remove props created to the internal list. A view
  // will access this list and add all props in it to its render list.
  virtual void AddProp ( vtkProp* iProp );
  virtual void RemoveProp ( vtkProp* iProp );
  virtual void RemoveAllProps ();

  // Description:
  // Subclasses should override these to add and remove their own
  // controls to the layers panel.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Subclasses should call this function if their info has changed.
  void InfoChanged ();

  // Description:
  // Subclasses can call this and we'll broadcast a PipelineChanged
  // message that the view will pickup, and then the view will redraw.
  void PipelineChanged ();

  //BTX

  // Pointer to this layer's shared properties. The layer does not own
  // this object.
  ScubaCollectionProperties const* mProperties;

  // Pointer to this layer's view properties. The layer does not own
  // this object.
  ScubaViewProperties const* mViewProperties;

  // The prop collection containing the props that should be drawn for
  // this layer.
  vtkPropCollection* mProps;

  // Flag for when the info this layer should display in the info area
  // has changed. Polled by view to see when to request new info.
  bool mbInfoChanged;

  // A map for associating props that the tool will hit to layers
  // owning those props. This is a static map and holds references too
  // all the layers' props. Entries are automatically made and deleted
  // when AddProp() and Remove*Prop[s]() are called.
  static std::map<vtkProp*,vtkKWScubaLayer*> mPropToLayerMap;
  //ETX
};

#endif
