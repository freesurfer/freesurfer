/**
 * @file  vtkKWScubaLayer2DMRIS.h
 * @brief A vtkKWScubaLayer that displays surface slices
 *
 * Displayes surface (from a vtkFSSurfaceSource) in a slice.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/07/25 19:53:47 $
 *    $Revision: 1.2 $
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

#ifndef vtkKWScubaLayer2DMRIS_h
#define vtkKWScubaLayer2DMRIS_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "vtkSmartPointer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesMRIS;
//ETX
class vtkPlane;
class vtkPolyDataMapper;
class vtkActor;

class vtkKWScubaLayer2DMRIS : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer2DMRIS* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer2DMRIS, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetMRISProperties ( ScubaCollectionPropertiesMRIS* const iProperties );
  //ETX

  // Description:
  // Listens for:
  // OpacityChanged - Changes the mesh's property's opacity.
  // FastModeChanged - Dis/connects the decimated mapper to the actor.
  // Layer2DInfoChanged - Changes the slice plane origin and normal,
  // and adjusts the actor's position.  
  //BTX
  void DoListenToMessage ( std::string const isMessage,
			   void* const iData );
  //ETX

  // Description:
  void Create ();

  // Description:
  // Load the surface and set up.
  void LoadDataFromProperties ();
  
  // Returns the bounds of the surface.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Returns the shortest edge distance.
  virtual void Get2DRASZIncrementHint ( float ioHint[3]) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX

protected:

  vtkKWScubaLayer2DMRIS ();
  virtual ~vtkKWScubaLayer2DMRIS ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  //BTX
  ScubaCollectionPropertiesMRIS const* mMRISProperties;

  // Pipeline -------------------------------------------------------------
  vtkSmartPointer<vtkPlane> mSlicePlane;
  vtkSmartPointer<vtkPolyDataMapper> mNormalMapper;
  vtkSmartPointer<vtkPolyDataMapper> mFastMapper;
  vtkSmartPointer<vtkActor> mActor;
  // ----------------------------------------------------------------------

  // Surface info.
  float mRASCenter[3];

  //ETX
};

#endif
