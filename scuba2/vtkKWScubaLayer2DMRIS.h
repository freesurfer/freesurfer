/**
 * @file  vtkKWScubaLayer2DMRIS.h
 * @brief A vtkKWScubaLayer that displays surface slices
 *
 * Displayes surface (from a vtkFSSurfaceSource) in a slice.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.5 $
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
  void SetMRISProperties ( ScubaCollectionPropertiesMRIS const* iProperties );
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

  //ETX
};

#endif
