/**
 * @file  vtkKWScubaLayer3DMRIS.h
 * @brief Displays 3D surfaces
 *
 * Displays surfaces as a 3D model.
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

#ifndef vtkKWScubaLayer3DMRIS_h
#define vtkKWScubaLayer3DMRIS_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "vtkSmartPointer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesMRIS;
//ETX
class vtkActor;
class vtkPlane;
class vtkPolyDataMapper;

class vtkKWScubaLayer3DMRIS : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer3DMRIS* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer3DMRIS, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetMRISProperties ( ScubaCollectionPropertiesMRIS const* iProperties );
  //ETX

  // Description:
  // Listens for:
  // OpacityChanged - Changes the mesh's property's opacity.
  // FastModeChanged - Dis/connects the decimated mapper to the actor.
  // Layer3DInfoChanged - Change the cut planes
  //BTX
  void DoListenToMessage ( std::string const isMessage, 
			   void* const iData );
  //ETX

  // Description:
  // Load the surface and set up.
  void Create ();

  // Description
  // Make and set up our pipeline.
  void LoadDataFromProperties ();

  // Returns the bounds of the surface.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX

  // Description:
  // Show or hide the 3D surface.
  void SetShowSurface ( int ibShow );

protected:

  vtkKWScubaLayer3DMRIS ();
  virtual ~vtkKWScubaLayer3DMRIS ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  //BTX
  ScubaCollectionPropertiesMRIS const* mMRISProperties;

  // Pipeline -------------------------------------------------------------
  vtkSmartPointer<vtkPolyDataMapper> m3DNormalMapper;
  vtkSmartPointer<vtkPolyDataMapper> m3DFastMapper;
  vtkSmartPointer<vtkActor> m3DActor;
  vtkSmartPointer<vtkPlane> m2DSlicePlane[3];
  vtkSmartPointer<vtkPolyDataMapper> m2DNormalMapper[3];
  vtkSmartPointer<vtkPolyDataMapper> m2DFastMapper[3];
  vtkSmartPointer<vtkActor> m2DActor[3];
  // ----------------------------------------------------------------------

  // Show the actual 3D surface.
  bool mbShowSurface;

  //ETX
};

#endif
