/**
 * @file  vtkKWScubaLayer3DMRIS.h
 * @brief Displays 3D surfaces
 *
 * Displays surfaces as a 3D model.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
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

#ifndef vtkKWScubaLayer3DMRIS_h
#define vtkKWScubaLayer3DMRIS_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesMRIS;
//ETX
class vtkPolyDataMapper;
class vtkActor;

class vtkKWScubaLayer3DMRIS : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer3DMRIS* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer3DMRIS, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetMRISProperties ( ScubaCollectionPropertiesMRIS* const iProperties );
  //ETX

  // Description:
  // Listens for:
  // OpacityChanged - Changes the mesh's property's opacity.
  // FastModeChanged - Dis/connects the decimated mapper to the actor.
  //BTX
  void DoListenToMessage ( std::string const isMessage, 
			   void* const iData );
  //ETX

  // Description:
  // Load the surface and set up.
  void Create ();

  // Returns the bounds of the surface.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX

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
  vtkPolyDataMapper* mNormalMapper;
  vtkPolyDataMapper* mFastMapper;
  vtkActor* mActor;
  // ----------------------------------------------------------------------

  // Surface info.
  float mRASCenter[3];

  //ETX
};

#endif
