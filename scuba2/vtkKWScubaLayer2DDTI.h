/**
 * @file  vtkKWScubaLayer2DDTI.h
 * @brief A vtkKWScubaLayer that displays DTI slices
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: dsjen $
 *    $Date: 2007/04/13 21:46:31 $
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

#ifndef vtkKWScubaLayer2DDTI_h
#define vtkKWScubaLayer2DDTI_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesDTI;
//ETX
class vtkActor;
class vtkLODActor;
class vtkFDTensorGlyph;
class vtkImageClip;
class vtkImageReslice;
class vtkImageShrink3D;
class vtkKWCheckButton;
class vtkPolyDataMapper;

class vtkKWScubaLayer2DDTI : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer2DDTI* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer2DDTI, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetDTIProperties ( ScubaCollectionPropertiesDTI* const iProperties );
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
  // Load the surface and set up.
  void Create ();

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

  vtkKWScubaLayer2DDTI ();
  virtual ~vtkKWScubaLayer2DDTI ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Update itself with various settings from the properties.
  void UpdateOpacity ();
  void UpdateDetail();
  void UpdateGlyphScaling();
  void Update2DInfo ();

  // Description:
  // Return a value corresponding to a shrinkage type.
  int GetCurrentShrinkageValue ();

  // Description:
  // Updates the tensor edge display.
  void UpdateEdges();

  //BTX
  ScubaCollectionPropertiesDTI const* mDTIProperties;

  // Pipeline -------------------------------------------------------------
  vtkImageReslice* mVolumeToRAS;
  vtkImageShrink3D* mReducedVolume;
  vtkImageClip* mClip;
  vtkFDTensorGlyph* mGlyph;
  vtkPolyDataMapper* mMapper;
  vtkLODActor* mActor;
  vtkLODActor* mEdgeActor;
  // ----------------------------------------------------------------------

  float mWorldCenter[3];
  float mWorldSize[3];

  //ETX
};

#endif
