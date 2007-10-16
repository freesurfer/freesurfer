/**
 * @file  vtkKWScubaLayer3DMRI.h
 * @brief A vtkKWScubaLayer that displayes 3DMRI slices
 *
 * Displayes MRI volumes (from a vtkFSVolumeSource) in three slices.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 16:01:50 $
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

#ifndef vtkKWScubaLayer3DMRI_h
#define vtkKWScubaLayer3DMRI_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesMRI;
//ETX
class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkKWScaleWithEntry;
class vtkFSVolumeSource;

class vtkKWScubaLayer3DMRI : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer3DMRI* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer3DMRI, vtkKWScubaLayer );

  // Description:
  // Set the properties item. Done by the collection.
  //BTX
  void SetMRIProperties ( ScubaCollectionPropertiesMRI const* iProperties );
  //ETX

  // Description:
  // Make and set up our pipeline.
  void Create();

  // Description:
  // Listens for:
  // OpacityChanged - 
  // ColorMapChanged - 
  // ResliceInterpolationChanged -
  // TextureSmoothingChanged -
  //BTX
  virtual void DoListenToMessage ( std::string const isMessage,
				   void* const iData );
  //ETX

  // Returns the bounds of the volume.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX

  // Description:
  // Get a pointer to the source for manipulation.
  vtkFSVolumeSource* GetSource () const;

protected:

  vtkKWScubaLayer3DMRI ();
  virtual ~vtkKWScubaLayer3DMRI ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Update itself with various settings from the properties.
  void UpdateOpacity ();
  void UpdateColorMap ();
  void UpdateResliceInterpolation ();
  void UpdateTextureSmoothing ();
  void UpdatePlanes ();

  //BTX
  ScubaCollectionPropertiesMRI const* mMRIProperties;

  // Pipeline ------------------------------------------------------------
  vtkImageReslice* mReslice[3];
  vtkImageMapToColors* mColorMap[3];
  vtkTransform* mPlaneTransform[3];
  vtkTexture* mTexture[3];
  vtkPolyDataMapper* mPlaneMapper[3];
  vtkActor* mPlaneActor[3];
  // ---------------------------------------------------------------------

  float mWorldCenter[3];
  float mWorldSize[3];
  //ETX
};

#endif
