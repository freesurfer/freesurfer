/**
 * @file  vtkKWScubaLayer3DMRI.h
 * @brief A vtkKWScubaLayer that displayes 3DMRI slices
 *
 * Displayes MRI volumes (from a vtkFSVolumeSource) in three slices.
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
