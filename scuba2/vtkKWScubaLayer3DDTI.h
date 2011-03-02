/**
 * @file  vtkKWScubaLayer3DDTI.h
 * @brief A vtkKWScubaLayer that displays DTI glyphs
 *
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

#ifndef vtkKWScubaLayer3DDTI_h
#define vtkKWScubaLayer3DDTI_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesDTI;
//ETX
class vtkFDTensorGlyph;
class vtkPolyDataMapper;
class vtkTransform;
class vtkImageReslice;
class vtkImageResample;
class vtkLODActor;
class vtkActor;
class vtkKWCheckButton;
class vtkMatrix4x4;

class vtkKWScubaLayer3DDTI : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer3DDTI* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer3DDTI, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetDTIProperties ( ScubaCollectionPropertiesDTI const* iProperties );
  //ETX

  // Description:
  // Listens for:
  // OpacityChanged - Changes the mesh's property's opacity.
  // FastModeChanged - Dis/connects the decimated mapper to the actor.
  // Layer3DInfoChanged - Changes the slice plane origin and normal,
  // and adjusts the actor's position.  
  //BTX
  void DoListenToMessage ( std::string const isMessage,
			   void* const iData );
  //ETX

  // Description:
  // Load the volume and set up.
  void Create ();

  // Returns the bounds of the volume.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Returns the shortest edge distance.
  virtual void Get3DRASZIncrementHint ( float ioHint[3]) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX
  
  // Description:
  // Changes whether the tensor on a plane is displayed.
  void SetIsPlaneXVisible( int ibIsVisible );
  void SetIsPlaneYVisible( int ibIsVisible );
  void SetIsPlaneZVisible( int ibIsVisible );
  
  // Description:
  // Callback for when displaying all tensors in the volume, rather than slice.
  void SetIsShowingAll( int ibIsShowingAll );
  

protected:

  vtkKWScubaLayer3DDTI ();
  virtual ~vtkKWScubaLayer3DDTI ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();
  
  void UpdatePlanes();
  void UpdatePlanes( bool hasDetailUpdated );
      
  //BTX
  ScubaCollectionPropertiesDTI const* mDTIProperties;

  // Pipeline -------------------------------------------------------------
  vtkImageReslice* mVolumeToRAS;
  vtkImageReslice* mVolumeToRASSlice[ 3 ];
  vtkTransform* mSliceTransform[3];
  vtkFDTensorGlyph* mGlyphs[ 3 ];
  vtkPolyDataMapper* mPlaneMappers[ 3 ];
  vtkLODActor* mPlaneActors[ 3 ];
  vtkLODActor* mGlyphEdgeActors[ 3 ];
  // ----------------------------------------------------------------------

  float mWorldCenter[3];
  float mWorldSize[3];

  vtkKWCheckButton* mIsPlaneXVisbleButton;
  vtkKWCheckButton* mIsPlaneYVisbleButton;
  vtkKWCheckButton* mIsPlaneZVisbleButton;
  
  vtkKWCheckButton* mIsShowingAllButton;
  
  enum {
    X_PLANE = 0,
    Y_PLANE = 1,
    Z_PLANE = 2
  };
  
  //ETX

private:

  void UpdateGlyphScaling () ;

  void UpdateDetail () ;
  
  int GetCurrentShrinkage ();
  
  void UpdateEdges ();
  
  bool mIsPlaneXVisible;
  bool mIsPlaneYVisible;
  bool mIsPlaneZVisible;
  bool mIsShowingAll;
};

#endif
