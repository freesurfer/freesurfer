/**
 * @file  vtkKWOrientMRIView3D.h
 * @brief Viewing pipeline for 3D volume and axes
 *
 * Pipeline objects for displaying the volume and axes in a 3D
 * view. This view is the one in which the user can rotate the
 * bounding box of the volume to reorient it freehand.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.4 $
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


#ifndef vtkKWOrientMRIView3D_h
#define vtkKWOrientMRIView3D_h

#include "vtkKWRenderWidget.h"
#include "vtkSmartPointer.h"

class vtkActor;
class vtkAxes;
class vtkCubeSource;
class vtkFollower;
class vtkFSVolumeSource;
class vtkImageData;
class vtkImagePlaneWidget;
class vtkOutlineFilter;
class vtkLookupTable;
class vtkMatrix4x4;
class vtkVectorText;

class vtkKWOrientMRIView3D : public vtkKWRenderWidget {

public:

  static vtkKWOrientMRIView3D* New ();
  vtkTypeRevisionMacro( vtkKWOrientMRIView3D, vtkKWRenderWidget );

  // Description:
  // Create our axes actor.
  void Create ();

  // Description:
  // Set the fs volume, image and color function. Setting the volume
  // will cause it to be drawn in the view.
  void SetFSVolumeAndImage ( vtkFSVolumeSource& iVolume,vtkImageData& iImage );
  void SetImageColors ( vtkLookupTable& iColors );

  // Description:
  // Go back to the original camera position.
  void RestoreView ();

  // Description:
  // Zoom in the view.
  void ZoomBy ( float iFactor );

  // Description: 
  // Rotate the cube (modifying the user transform) a number of
  // degrees around the axis. Does a cool little animation of rotating
  // the cube and updating all the views. Use 0 for X, 1 for Y, 2 for
  // Z.
  //BTX
  enum { kcRotateStepsToAnimate = 10 };
  //ETX
  void RotateUserTransform ( int iAxis, double iDegrees );

  // Description:
  // Get the user transform, maintained by the view here.
  vtkMatrix4x4& GetUserTransform ();

  //BTX
  // Description:
  // Compose the matrix with the current user transform.
  void ComposeWithUserTransform ( vtkMatrix4x4 const& iMatrix );
  //ETX
  
  // Description:
  // Called when our volume actor changes.
  void ActorPositionChanged ();

  // Description: 
  // Called when the volume to RAS transform changes. We
  // rebuild our cube source with the new boundaries and reset our
  // cube actor's rotation (which updates the user transform).
  void VolumeToRASTransformChanged ();

  // Description:
  // Our callback function for actor modified events.
  static void ActorModifiedCallback ( vtkObject* iCaller,
				      unsigned long iEventId,
				      void* iClientData, void* iCallData );

  // Description:
  // Our callback function for actor modified events.
  static void VolumeToRASTransformChangedCallback ( vtkObject* iCaller,
						    unsigned long iEventId,
				       void* iClientData, void* iCallData );

protected:

  vtkKWOrientMRIView3D ();
  virtual ~vtkKWOrientMRIView3D ();

  // Description:
  // Scale our axes and axes labels.
  void ScaleAxesAndLabels ( double iScale );

  //BTX
  // Set by the window; we don't own these.
  vtkSmartPointer<vtkFSVolumeSource> mFSVolume;
  vtkSmartPointer<vtkImageData> mImage;
  vtkSmartPointer<vtkLookupTable> mImageColors;

  // Our pipeline objects.
  vtkSmartPointer<vtkAxes> mAxes;
  vtkSmartPointer<vtkOutlineFilter> mOutlineFilter;
  vtkSmartPointer<vtkActor> mOutlineActor;
  vtkSmartPointer<vtkCubeSource> mCubeSource;
  vtkSmartPointer<vtkActor> mVolumeActor;
  vtkSmartPointer<vtkVectorText> mTextSource[3];
  vtkSmartPointer<vtkFollower> mTextActor[3];

  // Our user transform.
  vtkSmartPointer<vtkMatrix4x4> mUserTransform;

  // Default camera location. We go back here when RestoreView is
  // called.
  double mDefaultPosition[3];
  double mDefaultFocalPoint[3];
  double mDefaultViewUp[3];
  //ETX
};

#endif
