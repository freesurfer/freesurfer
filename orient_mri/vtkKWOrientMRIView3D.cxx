/**
 * @file  vtkKWOrientMRIView3D.cxx
 * @brief Viewing pipeline for 3D volume and axes
 *
 * Pipeline objects for displaying the volume and axes in a 3D
 * view. This view is the one in which the user can rotate the
 * bounding box of the volume to reorient it freehand.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/09/18 17:10:22 $
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


#include <sstream>
#include <stdexcept>
#include "vtkKWOrientMRIView3D.h"

#include "OrientMRIEvents.h"
#include "vtkActor.h"
#include "vtkAxes.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCornerAnnotation.h"
#include "vtkCubeSource.h"
#include "vtkFSVolumeSource.h"
// These are defined in the mri headers somewhere.
#ifdef X
#undef X
#endif
#ifdef Y
#undef Y
#endif
#include "vtkImageMapToColors.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkOrientMRIInteractorStyleView3D.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkTubeFilter.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIView3D );
vtkCxxRevisionMacro( vtkKWOrientMRIView3D, "$Revision: 1.2 $" );


vtkKWOrientMRIView3D::vtkKWOrientMRIView3D () {
  mUserTransform = vtkSmartPointer<vtkMatrix4x4>::New();
}

vtkKWOrientMRIView3D::~vtkKWOrientMRIView3D () {
  
}

void
vtkKWOrientMRIView3D::Create () {

  this->Superclass::Create();

  // Create our axes pipeline objects. We'll adjust the size when we
  // get a volume to display.
  mAxes = vtkSmartPointer<vtkAxes>::New();
  mAxes->SetOrigin( 0, 0, 0 );
  mAxes->SetScaleFactor( 100 );

  vtkSmartPointer<vtkTubeFilter> axesTubes =
    vtkSmartPointer<vtkTubeFilter>::New();
  axesTubes->SetInputConnection( mAxes->GetOutputPort() );
  axesTubes->SetRadius( 1 );
  axesTubes->SetNumberOfSides( 6 );

  vtkSmartPointer<vtkPolyDataMapper> axesMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  axesMapper->SetInputConnection( axesTubes->GetOutputPort() );

  vtkSmartPointer<vtkActor> axesActor =
    vtkSmartPointer<vtkActor>::New();
  axesActor->SetMapper( axesMapper );
  axesActor->SetPickable( false );

  this->GetRenderer()->AddActor( (vtkActor*)axesActor );

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetRenderWindow()->Render();

  // Set our interactor style.
  vtkSmartPointer<vtkOrientMRIInteractorStyleView3D> style = 
    vtkSmartPointer<vtkOrientMRIInteractorStyleView3D>::New();
  this->GetRenderWindow()->GetInteractor()->SetInteractorStyle( style );

}

void
vtkKWOrientMRIView3D::SetFSVolumeAndImage ( vtkFSVolumeSource& iVolume,
					    vtkImageData& iImage ) {

  // Save a reference.
  mFSVolume = &iVolume;
  mImage = &iImage;

  // Remove our props if we have them.
  if( mOutlineActor.GetPointer() )
    this->RemoveViewProp( mOutlineActor );
  if( mVolumeActor.GetPointer() )
    this->RemoveViewProp( mVolumeActor );

  // start with ID user transform.
  vtkSmartPointer<vtkMatrix4x4> identity =
    vtkSmartPointer<vtkMatrix4x4>::New();
  mUserTransform->DeepCopy( identity );

  // Volume outline.
  mOutlineFilter = vtkSmartPointer<vtkOutlineFilter>::New();
  mOutlineFilter->SetInput( mImage );

  vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInputConnection( mOutlineFilter->GetOutputPort() );

  mOutlineActor = vtkSmartPointer<vtkActor>::New();
  mOutlineActor->SetMapper( outlineMapper );

  this->GetRenderer()->AddActor( mOutlineActor );

  // Translucent volume.
  mCubeSource = vtkSmartPointer<vtkCubeSource>::New();
  mCubeSource->SetBounds( mImage->GetBounds() );

  vtkSmartPointer<vtkPolyDataMapper> cubeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  cubeMapper->SetInputConnection( mCubeSource->GetOutputPort() );

  mVolumeActor = vtkSmartPointer<vtkActor>::New();
  mVolumeActor->SetMapper( cubeMapper );
  mVolumeActor->GetProperty()->SetOpacity( 0.3 );
  mVolumeActor->GetProperty()->SetColor( 0, 0, 1 );

  this->GetRenderer()->AddActor( mVolumeActor );

  // Make a callback and observe our outline actor for Modified events.
  vtkSmartPointer<vtkCallbackCommand> callback =
    vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetCallback( ActorModifiedCallback );
  callback->SetClientData( this );
  mVolumeActor->AddObserver( vtkCommand::ModifiedEvent, callback );

  // That volume is our interactor style prop.
  vtkOrientMRIInteractorStyleView3D* style =
    dynamic_cast<vtkOrientMRIInteractorStyleView3D*>( this->GetRenderWindow()->GetInteractor()->GetInteractorStyle() );
  if( style ) {
    style->SetInteractionProp( mVolumeActor );
  }
  
  // Scale our axes to the volume. 0.6 because we only need half,
  // starting at the middle of the volume, but we want it to extend a
  // bit past that.
  int dimensions[3];
  mImage->GetDimensions( dimensions );
  mAxes->SetScaleFactor(MAX(dimensions[0]*0.6, MAX(dimensions[1]*0.6,
						   dimensions[2]*0.6)));

  // Initial camera position.
  this->GetRenderer()->GetActiveCamera()->SetPosition( 0, 1, 0 );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( 0, 0, 0 );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( 0, 0, 1 );

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetRenderWindow()->Render();

  // Grab the starting camera information so we can restore it later.
  this->GetRenderer()->GetActiveCamera()->GetPosition( mDefaultPosition );
  this->GetRenderer()->GetActiveCamera()->GetFocalPoint( mDefaultFocalPoint );
  this->GetRenderer()->GetActiveCamera()->GetViewUp( mDefaultViewUp );

}

void
vtkKWOrientMRIView3D::SetImageColors ( vtkLookupTable& iColors ) {

  // Save a pointer.
  mImageColors = &iColors;
}

void
vtkKWOrientMRIView3D::RestoreView () {

  // Set the camera position to our original settings.
  this->GetRenderer()->GetActiveCamera()->SetPosition( mDefaultPosition );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( mDefaultFocalPoint );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( mDefaultViewUp );

  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetRenderWindow()->Render();

}

void
vtkKWOrientMRIView3D::ZoomBy ( float iFactor ) {
  this->GetRenderer()->GetActiveCamera()->Zoom( iFactor );
  this->GetRenderer()->GetRenderWindow()->Render();
}

void
vtkKWOrientMRIView3D::RotateUserTransform ( int iAxis, double iDegrees ) {

  assert( mVolumeActor.GetPointer() );
  
  for( int nStep = 0; nStep < kcRotateStepsToAnimate; nStep++ ) {
    switch( iAxis ) {
    case 0:
      mVolumeActor->
	RotateWXYZ( iDegrees/(double)kcRotateStepsToAnimate, 1,0,0 );
      break;
    case 1:
      mVolumeActor->
	RotateWXYZ( iDegrees/(double)kcRotateStepsToAnimate, 0,1,0 );
      break;
    case 2:
      mVolumeActor->
	RotateWXYZ( iDegrees/(double)kcRotateStepsToAnimate, 0,0,1 );
      break;
    }

    this->Render();
  }

}

vtkMatrix4x4&
vtkKWOrientMRIView3D::GetUserTransform () {

  assert( mUserTransform.GetPointer() );

  return *mUserTransform;
}

void 
vtkKWOrientMRIView3D::ActorPositionChanged () {

  assert( mVolumeActor.GetPointer() );
  assert( mVolumeActor->GetMatrix() );

  // If we rotate it enough, we can introduce a translation, which is
  // not what we want, so get rid of that now.
  mVolumeActor->SetPosition( 0, 0, 0 );
    
  // Copy the user transform from the outline actor.
  mUserTransform->DeepCopy( mVolumeActor->GetMatrix() );

  // Copy the transform and invert it; this will be how we'll apply it
  // to the existing transform. Then get the orientation, and display
  // that in the corner annotation.
  vtkSmartPointer<vtkTransform> inverseTransform =
    vtkSmartPointer<vtkTransform>::New();
  inverseTransform->SetMatrix( mUserTransform );
  double orientation[3];
  inverseTransform->GetOrientation( orientation );
  stringstream ssLabel;
  ssLabel << setiosflags(ios::fixed) << setprecision(2)
	  << "Degrees around X: " << orientation[0]
	  << " Y: " << orientation[1]
	  << " Z: " << orientation[2];
  this->GetCornerAnnotation()->SetText( 0, ssLabel.str().c_str() );
  
  // Notify our observers.
  this->InvokeEvent( OrientMRIEvents::UserTransformChanged );
}

void
vtkKWOrientMRIView3D::VolumeToRASTransformChanged () {

  // I used to set the cube source's bounds to the new bounds of the
  // FS volume, but I'm not sure if this is useful. Forget it for now.
#if 0
  if( !mCubeSource.GetPointer() )
    return;

  assert( mCubeSource.GetPointer() );
  assert( mFSVolume.GetPointer() );

  // Get the bounds from the FS volume again.
  float bounds[6];
  mFSVolume->GetRASBounds( bounds );
    float offset[3];
  offset[0] = mFSVolume->GetRASCenterX();
  offset[1] = mFSVolume->GetRASCenterY();
  offset[2] = mFSVolume->GetRASCenterZ();

  cerr << "bounds " << bounds[0] << " " << bounds[1] << " "
       << bounds[2] << " " << bounds[3] << " "
       << bounds[4] << " " << bounds[5] << endl;
  cerr << "offset " << offset[0] << " " << offset[1] << " " << offset[2] << endl;
  cerr << "spacing " << mFSVolume->GetPixelSizeX() << " "
       << mFSVolume->GetPixelSizeY() << " "
       << mFSVolume->GetPixelSizeZ() << endl;

  // Set the cube bounds.
  mCubeSource->SetBounds( bounds[0]-offset[0], bounds[1]-offset[0],
 			  bounds[2]-offset[1], bounds[3]-offset[1],
 			  bounds[4]-offset[2], bounds[5]-offset[2] );
#endif

  // Reset the volume matrix. This will trigger a UserTransformChanged
  // event and notify our 2D views.
  assert( mVolumeActor.GetPointer() );
  mVolumeActor->SetOrientation( 0, 0, 0 );
}

void
vtkKWOrientMRIView3D::ActorModifiedCallback ( vtkObject* iCaller,
					      unsigned long iEventId,
					      void* iClientData,
					      void* iCallData ) {

  // We're listening for the actor's Modified event. When we get it,
  // it means the user transform has been changed, so we want to
  // notify the view.
  if( vtkCommand::ModifiedEvent == iEventId ) {

    // Get our view pointer from the client data.
    assert( iClientData );
    try {
      vtkKWOrientMRIView3D* view = 
	static_cast<vtkKWOrientMRIView3D*>( iClientData );

      if( view ) view->ActorPositionChanged();
    }
    catch(...) {
      cerr << "Invalid client data in ActorModified callback" << endl;
    }
  }

}

void
vtkKWOrientMRIView3D::VolumeToRASTransformChangedCallback ( vtkObject* iCaller,
						       unsigned long iEventId,
							    void* iClientData,
							    void* iCallData ) {

  // We're listening for the window's VolumeToRASTransformChanged
  // event. When we get it, it means the volume transform has been
  // changed, so we want to reset our user transform.
  if( OrientMRIEvents::VolumeToRASTransformChanged == iEventId ) {

    // Get our view pointer from the client data.
    assert( iClientData );
    try {
      vtkKWOrientMRIView3D* view = 
	static_cast<vtkKWOrientMRIView3D*>( iClientData );
      
      if( view ) view->VolumeToRASTransformChanged();
    } 
    catch(...) {
      cerr << "Invalid client data in VolumeToRASTransformChanged callback" 
	   << endl;
    }
  }

}
