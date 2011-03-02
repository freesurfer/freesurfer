/**
 * @file  vtkKWOrientMRIView2D.cxx
 * @brief Viewing pipeline for volume slice
 *
 * Pipeline objects for displaying a volume slice..
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.8 $
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


#include <sstream>
#include <stdexcept>

#include <assert.h>

#include "vtkKWOrientMRIView2D.h"

// These are defined in the mri headers somewhere.
#ifdef X
#undef X
#endif
#ifdef Y
#undef Y
#endif

#include "OrientMRIEvents.h"
#include "vtkActor.h"
#include "vtkArrowPipeline.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkColorTransferFunction.h"
#include "vtkCornerAnnotation.h"
#include "vtkFSVolumeSource.h"
#include "vtkImageData.h"
#include "vtkImageFlip.h"
#include "vtkImageMapToColors.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageReslice.h"
#include "vtkOrientMRIInteractorStyleView2D.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkOutlineFilter.h"
#include "vtkPlaneSource.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTextProperty.h"
#include "vtkTexture.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkWindowLevelLookupTable.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIView2D );
vtkCxxRevisionMacro( vtkKWOrientMRIView2D, "$Revision: 1.8 $" );

double vtkKWOrientMRIView2D::sReorthoPoints[3][2][3] = 
  { {{0,0,0}, {0,0,0}}, {{0,0,0}, {0,0,0}}, {{0,0,0}, {0,0,0}} };
vector<vtkKWOrientMRIView2D*> vtkKWOrientMRIView2D::slViews;

vtkKWOrientMRIView2D::vtkKWOrientMRIView2D () :
  mPlaneOrientation( UninitedPlaneOrientation ),
  mCurrentCameraOrientation( UninitedPlaneOrientation ),
  mThroughPlane( 0.0 ) {

  // Insert us in the list of instances.
  slViews.push_back( this );
}

vtkKWOrientMRIView2D::~vtkKWOrientMRIView2D () {

}

void
vtkKWOrientMRIView2D::Create () {

  this->Superclass::Create();

  // Set our interactor style.
  vtkSmartPointer<vtkOrientMRIInteractorStyleView2D> style = 
    vtkSmartPointer<vtkOrientMRIInteractorStyleView2D>::New();
  this->GetRenderWindow()->GetInteractor()->SetInteractorStyle( style );

  // Add the pipelines for our reortho lines.
  for( int nOrientation = 0; nOrientation < 3; nOrientation++ ) {

    // Make a pipeline and set our start and end poitns.
    mReorthoArrow[nOrientation] = vtkSmartPointer<vtkArrowPipeline>::New();
    mReorthoArrow[nOrientation]->
      SetStartAndEndPoint( sReorthoPoints[nOrientation][0],
			   sReorthoPoints[nOrientation][1] );

    // Get the actor and set its color and pickability.
    vtkSmartPointer<vtkActor> reorthoActor = 
      mReorthoArrow[nOrientation]->GetActor();
    reorthoActor->GetProperty()->
      SetColor( const_cast<double*>
		(this->GetColorForOrientation(nOrientation)) );
    reorthoActor->PickableOff();

    // Add it to the renderer.
    this->GetRenderer()->AddActor( reorthoActor );

  }

  // Tell our style to edit our points.
  style->SetPointsAndLines( sReorthoPoints[0][0], sReorthoPoints[0][1],
			    sReorthoPoints[1][0], sReorthoPoints[1][1],
			    sReorthoPoints[2][0], sReorthoPoints[2][1] );

  // Make a callback and observe the style for OrthoLineChanged events.
  vtkSmartPointer<vtkCallbackCommand> callback =
    vtkSmartPointer<vtkCallbackCommand>::New();
  callback->SetCallback( OrthoLineChanged );
  callback->SetClientData( this );
  style->AddObserver( OrientMRIEvents::OrthoLineChanged, callback );

  // Make our corner annotation bold.
  this->GetCornerAnnotation()->GetTextProperty()->BoldOn();
}

void
vtkKWOrientMRIView2D::SetFSVolumeAndImage ( vtkFSVolumeSource& iFSVolume,
					    vtkImageData& iImage ) {

  // Save a pointer.
  mFSVolume = &iFSVolume;
  mImage = &iImage;

  //
  // The reslice object just takes a slice out of the volume.
  //
  mReslice = vtkSmartPointer<vtkImageReslice>::New();
  mReslice->SetInput( mImage );
  //  mReslice->BorderOff();

  // This sets us to extract slices.
  mReslice->SetOutputDimensionality( 2 );

  // This will change depending what orienation we're in.
  mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					    0, 1, 0,
					    0, 0, 1 );

  // This will change to select a different slice.
  mReslice->SetResliceAxesOrigin( 0, 0, 0 );

  //
  // Flip over the x axis (left/right). This get us into neurological
  // view.
  //
  vtkSmartPointer<vtkImageFlip> imageFlip = 
    vtkSmartPointer<vtkImageFlip>::New();
  imageFlip->SetInputConnection( mReslice->GetOutputPort() );
  imageFlip->SetFilteredAxis( 0 ); // x axis


  //
  // A quick color table.
  //
  mColorTable = vtkSmartPointer<vtkWindowLevelLookupTable>::New();
  double range[2];
  mImage->GetScalarRange( range );
  mColorTable->SetWindow( range[1]-range[0] );
  mColorTable->SetLevel( (range[1]-range[0]) / 2.0 );

  //
  // Image to colors using color table.
  //
  mColorMap = vtkSmartPointer<vtkImageMapToColors>::New();
  mColorMap->SetInputConnection( imageFlip->GetOutputPort() );
  mColorMap->SetLookupTable( mColorTable );

  //
  // Colors to texture.
  //
  mTexture = vtkSmartPointer<vtkTexture>::New();
  mTexture->SetInputConnection( mColorMap->GetOutputPort() );
  mTexture->RepeatOff();
  mTexture->InterpolateOff();

  //
  // Plane mesh object.
  //
  vtkSmartPointer<vtkPlaneSource> plane = 
    vtkSmartPointer<vtkPlaneSource>::New();

  //
  // Plane mapper transform.
  //
  mPlaneTransform = vtkSmartPointer<vtkTransform>::New();

  //
  // Poly data from plane and plane transform.
  //
  vtkSmartPointer<vtkTransformPolyDataFilter> planePDF = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  planePDF->SetInput( plane->GetOutput() );
  planePDF->SetTransform( mPlaneTransform );

  //
  // Mapper for plane.
  //
  mPlaneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mPlaneMapper->SetInputConnection( planePDF->GetOutputPort() );


  //
  // Prop in scene with plane mesh and texture.
  //
  mPlaneActor = vtkSmartPointer<vtkActor>::New();
  mPlaneActor->SetMapper( mPlaneMapper );
  mPlaneActor->SetTexture( mTexture );

  // Add this actor to the view.
  this->GetRenderer()->AddActor( mPlaneActor );


  //
  // Plane outline.
  //
  vtkSmartPointer<vtkOutlineFilter> outline =
    vtkSmartPointer<vtkOutlineFilter>::New();
  outline->SetInput( planePDF->GetOutput() );

  vtkSmartPointer<vtkPolyDataMapper> outlineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInputConnection( outline->GetOutputPort() );

  mOutlineActor = vtkSmartPointer<vtkActor>::New();
  mOutlineActor->SetMapper( outlineMapper );

  // Add this actor to the view.
  this->GetRenderer()->AddActor( mOutlineActor );


  // Initialze our view.
  this->SetUpView();

  // That lookup table is controlled by our interactor style.
  vtkOrientMRIInteractorStyleView2D* style =
    dynamic_cast<vtkOrientMRIInteractorStyleView2D*>( this->GetRenderWindow()->GetInteractor()->GetInteractorStyle() );
  if( style ) {
    style->SetWindowLevelTable( mColorTable );
  }


#if 0
  // Initial camera position.
  this->GetRenderer()->GetActiveCamera()->SetPosition( 0, 1, 0 );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( 0, 0, 0 );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( 0, 0, 1 );

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetRenderWindow()->Render();
#endif


}

void
vtkKWOrientMRIView2D::SetImageColors ( vtkLookupTable& iColors ) {

  // Save a pointer.
  mImageColors = &iColors;

  // Set it in our plane if exist.
  if( mColorMap.GetPointer() )
    mColorMap->SetLookupTable( mImageColors );
}

void
vtkKWOrientMRIView2D::SetOrientation ( int iOrientation ) {
  
  // Set the orientation. This will determing the plane transform.
  mPlaneOrientation = iOrientation;

  // Set our annotation color.
  this->SetCornerAnnotationColor( const_cast<double*>
			  (this->GetColorForOrientation(mPlaneOrientation)) );

  if( mPlaneTransform.GetPointer() )
    this->SetUpView();
}

void
vtkKWOrientMRIView2D::SetThroughPlane ( float iZ ) {

  // Set the through plane. This will determine the reslice axes.
  mThroughPlane = iZ;

  if( mReslice.GetPointer() )
    this->SetUpView();
}

void
vtkKWOrientMRIView2D::UpdateOrthoLine ( int inOrientation ) {
  
  assert( inOrientation == 0 ||
	  inOrientation == 1 ||
	  inOrientation == 2 );

  // Tell the arrow to update.
  mReorthoArrow[inOrientation]->
    SetStartAndEndPoint( sReorthoPoints[inOrientation][0],
			 sReorthoPoints[inOrientation][1] );

  this->Render();
}

double const*
vtkKWOrientMRIView2D::GetColorForOrientation ( int inOrientation ) const {

  assert( inOrientation == 0 ||
	  inOrientation == 1 ||
	  inOrientation == 2 );

  static double sColors[3][3] = { { 1, 0, 0 }, { 1, 1, 0 }, {0, 1, 0} };

  return sColors[inOrientation];
}

void
vtkKWOrientMRIView2D::SetUpView () {
  
  if( UninitedPlaneOrientation == mPlaneOrientation ) 
    return;

  assert( mPlaneTransform.GetPointer() );
  assert( mReslice.GetPointer() );
  assert( mOutlineActor.GetPointer() );

  // Get the offset so our through plane coords are correct.
  float offset[3];
  offset[0] = mFSVolume->GetRASCenterX();
  offset[1] = mFSVolume->GetRASCenterY();
  offset[2] = mFSVolume->GetRASCenterZ();

  // Calcualte the size of the planes.
  float RASBounds[6];
  mFSVolume->GetRASBounds( RASBounds );
  float size[3];
  size[0] = RASBounds[1] - RASBounds[0];
  size[1] = RASBounds[3] - RASBounds[2];
  size[2] = RASBounds[5] - RASBounds[4];

  // If our camera is looking at a different orientation, set it up
  // now.
  if( mCurrentCameraOrientation != mPlaneOrientation ) {
    
    // Get a pointer to the camera.
    vtkSmartPointer<vtkCamera> camera = 
      this->GetRenderer()->GetActiveCamera();
    
    // Switch on the orientation and set our camera location
    // appropriately.
    switch ( mPlaneOrientation ) {
    case 0:
      camera->SetFocalPoint( 0, 0, 0 );
      camera->SetViewUp( 0, 0, 1 );
      camera->SetPosition( 10, 0, 0 );
      break;
    case 1:
      camera->SetFocalPoint( 0, 0, 0 );
      camera->SetViewUp( 0, 0, 1 );
      camera->SetPosition( 0, 10, 0 );
      break;
    case 2:
      camera->SetFocalPoint( 0, 0, 0 );
      camera->SetViewUp( 0, 1, 0 );
      camera->SetPosition( 0, 0, 10 );
      break;
    default:
      break;
    }
    
    // Reset the camera to show all the actors.
    this->GetRenderer()->ResetCamera();

    // Camera is now set for this orientation.
    mCurrentCameraOrientation = mPlaneOrientation;
  }

  // Switch on the orientation and set our plane rotation and reslice
  // axes appropriately.
  switch ( mPlaneOrientation ) {
  case 0:
    mPlaneTransform->Identity();
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 90 );
    //mPlaneTransform->Scale( size[1], size[2], 1 );
    mReslice->SetResliceAxesDirectionCosines( 0, -1, 0,  0, 0, 1,  1, 0, 0 );
    mReslice->SetResliceAxesOrigin( (int)(mThroughPlane - offset[0]), 0, 0 );
    break;
  case 1:
    mPlaneTransform->Identity();
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    //mPlaneTransform->Scale( size[0], size[2], 1 );
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,  0, 0, 1,  0, 1, 0 );
    mReslice->SetResliceAxesOrigin( 0, (int)(mThroughPlane - offset[1]), 0 );
    break;
  case 2:
    mPlaneTransform->Identity();
    //mPlaneTransform->Scale( size[0], size[1], 1 );
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,  0, 1, 0,  0, 0, -1 );
    mReslice->SetResliceAxesOrigin( 0, 0, (int)(mThroughPlane - offset[2]) );
    break;
  default:
    break;
  }

  // Set a color for the plane outline. This should match the axes
  // color in the 3D view.
  mOutlineActor->GetProperty()->
    SetColor( const_cast<double*>
	      (this->GetColorForOrientation(mPlaneOrientation)) );

  // Set the annotation to print our through plane label and number.
  stringstream ssLabel;
  if( 0 == mPlaneOrientation ) ssLabel << "X: ";
  else if( 1 == mPlaneOrientation ) ssLabel << "Y: ";
  else if( 2 == mPlaneOrientation ) ssLabel << "Z: ";
  ssLabel << mThroughPlane;
  this->GetCornerAnnotation()->SetText( 0, ssLabel.str().c_str() );

  this->GetRenderer()->GetRenderWindow()->Render();
}

void
vtkKWOrientMRIView2D::ResetOrthoPoints () {

  // Put all points back at 0,0,0.
  for( int nOrientation = 0; nOrientation < 3; ++nOrientation ) {

    for( int nPoint = 0; nPoint < 2; ++nPoint )
      for( int nCoord = 0; nCoord < 3; ++nCoord )
	sReorthoPoints[nOrientation][nPoint][nCoord] = 0;
    
    // Go through our static list of views and call their
    // UpdateOrthoLines functions.
    vector<vtkKWOrientMRIView2D*>::iterator tView;
    for( tView = slViews.begin(); tView != slViews.end(); ++tView )
      (*tView)->UpdateOrthoLine( nOrientation );
  }
}

void
vtkKWOrientMRIView2D::OrthoLineChanged ( vtkObject* iCaller, 
					 unsigned long iEventId,
					 void* iClientData,
					 void* iCallData ) {
  
  // We're listening for OrthoLineChanged events. When we get one,
  // call the views' UpdateOrthoLines function.
  if( OrientMRIEvents::OrthoLineChanged == iEventId ) {

    assert( iClientData );
    assert( iCallData );
    
    // Get the orientation.
    int nOrientation = *static_cast<int*>( iCallData );
    
    // Go through our static list of views and call their
    // UpdateOrthoLines functions.
    vector<vtkKWOrientMRIView2D*>::iterator tView;
    for( tView = slViews.begin(); tView != slViews.end(); ++tView )
      (*tView)->UpdateOrthoLine( nOrientation );
  }

}

void
vtkKWOrientMRIView2D::GetReorthoPoints ( double const*& oXStart, 
					 double const*& oXEnd,
					 double const*& oYStart,
					 double const*& oYEnd,
					 double const*& oZStart,
					 double const*& oZEnd ) {

  oXStart = sReorthoPoints[0][0];
  oXEnd = sReorthoPoints[0][1];
  oYStart = sReorthoPoints[1][0];
  oYEnd = sReorthoPoints[1][1];
  oZStart = sReorthoPoints[2][0];
  oZEnd = sReorthoPoints[2][1];
}
