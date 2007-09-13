/**
 * @file  vtkKWOrientMRIView2D.cxx
 * @brief Viewing pipeline for volume slice
 *
 * Pipeline objects for displaying a volume slice..
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/09/13 20:58:21 $
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


#include <sstream>
#include <stdexcept>
#include "vtkKWOrientMRIView2D.h"

// These are defined in the mri headers somewhere.
#ifdef X
#undef X
#endif
#ifdef Y
#undef Y
#endif

#include "vtkActor.h"
#include "vtkCamera.h"
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
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTexture.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkWindowLevelLookupTable.h"

using namespace std;

vtkStandardNewMacro( vtkKWOrientMRIView2D );
vtkCxxRevisionMacro( vtkKWOrientMRIView2D, "$Revision: 1.1 $" );


vtkKWOrientMRIView2D::vtkKWOrientMRIView2D () :
  mPlaneOrientation( UninitedPlaneOrientation ),
  mThroughPlane( 0.0 ) {

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
  //  mPlaneMapper->ImmediateModeRenderingOn();
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

  // Initial camera position.
  this->GetRenderer()->GetActiveCamera()->SetPosition( 0, 1, 0 );
  this->GetRenderer()->GetActiveCamera()->SetFocalPoint( 0, 0, 0 );
  this->GetRenderer()->GetActiveCamera()->SetViewUp( 0, 0, 1 );

  // Reset the camera to show all the actors.
  this->GetRenderer()->ResetCamera();
  this->GetRenderer()->GetRenderWindow()->Render();

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
vtkKWOrientMRIView2D::SetUpView () {
  
  if( UninitedPlaneOrientation == mPlaneOrientation ) 
    return;

  assert( mPlaneTransform.GetPointer() );
  assert( mReslice.GetPointer() );
  assert( mOutlineActor.GetPointer() );

  float offset[3];
  offset[0] = mFSVolume->GetRASCenterX();
  offset[1] = mFSVolume->GetRASCenterY();
  offset[2] = mFSVolume->GetRASCenterZ();

  // Set the orientation.
  switch ( mPlaneOrientation ) {
  case 0:
    mPlaneTransform->Identity();
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    mReslice->SetResliceAxesDirectionCosines( 0, -1, 0,
					      0, 0, 1,
					      1, 0, 0 );
    mReslice->SetResliceAxesOrigin( (int)(mThroughPlane - offset[0]), 0, 0 );
    break;
  case 1:
    mPlaneTransform->Identity();
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					      0, 0, 1,
					      0, 1, 0 );
    mReslice->SetResliceAxesOrigin( 0, (int)(mThroughPlane - offset[1]), 0 );
    break;
  case 2:
    mPlaneTransform->Identity();
    mPlaneTransform->RotateX( 90 );
    mPlaneTransform->RotateY( 180 );
    mReslice->SetResliceAxesDirectionCosines( 1, 0, 0,
					      0, 1, 0,
					      0, 0, -1 );
    mReslice->SetResliceAxesOrigin( 0, 0, (int)(mThroughPlane - offset[2]) );
    break;
  }

  // Set a color for the plane outline. This should match the axes
  // color in the 3D view.
  double color[3] = {0,0,0};
  switch ( mPlaneOrientation ) {
  case 0:
    color[0] = 1; color[1] = 0; color[2] = 0;
    break;
  case 1:
    color[0] = 1; color[1] = 1; color[2] = 0;
    break;
  case 2:
    color[0] = 0; color[1] = 1; color[2] = 0;
    break;
  }
  mOutlineActor->GetProperty()->SetColor( color );

  // Set the annotation to print our through plane label and number.
  stringstream ssLabel;
  if( 0 == mPlaneOrientation ) ssLabel << "X: ";
  else if( 1 == mPlaneOrientation ) ssLabel << "Y: ";
  else if( 2 == mPlaneOrientation ) ssLabel << "Z: ";
  ssLabel << mThroughPlane;
  this->GetCornerAnnotation()->SetText( 0, ssLabel.str().c_str() );

  this->GetRenderer()->GetRenderWindow()->Render();
}

