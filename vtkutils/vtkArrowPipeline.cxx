/**
 * @brief A simple wrapper for a array source, transform, mapper, and actor.
 *
 * This lets you easily instantiate and place a 3D arrow. You can set
 * the start and end points and get a pointer to the actor.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <limits>
#include <assert.h>

#include "vtkArrowPipeline.h"

#include "vtkActor.h"
#include "vtkArrowSource.h"
#include "vtkGeneralTransform.h"
#include "vtkLandmarkTransform.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkTransformPolyDataFilter.h"


vtkStandardNewMacro( vtkArrowPipeline );

vtkArrowPipeline::vtkArrowPipeline () {

  // Arrow source object.
  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetShaftResolution( 20 );
  arrowSource->SetTipResolution( 20 );

  // Initial points.
  mStartPoint[0] = mStartPoint[1] = mStartPoint[2] = 0;
  mEndPoint[0] = 1;
  mEndPoint[1] = mEndPoint[2] = 0;

  // This will be how we transform the source to go from point a to
  // b. Create a simple landmark transform that goes from 0,0,0 to
  // 1,0,0 at the source (the dimensions of the raw arrows source) to
  // our start and end point at the destination.
  mTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
  vtkSmartPointer<vtkPoints> sourcePoints =
    vtkSmartPointer<vtkPoints>::New();
  sourcePoints->SetNumberOfPoints( 2 );
  sourcePoints->SetPoint( 0, 0, 0, 0 );
  sourcePoints->SetPoint( 1, 1, 0, 0 );
  mTransform->SetSourceLandmarks( sourcePoints );
  vtkSmartPointer<vtkPoints> destPoints =
    vtkSmartPointer<vtkPoints>::New();
  destPoints->SetNumberOfPoints( 2 );
  destPoints->SetPoint( 0, mStartPoint );
  destPoints->SetPoint( 1, mEndPoint );
  mTransform->SetTargetLandmarks( destPoints );

  // Landmark transform goes into a general transform.
  vtkSmartPointer<vtkGeneralTransform> arrowTransform =
    vtkSmartPointer<vtkGeneralTransform>::New();
  arrowTransform->SetInput( mTransform );

  // General transform goes into a TPDF.
  vtkSmartPointer<vtkTransformPolyDataFilter> arrowTPDF =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  arrowTPDF->SetTransform( arrowTransform );
  arrowTPDF->SetInputConnection( arrowSource->GetOutputPort() );

  // Now generate normals.
  vtkSmartPointer<vtkPolyDataNormals> arrowNormals =
    vtkSmartPointer<vtkPolyDataNormals>::New();
  arrowNormals->SetInputConnection( arrowTPDF->GetOutputPort() );
  
  // Mapper.
  vtkSmartPointer<vtkPolyDataMapper> arrowMapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  arrowMapper->SetInputConnection( arrowNormals->GetOutputPort() );

  // Actor.
  mActor =  vtkSmartPointer<vtkActor>::New();
  mActor->SetMapper( arrowMapper );

}

vtkArrowPipeline::~vtkArrowPipeline () {

}


void 
vtkArrowPipeline::SetStartPoint ( double const* iPoint ) {

  assert( iPoint );
  assert( mTransform.GetPointer() );
  assert( mActor );

  // Set the first dest landmark point.
  mStartPoint[0] = iPoint[0];
  mStartPoint[1] = iPoint[1];
  mStartPoint[2] = iPoint[2];

  // Copy the new points in. We do it this way to force an update.
  vtkSmartPointer<vtkPoints> destPoints =
    vtkSmartPointer<vtkPoints>::New();
  destPoints->SetNumberOfPoints( 2 );
  destPoints->SetPoint( 0, mStartPoint );
  destPoints->SetPoint( 1, mEndPoint );
  mTransform->SetTargetLandmarks( destPoints );

  // If the points are the same, set visibility off. This is to
  // prevent a case on 64 bit machines where the landmark transform
  // will be way crazy and make the arrow infinitely large.
  if( fabs(mStartPoint[0] - mEndPoint[0]) < 
      std::numeric_limits<double>::epsilon() &&
      fabs(mStartPoint[1] - mEndPoint[2]) < 
      std::numeric_limits<double>::epsilon() &&
      fabs(mStartPoint[2] - mEndPoint[1]) < 
      std::numeric_limits<double>::epsilon() ) {
    mActor->VisibilityOff();
  } else {
    mActor->VisibilityOn();
  }
}

void 
vtkArrowPipeline::SetEndPoint ( double const* iPoint ) {

  assert( iPoint );
  assert( mTransform.GetPointer() );
  assert( mActor );

  // Set the first dest landmark point.
  mEndPoint[0] = iPoint[0];
  mEndPoint[1] = iPoint[1];
  mEndPoint[2] = iPoint[2];

  // Copy the new points in. We do it this way to force an update.
  vtkSmartPointer<vtkPoints> destPoints =
    vtkSmartPointer<vtkPoints>::New();
  destPoints->SetNumberOfPoints( 2 );
  destPoints->SetPoint( 0, mStartPoint );
  destPoints->SetPoint( 1, mEndPoint );
  mTransform->SetTargetLandmarks( destPoints );

  // If the points are the same, set visibility off. This is to
  // prevent a case on 64 bit machines where the landmark transform
  // will be way crazy and make the arrow infinitely large.
  if( fabs(mStartPoint[0] - mEndPoint[0]) < 
      std::numeric_limits<double>::epsilon() &&
      fabs(mStartPoint[1] - mEndPoint[2]) < 
      std::numeric_limits<double>::epsilon() &&
      fabs(mStartPoint[2] - mEndPoint[1]) < 
      std::numeric_limits<double>::epsilon() ) {
    mActor->VisibilityOff();
  } else {
    mActor->VisibilityOn();
  }
}

void 
vtkArrowPipeline::SetStartAndEndPoint ( double const* iStartPoint, 
					double const* iEndPoint ) {
  
  assert( iStartPoint );
  assert( iEndPoint );

  this->SetStartPoint( iStartPoint );
  this->SetEndPoint( iEndPoint );
}

vtkActor* 
vtkArrowPipeline::GetActor () const {

  return mActor.GetPointer();
}

