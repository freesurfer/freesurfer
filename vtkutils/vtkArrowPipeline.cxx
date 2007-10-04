/**
 * @file  vtkArrowPipeline.cxx
 * @brief A simple wrapper for a array source, transform, mapper, and actor.
 *
 * This lets you easily instantiate and place a 3D arrow. You can set
 * the start and end points and get a pointer to the actor.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/04 21:52:47 $
 *    $Revision: 1.4 $
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
vtkCxxRevisionMacro( vtkArrowPipeline, "$Revision: 1.4 $" );

vtkArrowPipeline::vtkArrowPipeline () {

  // Arrow source object.
  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();
  arrowSource->SetShaftResolution( 20 );
  arrowSource->SetTipResolution( 20 );

  // Initial points.
  mStartPoint[0] = mStartPoint[1] = mStartPoint[2] = 0;
  mEndPoint[0] = mEndPoint[1] = mEndPoint[2] = 1;

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
}

void 
vtkArrowPipeline::SetEndPoint ( double const* iPoint ) {

  assert( iPoint );
  assert( mTransform.GetPointer() );

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

