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
 *    $Date: 2007/10/04 18:23:58 $
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

#include "vtkArrowPipeline.h"

#include "vtkActor.h"
#include "vtkArrowSource.h"
#include "vtkGeneralTransform.h"
#include "vtkLandmarkTransform.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkTransformPolyDataFilter.h"


vtkStandardNewMacro( vtkArrowPipeline );
vtkCxxRevisionMacro( vtkArrowPipeline, "$Revision: 1.1 $" );

vtkArrowPipeline::vtkArrowPipeline () {

  // Arrow source object.
  vtkSmartPointer<vtkArrowSource> arrowSource =
    vtkSmartPointer<vtkArrowSource>::New();

  // This will be how we transform the source to go from point a to
  // b. Create a simple landmark transform with just one point in the
  // source and dest landmarks. Start tham out at 0,0,0.
  mTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
  vtkSmartPointer<vtkPoints> sourcePoints =
    vtkSmartPointer<vtkPoints>::New();
  sourcePoints->SetNumberOfPoints( 1 );
  sourcePoints->SetPoint( 0, 0, 0, 0 );
  mTransform->SetSourceLandmarks( sourcePoints );
  vtkSmartPointer<vtkPoints> destPoints =
    vtkSmartPointer<vtkPoints>::New();
  destPoints->SetNumberOfPoints( 1 );
  destPoints->SetPoint( 0, 0, 0, 0 );
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

  // Mapper.
  vtkSmartPointer<vtkPolyDataMapper> arrowMapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  arrowMapper->SetInputConnection( arrowTPDF->GetOutputPort() );

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

  // Just set the source landmark point.
  vtkSmartPointer<vtkPoints> sourcePoints =
    vtkSmartPointer<vtkPoints>::New();
  sourcePoints->SetNumberOfPoints( 1 );
  sourcePoints->SetPoint( 0, iPoint );
  mTransform->SetSourceLandmarks( sourcePoints );

}

void 
vtkArrowPipeline::SetEndPoint ( double const* iPoint ) {

  assert( iPoint );
  assert( mTransform.GetPointer() );

  // Just set the dest landmark point.
  vtkSmartPointer<vtkPoints> destPoints =
    vtkSmartPointer<vtkPoints>::New();
  destPoints->SetNumberOfPoints( 1 );
  destPoints->SetPoint( 0, iPoint );
  mTransform->SetTargetLandmarks( destPoints );
}

void 
vtkArrowPipeline::SetStartAndEndPoint ( double const* iStartPoint, 
					double const* iEndPoint ) {
  
  assert( iStartPoint );
  assert( iEndPoint );

  this->SetStartPoint( iStartPoint );
  this->SetStartPoint( iEndPoint );
}

vtkActor* 
vtkArrowPipeline::GetActor () const {

  return mActor.GetPointer();
}

