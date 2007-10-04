/**
 * @file  vtkArrowPipeline.h
 * @brief A simple wrapper for a array source, transform, mapper, and actor.
 *
 * This lets you easily instantiate and place a 3D arrow. You can set
 * the start and end points and get a pointer to the actor.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/04 20:49:39 $
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

#ifndef vtkArrowPipeline_h
#define vtkArrowPipeline_h

#include "vtkObject.h"

#include "vtkSmartPointer.h"

class vtkActor;
class vtkLandmarkTransform;
class vtkTransformPolyDataFilter;

class vtkArrowPipeline : public vtkObject {

 public:

  static vtkArrowPipeline* New();
  vtkTypeRevisionMacro( vtkArrowPipeline, vtkObject );

  // Description:
  // Set the start and end points of the arrow.
  void SetStartPoint ( double const* iPoint );
  void SetEndPoint ( double const* iPoint );
  void SetStartAndEndPoint ( double const* iStartPoint, 
			     double const* iEndPoint );
  // Description:
  // Get a pointer to the actor.
  vtkActor* GetActor () const;

 protected:

  vtkArrowPipeline ();
  virtual ~vtkArrowPipeline ();

  vtkSmartPointer<vtkLandmarkTransform> mTransform;
  vtkSmartPointer<vtkActor> mActor;

  double mStartPoint[3];
  double mEndPoint[3];

 private:
  vtkArrowPipeline ( const vtkArrowPipeline& );
  void operator= ( const vtkArrowPipeline& );
};


#endif
