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
  vtkTypeMacro( vtkArrowPipeline, vtkObject );

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
