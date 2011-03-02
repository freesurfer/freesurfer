/**
 * @file  vtkOrientMRIInteractorStyleView3D.h
 * @brief Interactor style for the 3D view
 *
 * This interactor style allows camera trackball rotation with left
 * button, volume trackball rotation with middle button down, and
 * zooming with right button. Most code is yanked from VTK's
 * vtkInteractorStyle* classes.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:35 $
 *    $Revision: 1.2 $
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

#ifndef vtkOrientMRIInteractorStyleView3D_h
#define vtkOrientMRIInteractorStyleView3D_h

#include "vtkInteractorStyle.h"
#include "vtkSmartPointer.h"

class vtkOrientMRIInteractorStyleView3D : public vtkInteractorStyle {

public:
  static vtkOrientMRIInteractorStyleView3D *New();
  vtkTypeRevisionMacro(vtkOrientMRIInteractorStyleView3D,vtkInteractorStyle);

  // Description:
  // Sets the prop to rotate. This interactor style does not use a picker
  // to find an prop to rotate, it just uses this one.
  void SetInteractionProp ( vtkProp3D* iProp );

  // Description:
  // Event bindings controlling the effects of pressing mouse buttons
  // or moving the mouse.
  virtual void OnMouseMove ();
  virtual void OnLeftButtonDown ();
  virtual void OnLeftButtonUp ();
  virtual void OnMiddleButtonDown ();
  virtual void OnMiddleButtonUp ();
  virtual void OnRightButtonDown ();
  virtual void OnRightButtonUp ();

  // Description:
  // Methods for manipulating the prop. These are copied from
  // vtkInteractorStyleTrackballActor.
  virtual void RotateProp ();
  virtual void SpinProp ();

  // Description:
  // Methods for manipulating the camera. These are copied from
  // vtkInteractorStyleTrackballCamera.
  virtual void RotateCamera ();
  virtual void SpinCamera ();
  virtual void DollyCamera ();
  virtual void DollyCamera ( double factor );


protected:
  vtkOrientMRIInteractorStyleView3D();
  ~vtkOrientMRIInteractorStyleView3D();

  // Description:
  // Method for manipulating the prop. This is copied from
  // vtkInteractorStyleTrackballActor.
  void Prop3DTransform(vtkProp3D *prop3D,
                       double *boxCenter,
                       int NumRotation,
                       double **rotate,
                       double *scale);

  // A constant for camera movement.
  double MotionFactor;

  // The prop we'll manipulate.
  vtkSmartPointer<vtkProp3D> InteractionProp;

  // Our custom states. We can't use the normal VTKIS_* states because
  // they don't differentiate between the actor and the camera.
  enum { IS_RotateProp = 1000, IS_SpinProp, 
	 IS_RotateCamera, IS_SpinCamera, IS_DollyCamera };

  // Start and end our states.
  void StartRotateProp ();
  void EndRotateProp ();
  void StartSpinProp ();
  void EndSpinProp ();
  void StartRotateCamera ();
  void EndRotateCamera ();
  void StartSpinCamera ();
  void EndSpinCamera ();
  void StartDollyCamera ();
  void EndDollyCamera ();

private:
  vtkOrientMRIInteractorStyleView3D(const vtkOrientMRIInteractorStyleView3D&);
  void operator=(const vtkOrientMRIInteractorStyleView3D&);
};

#endif
