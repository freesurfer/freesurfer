/**
 * @file  vtkOrientMRIInteractorStyleView3D.cxx
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

#include "vtkOrientMRIInteractorStyleView3D.h"

#include "vtkCamera.h"
#include "vtkCallbackCommand.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkProp3D.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkTransform.h"

vtkCxxRevisionMacro(vtkOrientMRIInteractorStyleView3D, "$Revision: 1.4 $");
vtkStandardNewMacro(vtkOrientMRIInteractorStyleView3D);

vtkOrientMRIInteractorStyleView3D::vtkOrientMRIInteractorStyleView3D () :
  MotionFactor( 10.0 ) {
}

vtkOrientMRIInteractorStyleView3D::~vtkOrientMRIInteractorStyleView3D () {
}

void
vtkOrientMRIInteractorStyleView3D::SetInteractionProp ( vtkProp3D* iProp ) {

  // Save a pointer to our prop.
  this->InteractionProp = iProp;
}

void
vtkOrientMRIInteractorStyleView3D::OnMouseMove () {
  
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  // If we're in one of our states, find the renderer and call the
  // function to do the action.
  switch( this->State ) {
  
  case VTKIS_NONE:
    return;
    
  case IS_RotateProp:
    this->FindPokedRenderer( x, y );
    this->RotateProp();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    break;
    
  case IS_SpinProp:
    this->FindPokedRenderer( x, y );
    this->SpinProp();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    break;
    
  case IS_RotateCamera:
    this->FindPokedRenderer( x, y );
    this->RotateCamera();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    break;
    
  case IS_SpinCamera:
    this->FindPokedRenderer( x, y );
    this->SpinCamera();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    break;
    
  case IS_DollyCamera:
    this->FindPokedRenderer( x, y );
    this->DollyCamera();
    this->InvokeEvent( vtkCommand::InteractionEvent, NULL );
    break;
    
  default:
    vtkErrorMacro( << "Unrecognized style state" );
  }
  
}

void
vtkOrientMRIInteractorStyleView3D::OnLeftButtonDown () {
  
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];
  
  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer)
    return;

  // Grab the focus and work on the prop: if ctrl is down, spin it,
  // othewise rotate it.
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  this->GrabFocus( this->EventCallbackCommand );
#endif
  if( this->Interactor->GetControlKey() )
    this->StartSpinProp();
  else
    this->StartRotateProp();
}

void
vtkOrientMRIInteractorStyleView3D::OnLeftButtonUp () {

  // If we're in one of the states that we enter on left button down,
  // stop doing that action.
  switch( this->State ) {
    
  case VTKIS_NONE:
    return;
    
  case IS_SpinProp:
    this->EndSpinProp();
    break;
    
  case IS_RotateProp:
    this->EndRotateProp();
      break;
      
  default:
    vtkErrorMacro( << "Unrecognized style state" );
  }
  
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  if( this->Interactor )
    this->ReleaseFocus();
#endif
}

void
vtkOrientMRIInteractorStyleView3D::OnMiddleButtonDown () {

  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer( x, y );
  if( NULL == this->CurrentRenderer )
    return;
  
  // Grab the focus and work on the camera: if ctrl is down, spin it,
  // othewise rotate it.
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  this->GrabFocus( this->EventCallbackCommand );
#endif
  if( this->Interactor->GetControlKey() )
    this->StartSpinCamera();
  else
    this->StartRotateCamera();
}

void 
vtkOrientMRIInteractorStyleView3D::OnMiddleButtonUp () {

  // If we're in one of the states that we enter on middle button
  // down, stop doing that action.
  switch( this->State ) {

  case VTKIS_NONE:
    return;
    
  case IS_SpinCamera:
    this->EndSpinCamera();
    break;
    
  case IS_RotateCamera:
    this->EndRotateCamera();
    break;
    
  default:
    vtkErrorMacro( << "Unrecognized style state" );
  }
  
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  if( this->Interactor )
    this->ReleaseFocus();
#endif
}

void
vtkOrientMRIInteractorStyleView3D::OnRightButtonDown () {

  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer(x, y);
  if( NULL == this->CurrentRenderer )
    return;
  
  // Grab the focus and start dollying the camera.
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  this->GrabFocus(this->EventCallbackCommand);
#endif
  this->StartDollyCamera();
}

void
vtkOrientMRIInteractorStyleView3D::OnRightButtonUp () {

  // If we were in dolly state, end it.
  switch( this->State ) {
    
  case VTKIS_NONE:
    return;
    
  case IS_DollyCamera:
    this->EndDollyCamera();
    break;
    
  default:
    vtkErrorMacro( << "Unrecognized style state" );
  }
  
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 0))
  if( this->Interactor )
    this->ReleaseFocus();
#endif
}

void
vtkOrientMRIInteractorStyleView3D::RotateProp () {

  if( this->CurrentRenderer == NULL || 
      this->InteractionProp.GetPointer() == NULL )
    return;

  // Copied from vtkInteractorStyleTrackballCamera.
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();

  // First get the origin of the assembly
  double *obj_center = this->InteractionProp->GetCenter();
  
  // GetLength gets the length of the diagonal of the bounding box
  double boundRadius = this->InteractionProp->GetLength() * 0.5;

  // Get the view up and view right vectors
  double view_up[3], view_look[3], view_right[3];

  cam->OrthogonalizeViewUp();
  cam->ComputeViewPlaneNormal();
  cam->GetViewUp(view_up);
  vtkMath::Normalize(view_up);
  cam->GetViewPlaneNormal(view_look);
  vtkMath::Cross(view_up, view_look, view_right);
  vtkMath::Normalize(view_right);
  
  // Get the furtherest point from object position+origin
  double outsidept[3];

  outsidept[0] = obj_center[0] + view_right[0] * boundRadius;
  outsidept[1] = obj_center[1] + view_right[1] * boundRadius;
  outsidept[2] = obj_center[2] + view_right[2] * boundRadius;
  
  // Convert them to display coord
  double disp_obj_center[3];

  this->ComputeWorldToDisplay(obj_center[0], obj_center[1], obj_center[2], 
                              disp_obj_center);

  this->ComputeWorldToDisplay(outsidept[0], outsidept[1], outsidept[2], 
                              outsidept);
  
  double radius = sqrt(vtkMath::Distance2BetweenPoints(disp_obj_center,
                                                       outsidept));
  double nxf = 
    ((double)rwi->GetEventPosition()[0] - (double)disp_obj_center[0]) / radius;
  
  double nyf = 
    ((double)rwi->GetEventPosition()[1] - (double)disp_obj_center[1]) / radius;

  double oxf = 
    ((double)rwi->GetLastEventPosition()[0] - (double)disp_obj_center[0]) / radius;

  double oyf = 
    ((double)rwi->GetLastEventPosition()[1] - (double)disp_obj_center[1]) / radius;

  if (((nxf * nxf + nyf * nyf) <= 1.0) &&
      ((oxf * oxf + oyf * oyf) <= 1.0))
    {
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 3))
    double newXAngle = asin(nxf) * vtkMath::DegreesFromRadians(1.f);
    double newYAngle = asin(nyf) * vtkMath::DegreesFromRadians(1.f);
    double oldXAngle = asin(oxf) * vtkMath::DegreesFromRadians(1.f);
    double oldYAngle = asin(oyf) * vtkMath::DegreesFromRadians(1.f);
#else
    double newXAngle = asin(nxf) * vtkMath::RadiansToDegrees();
    double newYAngle = asin(nyf) * vtkMath::RadiansToDegrees();
    double oldXAngle = asin(oxf) * vtkMath::RadiansToDegrees();
    double oldYAngle = asin(oyf) * vtkMath::RadiansToDegrees();
#endif    
    double scale[3];
    scale[0] = scale[1] = scale[2] = 1.0;

    double **rotate = new double*[2];

    rotate[0] = new double[4];
    rotate[1] = new double[4];
    
    rotate[0][0] = newXAngle - oldXAngle;
    rotate[0][1] = view_up[0];
    rotate[0][2] = view_up[1];
    rotate[0][3] = view_up[2];
    
    rotate[1][0] = oldYAngle - newYAngle;
    rotate[1][1] = view_right[0];
    rotate[1][2] = view_right[1];
    rotate[1][3] = view_right[2];
    
    
    this->Prop3DTransform(this->InteractionProp,
                          obj_center,
                          2, 
                          rotate, 
                          scale);
    
    delete [] rotate[0];
    delete [] rotate[1];
    delete [] rotate;
    
    if (this->AutoAdjustCameraClippingRange)
      {
      this->CurrentRenderer->ResetCameraClippingRange();
      }
    
    rwi->Render();
    }
}
  
void
vtkOrientMRIInteractorStyleView3D::SpinProp() {

  if( this->CurrentRenderer == NULL || 
      this->InteractionProp.GetPointer() == NULL )
    return;

  // Copied from vtkInteractorStyleTrackballCamera.

  vtkRenderWindowInteractor *rwi = this->Interactor;
  vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();
  
  // Get the axis to rotate around = vector from eye to origin

  double *obj_center = this->InteractionProp->GetCenter();

  double motion_vector[3];
  double view_point[3];

  if (cam->GetParallelProjection())
    {
    // If parallel projection, want to get the view plane normal...
    cam->ComputeViewPlaneNormal();
    cam->GetViewPlaneNormal(motion_vector);
    }
  else
    {   
    // Perspective projection, get vector from eye to center of actor
    cam->GetPosition(view_point);
    motion_vector[0] = view_point[0] - obj_center[0];
    motion_vector[1] = view_point[1] - obj_center[1];
    motion_vector[2] = view_point[2] - obj_center[2];
    vtkMath::Normalize(motion_vector);
    }
  
  double disp_obj_center[3];
  
  this->ComputeWorldToDisplay(obj_center[0], obj_center[1], obj_center[2], 
                              disp_obj_center);
  
  double newAngle = 
    atan2((double)rwi->GetEventPosition()[1] - (double)disp_obj_center[1],
          (double)rwi->GetEventPosition()[0] - (double)disp_obj_center[0]);

  double oldAngle = 
    atan2((double)rwi->GetLastEventPosition()[1] - (double)disp_obj_center[1],
          (double)rwi->GetLastEventPosition()[0] - (double)disp_obj_center[0]);

#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 3))
  newAngle *= vtkMath::DegreesFromRadians(1.f);
  oldAngle *= vtkMath::DegreesFromRadians(1.f);
#else
  newAngle *= vtkMath::RadiansToDegrees();
  oldAngle *= vtkMath::RadiansToDegrees();
#endif
  
  double scale[3];
  scale[0] = scale[1] = scale[2] = 1.0;

  double **rotate = new double*[1];
  rotate[0] = new double[4];
  
  rotate[0][0] = newAngle - oldAngle;
  rotate[0][1] = motion_vector[0];
  rotate[0][2] = motion_vector[1];
  rotate[0][3] = motion_vector[2];
  
  this->Prop3DTransform(this->InteractionProp,
                        obj_center,
                        1, 
                        rotate, 
                        scale);
  
  delete [] rotate[0];
  delete [] rotate;
  
  if (this->AutoAdjustCameraClippingRange)
    {
    this->CurrentRenderer->ResetCameraClippingRange();
    }

  rwi->Render();
}

void
vtkOrientMRIInteractorStyleView3D::RotateCamera () {

  // Copied from vtkInteractorStyleTrackballCamera.
  
  if (this->CurrentRenderer == NULL)
    {
    return;
    }

  vtkRenderWindowInteractor *rwi = this->Interactor;

  int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  
  int *size = this->CurrentRenderer->GetRenderWindow()->GetSize();

  double delta_elevation = -20.0 / size[1];
  double delta_azimuth = -20.0 / size[0];
  
  double rxf = (double)dx * delta_azimuth * this->MotionFactor;
  double ryf = (double)dy * delta_elevation * this->MotionFactor;
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->Azimuth(rxf);
  camera->Elevation(ryf);
  camera->OrthogonalizeViewUp();

  if (this->AutoAdjustCameraClippingRange)
    {
    this->CurrentRenderer->ResetCameraClippingRange();
    }

  if (rwi->GetLightFollowCamera()) 
    {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
    }

  rwi->Render();
}

void
vtkOrientMRIInteractorStyleView3D::SpinCamera () {

  // Copied from vtkInteractorStyleTrackballCamera.
  
  if (this->CurrentRenderer == NULL)
    {
    return;
    }

  vtkRenderWindowInteractor *rwi = this->Interactor;

  double *center = this->CurrentRenderer->GetCenter();

  double newAngle = 
    atan2((double)rwi->GetEventPosition()[1] - (double)center[1],
          (double)rwi->GetEventPosition()[0] - (double)center[0]);

  double oldAngle = 
    atan2((double)rwi->GetLastEventPosition()[1] - (double)center[1],
          (double)rwi->GetLastEventPosition()[0] - (double)center[0]);
  
#if ((VTK_MAJOR_VERSION == 5)&&(VTK_MINOR_VERSION > 3))
  newAngle *= vtkMath::DegreesFromRadians(1.f);
  oldAngle *= vtkMath::DegreesFromRadians(1.f);
#else
  newAngle *= vtkMath::RadiansToDegrees();
  oldAngle *= vtkMath::RadiansToDegrees();
#endif

  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->Roll(newAngle - oldAngle);
  camera->OrthogonalizeViewUp();
      
  rwi->Render();
}

void 
vtkOrientMRIInteractorStyleView3D::DollyCamera () {

  // Copied from vtkInteractorStyleTrackballCamera.
  
  if (this->CurrentRenderer == NULL)
    {
    return;
    }
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  double *center = this->CurrentRenderer->GetCenter();
  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  double dyf = this->MotionFactor * (double)(dy) / (double)(center[1]);
  this->DollyCamera(pow((double)1.1, dyf));
}

void
vtkOrientMRIInteractorStyleView3D::DollyCamera ( double factor ) {

  // Copied from vtkInteractorStyleTrackballCamera.
  
  if (this->CurrentRenderer == NULL)
    {
    return;
    }
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  if (camera->GetParallelProjection())
    {
    camera->SetParallelScale(camera->GetParallelScale() / factor);
    }
  else
    {
    camera->Dolly(factor);
    if (this->AutoAdjustCameraClippingRange)
      {
      this->CurrentRenderer->ResetCameraClippingRange();
      }
    }
  
  if (this->Interactor->GetLightFollowCamera()) 
    {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
    }
  
  this->Interactor->Render();
}

void
vtkOrientMRIInteractorStyleView3D::Prop3DTransform ( vtkProp3D *prop3D,
						     double *boxCenter,
						     int numRotation,
						     double **rotate,
						     double *scale ) {
  // Copied from vtkInteractorStyleTrackballActor.
  
  vtkMatrix4x4 *oldMatrix = vtkMatrix4x4::New();
  prop3D->GetMatrix(oldMatrix);
  
  double orig[3];
  prop3D->GetOrigin(orig);
  
  vtkTransform *newTransform = vtkTransform::New();
  newTransform->PostMultiply();
  if (prop3D->GetUserMatrix() != NULL) 
    {
    newTransform->SetMatrix(prop3D->GetUserMatrix());
    }
  else 
    {
    newTransform->SetMatrix(oldMatrix);
    }
  
  newTransform->Translate(-(boxCenter[0]), -(boxCenter[1]), -(boxCenter[2]));
  
  for (int i = 0; i < numRotation; i++) 
    {
    newTransform->RotateWXYZ(rotate[i][0], rotate[i][1],
                             rotate[i][2], rotate[i][3]);
    }
  
  if ((scale[0] * scale[1] * scale[2]) != 0.0) 
    {
    newTransform->Scale(scale[0], scale[1], scale[2]);
    }
  
  newTransform->Translate(boxCenter[0], boxCenter[1], boxCenter[2]);
  
  // now try to get the composit of translate, rotate, and scale
  newTransform->Translate(-(orig[0]), -(orig[1]), -(orig[2]));
  newTransform->PreMultiply();
  newTransform->Translate(orig[0], orig[1], orig[2]);
  
  if (prop3D->GetUserMatrix() != NULL) 
    {
    newTransform->GetMatrix(prop3D->GetUserMatrix());
    }
  else 
    {
    prop3D->SetPosition(newTransform->GetPosition());
    prop3D->SetScale(newTransform->GetScale());
    prop3D->SetOrientation(newTransform->GetOrientation());
    }
  oldMatrix->Delete();
  newTransform->Delete();
}

void
vtkOrientMRIInteractorStyleView3D::StartRotateProp () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_RotateProp );
}

void
vtkOrientMRIInteractorStyleView3D::EndRotateProp () {

  if( this->State != IS_RotateProp )
    return;

  this->StopState();
}

void
vtkOrientMRIInteractorStyleView3D::StartSpinProp () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_SpinProp );
}

void
vtkOrientMRIInteractorStyleView3D::EndSpinProp () {

  if( this->State != IS_SpinProp )
    return;

  this->StopState();
}

void
vtkOrientMRIInteractorStyleView3D::StartRotateCamera () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_RotateCamera );
}

void
vtkOrientMRIInteractorStyleView3D::EndRotateCamera () {

  if( this->State != IS_RotateCamera )
    return;

  this->StopState();
}

void
vtkOrientMRIInteractorStyleView3D::StartSpinCamera () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_SpinCamera );
}

void
vtkOrientMRIInteractorStyleView3D::EndSpinCamera () {

  if( this->State != IS_SpinCamera )
    return;

  this->StopState();
}

void
vtkOrientMRIInteractorStyleView3D::StartDollyCamera () {

  if( this->State != VTKIS_NONE )
    return;

  this->StartState( IS_DollyCamera );
}

void
vtkOrientMRIInteractorStyleView3D::EndDollyCamera () {

  if( this->State != IS_DollyCamera )
    return;

  this->StopState();
}
