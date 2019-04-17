/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkInteractorStyleMyTrackballCamera.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkInteractorStyleMyTrackballCamera.h"

#include "vtkCamera.h"
#include "vtkCallbackCommand.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"

vtkStandardNewMacro(vtkInteractorStyleMyTrackballCamera);

//----------------------------------------------------------------------------
vtkInteractorStyleMyTrackballCamera::vtkInteractorStyleMyTrackballCamera()
{
  m_bRotateAroundPoint = false;
}

vtkInteractorStyleMyTrackballCamera::~vtkInteractorStyleMyTrackballCamera () {}

void vtkInteractorStyleMyTrackballCamera::SetRotateByPoint(bool b, double *dPos)
{
  m_bRotateAroundPoint = b;
  if (dPos)
  {
    m_dCenterPoint[0] = dPos[0];
    m_dCenterPoint[1] = dPos[1];
    m_dCenterPoint[2] = dPos[2];
  }
}

//----------------------------------------------------------------------------
void vtkInteractorStyleMyTrackballCamera::Rotate()
{
  if (this->CurrentRenderer == nullptr)
  {
    return;
  }

  //  vtkRenderWindowInteractor *rwi = this->Interactor;

  //  int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
  //  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];

  //  int *size = this->CurrentRenderer->GetRenderWindow()->GetSize();

  //  double delta_elevation = -20.0 / size[1];
  //  double delta_azimuth = -20.0 / size[0];

  //  double rxf = dx * delta_azimuth * this->MotionFactor;
  //  double ryf = dy * delta_elevation * this->MotionFactor;

  //  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  //  camera->Azimuth(rxf);
  //  camera->Elevation(ryf);
  //  camera->OrthogonalizeViewUp();

  //  if (this->AutoAdjustCameraClippingRange)
  //  {
  //    this->CurrentRenderer->ResetCameraClippingRange();
  //  }

  //  if (rwi->GetLightFollowCamera())
  //  {
  //    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  //  }

  //  rwi->Render();

  if (m_bRotateAroundPoint)
  {
    vtkRenderWindowInteractor *rwi = this->Interactor;;
    vtkRenderer *renderer = this->CurrentRenderer;
    vtkCamera *camera = renderer->GetActiveCamera();

    int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
    int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];

    int *size = this->CurrentRenderer->GetRenderWindow()->GetSize();

    double delta_elevation = -20.0 / size[1];
    double delta_azimuth = -20.0 / size[0];

    double rxf = dx * delta_azimuth * this->MotionFactor;
    double ryf = dy * delta_elevation * this->MotionFactor;

    // Camera Parameters ///////////////////////////////////////////////////
    double *focalPoint = camera->GetFocalPoint();
    double *viewUp = camera->GetViewUp();
    double *position = camera->GetPosition();
    double axis[3];
    axis[0] = -camera->GetViewTransformMatrix()->GetElement(0,0);
    axis[1] = -camera->GetViewTransformMatrix()->GetElement(0,1);
    axis[2] = -camera->GetViewTransformMatrix()->GetElement(0,2);

    // Build The transformatio /////////////////////////////////////////////////
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();

    transform->Translate(m_dCenterPoint[0], m_dCenterPoint[1], m_dCenterPoint[2]);
    transform->RotateWXYZ(rxf, viewUp); // Azimuth
    transform->RotateWXYZ(ryf, axis);   // Elevation
    transform->Translate(-m_dCenterPoint[0], -m_dCenterPoint[1], -m_dCenterPoint[2]);

    double newPosition[3];
    transform->TransformPoint(position,newPosition); // Transform Position
    double newFocalPoint[3];
    transform->TransformPoint(focalPoint, newFocalPoint); // Transform Focal Point

    camera->SetPosition(newPosition);
    camera->SetFocalPoint(newFocalPoint);

    // Orhthogonalize View Up //////////////////////////////////////////////////
    camera->OrthogonalizeViewUp();
    //    renderer->ResetCameraClippingRange();

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
  else
    vtkInteractorStyleTrackballCamera::Rotate();
}

