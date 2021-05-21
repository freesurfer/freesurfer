/**
 * @brief Layer class for structural landmarks.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#include "LayerLandmarks.h"
#include <QDebug>
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSphereSource.h"
#include "vtkRenderer.h"
#include "vtkCutter.h"
#include "vtkStripper.h"
#include "vtkTriangleFilter.h"
#include "vtkPlane.h"
#include "LayerMRI.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

Landmark::Landmark()
{
  actorSphere = vtkSmartPointer<vtkActor>::New();
  for (int i = 0; i < 3; i++)
  {
    actorSlice[i] = vtkSmartPointer<vtkActor>::New();
    actorSlice[i]->GetProperty()->SetInterpolationToFlat();
    actorSlice[i]->GetProperty()->SetAmbient( 1 );
    actorSlice[i]->GetProperty()->SetDiffuse( 0 );
    actorSlice[i]->GetProperty()->SetOpacity(0.7);
    actorSlice[i]->VisibilityOff();
  }
}

LayerLandmarks::LayerLandmarks(QObject *parent) :
  LayerEditable(parent),
  m_dRadius(3.0)
{
  this->m_strTypeNames << "Supplement" << "Landmarks";
}

LayerLandmarks::~LayerLandmarks()
{

}

void LayerLandmarks::SetMRIRef(LayerMRI *mri)
{
  m_mriRef = mri;
}

void LayerLandmarks::Append2DProps(vtkRenderer *renderer, int nPlane)
{
  if ( nPlane < 3 && nPlane >= 0 )
  {
    for (int i = 0; i < m_landmarks.size(); i++)
      renderer->AddViewProp( m_landmarks[i].actorSlice[nPlane] );
  }
}

void LayerLandmarks::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  Q_UNUSED(bPlaneVisibility);
  for (int i = 0; i < m_landmarks.size(); i++)
    renderer->AddViewProp(m_landmarks[i].actorSphere);
}

bool LayerLandmarks::HasProp(vtkProp *prop)
{
  for (int i = 0; i < m_landmarks.size(); i++)
  {
    if ( (vtkProp*)m_landmarks[i].actorSphere == prop)
      return true;
  }
  return false;
}

bool LayerLandmarks::IsVisible()
{
  if (m_landmarks.isEmpty())
    return false;
  else
    return m_landmarks[0].actorSphere->GetVisibility() > 0;
}

void LayerLandmarks::SetVisible(bool bVisible)
{
  for (int i = 0; i < m_landmarks.size(); i++)
  {
    m_landmarks[i].actorSlice[0]->SetVisibility(bVisible);
    m_landmarks[i].actorSlice[1]->SetVisibility(bVisible);
    m_landmarks[i].actorSlice[2]->SetVisibility(bVisible);
    m_landmarks[i].actorSphere->SetVisibility(bVisible);
  }
  LayerEditable::SetVisible(bVisible);
}


void LayerLandmarks::OnSlicePositionChanged(int nPlane)
{
  Q_UNUSED(nPlane);
  UpdateActors(false); // no need to rebuild 3D actors
}

void LayerLandmarks::SetLandmarkColor(int n, const QColor &color)
{
  MakeSureLandmarkExist(n);
  m_landmarks[n].color = color;
  m_landmarks[n].actorSphere->GetProperty()->SetColor(color.redF(), color.greenF(), color.blueF());
  for (int i = 0; i < 3; i++)
    m_landmarks[n].actorSlice[i]->GetProperty()->SetColor(color.redF(), color.greenF(), color.blueF());
  emit ActorUpdated();
}

Landmark& LayerLandmarks::GetLandmark(int n)
{
  MakeSureLandmarkExist(n);
  return m_landmarks[n];
}

void LayerLandmarks::SetLandmarkPosition(int n, double *pos)
{
  SetLandmarkPosition(n, pos[0], pos[1], pos[2]);
}

void LayerLandmarks::SetLandmarkPosition(int n, double x, double y, double z)
{
  MakeSureLandmarkExist(n);
  m_landmarks[n].pos[0] = x;
  m_landmarks[n].pos[1] = y;
  m_landmarks[n].pos[2] = z;
  m_landmarks[n].valid = true;
  UpdateActors();
}

void LayerLandmarks::SetRadius(double dRadius)
{
  m_dRadius = dRadius;
  UpdateActors();
}

bool LayerLandmarks::MakeSureLandmarkExist(int n)
{
  bool expanded = false;
  for (int i = m_landmarks.size()-1; i < n; i++)
  {
    m_landmarks << Landmark();
    expanded = true;
  }
  if (expanded)
    emit LandmarkAdded();
  return expanded;
}

void LayerLandmarks::UpdateActors(bool bBuild3D)
{
  double scale = 1;
  double voxel_size[3] = {1, 1, 1};
  if (m_mriRef)
  {
    m_mriRef->GetImageData()->GetSpacing(voxel_size);
    scale = qMin(voxel_size[0], qMin(voxel_size[1], voxel_size[2]));
  }
  for (int i = 0; i < m_landmarks.size(); i++)
  {
    if (bBuild3D)
    {
      vtkSmartPointer<vtkSphereSource> ball = vtkSmartPointer<vtkSphereSource>::New();
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      ball->SetRadius(m_dRadius);
      mapper->SetInputConnection(ball->GetOutputPort());
      m_landmarks[i].actorSphere->SetMapper(mapper);
      m_landmarks[i].actorSphere->SetPosition(m_landmarks[i].pos);
    }
    for ( int j = 0; j < 3; j++ )
    {
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      if ( m_landmarks[i].valid && fabs( m_dSlicePosition[j] - m_landmarks[i].pos[j] ) < voxel_size[j]/2 )
      {
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        double point[3] = { m_landmarks[i].pos[0],
                            m_landmarks[i].pos[1],
                            m_landmarks[i].pos[2] };
        point[j] = m_dSlicePosition[j];
        sphere->SetCenter( point );
        sphere->SetRadius( m_dRadius * scale );
        sphere->SetThetaResolution( 12 );
        sphere->SetPhiResolution( 24 );

        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetOrigin( point );
        plane->SetNormal( j==0?1:0, j==1?1:0, j==2?1:0 );

        vtkSmartPointer<vtkCutter> cutter =
            vtkSmartPointer<vtkCutter>::New();
        cutter->SetInputConnection( sphere->GetOutputPort() );
        cutter->SetCutFunction( plane );

        vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
        stripper->SetInputConnection( cutter->GetOutputPort() );
        stripper->Update();

        vtkSmartPointer<vtkPolyData> cutpoly = vtkSmartPointer<vtkPolyData>::New();
        cutpoly->SetPoints( stripper->GetOutput()->GetPoints() );
        cutpoly->SetPolys( stripper->GetOutput()->GetLines() );

        vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
#if VTK_MAJOR_VERSION > 5
        triangleFilter->SetInputData( cutpoly );
#else
        triangleFilter->SetInput( cutpoly );
#endif
        mapper->SetInputConnection( triangleFilter->GetOutputPort() );
      }
      else
#if VTK_MAJOR_VERSION > 5
        mapper->SetInputData( vtkSmartPointer<vtkPolyData>::New() );
#else
        mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
#endif

      m_landmarks[i].actorSlice[j]->SetMapper( mapper );
    }
  }

  emit ActorUpdated();
}

void LayerLandmarks::DoTransform(double *mat, int sample_method)
{
  Q_UNUSED(sample_method);
  if (m_landmarksOriginal.isEmpty())
    m_landmarksOriginal = m_landmarks;

  for (int i = 0; i < m_landmarks.size(); i++)
  {
    double pos[4] = {m_landmarks[i].pos[0],
                     m_landmarks[i].pos[1],
                     m_landmarks[i].pos[2],
                     1};
    vtkMatrix4x4::MultiplyPoint(mat, pos, pos);
    m_landmarks[i].pos[0] = pos[0];
    m_landmarks[i].pos[1] = pos[1];
    m_landmarks[i].pos[2] = pos[2];
  }
  UpdateActors();
}

void LayerLandmarks::DoRestore()
{
  if (!m_landmarksOriginal.isEmpty())
  {
    m_landmarks = m_landmarksOriginal;
    UpdateActors();
  }
}
