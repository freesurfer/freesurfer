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

#include "LayerEditRef.h"
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
#include "vtkMath.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkAppendPolyData.h"
#include "vtkCubeSource.h"
#include "BrushProperty.h"

LayerEditRef::LayerEditRef(QObject *parent) :
  LayerEditable(parent), m_mriRef(NULL)
{
  this->m_strTypeNames << "Supplement" << "EditRef";
  m_points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < 3; i++)
  {
    m_actorSlice[i] = vtkSmartPointer<vtkActor>::New();
    m_actorSlice[i]->GetProperty()->SetInterpolationToFlat();
    m_actorSlice[i]->GetProperty()->SetAmbient( 1 );
    m_actorSlice[i]->GetProperty()->SetDiffuse( 0 );
    m_actorSlice[i]->GetProperty()->SetOpacity(0.7);
//    m_actorSlice[i]->VisibilityOff();
  }
  m_actor = vtkSmartPointer<vtkActor>::New();
//  m_actor->VisibilityOff();
}

LayerEditRef::~LayerEditRef()
{

}

void LayerEditRef::SetMRIRef(LayerMRI *mri)
{
  m_mriRef = mri;
}

void LayerEditRef::Append2DProps(vtkRenderer *renderer, int nPlane)
{
  if ( nPlane < 3 && nPlane >= 0 )
  {
    renderer->AddViewProp(m_actorSlice[nPlane]);
  }
}

void LayerEditRef::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  Q_UNUSED(bPlaneVisibility);
  renderer->AddViewProp(m_actor);
}

bool LayerEditRef::HasProp(vtkProp *prop)
{
  return ((vtkProp*)m_actor == prop);
}

bool LayerEditRef::IsVisible()
{
  return m_actor->GetVisibility() > 0;
}

void LayerEditRef::SetVisible(bool bVisible)
{
  for (int i = 0; i < 3; i++)
  {
    m_actorSlice[i]->SetVisibility(bVisible?1:0);
  }
  m_actor->SetVisibility(bVisible?1:0);
  LayerEditable::SetVisible(bVisible);
}

void LayerEditRef::OnSlicePositionChanged(int nPlane)
{
  Q_UNUSED(nPlane);
  UpdateActors(false); // no need to rebuild 3D actors
}

void LayerEditRef::SetColor(const QColor &color)
{
  m_actor->GetProperty()->SetColor(color.redF(), color.greenF(), color.blueF());
  for (int i = 0; i < 3; i++)
    m_actorSlice[i]->GetProperty()->SetColor(color.redF(), color.greenF(), color.blueF());
  emit ActorUpdated();
}

void LayerEditRef::SetStartPosition(double *pos)
{
  m_points->Reset();
  m_points->InsertNextPoint(pos);
  UpdateActors();
}

void LayerEditRef::SetEndPosition(double* pos)
{
  m_points->InsertNextPoint(pos);
  UpdateActors();
}

void LayerEditRef::UpdateActors(bool bBuild3D)
{
  int dim[3] = {0, 0, 0};
  double voxel_size[3] = {1, 1, 1};
  double origin[3] = {0};
  if (m_mriRef)
  {
    m_mriRef->GetImageData()->GetSpacing(voxel_size);
    m_mriRef->GetImageData()->GetOrigin(origin);
    m_mriRef->GetImageData()->GetDimensions(dim);
  }

  int n1[3];
  if (m_points->GetNumberOfPoints() == 0)
  {
    Reset();
    return;
  }

  m_voxels.clear();
  double pos[3], pos1[3];
  POINT pt;
  m_points->GetPoint(0, pos);
  if (m_points->GetNumberOfPoints() > 1)
  {
    m_points->GetPoint(1, pos1);
    for ( int i = 0; i < 3; i++ )
    {
      n1[i] = (int)((pos[i] - origin[i]) / voxel_size[i] + 0.5);
    }
  }
  else
  {
    pos1[0] = pos[0];
    pos1[1] = pos[1];
    pos1[2] = pos[2];
  }

  int nBrushSize = 1;
  if (m_mriRef)
    nBrushSize = m_mriRef->GetBrushProperty()->GetBrushSize();

  double dstep = qMin(voxel_size[0], qMin(voxel_size[1], voxel_size[2]))/4;
  double dist = sqrt(vtkMath::Distance2BetweenPoints(pos, pos1));
  double v[3];
  for (int i = 0; i < 3; i++)
    v[i] = pos1[i]-pos[i];

  double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if (norm > 0)
  {
    for (int i = 0; i < 3; i++)
      v[i] /= norm;
  }
  else
    v[0] = 1;

  while (dist >= 0)
  {
    for (int i = 0; i < 3; i++)
    {
      pt.n[i] = (int)((pos[i] - origin[i]) / voxel_size[i] + 0.5);
    }
    if (!m_voxels.contains(pt))
    {
      m_voxels << pt; // qDebug() << pt.n[0] << pt.n[1] << pt.n[2];
    }

    for (int i = 0; i < 3; i++)
    {
      pos[i] += v[i]*dstep;
    }
    dist -= dstep;
  }

  if (nBrushSize > 1)
  {
    QVector<POINT> old_vox = m_voxels;
    foreach (POINT pt_in, old_vox)
    {
      for (int i = -nBrushSize+1; i < nBrushSize; i++)
      {
        for (int j = -nBrushSize+1; j < nBrushSize; j++)
        {
          for (int k = -nBrushSize+1; k < nBrushSize; k++)
          {
            if (i*i+j*j+k*k <= (nBrushSize-1)*(nBrushSize-1))
            {
              POINT pt = pt_in;
              pt.n[0] += i;
              pt.n[1] += j;
              pt.n[2] += k;
              for (int ii = 0; ii < 3; ii++)
              {
                if (pt.n[ii] < 0)
                  pt.n[ii] = 0;
                else if (pt.n[ii] >= dim[ii])
                  pt.n[ii] = dim[ii]-1;

                if (!m_voxels.contains(pt))
                  m_voxels << pt;
              }
            }
          }
        }
      }
    }
  }
//  qDebug() << "n vox" << m_voxels.size() << pos1[0] << pos1[1] << pos1[2];

  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  for (int i = 0; i < m_voxels.size(); i++)
  {
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    cube->SetXLength(voxel_size[0]);
    cube->SetYLength(voxel_size[1]);
    cube->SetZLength(voxel_size[2]);
    for (int j = 0; j < 3; j++)
      pos[j] = m_voxels[i].n[j]*voxel_size[j] + origin[j];
    cube->SetCenter(pos);
    append->AddInputConnection(cube->GetOutputPort());
  }

  if (bBuild3D)
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(append->GetOutputPort());
    m_actor->SetMapper(mapper);
  }

  for ( int i = 0; i < 3; i++ )
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin( m_dSlicePosition );
    plane->SetNormal( i==0?1:0, i==1?1:0, i==2?1:0 );

    vtkSmartPointer<vtkCutter> cutter =
        vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputConnection(append->GetOutputPort());
    cutter->SetCutFunction( plane );

    vtkSmartPointer<vtkStripper> stripper = vtkSmartPointer<vtkStripper>::New();
    stripper->SetInputConnection(cutter->GetOutputPort());
    stripper->Update();

    vtkSmartPointer<vtkPolyData> cutpoly = vtkSmartPointer<vtkPolyData>::New();
    cutpoly->SetPoints(stripper->GetOutput()->GetPoints());
    cutpoly->SetPolys(stripper->GetOutput()->GetLines());

    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
#if VTK_MAJOR_VERSION > 5
    triangleFilter->SetInputData( cutpoly );
#else
    triangleFilter->SetInput( cutpoly );
#endif
//    mapper->SetInputConnection( triangleFilter->GetOutputPort() );
    mapper->SetInputData(cutpoly);
//    mapper->SetInputConnection( stripper->GetOutputPort() );

    m_actorSlice[i]->SetMapper( mapper );
    switch ( i )
    {
    case 0:
      m_actorSlice[i]->SetPosition(0.1, 0, 0);
      break;
    case 1:
      m_actorSlice[i]->SetPosition(0, m_dTinyOffset, 0);
      break;
    case 2:
      m_actorSlice[i]->SetPosition(0, 0, -m_dTinyOffset);
      break;
    }
  }

  emit ActorUpdated();
}

int LayerEditRef::GetNumberOfMarks()
{
  return m_points->GetNumberOfPoints();
}

void LayerEditRef::Reset()
{
  m_points->Reset();
  m_voxels.clear();
  for (int i = 0; i < 3; i++)
    m_actorSlice[i]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());

  m_actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
  emit ActorUpdated();
}

void LayerEditRef::ApplyToMRI(LayerMRI *mri_in)
{
  LayerMRI* mri = mri_in;
  if (!mri)
    mri = m_mriRef;
  if (!mri)
  {
    qDebug() << "Refrence volume not set";
    return;
  }

  mri->SaveForUndo();
  foreach (POINT vox, m_voxels)
  {
    mri->SetVoxelByIndex(vox.n, 0, true, true);
  }
  mri->MarkDataModified();
  if (mri->GetProperty()->GetShowAsContour() && mri->GetProperty()->GetShowVoxelizedContour())
  {
    mri->RebuildContour(mri->GetFillValue());
  }
  Reset();
}
