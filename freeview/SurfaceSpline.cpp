#include "SurfaceSpline.h"
#include "vtkActor.h"
#include "vtkSplineFilter.h"
#include "vtkTubeFilter.h"
#include "vtkPolyDataMapper.h"
#include "LayerSurface.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include <QDebug>
#include "FSSurface.h"

SurfaceSpline::SurfaceSpline(LayerSurface *parent) :
    QObject(parent),
    m_mri(NULL),
    m_nActiveVertex(-1)
{
  m_actor = vtkSmartPointer<vtkActor>::New();
  m_actor->GetProperty()->SetColor(1.0, 0.5, 0);
  for (int i = 0; i < 3; i++)
  {
    m_actor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_actor2D[i]->GetProperty()->SetColor(1.0, 0.5, 0);
  }

  MRIS* mris = parent->GetSourceSurface()->GetMRIS();
  m_mriSurf = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR, 1);
  useVolGeomToMRI(&mris->vg, m_mriSurf);
}

SurfaceSpline::~SurfaceSpline()
{
  if (m_mri)
    ::MRIfree(&m_mri);
  if (m_mriSurf)
    ::MRIfree(&m_mriSurf);
}

bool SurfaceSpline::Load(const QString &filename)
{
  m_mri = ::MRIread( filename.toAscii().data() );      // could be long process

  if ( m_mri == NULL )
  {
    return false;
  }

  return true;
}

vtkActor* SurfaceSpline::GetActor()
{
  return m_actor;
}

vtkActor* SurfaceSpline::GetActor2D(int nPlane)
{
  return m_actor2D[nPlane];
}

void SurfaceSpline::SetActiveVertex(int n)
{
  m_nActiveVertex = n;
  RebuildActors();
}

void SurfaceSpline::RebuildActors()
{
  if (!m_mri)
    return;
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  m_actor->SetMapper(mapper);
  vtkSmartPointer<vtkPolyDataMapper> mapper2d[3];
  for (int i = 0; i < 3; i++)
  {
    mapper2d[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_actor2D[i]->SetMapper(mapper2d[i]);
  }

  LayerSurface* surf = qobject_cast<LayerSurface*>(parent());
  double slice_pos[3];
  surf->GetSlicePosition(slice_pos);
  LayerMRI* mri = surf->GetRefVolume();
  if (m_nActiveVertex >= 0)
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(m_mri->height);
    qDebug() << m_mri->height;
    double pos[3];
    for (int i = 0; i < m_mri->height; i++)
    {
      Real x = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 0, 0);
      Real y = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 1, 0);
      Real z = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 2, 0);
      if (i == 0 && x == 0 && y == 0 && z == 0)
        return;
      qDebug() << x << y << z;
      ::MRIvoxelToWorld( m_mriSurf, x, y, z, &x, &y, &z );
      pos[0] = x; pos[1] = y; pos[2] = z;
      if (mri)
      {
        mri->RASToTarget(pos, pos);
      }
      points->InsertNextPoint(pos);
      lines->InsertNextCell(i);
    //  qDebug() << pos[0] << pos[1] << pos[2];
    }
    {
      vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(points);
      polydata->SetLines(lines);
      vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
      spline->SetInput(polydata);
      vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
      tube->SetInput(spline->GetOutput());
      tube->SetNumberOfSides(8);
      tube->SetRadius(0.25);
      mapper->SetInput(tube->GetOutput());
    }
    for (int n = 0; n < 3; n++)
    {
      vtkSmartPointer<vtkPoints> slice_points = vtkSmartPointer<vtkPoints>::New();
      for (int i = 0; i < points->GetNumberOfPoints(); i++)
      {
        double pt[3];
        points->GetPoint(i, pt);
        pt[n] = slice_pos[n];
        slice_points->InsertNextPoint(pt);
      }
      vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(slice_points);
      polydata->SetLines(lines);
      vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
      spline->SetInput(polydata);
      vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
      tube->SetInput(spline->GetOutput());
      tube->SetNumberOfSides(8);
      tube->SetRadius(0.25);
      mapper2d[n]->SetInput(tube->GetOutput());
    }
  }
}

void SurfaceSpline::SetVisible(bool visible)
{
  m_actor->SetVisibility(visible);
  for (int i = 0; i < 3; i++)
    m_actor2D[i]->SetVisibility(visible);
}
