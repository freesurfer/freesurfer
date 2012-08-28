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
#include <QFileInfo>
#include "vtkRenderer.h"
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>

SurfaceSpline::SurfaceSpline(LayerSurface *parent) :
    QObject(parent),
    m_mri(NULL),
    m_mriSurf(NULL),
    m_nActiveVertex(-1),
    m_bProjection(true)
{
  m_actor = vtkSmartPointer<vtkActor>::New();
  m_actorSpheres = vtkSmartPointer<vtkActor>::New();
  for (int i = 0; i < 3; i++)
  {
    m_actor2D[i] = vtkSmartPointer<vtkActor>::New();
    m_actor2DSpheres[i] = vtkSmartPointer<vtkActor>::New();
  }
  SetColor(QColor(255, 128, 0));
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
  if (m_mri)
    ::MRIfree(&m_mri);

  if (!m_mriSurf)
  {
    MRIS* mris = qobject_cast<LayerSurface*>(parent())->GetSourceSurface()->GetMRIS();
    m_mriSurf = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR, 1);
    useVolGeomToMRI(&mris->vg, m_mriSurf);
  }

  m_mri = ::MRIread( filename.toAscii().data() );      // could be long process

  if ( m_mri == NULL )
  {
    return false;
  }

  m_strName = QFileInfo(filename).completeBaseName();
  RebuildActors();

  return true;
}

void SurfaceSpline::AppendProp2D(vtkRenderer *ren, int nPlane)
{
  ren->AddViewProp(m_actor2D[nPlane]);
  ren->AddViewProp(m_actor2DSpheres[nPlane]);
}

void SurfaceSpline::AppendProp3D(vtkRenderer *ren)
{
  ren->AddViewProp(m_actor);
  ren->AddViewProp(m_actorSpheres);
}

void SurfaceSpline::SetActiveVertex(int n)
{
  m_nActiveVertex = n;
  RebuildActors();
}

void SurfaceSpline::SetProjection(bool bProjection)
{
  m_bProjection = bProjection;
  RebuildActors();
}

void SurfaceSpline::BuildSphereActor(vtkActor* actor, vtkPoints* pts)
{
  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  for ( int i = 0; i < pts->GetNumberOfPoints(); i++ )
  {
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
    double pt[3];
    pts->GetPoint(i, pt);
    sphere->SetCenter( pt );
    sphere->SetRadius( 0.625 );
    sphere->SetThetaResolution( 10 );
    sphere->SetPhiResolution( 20 );
    append->AddInput( sphere->GetOutput() );
  }
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( append->GetOutput() );
  actor->SetMapper(mapper);
}


void SurfaceSpline::RebuildActors()
{
  if (!m_mri)
    return;
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  m_actor->SetMapper(mapper);
  m_actorSpheres->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
  vtkSmartPointer<vtkPolyDataMapper> mapper2d[3];
  for (int i = 0; i < 3; i++)
  {
    mapper2d[i] = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_actor2D[i]->SetMapper(mapper2d[i]);
    m_actor2DSpheres[i]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
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
    double pos[3];
    for (int i = 0; i < m_mri->height; i++)
    {
      Real x = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 0, 0);
      Real y = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 1, 0);
      Real z = ::MRIgetVoxVal(m_mri, m_nActiveVertex, i, 2, 0);
      if (i == 0 && x == 0 && y == 0 && z == 0)
      {
        emit SplineChanged();
        return;
      }
    //  qDebug() << x << y << z;
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
      BuildSphereActor(m_actorSpheres, points);
    }
    for (int n = 0; n < 3; n++)
    {
      vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(points);
      polydata->SetLines(lines);
      vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
      spline->SetInput(polydata);
      spline->Update();
      vtkPolyData* spline_poly = spline->GetOutput();
      vtkPoints* spline_points = spline_poly->GetPoints();
      vtkSmartPointer<vtkPoints> ctrl_points = vtkSmartPointer<vtkPoints>::New();
      for (int i = 0; i < points->GetNumberOfPoints(); i++)
      {
        double pt[3];
        points->GetPoint(i, pt);
        if (m_bProjection)
          pt[n] = slice_pos[n];
        ctrl_points->InsertNextPoint(pt);
      }
      if (m_bProjection)
      {
        for (int i = 0; i < spline_points->GetNumberOfPoints(); i++)
        {
          double pt[3];
          spline_points->GetPoint(i, pt);
          pt[n] = slice_pos[n];
          spline_points->SetPoint(i, pt);
        }
      }

      vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
      tube->SetInput(spline_poly);
      tube->SetNumberOfSides(8);
      tube->SetRadius(0.25);
      mapper2d[n]->SetInput(tube->GetOutput());
      BuildSphereActor(m_actor2DSpheres[n], ctrl_points);
    }
  }

  emit SplineChanged();
}

void SurfaceSpline::SetVisible(bool visible)
{
  m_actor->SetVisibility(visible);
  m_actorSpheres->SetVisibility(visible);
  for (int i = 0; i < 3; i++)
  {
    m_actor2D[i]->SetVisibility(visible);
    m_actor2DSpheres[i]->SetVisibility(visible);
  }

  emit SplineChanged();
}

bool SurfaceSpline::IsVisible()
{
  return m_actor->GetVisibility();
}

void SurfaceSpline::SetColor(const QColor &c)
{
  m_color = c;
  QColor c2 = c.lighter();
  m_actor->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  m_actorSpheres->GetProperty()->SetColor(c2.redF(), c2.greenF(), c2.blueF());
  for (int i = 0; i < 3; i++)
  {
    m_actor2D[i]->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
    m_actor2DSpheres[i]->GetProperty()->SetColor(c2.redF(), c2.greenF(), c2.blueF());
  }
  emit SplineChanged();
}
