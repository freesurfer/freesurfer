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
 */
#include "LayerTrack.h"
#include "FSTrack.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "LayerPropertyTrack.h"
#include <QFileInfo>
#include <QDir>
#include <QDebug>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkLODActor.h>
#include "MyUtils.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkDataArray.h"

LayerTrack::LayerTrack(LayerMRI* ref, QObject* parent, bool bCluster) : Layer(parent),
  m_trackData(0),
  m_layerMRIRef(ref)
{
  this->m_strTypeNames << "Tract";

  mProperty = new LayerPropertyTrack( this );
  connect(mProperty, SIGNAL(ColorCodeChanged(int)), this, SLOT(RebuildActors()));
  connect(mProperty, SIGNAL(DirectionSchemeChanged(int)), this, SLOT(RebuildActors()));
  connect(mProperty, SIGNAL(DirectionMappingChanged(int)), this, SLOT(RebuildActors()));
  connect(mProperty, SIGNAL(SolidColorChanged(QColor)), this, SLOT(UpdateColor()));
  connect(mProperty, SIGNAL(ScalarColorMapChanged(int)), this, SLOT(UpdateColor()));
  connect(mProperty, SIGNAL(ScalarIndexChanged(int)), this, SLOT(UpdateColor()));
  connect(mProperty, SIGNAL(ScalarThresholdChanged(double, double)), this, SLOT(UpdateColor()));
  connect(mProperty, SIGNAL(RenderRepChanged()), this, SLOT(RebuildActors()));
  connect(mProperty, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity(double)));

  m_colorTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
}

LayerTrack::~LayerTrack()
{
  if (m_trackData)
  {
    delete m_trackData;
  }
  foreach (vtkActor* actor, m_actors)
    actor->Delete();
  m_actors.clear();
}

QStringList LayerTrack::GetScalarNames()
{
  return m_trackData?m_trackData->m_scalarNames:QStringList();
}

bool LayerTrack::LoadTrackFromFiles()
{
  if (this->m_sFilename.isEmpty())
  {
    return false;
  }
  FSVolume* refVol = NULL;
  if (m_layerMRIRef)
  {
    refVol = m_layerMRIRef->GetSourceVolume();
  }
  m_trackData = new FSTrack(refVol);
  connect(m_trackData, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)));
  if (!m_trackData->LoadFromFiles(m_listFilenames))
  {
    delete m_trackData;
    m_trackData = 0;
    cerr << "Failed to load from file " << qPrintable(m_sFilename) << endl;
    return false;
  }
  if (IsCluster())
    SetName(QFileInfo(m_sFilename).dir().dirName());
  else
    SetName(QFileInfo(m_sFilename).completeBaseName());

  GetProperty()->InitializeScalarThreshold(m_trackData->m_rangeScalar);

  if (m_trackData->HasEmbeddedColor())
  {
    GetProperty()->blockSignals(true);
    GetProperty()->SetColorCode(LayerPropertyTrack::EmbeddedColor);
    GetProperty()->blockSignals(false);
  }

  RebuildActors();

  double dval[6];
  m_trackData->GetRASBounds(dval);
  m_dWorldOrigin[0] = dval[0];
  m_dWorldOrigin[1] = dval[2];
  m_dWorldOrigin[2] = dval[4];
  m_dWorldSize[0] = dval[1]-dval[0];
  m_dWorldSize[1] = dval[3]-dval[2];
  m_dWorldSize[2] = dval[5]-dval[4];
  m_dWorldVoxelSize[0] = m_trackData->m_dVoxelSize[0];
  m_dWorldVoxelSize[1] = m_trackData->m_dVoxelSize[1];
  m_dWorldVoxelSize[2] = m_trackData->m_dVoxelSize[2];

  return true;
}

void LayerTrack::Append2DProps(vtkRenderer *renderer, int nPlane)
{
  Q_UNUSED(renderer);
  Q_UNUSED(nPlane);
}

void LayerTrack::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  Q_UNUSED(bPlaneVisibility);
  foreach (vtkActor* actor, m_actors)
    renderer->AddViewProp(actor);
}

bool LayerTrack::HasProp(vtkProp *prop)
{
  foreach(vtkActor* actor, m_actors)
  {
    if (actor == prop)
      return true;
  }

  return false;
}

void LayerTrack::OnSlicePositionChanged(int nPlane)
{
  Q_UNUSED(nPlane);
}

void LayerTrack::RebuildActors()
{
  foreach (vtkActor* actor, m_actors)
    actor->Delete();
  m_actors.clear();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkFloatArray> scalar_array[100];
  scalars->SetNumberOfComponents(4);
  QList<vtkDataArray*> scalarList;
  QStringList scalarNames = GetScalarNames();
  if (GetProperty()->GetColorCode() == LayerPropertyTrack::Scalar)
  {
    int m = 0;
    foreach (QString name, scalarNames)
    {
      vtkSmartPointer<vtkFloatArray> fa = vtkSmartPointer<vtkFloatArray>::New();
      fa->SetName(qPrintable(name));
      scalarList << fa;
      scalar_array[m++] = fa;
    }
  }
  else
  {
    scalarList << scalars;
  }
  int nLimit = 1000000;
  int nCount = 0;
  float vals[4] = { 0,0,0,255 };
  LayerPropertyTrack* p = GetProperty();
  for (int i = 0; i < m_trackData->m_tracks.size(); i++)
  {
    Track& t = m_trackData->m_tracks[i];
    lines->InsertNextCell(t.nNum);    
    if (p->GetColorCode() == LayerPropertyTrack::Scalar)
    {
      for (int n = 0; n < t.nNum; n++)
      {
        points->InsertNextPoint(t.fPts + n*3);
        lines->InsertCellPoint(nCount);
        for (int j = 0; j < scalarList.size(); j++)
        {
            scalar_array[j]->InsertNextValue(t.fScalars[j][n]);
        }
        nCount++;
      }
    }
    else
    {
      if (p->GetColorCode() == LayerPropertyTrack::EmbeddedColor)
      {
        vals[0] = t.charColor[0];
        vals[1] = t.charColor[1];
        vals[2] = t.charColor[2];
      }
      else if (p->GetDirectionScheme() == LayerPropertyTrack::EndPoints)
        VectorToColor(t.fPts, t.fPts + (t.nNum-1)*3, vals, p->GetDirectionMapping());
      else if (p->GetDirectionScheme() == LayerPropertyTrack::MidSegment)
        VectorToColor(t.fPts+t.nNum/2*3, t.fPts+(t.nNum/2-1)*3, vals, p->GetDirectionMapping());
      for (int n = 0; n < t.nNum; n++)
      {
        points->InsertNextPoint(t.fPts + n*3);
        lines->InsertCellPoint(nCount);
        if (p->GetDirectionScheme() == LayerPropertyTrack::EverySegment)
        {
          if (n == 0)
            VectorToColor(t.fPts+n*3, t.fPts+n*3+3, vals, p->GetDirectionMapping());
          else
            VectorToColor(t.fPts+n*3, t.fPts+n*3-3, vals, p->GetDirectionMapping());
        }

        scalars->InsertNextTuple(vals);
        nCount++;
      }
    }
    if (nCount > nLimit)
    {
      vtkActor* actor = ConstructActor(points, lines, scalarList);
      m_actors << actor;
      points = vtkSmartPointer<vtkPoints>::New();
      lines = vtkSmartPointer<vtkCellArray>::New();
      scalarList.clear();
      if (GetProperty()->GetColorCode() == LayerPropertyTrack::Scalar)
      {
        int m = 0;
        foreach (QString name, scalarNames)
        {
          vtkSmartPointer<vtkFloatArray> fa = vtkSmartPointer<vtkFloatArray>::New();
          fa->SetName(qPrintable(name));
          scalarList << fa;
          scalar_array[m++] = fa;
        }
      }
      else
      {
        scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
        scalars->SetNumberOfComponents(4);
        scalarList << scalars;
      }
      nCount = 0;
    }
  }
  if (nCount > 0)
  {
    vtkActor* actor = ConstructActor(points, lines, scalarList);
    m_actors << actor;
  }

  UpdateColor(false);

  emit ActorChanged();
}

void LayerTrack::UpdateColor(bool emitSignal)
{
  int nCode = GetProperty()->GetColorCode();
  foreach (vtkActor* actor, m_actors)
    actor->GetMapper()->SetScalarVisibility(nCode != LayerPropertyTrack::SolidColor);
  if (nCode == LayerPropertyTrack::SolidColor)
  {
    QColor qc = GetProperty()->GetSolidColor();
    double c[3] = {qc.redF(), qc.greenF(), qc.blueF()};
    foreach (vtkActor* actor, m_actors)
      actor->GetProperty()->SetColor(c);
  }
  else if (nCode == LayerPropertyTrack::Scalar)
  {
    int nMap = GetProperty()->GetScalarColorMap();
    double th[2];
    GetProperty()->GetScalarThreshold(th);
    if (th[1] <= th[0])
      th[1] = th[0]+1;
    m_colorTable->RemoveAllPoints();
    switch (nMap)
    {
    case LayerPropertyTrack::Heatscale:
      m_colorTable->AddRGBAPoint( th[0], 1, 0, 0, 1 );
      m_colorTable->AddRGBAPoint( th[1], 1, 1, 0, 1 );
      m_colorTable->Build();
      break;
    case LayerPropertyTrack::Jet:
      m_colorTable->AddRGBAPoint( qMin( 0.0, th[0]), 0, 0, 0, 0 );
      m_colorTable->AddRGBAPoint( th[0], 0, 0, 1, 1 );
      m_colorTable->AddRGBAPoint( th[0] + (th[1] - th[0]) / 4, 0, 1, 1, 1 );
      m_colorTable->AddRGBAPoint( (th[0] + th[1]) / 2, 0, 1, 0, 1 );
      m_colorTable->AddRGBAPoint( th[1] - (th[1] - th[0]) / 4, 1, 1, 0, 1 );
      m_colorTable->AddRGBAPoint( th[1], 1, 0, 0, 1 );
      m_colorTable->Build();
      break;
    }

    QStringList names = GetScalarNames();
    int n = GetProperty()->GetScalarIndex();
    foreach (vtkActor* actor, m_actors)
    {
      actor->GetMapper()->SetLookupTable(m_colorTable);
      vtkPolyData* poly = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInput());
      if (poly)
        poly->GetPointData()->SetActiveScalars(qPrintable(names[n]));
    }
  }

  if (emitSignal)
    emit ActorUpdated();
}

void LayerTrack::VectorToColor(float *pt1, float *pt2, float *c_out, int nMappingType)
{
  double v[3] = { 1, 0 ,0 }, t;
  MyUtils::GetVector(pt1, pt2, v);
  switch (nMappingType)
  {
  case LayerPropertyTrack::RBG:
    t = v[1];
    v[1] = v[2];
    v[2] = t;
    break;
  case LayerPropertyTrack::GBR:
    t = v[0];
    v[0] = v[1];
    v[1] = v[2];
    v[2] = t;
    break;
  case LayerPropertyTrack::GRB:
    t = v[0];
    v[0] = v[1];
    v[1] = t;
    break;
  case LayerPropertyTrack::BGR:
    t = v[0];
    v[0] = v[2];
    v[2] = t;
    break;
  case LayerPropertyTrack::BRG:
    t = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = t;
    break;
  }

  c_out[0] = fabs(v[0])*255;
  c_out[1] = fabs(v[1])*255;
  c_out[2] = fabs(v[2])*255;
}

vtkActor* LayerTrack::ConstructActor(vtkPoints *points, vtkCellArray *lines, QList<vtkDataArray*> scalars)
{
  vtkActor* actor = vtkActor::New();
#if VTK_MAJOR_VERSION > 5
  actor->ForceOpaqueOn();
#endif
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetLines(lines);
  polydata->GetPointData()->SetScalars(scalars[0]);
  for (int i = 1; i < scalars.size(); i++)
    polydata->GetPointData()->AddArray(scalars[i]);
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  if (GetProperty()->GetRenderRep() == LayerPropertyTrack::Tube)
  {
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->CappingOn();
#if VTK_MAJOR_VERSION > 5
    tube->SetInputData(polydata);
#else
    tube->SetInput(polydata);
#endif
    tube->SetRadius(GetProperty()->GetTubeRadius());
    tube->SetNumberOfSides(GetProperty()->GetNumberOfSides());
    mapper->SetInputConnection(tube->GetOutputPort());
    actor->SetMapper(mapper);
    /*
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(polydata);
    actor->AddLODMapper(mapper);
    */
  }
  else
  {
#if VTK_MAJOR_VERSION > 5
    mapper->SetInputData(polydata);
#else
    mapper->SetInput(polydata);
#endif
    actor->SetMapper(mapper);
  }
  return actor;
}

void LayerTrack::SetVisible(bool bVisible)
{
  foreach(vtkActor* actor, m_actors)
    actor->SetVisibility(bVisible?1:0);
  Layer::SetVisible(bVisible);
}

bool LayerTrack::IsVisible()
{
  if (m_actors.isEmpty())
    return false;
  else
    return m_actors[0]->GetVisibility();
}

void LayerTrack::SetFileName(const QString &filename)
{
  Layer::SetFileName(filename);
  m_listFilenames.clear();
  m_listFilenames << filename;
}

void LayerTrack::SetFileNames(const QStringList &filenames)
{
  Layer::SetFileName(filenames.first());
  m_listFilenames = filenames;
}

void LayerTrack::SetClusterData(const QVariantMap &data)
{
  m_mapCluster = data;
  if (data.contains("filenames"))
    SetFileNames(data["filenames"].toStringList());
}

bool LayerTrack::IsCluster()
{
  return !m_mapCluster.isEmpty();
}

bool LayerTrack::HasEmbeddedColor()
{
  return (m_trackData && m_trackData->HasEmbeddedColor());
}

void LayerTrack::UpdateOpacity(double val)
{
  foreach(vtkActor* actor, m_actors)
    actor->GetProperty()->SetOpacity(val);
  emit ActorUpdated();
}

void LayerTrack::GetScalarRange(double *range, int nIndex)
{
  if (nIndex < 0)
    nIndex = GetProperty()->GetScalarIndex();

  QList< QPair<double, double> > ranges = m_trackData->m_rangeScalar;
  if (nIndex < ranges.size())
  {
    range[0] = ranges[nIndex].first;
    range[1] = ranges[nIndex].second;
  }
}

vtkRGBAColorTransferFunction* LayerTrack::GetColorTable() const
{
  return m_colorTable;
}
