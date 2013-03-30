#include "LayerConnectomeMatrix.h"
#include "LayerMRI.h"
#include "LayerPropertyConnectomeMatrix.h"
#include <QDebug>
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkSplineFilter.h"
#include "vtkTubeFilter.h"
#include "vtkPoints.h"
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include "MyVTKUtils.h"

LayerConnectomeMatrix::LayerConnectomeMatrix(LayerMRI* ref, QObject *parent) :
  Layer(parent),
  m_mriRef(ref),
  m_mriParcel(NULL),
  m_cmat(NULL),
  m_ctab(NULL),
  m_nFromLabelIndex(-1),
  m_nToLabelIndex(-1),
  m_bToAllLabels(false)
{
  this->m_strTypeNames << "CMAT";
  mProperty = new LayerPropertyConnectomeMatrix( this );

  m_actorSplines = vtkSmartPointer<vtkActor>::New();
}

LayerConnectomeMatrix::~LayerConnectomeMatrix()
{
  if (m_cmat)
    ::CMATfree(&m_cmat);

  if (m_mriParcel)
    delete m_mriParcel;
}

void LayerConnectomeMatrix::SetColorTable(COLOR_TABLE* ctab)
{
  m_ctab = ctab;
}

bool LayerConnectomeMatrix::LoadFromFile(const QString &fn_cmat, const QString &fn_parcel)
{
  if (!m_ctab)
  {
    cerr << "Did not set color table. Abort loading.";
    return false;
  }

  // read cmat file
  if (m_cmat)
    ::CMATfree(&m_cmat);
  m_cmat = ::CMATread(qPrintable(fn_cmat));
  if (!m_cmat)
  {
    cerr << "Could not load CMAT file " << qPrintable(fn_cmat) << endl;
    return false;
  }
  m_listLabels.clear();
  for (int i = 0; i < m_cmat->nlabels; i++)
    m_listLabels << m_cmat->labels[i];

  // read parcellation file
  if (m_mriParcel)
    delete m_mriParcel;
  m_mriParcel = new LayerMRI(m_mriRef);
  m_mriParcel->SetFileName(fn_parcel);
  if (!m_mriParcel->LoadVolumeFromFile())
  {
    cerr << "Could not load parcellation file " << qPrintable(fn_parcel) << endl;
    delete m_mriParcel;
    m_mriParcel = NULL;
    ::CMATfree(&m_cmat);
    return false;
  }

  SetWorldOrigin(m_mriParcel->GetWorldOrigin());
  SetWorldSize(m_mriParcel->GetWorldSize());
  SetWorldVoxelSize(m_mriParcel->GetWorldVoxelSize());

  // update label coords to target coords
  double pt[3];
  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    for (int j = 0; j < m_cmat->nlabels; j++)
    {
      LABEL* label = m_cmat->splines[i][j];
      if (label)
      {
        for (int n = 0; n < label->n_points; n++)
        {
          pt[0] = label->lv[n].x;
          pt[1] = label->lv[n].y;
          pt[2] = label->lv[n].z;
          m_mriParcel->RASToTarget(pt, pt);
          label->lv[n].x = pt[0];
          label->lv[n].y = pt[1];
          label->lv[n].z = pt[2];
        }
      }
    }
  }

  BuildLabelActors();

  return true;
}

void LayerConnectomeMatrix::BuildLabelActors()
{
  for (int i = 0; i < m_listLabels.size(); i++)
  {
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    actor->SetMapper(mapper);
    MyVTKUtils::BuildContourActor(m_mriParcel->GetImageData(),
                                  m_listLabels[i], m_listLabels[i],
                                  actor);
    QColor c = GetLabelColor(m_listLabels[i]);
    actor->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
    actor->VisibilityOff();
    actor->GetMapper()->SetScalarVisibility(0);
    m_actorLabels << actor;
  }
}

QList<int> LayerConnectomeMatrix::GetLabelList()
{
  return m_listLabels;
}

QString LayerConnectomeMatrix::GetLabelName(int n)
{
  int valid;
  CTABisEntryValid(m_ctab, n, &valid);
  if (valid)
  {
    char name[1000];
    CTABcopyName( m_ctab, n, name, 1000 );
    return name;
  }
  else
    return "";
}

QColor LayerConnectomeMatrix::GetLabelColor(int n)
{
  int valid;
  CTABisEntryValid(m_ctab, n, &valid);
  if (valid)
  {
    int nr, ng, nb;
    CTABrgbAtIndexi( m_ctab, n, &nr, &ng, &nb );
    return QColor( nr, ng, nb );
  }
  else
    return QColor(0, 0, 0);
}

bool LayerConnectomeMatrix::LoadFromFile()
{
  return LoadFromFile(m_sFilename, m_sFilenameParcel);
}

void LayerConnectomeMatrix::Append2DProps(vtkRenderer *renderer, int nPlane)
{

}

void LayerConnectomeMatrix::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  renderer->AddViewProp(m_actorSplines);
  for (int i = 0; i < m_listLabels.size(); i++)
    renderer->AddViewProp(m_actorLabels[i]);
}

void LayerConnectomeMatrix::UpdateLabelActors()
{
  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    m_actorLabels[i]->VisibilityOff();
  }

  if (m_nFromLabelIndex >= 0)
    m_actorLabels[m_nFromLabelIndex]->VisibilityOn();

  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    if (i != m_nFromLabelIndex && (m_bToAllLabels || i == m_nToLabelIndex))
    {
      if (m_cmat->splines[m_nFromLabelIndex][i])
        m_actorLabels[i]->VisibilityOn();
    }
  }
}

bool LayerConnectomeMatrix::HasProp(vtkProp *prop)
{

}

bool LayerConnectomeMatrix::IsVisible()
{

}

void LayerConnectomeMatrix::OnSlicePositionChanged(int nPlane)
{

}

void LayerConnectomeMatrix::SetFromLabelIndex(int n)
{
  m_nFromLabelIndex = n;
  UpdateLabelActors();
  RebuildSplineActors();
}

void LayerConnectomeMatrix::SetToLabelIndex(int n)
{
  m_nToLabelIndex = n;
  if (!m_bToAllLabels)
  {
    UpdateLabelActors();
    RebuildSplineActors();
  }
}

void LayerConnectomeMatrix::SetToAllLabels(bool bAll)
{
  m_bToAllLabels = bAll;
  UpdateLabelActors();
  RebuildSplineActors();
}

void LayerConnectomeMatrix::RebuildSplineActors()
{
  if (m_nFromLabelIndex < 0)
    return;

  if (!m_bToAllLabels && m_nToLabelIndex < 0)
    return;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  int nCount = 0;
  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    if (m_bToAllLabels || i == m_nToLabelIndex)
    {
      LABEL* label = m_cmat->splines[m_nFromLabelIndex][i];
      if (label)
      {
        lines->InsertNextCell(label->n_points);
        for (int n = 0; n < label->n_points; n++)
        {
          points->InsertNextPoint(label->lv[n].x, label->lv[n].y, label->lv[n].z);
          lines->InsertCellPoint(nCount);
          nCount++;
        }
      }
    }
  }
  double voxel_len = qMin(m_dWorldVoxelSize[0], qMin(m_dWorldVoxelSize[1], m_dWorldVoxelSize[2]));
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetLines(lines);
  vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
  spline->SetInput(polydata);
  spline->SetSubdivideToLength();
  spline->SetLength(2*voxel_len);
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
  tube->SetInput(spline->GetOutput());
  tube->SetRadius(GetProperty()->GetSplineRadius()*voxel_len);
  tube->SetNumberOfSides(GetProperty()->GetNumberOfSides());
  mapper->SetInput(tube->GetOutput());
  m_actorSplines->SetMapper(mapper);
  emit ActorUpdated();
}

bool LayerConnectomeMatrix::HasConnection(int i, int j)
{
  if (!m_cmat || i < 0 || i >= m_cmat->nlabels || j < 0 || j >= m_cmat->nlabels)
    return false;

  return m_cmat->splines[i][j];
}
