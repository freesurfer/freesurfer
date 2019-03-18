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
#include <vtkBox.h>
#include <vtkCutter.h>
#include "MyVTKUtils.h"

LayerConnectomeMatrix::LayerConnectomeMatrix(LayerMRI* ref, QObject *parent) :
  Layer(parent),
  m_mriRef(ref),
  m_mriParcel(NULL),
  m_cmat(NULL),
  m_ctab(NULL),
  m_bToAllLabels(false),
  m_bVisible(true)
{
  this->m_strTypeNames << "CMAT";
  mProperty = new LayerPropertyConnectomeMatrix( this );

  LayerPropertyConnectomeMatrix* p = GetProperty();
  connect(p, SIGNAL(OpacityChanged()), this, SLOT(UpdateOpacity()));
  connect(p, SIGNAL(SplineRadiusChanged()), this, SLOT(RebuildSplineActors()));
  connect(p, SIGNAL(SplineColorChanged()), this, SLOT(UpdateSplineColor()));

  m_actorSplines = vtkSmartPointer<vtkActor>::New();
  QColor c = p->GetSplineColor();
  m_actorSplines->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  for (int i = 0; i < 3; i++)
  {
    m_actorSlice[i] = vtkSmartPointer<vtkActor>::New();
    m_actorSlice[i]->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  }
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
  bool bVoxelCoords = (m_cmat->coords == LABEL_COORDS_VOXEL);
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
          if (bVoxelCoords)
          {
            m_mriParcel->OriginalVoxelToRAS(pt, pt);
          }
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
                                  actor, 0, 0, true);
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
  if (m_actorSlice[nPlane].GetPointer())
    renderer->AddViewProp(m_actorSlice[nPlane]);
}

void LayerConnectomeMatrix::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  Q_UNUSED(bPlaneVisibility);
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

  if (!m_bVisible)
    return;

  for (int i = 0; i < m_listFromLabelIndices.size(); i++)
    m_actorLabels[m_listFromLabelIndices[i]]->VisibilityOn();

  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    for (int j = 0; j < m_listFromLabelIndices.size(); j++)
    {
      int nFrom = m_listFromLabelIndices[j];
      if (i != nFrom && (m_bToAllLabels || m_listToLabelIndices.contains(i)))
      {
        if (m_cmat->splines[nFrom][i] || m_cmat->splines[i][nFrom])
          m_actorLabels[i]->VisibilityOn();
      }
    }
  }
}

bool LayerConnectomeMatrix::HasProp(vtkProp *prop)
{
  Q_UNUSED(prop);
  return false;
}

bool LayerConnectomeMatrix::IsVisible()
{
  return m_bVisible;
}

void LayerConnectomeMatrix::SetVisible(bool bVisible)
{
  m_bVisible = bVisible;
  m_actorSplines->SetVisibility(bVisible?1:0);
  for (int i = 0; i < 3; i++)
    m_actorSlice[i]->SetVisibility(bVisible?1:0);

  UpdateLabelActors();

  Layer::SetVisible(bVisible);
}

void LayerConnectomeMatrix::OnSlicePositionChanged(int nPlane)
{
  Q_UNUSED(nPlane);
  RebuildSplineActors();
}

void LayerConnectomeMatrix::SetFromLabelIndices(const QList<int> &indices)
{
  m_listFromLabelIndices = indices;
  UpdateLabelActors();
  UpdateOpacity();
  RebuildSplineActors();
}

void LayerConnectomeMatrix::SetToLabelIndices(const QList<int> &indices)
{
  m_listToLabelIndices = indices;
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

void LayerConnectomeMatrix::UpdateSplineColor()
{
  QColor c = GetProperty()->GetSplineColor();
  m_actorSplines->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  emit ActorUpdated();
}

void LayerConnectomeMatrix::RebuildSplineActors()
{
  if (m_listFromLabelIndices.isEmpty())
    return;

  if (!m_bToAllLabels && m_listToLabelIndices.isEmpty())
  {
    m_actorSplines->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    emit ActorUpdated();
    return;
  }

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  int nCount = 0;
  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    for (int j = 0; j < m_listFromLabelIndices.size(); j++)
    {
      int nFrom = m_listFromLabelIndices[j];
      if (m_bToAllLabels || m_listToLabelIndices.contains(i))
      {
        LABEL* label = m_cmat->splines[nFrom][i];
        if (!label)
          label = m_cmat->splines[i][nFrom];
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
  }
  double voxel_len = qMin(m_dWorldVoxelSize[0], qMin(m_dWorldVoxelSize[1], m_dWorldVoxelSize[2]));
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetLines(lines);
  vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION > 5
  spline->SetInputData(polydata);
#else
  spline->SetInput(polydata);
#endif
  spline->SetSubdivideToLength();
  spline->SetLength(2*voxel_len);
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
  tube->SetInputConnection(spline->GetOutputPort());
  tube->SetRadius(GetProperty()->GetSplineRadius()*voxel_len);
  tube->SetNumberOfSides(GetProperty()->GetNumberOfSides());
  mapper->SetInputConnection(tube->GetOutputPort());
  mapper->SetScalarVisibility(0);
  m_actorSplines->SetMapper(mapper);

  // 2D actors
  double wsize[3], worigin[3], vsize[3], pos[3];
  GetWorldOrigin(worigin);
  GetWorldSize(wsize);
  GetWorldVoxelSize(vsize);
  GetSlicePosition(pos);
  for (int i = 0; i < 3; i++)
  {
    vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();

    double bound[6] = { worigin[0], worigin[0]+wsize[0], worigin[1], worigin[1]+wsize[1],
                        worigin[2], worigin[2]+wsize[2] };
    bound[i*2] = pos[i] - vsize[i]/2;
    bound[i*2+1] = pos[i] + vsize[i]/2;
    box->SetBounds(bound);

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputConnection( spline->GetOutputPort() );
    cutter->SetCutFunction( box );

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInputConnection(cutter->GetOutputPort());
    tube->SetRadius(GetProperty()->GetSplineRadius()*voxel_len);
    tube->SetNumberOfSides(GetProperty()->GetNumberOfSides());
    mapper->SetInputConnection(tube->GetOutputPort());
    mapper->SetScalarVisibility(0);
    m_actorSlice[i]->SetMapper(mapper);
  }

  emit ActorUpdated();
}

bool LayerConnectomeMatrix::HasConnection(int i, int j)
{
  if (!m_cmat || i < 0 || i >= m_cmat->nlabels || j < 0 || j >= m_cmat->nlabels)
    return false;

  return (m_cmat->splines[i][j] || m_cmat->splines[j][i]);
}

void LayerConnectomeMatrix::UpdateOpacity()
{
  for (int i = 0; i < m_actorLabels.size(); i++)
  {
    if (m_listFromLabelIndices.contains(i))
      m_actorLabels[i]->GetProperty()->SetOpacity(GetProperty()->GetFromLabelOpacity());
    else
      m_actorLabels[i]->GetProperty()->SetOpacity(GetProperty()->GetToLabelOpacity());
  }
  emit ActorUpdated();
}
