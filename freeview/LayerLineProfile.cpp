#include "LayerLineProfile.h"
#include "LayerPointSet.h"
#include "LayerMRI.h"
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkTubeFilter.h>
#include <vtkSplineFilter.h>
#include <vtkCutter.h>
#include <vtkStripper.h>
#include <vtkTriangleFilter.h>
#include <vtkPlane.h>
#include <vtkProperty.h>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include "LayerPropertyLineProfile.h"
#include "FSPointSet.h"
#include "LayerPropertyPointSet.h"
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include "MyUtils.h"
#include "MigrationDefs.h"

LayerLineProfile::LayerLineProfile(int nPlane, QObject *parent, LayerPointSet *line1, LayerPointSet *line2) :
  Layer(parent),
  m_nPlane(nPlane),
  m_dResolution(2),
  m_dSpacing(2),
  m_nSamples(100),
  m_dOffset(10),
  m_nActiveLineId(-1),
  m_dReferenceSize(1.0)
{
  this->m_strTypeNames << "Supplement" << "LineProfile";
  SetSourceLayers(line1, line2);

  mProperty = new LayerPropertyLineProfile( this );

  m_endLines = vtkSmartPointer<vtkActor>::New();
  m_profileLines = vtkSmartPointer<vtkActor>::New();
  m_endLines->GetProperty()->SetColor(0.3, 0.3, 1.0);
  QColor c = GetProperty()->GetColor();
  m_profileLines->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  m_profileLines->GetProperty()->SetOpacity(GetProperty()->GetOpacity());

  m_activeLine = vtkSmartPointer<vtkActor>::New();
  m_activeLine->GetProperty()->SetColor(0.3, 0.3, 1.0);

  LayerPropertyLineProfile* p = GetProperty();
  connect(p, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()));
  connect(p, SIGNAL(ColorChanged(QColor)), this, SLOT(UpdateColor()));
  connect(p, SIGNAL(RadiusChanged(double)), this, SLOT(UpdateActors()));
}

LayerLineProfile::~LayerLineProfile()
{

}

void LayerLineProfile::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  if (nPlane == m_nPlane)
  {
    renderer->AddViewProp(m_profileLines);
    renderer->AddViewProp(m_endLines);
    renderer->AddViewProp(m_activeLine);
  }
}

bool LayerLineProfile::HasProp( vtkProp* prop )
{
  return (prop == m_profileLines.GetPointer() || prop == m_endLines.GetPointer());
}

bool LayerLineProfile::IsVisible()
{
  return (this->m_profileLines->GetVisibility() > 0);
}

void LayerLineProfile::SetVisible( bool bVisible )
{
  m_profileLines->SetVisibility(bVisible?1:0);
  m_endLines->SetVisibility(bVisible?1:0);
  m_activeLine->SetVisibility(bVisible?1:0);
}

void LayerLineProfile::SetSourceLayers(LayerPointSet *line1, LayerPointSet *line2)
{
  m_spline0 = line1;
  m_spline1 = line2;
  if (m_spline0)
    connect(m_spline0, SIGNAL(destroyed()), this, SLOT(OnSourceLineDestroyed()));
  if (m_spline1)
    connect(m_spline1, SIGNAL(destroyed()), this, SLOT(OnSourceLineDestroyed()));
}

void LayerLineProfile::OnSlicePositionChanged(int nPlane)
{
  Q_UNUSED(nPlane);
}

void LayerLineProfile::UpdateColor()
{
  QColor c = GetProperty()->GetColor();
  m_profileLines->GetProperty()->SetColor(c.redF(), c.greenF(), c.blueF());
  emit this->ActorUpdated();
}

void LayerLineProfile::UpdateOpacity()
{
  m_profileLines->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  m_endLines->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  // m_activeLine->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  emit ActorUpdated();
}

std::vector < std::vector < double > > LayerLineProfile::Points3DToSpline2D(std::vector<double> pts3d, double distance)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell(pts3d.size()/3);
  for (size_t i = 0; i < pts3d.size(); i+= 3)
  {
    points->InsertNextPoint(pts3d[i], pts3d[i+1], pts3d[i+2]);
    lines->InsertNextCell(i/3);
  }
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
  spline->SetLength(distance);
  spline->Update();
  polydata = spline->GetOutput();
  points = polydata->GetPoints();
  std::vector < std::vector < double > > points2d;
  for (int i = 0; i < points->GetNumberOfPoints(); i++)
  {
    std::vector<double> pt;
    double x[3];
    points->GetPoint(i, x);
    if (m_nPlane != 0)
      pt.push_back(x[0]);
    if (m_nPlane != 1)
      pt.push_back(x[1]);
    if (m_nPlane != 2)
      pt.push_back(x[2]);
    points2d.push_back(pt);
  }
  return points2d;
}

std::vector< std::vector<double> > LayerLineProfile::Points2DToSpline3D(std::vector < std::vector<double> > pts2d, int nSample)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell(pts2d.size());
  for (size_t i = 0; i < pts2d.size(); i++)
  {
    if (m_nPlane == 0)
      points->InsertNextPoint(m_dSliceLocation, pts2d[i][0], pts2d[i][1]);
    else if (m_nPlane == 1)
      points->InsertNextPoint(pts2d[i][0], m_dSliceLocation, pts2d[i][1]);
    else if (m_nPlane == 2)
      points->InsertNextPoint(pts2d[i][0], pts2d[i][1], m_dSliceLocation);
    lines->InsertNextCell(i);
  }
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetLines(lines);
  vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION > 5
  spline->SetInputData(polydata);
#else
  spline->SetInput(polydata);
#endif
  spline->SetNumberOfSubdivisions(nSample);
  spline->Update();
  polydata = spline->GetOutput();
  points = polydata->GetPoints();
  std::vector < std::vector < double > > points3d;
  for (int i = 0; i < points->GetNumberOfPoints(); i++)
  {
    std::vector<double> pt;
    double x[3];
    points->GetPoint(i, x);
    pt.push_back(x[0]);
    pt.push_back(x[1]);
    pt.push_back(x[2]);
    points3d.push_back(pt);
  }
  return points3d;
}

bool LayerLineProfile::Solve(double profileSpacing, double referenceSize, double laplaceResolution, double offset)
{
  //double laplaceResolution = ;
  //  double referenceSize     = voxel_length; // smallest voxel lenth
  double resolution        = laplaceResolution * referenceSize;
  double spacing           = profileSpacing * referenceSize;
  //double resolution        = spacing * 10.0;
  double distance          = referenceSize / 3.0;

  m_dReferenceSize = referenceSize;

  std::vector < std::vector < double > > points2d;
  std::vector < int > segment0;
  std::vector < int > segment1;
  std::vector < int > segmentL;
  std::vector < int > segmentR;

  std::vector<double> ctrl_pts0  = m_spline0->GetPoints();
  std::vector < std::vector < double > > pts0 = Points3DToSpline2D(ctrl_pts0, distance);
  for (size_t i = 0; i < pts0.size(); i++)
    segment0.push_back(i);

  std::vector<double> ctrl_pts1  = m_spline1->GetPoints();
  std::vector < std::vector < double > > pts1 = Points3DToSpline2D(ctrl_pts1, distance);
  for (size_t i = 0; i < pts1.size(); i++)
    segment1.push_back(pts0.size()+i);

  std::vector<double> ctrl_ptsL;
  ctrl_ptsL.push_back(ctrl_pts0[0]);
  ctrl_ptsL.push_back(ctrl_pts0[1]);
  ctrl_ptsL.push_back(ctrl_pts0[2]);
  ctrl_ptsL.push_back(ctrl_pts1[0]);
  ctrl_ptsL.push_back(ctrl_pts1[1]);
  ctrl_ptsL.push_back(ctrl_pts1[2]);
  std::vector < std::vector < double > > ptsL = Points3DToSpline2D(ctrl_ptsL, distance);
  for (size_t i = 0; i < ptsL.size(); i++)
    segmentL.push_back(pts0.size() + pts1.size() + i);

  std::vector<double> ctrl_ptsR;
  ctrl_ptsR.push_back(ctrl_pts0[ctrl_pts0.size()-3]);
  ctrl_ptsR.push_back(ctrl_pts0[ctrl_pts0.size()-2]);
  ctrl_ptsR.push_back(ctrl_pts0[ctrl_pts0.size()-1]);
  ctrl_ptsR.push_back(ctrl_pts1[ctrl_pts1.size()-3]);
  ctrl_ptsR.push_back(ctrl_pts1[ctrl_pts1.size()-2]);
  ctrl_ptsR.push_back(ctrl_pts1[ctrl_pts1.size()-1]);
  std::vector < std::vector < double > > ptsR = Points3DToSpline2D(ctrl_ptsR, distance);
  for (size_t i = 0; i < ptsR.size(); i++)
    segmentR.push_back(pts0.size() + pts1.size() + ptsL.size() + i);

  points2d = pts0;
  for (size_t i = 0; i < pts1.size(); i++)
    points2d.push_back(pts1[i]);
  for (size_t i = 0; i < ptsL.size(); i++)
    points2d.push_back(ptsL[i]);
  for (size_t i = 0; i < ptsR.size(); i++)
    points2d.push_back(ptsR[i]);

  m_dSliceLocation = ctrl_pts0[m_nPlane];

#if !defined(ARM64) && !defined(DISABLE_LINEPROF)
  LineProf LP(points2d, segment0, segment1, segmentL, segmentR);
  // Next we solve the Laplace on the domain
  int paddingL       = 10;
  int paddingR       = 10;
  int convergence    = 8;
  LP.solveLaplace(paddingL,paddingR,resolution,convergence);

  // And finally compute line profiles
  m_ptsProfile = LP.ComputeProfiles(offset*referenceSize, spacing);
#endif

  m_nActiveLineId = -1;
  m_activeLine->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());

  m_dResolution = laplaceResolution;
  m_dSpacing = profileSpacing;
  m_dOffset = offset;

  UpdateActors();
  return true;
}

void LayerLineProfile::MakeFlatTube(vtkPoints* points, vtkCellArray* lines, vtkActor* actor_in, double radius)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetLines(lines);
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION > 5
  tube->SetInputData(polydata);
#else
  tube->SetInput(polydata);
#endif
  tube->SetNumberOfSides(6);
  tube->SetRadius(radius);

  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  double origin[3] = {0, 0, 0};
  origin[m_nPlane] = m_dSliceLocation;
  plane->SetOrigin( origin );
  plane->SetNormal( m_nPlane==0?1:0, m_nPlane==1?1:0, m_nPlane==2?1:0 );

  vtkSmartPointer<vtkCutter> cutter =
      vtkSmartPointer<vtkCutter>::New();
  cutter->SetInputConnection( tube->GetOutputPort() );
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
  mapper->SetInputConnection(tube->GetOutputPort());
  /*
  mapper->SetScalarVisibility(true);
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfColors(2);
  QColor c = GetProperty()->GetColor();
  lut->SetTableValue(0, c.redF(), c.greenF(), c.blueF());
  lut->SetTableValue(1, 0, 0, 1.0);
  mapper->SetLookupTable(lut);
  */
  actor_in->SetMapper(mapper);
}

void LayerLineProfile::UpdateActiveLine()
{
  if (m_nActiveLineId >= 0)
  {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    int nCount = 0;
    for (int i = m_nActiveLineId; i < m_nActiveLineId+1; i++)
    {
      std::vector < std::vector < double > > line = m_ptsProfile[i];
      lines->InsertNextCell(line.size());
      for (size_t j = 0; j < line.size(); j++)
      {
        lines->InsertCellPoint(nCount++);
        double pt[3] = {m_dSliceLocation, m_dSliceLocation, m_dSliceLocation};
        if (m_nPlane == 0)
        {
          pt[1] = line[j][0];
          pt[2] = line[j][1];
        }
        else if (m_nPlane == 1)
        {
          pt[0] = line[j][0];
          pt[2] = line[j][1];
        }
        else
        {
          pt[0] = line[j][0];
          pt[1] = line[j][1];
        }
        points->InsertNextPoint(pt);
      }
    }

    MakeFlatTube(points, lines, m_activeLine, GetProperty()->GetRadius()*1.01);
    emit ActorUpdated();
  }
}

void LayerLineProfile::UpdateActors()
{
  if (!m_spline0 || !m_spline1)
    return;

  std::vector<double> pts1 = m_spline0->GetPoints(), pts2 = m_spline1->GetPoints();
  if (pts1.size() < 4 || pts2.size() < 4)
    return;

  // update boundary lines
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  points->InsertNextPoint(pts1[0], pts1[1], pts1[2]);
  points->InsertNextPoint(pts2[0], pts2[1], pts2[2]);
  points->InsertNextPoint(pts1[pts1.size()-3], pts1[pts1.size()-2], pts1[pts1.size()-1]);
  points->InsertNextPoint(pts2[pts2.size()-3], pts2[pts2.size()-2], pts2[pts2.size()-1]);
  lines->InsertNextCell(2);
  lines->InsertCellPoint(0);
  lines->InsertCellPoint(1);
  lines->InsertNextCell(2);
  lines->InsertCellPoint(2);
  lines->InsertCellPoint(3);
  MakeFlatTube(points, lines, m_endLines, GetProperty()->GetRadius());

  // update profile lines
  points = vtkSmartPointer<vtkPoints>::New();
  lines = vtkSmartPointer<vtkCellArray>::New();
  int nCount = 0;
  for (size_t i = 0; i < m_ptsProfile.size(); i++)
  {
    std::vector < std::vector < double > > line = m_ptsProfile[i];
    lines->InsertNextCell(line.size());
    for (size_t j = 0; j < line.size(); j++)
    {
      lines->InsertCellPoint(nCount++);
      double pt[3] = {m_dSliceLocation, m_dSliceLocation, m_dSliceLocation};
      if (m_nPlane == 0)
      {
        pt[1] = line[j][0];
        pt[2] = line[j][1];
      }
      else if (m_nPlane == 1)
      {
        pt[0] = line[j][0];
        pt[2] = line[j][1];
      }
      else
      {
        pt[0] = line[j][0];
        pt[1] = line[j][1];
      }
      points->InsertNextPoint(pt);
    }
  }

  MakeFlatTube(points, lines, m_profileLines, GetProperty()->GetRadius());

  UpdateActiveLine();
  emit this->ActorChanged();
}

void LayerLineProfile::OnSourceLineDestroyed()
{
  if (sender() == m_spline0)
    m_spline0 = NULL;
  if (sender() == m_spline1)
    m_spline1 = NULL;

  if (!m_spline0 || !m_spline1)
    this->Hide();
}

bool LayerLineProfile::Export(const QString &filename, LayerMRI *mri, int nSamples)
{
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return false;

  QTextStream out(&file);
  for (size_t i = 0; i < m_ptsProfile.size(); i++)
  {
    std::vector < std::vector <double> > line3d = Points2DToSpline3D(m_ptsProfile[i], nSamples);
    std::vector<double> vals = mri->GetMeanSegmentValues(line3d);
    double len = 0;
    for (size_t j = 0; j < line3d.size()-1; j++)
    {
      double pt0[3] = { line3d[j][0], line3d[j][1], line3d[j][2] };
      double pt1[3] = { line3d[j+1][0], line3d[j+1][1], line3d[j+1][2] };
      len += sqrt(vtkMath::Distance2BetweenPoints(pt0, pt1));
    }
    out << i << "," << len <<  "," << nSamples;
    for (size_t j = 0; j < vals.size(); j++)
    {
      out << "," << vals[j];
    }
    out << ",,";
    double pt[3];
    for (size_t n = 0; n < line3d.size(); n++)
    {
      pt[0] = line3d[n][0]; pt[1] = line3d[n][1]; pt[2] = line3d[n][2];
      mri->TargetToRAS(pt, pt);
      out << pt[0] << "," << pt[1] << "," << pt[2] << ",";
    }
    out << "\n";
  }
  m_nSamples = nSamples;

  return true;
}

QList<double> LayerLineProfile::ComputeThicknessAlongLineProfile(std::vector<std::vector<double> > &line2d_in,
                                                                 const QList<LayerPointSet *> &splines, int nSample)
{
  QList<double> thicknesses;
//  double x0 = line2d[0][0], y0 = line2d[0][1];
  double x, y, len = 0;
  int n0 = 0, n;
  std::vector< std::vector<double> > line3d = Points2DToSpline3D(line2d_in, nSample);
  std::vector<std::vector<double> > line2d;
  for (size_t i = 0; i < line3d.size(); i++)
  {
    std::vector<double> point2d;
    point2d.push_back(line3d[i][0]);
    point2d.push_back(line3d[i][1]);
    switch (m_nPlane)
    {
    case 0:
      point2d[0] = line3d[i][1];
      point2d[1] = line3d[i][2];
      break;
    case 1:
      point2d[1] = line3d[i][2];
      break;
    default:
      break;
    }
    line2d.push_back(point2d);
  }
  for (int i = 0; i < splines.size(); i++)
  {
    std::vector<double> ctrl_pts0  = splines[i]->GetPoints();
    std::vector < std::vector < double > > spline2d = Points3DToSpline2D(ctrl_pts0, m_dWorldVoxelSize[m_nPlane]/3.0);
    if (MyUtils::FindIntersection(line2d, spline2d, &x, &y, &n))
    {
      for (int j = n0; j < n; j++)
      {
        len += sqrt((line2d[j][0]-line2d[j+1][0])*(line2d[j][0]-line2d[j+1][0]) +
            (line2d[j][1]-line2d[j+1][1])*(line2d[j][1]-line2d[j+1][1]));
      }
      len += sqrt((line2d[n][0]-x)*(line2d[n][0]-x)+(line2d[n][1]-y)*(line2d[n][1]-y));
      thicknesses << len;
      if (n+1 < line2d.size())
      {
        len = sqrt((line2d[n+1][0]-x)*(line2d[n+1][0]-x)+(line2d[n+1][1]-y)*(line2d[n+1][1]-y));
        n0 = n+1;
      }
    }
    else
    {
      thicknesses << -1;
      cerr << "Failed to find intersection point";
    }
  }
  n = line2d.size()-1;
  for (int j = n0; j < n; j++)
  {
    len += sqrt((line2d[j][0]-line2d[j+1][0])*(line2d[j][0]-line2d[j+1][0]) +
        (line2d[j][1]-line2d[j+1][1])*(line2d[j][1]-line2d[j+1][1]));
  }
  thicknesses << len;
  return thicknesses;
}

bool LayerLineProfile::ExportThickness(const QString &filename, const QList<LayerPointSet *> &splines, int nSample)
{
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return false;

  QTextStream out(&file);
  out << "Laplace line";
  for (int i = 0; i <= splines.size(); i++)
    out << "," << QString("Layer %1").arg(i+1);
  out << "\n";
  for (size_t i = 0; i < m_ptsProfile.size(); i++)
  {
    QList<double> vals = ComputeThicknessAlongLineProfile(m_ptsProfile[i], splines, nSample);
    out << i;
    for (int j = 0; j < vals.size(); j++)
    {
      out << "," << vals[j];
    }
    out << "\n";
  }

  return true;
}

bool LayerLineProfile::Save(const QString &filename)
{
  QFile file(filename);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text) || !m_spline0 || !m_spline1)
    return false;

  QTextStream out(&file);
  m_spline0->UpdateLabelData();
  out << m_spline0->GetPointSetData()->WriteAsControlPointsToString();
  out << "\n---\n";
  m_spline1->UpdateLabelData();
  out << m_spline1->GetPointSetData()->WriteAsControlPointsToString();
  out << "\n---\n";

  out << QString("\nViewport %1\nResolution %2\nSample %3\nOffset %4\n").arg(m_nPlane).arg(m_dResolution).arg(m_nSamples)
         .arg(m_dOffset);
  return true;
}

LayerLineProfile* LayerLineProfile::Load(const QString &filename, LayerMRI* ref)
{
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    return NULL;

  QTextStream in (&file);
  QString content = in.readAll();
  QStringList slist = content.split("\n---\n", MD_SkipEmptyParts);
  if (slist.size() < 3)
    return NULL;

  LayerPointSet* ptset0 = new LayerPointSet(ref, LayerPropertyPointSet::ControlPoint);
  ptset0->SetName("Spline 1");
  ptset0->GetProperty()->SetSnapToVoxelCenter(false);
  ptset0->GetProperty()->SetShowSpline(true);
  if (!ptset0->LoadFromString(slist[0]))
    return NULL;

  LayerPointSet* ptset1 = new LayerPointSet(ref, LayerPropertyPointSet::ControlPoint);
  ptset1->SetName("Spline 2");
  ptset1->GetProperty()->SetSnapToVoxelCenter(false);
  ptset1->GetProperty()->SetShowSpline(true);
  if (!ptset1->LoadFromString(slist[1]))
    return NULL;

  QStringList ar = slist[2].split("\n", MD_SkipEmptyParts);
  if (ar.size() < 3)
  {
    delete ptset0;
    delete ptset1;
    return NULL;
  }

  QStringList sublist = ar[0].split(" ", MD_SkipEmptyParts);
  if (sublist.size() < 2)
    return NULL;
  int nPlane = sublist[1].toInt();

  sublist = ar[1].split(" ", MD_SkipEmptyParts);
  if (sublist.size() < 2)
    return NULL;
  double dResolution = sublist[1].toDouble();

  sublist = ar[2].split(" ", MD_SkipEmptyParts);
  if (sublist.size() < 2)
    return NULL;
  int nSamples = sublist[1].toInt();

  LayerLineProfile* lp = new LayerLineProfile(nPlane, NULL, ptset0, ptset1);
  lp->m_dResolution = dResolution;
  lp->m_nSamples = nSamples;

  if (ar.size() > 3)
  {
    sublist = ar[3].split(" ", MD_SkipEmptyParts);
    if (sublist.size() > 1)
    {
      lp->m_dOffset = sublist[1].toDouble();
    }
  }

  return lp;
}

vtkActor* LayerLineProfile::GetLineProfileActor()
{
  return m_profileLines;
}

void LayerLineProfile::SetActiveLineId(int nId)
{
  m_nActiveLineId = nId;
  UpdateActiveLine();
}
