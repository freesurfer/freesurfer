#include "Region3D.h"
#include "vtkRenderer.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "MainWindow.h"
#include "RenderView3D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "vtkSplineFilter.h"
#include <vtkKdTreePointLocator.h>
#include "vtkCleanPolyData.h"
#include <QDebug>

Region3D::Region3D(LayerMRI* owner ) : QObject(owner), m_mri(owner)
{
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif

  m_actor = vtkSmartPointer<vtkActor>::New();
  m_actor->GetProperty()->SetColor( 0, 1, 0 );
  m_actor->GetProperty()->SetLineWidth( 5*ratio );

  m_points = vtkSmartPointer<vtkPoints>::New();
  m_interpolatedPoints = vtkSmartPointer<vtkPoints>::New();
  m_color = Qt::blue;

  m_locator = vtkSmartPointer<vtkKdTreePointLocator>::New();
  QList<vtkActor*> list = m_mri->GetContourActors(true);
  if (!list.isEmpty())
  {
    m_locator->SetDataSet(vtkPolyData::SafeDownCast( list.first()->GetMapper()->GetInput() ));
    m_locator->BuildLocator();
  }
}

Region3D::~Region3D()
{}

void Region3D::AddPoint( double* pt )
{
  m_points->InsertNextPoint( pt );
  RebuildOutline( false );
}

void Region3D::RebuildOutline( bool bInterpolate )
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//  if ( bClose && m_points->GetNumberOfPoints() > 0 )
//  {
//    m_points->InsertNextPoint( m_points->GetPoint( 0 ) );
//  }
  lines->InsertNextCell( m_points->GetNumberOfPoints() );
  for ( int i = 0; i < m_points->GetNumberOfPoints(); i++ )
  {
    lines->InsertCellPoint( i );
  }
  polydata->SetPoints( m_points );
  polydata->SetLines( lines );
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  if (bInterpolate)
  {
    vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
    spline->SetInputData(polydata);
    int dim[3];
    double vs[3];
    m_mri->GetVolumeInfo(dim, vs);
    spline->SetSubdivideToLength();
    spline->SetLength(qMin(qMin(vs[0],vs[1]), vs[2])/3);
    spline->Update();
    m_interpolatedPoints = spline->GetOutput()->GetPoints();
    polydata = vtkSmartPointer<vtkPolyData>::New();
    lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell( m_interpolatedPoints->GetNumberOfPoints() );
    for ( int i = 0; i < m_interpolatedPoints->GetNumberOfPoints(); i++ )
    {
      lines->InsertCellPoint( i );
    }
    polydata->SetPoints( m_interpolatedPoints );
    polydata->SetLines( lines );
    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(polydata);
    mapper->SetInputConnection(clean->GetOutputPort());
  }
  else
  {
  #if VTK_MAJOR_VERSION > 5
    mapper->SetInputData( polydata );
  #else
    mapper->SetInput( polydata );
  #endif
  }

  m_actor->SetMapper( mapper );
}

bool Region3D::Close()
{
  RebuildOutline( true );
  return true;
}

void Region3D::Update()
{}

void Region3D::AppendProps( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actor );
}

void Region3D::Highlight( bool bHighlight )
{
  if ( bHighlight )
  {
    m_actor->GetProperty()->SetColor( 1, 1, 1 );
  }
  else
  {
    m_actor->GetProperty()->SetColor( m_color.redF(), m_color.greenF(), m_color.blueF() );
  }
}

void Region3D::Show( bool bShow )
{
  m_actor->SetVisibility( bShow?1:0 );
}

void Region3D::SetColor( const QColor& color )
{
  m_color = color;
  m_actor->GetProperty()->SetColor( color.redF(), color.greenF(), color.blueF() );
  emit ColorChanged( color );
}

QColor Region3D::GetColor()
{
  return m_color;
}

bool Region3D::HasPoint(double *pt, double dist)
{
  double th = dist*dist;
  vtkPoints* pts = m_points;
  if (m_interpolatedPoints.GetPointer())
    pts = m_interpolatedPoints;
  for (int i = 0; i < pts->GetNumberOfPoints(); i++)
  {
    if (vtkMath::Distance2BetweenPoints(pt, pts->GetPoint(i)) < th)
      return true;
  }
  return false;
}

bool Region3D::WriteHeader( FILE* fp, LayerMRI* mri, int nNum )
{
  QString strg = QString( "VOLUME_PATH %1\nVOLUME_THRESHOLD %2 %3\nNUM_OF_REGIONS %4\n" )
      .arg( mri->GetFileName() )
      .arg( mri->GetProperty()->GetContourMinThreshold() )
      .arg( mri->GetProperty()->GetContourMaxThreshold() )
      .arg( nNum );
  QFile file;
  file.open( fp, QIODevice::Append );
  QByteArray ba = strg.toLatin1();
  int nsize = file.write( ba );
  return nsize == ba.size();
}

bool Region3D::WriteBody( FILE* fp )
{
  QString strg = QString( "REGION\nPOINTS %1\n" )
      .arg(m_points->GetNumberOfPoints());
  double pt[3];
  for ( vtkIdType i = 0; i < m_points->GetNumberOfPoints(); i++ )
  {
    m_points->GetPoint( i, pt );
    m_mri->TargetToRAS( pt, pt );
    strg += QString("%1 %2 %3\n").arg( pt[0] ).arg( pt[1] ).arg( pt[2] );
  }
  strg += QString("COLOR %1 %2 %3\n").arg( m_color.red() ).arg( m_color.green() ).arg( m_color.blue() );
  QFile file;
  file.open( fp, QIODevice::Append );
  QByteArray ba = strg.toLatin1();
  int nsize = file.write( ba );
  return nsize == ba.size();
}

bool Region3D::Load( FILE* fp )
{
  char tmp_strg[1000];
  QString id_strg = "REGION";
  while ( fscanf( fp, "%s\n", tmp_strg ) != EOF && id_strg != tmp_strg )
  {
    ;
  }
  if ( id_strg != tmp_strg )
  {
    return false;
  }

  int nPts = 0;
  float x, y, z;
  fscanf( fp, "POINTS %d", &nPts );
  double pt[3];
  for ( int i = 0; i < nPts; i++ )
  {
    fscanf( fp, "%f %f %f", &x, &y, &z );
    pt[0] = x;
    pt[1] = y;
    pt[2] = z;
    m_mri->RASToTarget( pt, pt );
    m_points->InsertNextPoint( pt );
  }

  int r, g, b;
  if ( fscanf( fp, "\nCOLOR %d %d %d", &r, &g, &b ) == EOF)
  {
    return false;
  }
  m_color = QColor(r, g, b);

  RebuildOutline( true );

  return true;
}
