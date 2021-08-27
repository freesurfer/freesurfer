/**
 * @brief Cursor for 2D view.
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

#include "Cursor2D.h"
#include "vtkRenderer.h"
#include "vtkActor2D.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkPolyData.h"
#include "MainWindow.h"
#include "RenderView2D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include <vtkProperty2D.h>
#include <vtkCursor2D.h>
#include "MyUtils.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include <QDebug>

Cursor2D::Cursor2D( RenderView2D* view ) : QObject( view ),
  m_view( view ),
  m_nSize( 5 )
{
  m_actorCursor = vtkSmartPointer<vtkActor2D>::New();
  m_actorCursor->GetProperty()->SetColor( 1, 0, 0 );

  Update();
}

Cursor2D::~Cursor2D()
{}

void Cursor2D::SetPosition2( double* pos)
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dPosition2[i] = pos[i];
  }
}

void Cursor2D::SetPosition( double* pos, bool bConnectPrevious  )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dPosition[i] = pos[i];
  }

  Update( bConnectPrevious );
}

void Cursor2D::SetSize(int nSize)
{
  m_nSize = qMax(5, nSize);
  Update();
}

void Cursor2D::SetThickness(int nThickness)
{
  m_nThickness = qMax(1, nThickness);
  Update();
}

void Cursor2D::Update( bool bConnectPrevious )
{
  // vtkRenderer* renderer = m_view->GetRenderer();

  double dLen = ( m_nSize == 100 ? 100000 : m_nSize );
#if VTK_MAJOR_VERSION > 7
  dLen *= m_view->devicePixelRatio();
#endif

  int n = 0;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  double pos[3];
  if ( bConnectPrevious )
  {
    m_view->WorldToViewport( m_dPosition2[0], m_dPosition2[1], m_dPosition2[2], pos[0], pos[1], pos[2] );
    points->InsertNextPoint( pos );
    lines->InsertNextCell( 2 + m_dInterpolationPoints.size() / 3 );
    lines->InsertCellPoint( n++ );
    for ( size_t i = 0; i < m_dInterpolationPoints.size(); i+=3 )
    {
      m_view->WorldToViewport( m_dInterpolationPoints[i],
                               m_dInterpolationPoints[i+1],
          m_dInterpolationPoints[i+2],
          pos[0],
          pos[1],
          pos[2] );
      points->InsertNextPoint( pos[0],
          pos[1],
          pos[2] );
      lines->InsertCellPoint( n++ );
    }
    m_view->WorldToViewport( m_dPosition[0], m_dPosition[1], m_dPosition[2], pos[0], pos[1], pos[2] );
    points->InsertNextPoint( pos );
    lines->InsertCellPoint( n++ );
  }
  else
  {
    m_view->WorldToViewport( m_dPosition[0], m_dPosition[1], m_dPosition[2], pos[0], pos[1], pos[2] );
  }

  points->InsertNextPoint( pos[0] + dLen, pos[1], pos[2] );
  points->InsertNextPoint( pos[0] - dLen, pos[1], pos[2] );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  points->InsertNextPoint( pos[0], pos[1] + dLen, pos[2] );
  points->InsertNextPoint( pos[0], pos[1] - dLen, pos[2] );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  if ( m_nSize < 100 )
  {
    int w, h;
    w = m_view->rect().width();
    h = m_view->rect().height();
    int nd = 9;
#if VTK_MAJOR_VERSION > 7
  if (m_view->devicePixelRatio() > 1)
  {
      w *= m_view->devicePixelRatio();
      h *= m_view->devicePixelRatio();
      nd *= m_view->devicePixelRatio();
  }
#endif
    points->InsertNextPoint( 0, pos[1], pos[2] );
    points->InsertNextPoint( nd, pos[1], pos[2] );
    lines->InsertNextCell( 2 );
    lines->InsertCellPoint( n++ );
    lines->InsertCellPoint( n++ );
    points->InsertNextPoint( w-nd, pos[1], pos[2] );
    points->InsertNextPoint( w, pos[1], pos[2] );
    lines->InsertNextCell( 2 );
    lines->InsertCellPoint( n++ );
    lines->InsertCellPoint( n++ );
    points->InsertNextPoint( pos[0], 0, pos[2] );
    points->InsertNextPoint( pos[0], nd, pos[2] );
    lines->InsertNextCell( 2 );
    lines->InsertCellPoint( n++ );
    lines->InsertCellPoint( n++ );
    points->InsertNextPoint( pos[0], h-nd, pos[2] );
    points->InsertNextPoint( pos[0], h, pos[2] );
    lines->InsertNextCell( 2 );
    lines->InsertCellPoint( n++ );
    lines->InsertCellPoint( n++ );
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData( polydata );
#else
  mapper->SetInput(polydata);
#endif
  vtkSmartPointer<vtkCoordinate> coords = vtkSmartPointer<vtkCoordinate>::New();
  coords->SetCoordinateSystemToViewport();
  mapper->SetTransformCoordinate( coords );

  m_actorCursor->SetMapper( mapper );
#if VTK_MAJOR_VERSION > 7
  m_actorCursor->GetProperty()->SetLineWidth(m_nThickness*m_view->devicePixelRatio());
#else
  m_actorCursor->GetProperty()->SetLineWidth(m_nThickness);
#endif
  emit Updated();
}

/*
void Cursor2D::Update( bool bConnectPrevious )
{
// vtkRenderer* renderer = m_view->GetRenderer();

  double dLen = m_nRadius;

  int n = 0;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  if ( bConnectPrevious )
  {
    points->InsertNextPoint( m_dPosition2 );
    lines->InsertNextCell( 2 + m_dInterpolationPoints.size() / 3 );
    lines->InsertCellPoint( n++ );
    for ( size_t i = 0; i < m_dInterpolationPoints.size(); i+=3 )
    {
      points->InsertNextPoint( m_dInterpolationPoints[i],
                               m_dInterpolationPoints[i+1],
                               m_dInterpolationPoints[i+2] );
      lines->InsertCellPoint( n++ );
    }
    points->InsertNextPoint( m_dPosition );
    lines->InsertCellPoint( n++ );
  }
  points->InsertNextPoint( m_dPosition[0] + dLen, m_dPosition[1], m_dPosition[2] );
  points->InsertNextPoint( m_dPosition[0] - dLen, m_dPosition[1], m_dPosition[2] );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  points->InsertNextPoint( m_dPosition[0], m_dPosition[1] + dLen, m_dPosition[2] );
  points->InsertNextPoint( m_dPosition[0], m_dPosition[1] - dLen, m_dPosition[2] );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  points->InsertNextPoint( m_dPosition[0], m_dPosition[1], m_dPosition[2] + dLen );
  points->InsertNextPoint( m_dPosition[0], m_dPosition[1], m_dPosition[2] - dLen );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
  mapper->SetInput( polydata );

  vtkSmartPointer<vtkCoordinate> coords = vtkSmartPointer<vtkCoordinate>::New();
  coords->SetCoordinateSystemToWorld();
  mapper->SetTransformCoordinate( coords );

  m_actorCursor->SetMapper( mapper );
}
*/

void Cursor2D::AppendActor( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorCursor );
}

void Cursor2D::SetColor( double r, double g, double b )
{
  m_actorCursor->GetProperty()->SetColor( r, g, b );
  emit Updated();
}

void Cursor2D::SetColor( const QColor& color )
{
  SetColor( color.redF(), color.greenF(), color.blueF() );
}

void Cursor2D::GetColor( double* rgb )
{
  m_actorCursor->GetProperty()->GetColor( rgb );
}

QColor Cursor2D::GetColor()
{
  double* rgb = m_actorCursor->GetProperty()->GetColor();
  return QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) );
}

int Cursor2D::GetSize()
{
  return m_nSize;
}

double* Cursor2D::GetPosition()
{
  return m_dPosition;
}

void Cursor2D::GetPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = m_dPosition[i];
  }
}

void Cursor2D::SetInterpolationPoints( std::vector<double> pts )
{
  m_dInterpolationPoints = pts;
  // Update( true );
}

void Cursor2D::Show( bool bShow )
{
  m_actorCursor->SetVisibility( bShow?1:0 );
}

bool Cursor2D::IsShown()
{
  return m_actorCursor->GetVisibility() > 0;
}
