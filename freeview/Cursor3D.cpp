/**
 * @brief Cursor for 3D view.
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

#include "Cursor3D.h"
#include "vtkRenderer.h"
#include "vtkActor2D.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include <vtkTubeFilter.h>
#include "vtkPolyData.h"
#include "MainWindow.h"
#include "RenderView3D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include <vtkProperty.h>
#include "MyUtils.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"


Cursor3D::Cursor3D( RenderView3D* view ) : QObject( view ),
  m_view( view ), m_nSize(5), m_nThickness(1), m_dScale(1.0)
{
  m_actorCursor = vtkSmartPointer<vtkActor>::New();
  m_actorCursor->GetProperty()->SetColor( 1, 0, 0 );
  m_actorCursor->PickableOff();

  RebuildActor();
}

Cursor3D::~Cursor3D()
{}

void Cursor3D::RebuildActor(double scale)
{
  // vtkRenderer* renderer = m_view->GetRenderer();

  if (scale > 0)
    m_dScale = scale;

  double dLen = 1.5*m_dScale * (1 + (m_nSize-1.0)/5.0);
  int n = 0;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  points->InsertNextPoint( dLen, 0, 0 );
  points->InsertNextPoint( -dLen, 0, 0 );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  points->InsertNextPoint( 0, dLen, 0 );
  points->InsertNextPoint( 0, -dLen, 0 );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );
  points->InsertNextPoint( 0, 0, dLen);
  points->InsertNextPoint( 0, 0, -dLen );
  lines->InsertNextCell( 2 );
  lines->InsertCellPoint( n++ );
  lines->InsertCellPoint( n++ );

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );

  vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION > 5
  tube->SetInputData( polydata );
#else
  tube->SetInput( polydata );
#endif
  tube->SetNumberOfSides( 12 );
  tube->SetRadius( 0.15*(1+(m_nThickness-1)/3.0) );
  tube->CappingOn();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( tube->GetOutputPort() );

  m_actorCursor->SetMapper( mapper );
  emit Updated();
}

void Cursor3D::SetSize(int nSize)
{
  m_nSize = qMax(1, nSize);
  RebuildActor();
}

void Cursor3D::SetThickness(int nThickness)
{
  m_nThickness = qMax(1, nThickness);
  RebuildActor();
}

void Cursor3D::SetPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dPosition[i] = pos[i];
  }

  m_actorCursor->SetPosition( pos );
}

void Cursor3D::Update()
{}

void Cursor3D::AppendActor( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorCursor );
}

void Cursor3D::SetColor( double r, double g, double b )
{
  m_actorCursor->GetProperty()->SetColor( r, g, b );
  emit Updated();
}

void Cursor3D::SetColor( const QColor& color )
{
  SetColor( color.redF(), color.greenF(), color.blueF() );
}

void Cursor3D::GetColor( double* rgb )
{
  m_actorCursor->GetProperty()->GetColor( rgb );
}

QColor Cursor3D::GetColor()
{
  double* rgb = m_actorCursor->GetProperty()->GetColor();
  return QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) );
}

double* Cursor3D::GetPosition()
{
  return m_dPosition;
}

void Cursor3D::GetPosition( double* pos )
{
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = m_dPosition[i];
  }
}


void Cursor3D::Show( bool bShow )
{
  m_actorCursor->SetVisibility( bShow?1:0 );
}

bool Cursor3D::IsShown()
{
  return m_actorCursor->GetVisibility() > 0;
}
