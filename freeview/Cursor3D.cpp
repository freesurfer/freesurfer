/**
 * @file  Cursor3D.cpp
 * @brief Cursor for 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:57 $
 *    $Revision: 1.12 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
    m_view( view )
{
  m_actorCursor = vtkSmartPointer<vtkActor>::New();
  m_actorCursor->GetProperty()->SetColor( 1, 0, 0 );
  m_actorCursor->PickableOff();

  RebuildActor();
}

Cursor3D::~Cursor3D()
{}

void Cursor3D::RebuildActor()
{
// vtkRenderer* renderer = m_view->GetRenderer();

  double dLen = 1.5;
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
  tube->SetInput( polydata );
  tube->SetNumberOfSides( 12 );
  tube->SetRadius( 0.1 );
  tube->CappingOn();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( tube->GetOutputPort() );

  m_actorCursor->SetMapper( mapper );
  emit Updated();
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
    pos[i] = m_dPosition[i];
}


void Cursor3D::Show( bool bShow )
{
  m_actorCursor->SetVisibility( bShow?1:0 );
}

bool Cursor3D::IsShown()
{
  return m_actorCursor->GetVisibility() > 0;
}
