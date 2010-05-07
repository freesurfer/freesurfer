/**
 * @file  SurfaceRegion.cpp
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/05/07 20:06:30 $
 *    $Revision: 1.2 $
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

#include "SurfaceRegion.h"
#include "vtkRenderer.h"
#include "vtkActor2D.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "MainWindow.h"
#include "RenderView3D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkSelectPolyData.h"
#include "vtkProperty.h"
#include "vtkClipPolyData.h"
#include "MyUtils.h"

SurfaceRegion::SurfaceRegion()
{
  m_actorMesh = vtkSmartPointer<vtkActor>::New();
  m_actorMesh->GetProperty()->SetColor( 0, 0, 1 );
  m_actorMesh->GetProperty()->SetRepresentationToWireframe();
  m_actorMesh->GetProperty()->SetLineWidth( 2 );

  m_actorOutline = vtkSmartPointer<vtkActor>::New();
  m_actorOutline->GetProperty()->SetColor( 0, 0, 1 );
  m_actorOutline->GetProperty()->SetLineWidth( 3 );
  
  m_points = vtkSmartPointer<vtkPoints>::New();
  m_selector = vtkSmartPointer<vtkSelectPolyData>::New();
  m_selector->SetSelectionModeToSmallestRegion();
}

SurfaceRegion::~SurfaceRegion()
{}

void SurfaceRegion::RebuildOutline( bool bClose )
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell( m_points->GetNumberOfPoints() + (bClose?1:0) );
  for ( int i = 0; i < m_points->GetNumberOfPoints(); i++ )
    lines->InsertCellPoint( i );
  if ( bClose )
    lines->InsertCellPoint( 0 );
  polydata->SetPoints( m_points );
  polydata->SetLines( lines );
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput( polydata );
  m_actorOutline->SetMapper( mapper );
}

void SurfaceRegion::SetInput( vtkPolyData* polydata )
{
  m_selector->SetInput( polydata );
}

void SurfaceRegion::AddPoint( double* pt )
{
  m_points->InsertNextPoint( pt );
  RebuildOutline( false );
}

void SurfaceRegion::Close()
{
  RebuildOutline( true );
  if ( m_points->GetNumberOfPoints() > 3 )
  {
    m_selector->SetLoop( m_points );
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_selector->GetOutputPort() );
    m_actorMesh->SetMapper( mapper );
    m_points->Modified();
  }
}

void SurfaceRegion::Update()
{}

void SurfaceRegion::AppendActor( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorMesh );
  renderer->AddViewProp( m_actorOutline );
}

void SurfaceRegion::SetColor( const wxColour& color )
{
  m_actorMesh->GetProperty()->SetColor( color.Red()/255.0, color.Green()/255.0, color.Blue()/255.0 );
}

wxColour SurfaceRegion::GetColor()
{
  double* rgb = m_actorMesh->GetProperty()->GetColor();
  return wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) );
}

void SurfaceRegion::Show( bool bShow )
{
  m_actorMesh->SetVisibility( bShow?1:0 );
}

