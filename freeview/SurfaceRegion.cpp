/**
 * @file  SurfaceRegion.cpp
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/05/17 20:06:22 $
 *    $Revision: 1.4 $
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
#include "vtkBox.h"
#include "vtkMath.h"
#include "MyUtils.h"

SurfaceRegion::SurfaceRegion()
{
  m_actorMesh = vtkSmartPointer<vtkActor>::New();
  m_actorMesh->GetProperty()->SetColor( 0, 1, 0 );
  m_actorMesh->GetProperty()->SetRepresentationToWireframe();
  m_actorMesh->GetProperty()->SetLineWidth( 2 );

  m_actorOutline = vtkSmartPointer<vtkActor>::New();
  m_actorOutline->GetProperty()->SetColor( 0, 1, 0 );
  m_actorOutline->GetProperty()->SetLineWidth( 3 );
  
  m_points = vtkSmartPointer<vtkPoints>::New();
  m_selector = vtkSmartPointer<vtkSelectPolyData>::New();
  m_selector->SetSelectionModeToSmallestRegion();
  
  m_clipbox = vtkSmartPointer<vtkBox>::New();
  m_clipper = vtkSmartPointer<vtkClipPolyData>::New();
  m_clipper->SetClipFunction( m_clipbox );
//  m_clipper->GenerateClippedOutputOn();
  m_clipper->InsideOutOn();
  m_selector->SetInputConnection( m_clipper->GetOutputPort() );
}

SurfaceRegion::~SurfaceRegion()
{}

vtkActor* SurfaceRegion::GetMeshActor()
{
  return m_actorMesh;
}

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
  m_actorOutline->SetVisibility( bClose ? 0 : 1 );
}

void SurfaceRegion::SetInput( vtkPolyData* polydata )
{
  m_clipper->SetInput( polydata );
  m_clipbox->SetBounds( polydata->GetBounds() );
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
    double bounds[6], cpt[3], len = 0;
    vtkPolyDataMapper::SafeDownCast( m_actorOutline->GetMapper() )->GetInput()->GetBounds( bounds );
    for ( int i = 0; i < 3; i++ )
    {
      cpt[i] = (bounds[i*2+1] + bounds[i*2]) / 2.0;
      len += (bounds[i*2+1] - bounds[i*2]);
    }
    for ( int i = 0; i < 3; i++ )
    {
      bounds[i*2] = cpt[i] - len/2.0;
      bounds[i*2+1] = cpt[i] + len/2.0;
    }
    m_clipbox->SetBounds( bounds );
    m_selector->SetLoop( m_points );
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_selector->GetOutputPort() );
    mapper->ScalarVisibilityOff();
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

bool SurfaceRegion::HasPoint( double* pos )
{
  double delta[3] = { 0, 0, 0 };
  return vtkMath::PointIsWithinBounds( pos, m_actorMesh->GetBounds(), delta );
}
  
void SurfaceRegion::Highlight( bool bHighlight )
{
  if ( bHighlight )
    m_actorMesh->GetProperty()->SetColor( 0, 1, 0 );
  else 
    m_actorMesh->GetProperty()->SetColor( 0, 0, 1 );
}

