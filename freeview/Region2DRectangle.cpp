/**
 * @file  Region2DRectangle.cpp
 * @brief Region2DRectangle.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/08/03 20:29:27 $
 *    $Revision: 1.1 $
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

#include "Region2DRectangle.h"
#include "RenderView2D.h"
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCoordinate.h>
#include <vtkPolyData.h>

Region2DRectangle::Region2DRectangle( RenderView2D* view ) :
  Region2D( view )
{
  m_actorRect = vtkSmartPointer<vtkActor2D>::New();
  m_actorRect->GetProperty()->SetColor( 1, 0, 0 );
  m_actorRect->GetProperty()->SetOpacity( 0.5 );
  m_actorRect->VisibilityOff();
}

Region2DRectangle::~Region2DRectangle()
{}


void Region2DRectangle::SetRect( int left, int top, int w, int h )
{
  m_nX1 = left;
  m_nY1 = top;
  m_nX2 = left+w;
  m_nY2 = top+h;  
  
  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::UpdateWorldCoords()
{
  m_view->MousePositionToRAS( m_nX1, m_nY1, m_dPt0 ); 
  m_view->MousePositionToRAS( m_nX1, m_nY2, m_dPt1 ); 
  m_view->MousePositionToRAS( m_nX2, m_nY2, m_dPt2 ); 
  m_view->MousePositionToRAS( m_nX2, m_nY1, m_dPt3 ); 
}

void Region2DRectangle::SetTopLeft( int left, int top )
{
  m_nX1 = left;
  m_nY1 = top;
  
  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::SetBottomRight( int right, int bottom )
{
  m_nX2 = right;
  m_nY2 = bottom;
  
  UpdateWorldCoords();
  Update();
}

void Region2DRectangle::Update()
{
  if ( !m_actorRect->GetVisibility() )
    return;
  
  double pt0[3], pt1[3], pt2[3], pt3[3];
 
  m_view->WorldToViewport( m_dPt0[0], m_dPt0[1], m_dPt0[2], pt0[0], pt0[1], pt0[2] );
  m_view->WorldToViewport( m_dPt1[0], m_dPt1[1], m_dPt1[2], pt1[0], pt1[1], pt1[2] );
  m_view->WorldToViewport( m_dPt2[0], m_dPt2[1], m_dPt2[2], pt2[0], pt2[1], pt2[2] );
  m_view->WorldToViewport( m_dPt3[0], m_dPt3[1], m_dPt3[2], pt3[0], pt3[1], pt3[2] );
  
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> poly = vtkSmartPointer<vtkCellArray>::New();
  pts->InsertNextPoint( pt0 );
  pts->InsertNextPoint( pt1 );
  pts->InsertNextPoint( pt2 );
  pts->InsertNextPoint( pt3 );
  vtkIdType face[4] = { 0, 1, 2, 3 };
  poly->InsertNextCell( 4, face );
  
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( pts );
  polydata->SetPolys( poly );
  vtkSmartPointer<vtkPolyDataMapper2D> mapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
  mapper->SetInput( polydata );
  
  vtkSmartPointer<vtkCoordinate> coords = vtkSmartPointer<vtkCoordinate>::New();
  coords->SetCoordinateSystemToViewport();
  mapper->SetTransformCoordinate( coords );
  m_actorRect->SetMapper( mapper );
}

void Region2DRectangle::AppendProp( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorRect );
}

void Region2DRectangle::Show( bool bShow )
{
  m_actorRect->SetVisibility( bShow?1:0 );
  Update();
}

void Region2DRectangle::GetWorldPoint1( double* pt )
{
  pt[0] = m_dPt0[0];
  pt[1] = m_dPt0[1];
  pt[2] = m_dPt0[2];
}
  
void Region2DRectangle::GetWorldPoint2( double* pt )
{
  pt[0] = m_dPt2[0];
  pt[1] = m_dPt2[1];
  pt[2] = m_dPt2[2];
}

bool Region2DRectangle::Contains( int nX, int nY, int* nIndexOut )
{
  return false;
}

void Region2DRectangle::Offset( int nX, int nY )
{
}
  
void Region2DRectangle::UpdatePoint( int nIndex, int nX, int nY )
{
}
