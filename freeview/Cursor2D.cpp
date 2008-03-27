/**
 * @file  Cursor2D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:14 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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
 
#include <wx/xrc/xmlres.h> 
#include "Cursor2D.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "MainWindow.h"
#include "RenderView2D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "MyUtils.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"


Cursor2D::Cursor2D( RenderView2D* view )
{
	m_actorCursor = vtkSmartPointer<vtkActor>::New();
	m_actorCursor->GetProperty()->SetColor( 1, 0, 0 );
	m_nRadius = 5;
	m_view = view;
}

Cursor2D::~Cursor2D()
{
}

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

void Cursor2D::Update( bool bConnectPrevious )
{
	vtkRenderer* renderer = m_view->GetRenderer();
	double pos1[3] = { 0, 0, 0 }, pos2[3] = { m_nRadius, 0, 0 };
	renderer->ViewportToNormalizedViewport( pos1[0], pos1[1] );
	renderer->ViewportToNormalizedViewport( pos2[0], pos2[1] );	
	renderer->NormalizedViewportToView( pos1[0], pos1[1], pos1[2] );
	renderer->NormalizedViewportToView( pos2[0], pos2[1], pos2[2] );
	renderer->ViewToWorld( pos1[0], pos1[1], pos1[2] );
	renderer->ViewToWorld( pos2[0], pos2[1], pos2[2] );
	
	double dLen = MyUtils::GetDistance<double>( pos1, pos2 );
	double pos[3] = { m_dPosition[0], m_dPosition[1], m_dPosition[2] };
	double prev_pos[3] = { m_dPosition2[0], m_dPosition2[1], m_dPosition2[2] };
	int nPlane = m_view->GetViewPlane();
	pos[nPlane] += ( nPlane == 2 ? -0.001 : 0.001 );
	prev_pos[nPlane] += ( nPlane == 2 ? -0.001 : 0.001 );
	int n = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	if ( bConnectPrevious )
	{
		points->InsertNextPoint( prev_pos );
		points->InsertNextPoint( pos );
		lines->InsertNextCell( 2 );
		lines->InsertCellPoint( n++ );
		lines->InsertCellPoint( n++ );
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
	points->InsertNextPoint( pos[0], pos[1], pos[2] + dLen );
	points->InsertNextPoint( pos[0], pos[1], pos[2] - dLen );
	lines->InsertNextCell( 2 );
	lines->InsertCellPoint( n++ );	
	lines->InsertCellPoint( n++ );
	
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints( points );
	polydata->SetLines( lines );
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput( polydata );
	m_actorCursor->SetMapper( mapper );	
}

void Cursor2D::AppendCursor( vtkRenderer* renderer )
{
	renderer->AddViewProp( m_actorCursor );
}

void Cursor2D::SetColor( double r, double g, double b )
{
	m_actorCursor->GetProperty()->SetColor( r, g, b );
}

void Cursor2D::GetColor( double* rgb )
{
	m_actorCursor->GetProperty()->GetColor( rgb );
}

int Cursor2D::GetRadius()
{
	return m_nRadius;
}

void Cursor2D::SetRadius( int nRadius )
{
	if ( nRadius > 1 )
	{
		m_nRadius = nRadius;
		Update();
	}
}

double* Cursor2D::GetPosition()
{
	return m_dPosition;
}

void Cursor2D::GetPosition( double* pos )
{
	for ( int i = 0; i < 3; i++ )
		pos[i] = m_dPosition[i];
}
