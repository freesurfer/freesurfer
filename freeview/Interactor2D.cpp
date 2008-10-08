/**
 * @file  Interactor2D.cpp
 * @brief Interactor2D to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2002-2009,
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

#include "Interactor2D.h"
#include "RenderView2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>

Interactor2D::Interactor2D() : Interactor(),
		m_nMousePosX( -1 ),
		m_nMousePosY( -1 ),
		m_bWindowLevel( false ),
		m_bChangeSlice( false ),
		m_bMovingCursor( false )
{
}

Interactor2D::~Interactor2D()
{
}

bool Interactor2D::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
	
	m_nMousePosX = event.GetX();
	m_nMousePosY = event.GetY();
	
	view->UpdateAnnotation();
	
	if ( event.ControlDown() && !event.ShiftDown() )
	{
		if ( event.LeftDown() )
		{
			view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, true );
			return false;
		}
		else if ( event.RightDown() )
		{
			view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, false );
			return false;
		}
	}
	
	if ( event.LeftDown() )
	{
		m_nDownPosX = m_nMousePosX;
		m_nDownPosY = m_nMousePosY;
		
		if ( event.ShiftDown() && !event.ControlDown() )
			m_bWindowLevel = true;
		else 
		{
			m_bMovingCursor = true;
			view->UpdateCursorRASPosition( m_nMousePosX, m_nMousePosY );
			view->NeedRedraw();
		}
	}
	else if ( event.MiddleDown() && event.ControlDown() )
	{
		m_bWindowLevel = true;	
	}	
	else
	{				
		return Interactor::ProcessMouseDownEvent( event, renderview );	// pass down the event
	}	
		
	return false;	// do not pass down the event
}

bool Interactor2D::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
	
	m_nMousePosX = event.GetX();
	m_nMousePosY = event.GetY();
	m_bWindowLevel = false;
	m_bChangeSlice = false;
	m_bMovingCursor = false;
	
	view->UpdateAnnotation();
	view->UpdateCursor2D();
	
	if ( event.LeftUp() )
	{
		return false;	
	}
	else
	{		
		return Interactor::ProcessMouseUpEvent( event, renderview );
	}
}

bool Interactor2D::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
		
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	if ( !lcm->HasAnyLayer() )
	{
		return Interactor::ProcessMouseMoveEvent( event, renderview );
	}
	
	int posX = event.GetX();
	int posY = event.GetY();
				
	if ( m_bChangeSlice )
	{
		double* voxelSize =	lcm->GetLayerCollection( "MRI" )->GetWorldVoxelSize();
		int nPlane = view->GetViewPlane();
		double dPixelPer = -0.20;
		double dPosDiff =  ( ( (int)( dPixelPer * ( posY - m_nDownPosY ) ) ) / dPixelPer - 
					( (int)( dPixelPer * ( m_nMousePosY - m_nDownPosY ) ) ) / dPixelPer )
					* dPixelPer * voxelSize[nPlane];
		if ( lcm->OffsetSlicePosition( nPlane, dPosDiff ) )
		{
			m_nMousePosX = posX;
			m_nMousePosY = posY;	
		}
	}
	else if ( m_bMovingCursor )
	{
		view->UpdateCursorRASPosition( posX, posY );
		view->NeedRedraw();
	}
	else if ( m_bWindowLevel )
	{
		LayerMRI* layer = ( LayerMRI* )(lcm->GetLayerCollection( "MRI" )->GetActiveLayer());
		if ( layer && layer->IsVisible() && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::Grayscale )
		{
			double scaleX = 0.002;
			double scaleY = 0.002;
			double w = ( posX - m_nMousePosX ) * scaleX;
			double l = ( posY - m_nMousePosY ) * scaleY;
			double scaleOverall = layer->GetProperties()->GetMaxValue() - 
						layer->GetProperties()->GetMinValue();
			w *= scaleOverall;
			l *= scaleOverall;
			w += layer->GetProperties()->GetWindow();
			l += layer->GetProperties()->GetLevel();
			if ( w < 0 )
				w = 0;
			layer->GetProperties()->SetWindowLevel( w, l );
		}
		m_nMousePosX = posX;
		m_nMousePosY = posY;	
	}
	else
	{	
		if ( event.MiddleIsDown() || event.RightIsDown() )
		{
			view->UpdateAnnotation();
			view->UpdateCursor2D();
			if ( event.RightIsDown() )
				view->SendBroadcast( "Zooming", view );
		}
		else
			view->UpdateMouseRASPosition( posX, posY );

		return Interactor::ProcessMouseMoveEvent( event, renderview );
	}
	
	return false;
}

void Interactor2D::ProcessPostMouseWheelEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
	view->UpdateAnnotation();
	view->UpdateCursor2D();
	view->NeedRedraw();
	view->SendBroadcast( "Zooming", view );
	
	Interactor::ProcessPostMouseWheelEvent( event, renderview );
}

void Interactor2D::ProcessPostMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
	if ( event.RightIsDown() )
	{
		view->UpdateCursor2D();
		view->NeedRedraw();
	}
	
	Interactor::ProcessPostMouseMoveEvent( event, renderview );
}

bool Interactor2D::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
		
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	if ( !lcm->HasAnyLayer() )
	{
		return Interactor::ProcessKeyDownEvent( event, renderview );
	}
		
	int nKeyCode = event.GetKeyCode();
	if ( nKeyCode == WXK_PAGEUP )
	{
		double* voxelSize = lcm->GetLayerCollection( "MRI" )->GetWorldVoxelSize();
		int nPlane = view->GetViewPlane();
		lcm->OffsetSlicePosition( nPlane, voxelSize[nPlane] );
	}
	else if ( nKeyCode == WXK_PAGEDOWN)
	{
		double* voxelSize = lcm->GetLayerCollection( "MRI" )->GetWorldVoxelSize();
		int nPlane = view->GetViewPlane();
		lcm->OffsetSlicePosition( nPlane, -voxelSize[nPlane] );
	}
	else if ( nKeyCode == WXK_UP )
	{
		view->MoveUp();
	}
	else if ( nKeyCode == WXK_DOWN )
	{
		view->MoveDown();
	}
	else if ( nKeyCode == WXK_LEFT )
	{
		view->MoveLeft();
	}
	else if ( nKeyCode == WXK_RIGHT )
	{
		view->MoveRight();
	}
	else if ( nKeyCode == '3' || nKeyCode == 'W' || nKeyCode == 'S' || nKeyCode == 'R' || nKeyCode == 'F' )
	{
		// do nothing, just intercept these keycodes
	}
	else
		return Interactor::ProcessKeyDownEvent( event, view );
	
	return false;
}
