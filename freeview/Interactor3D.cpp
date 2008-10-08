/**
 * @file  Interactor3D.cpp
 * @brief Interactor3D to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.7 $
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

#include "Interactor3D.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>

Interactor3D::Interactor3D() : 
	Interactor(),
	m_nMousePosX( -1 ),
	m_nMousePosY( -1 ),
	m_bWindowLevel( false )
{
}

Interactor3D::~Interactor3D()
{
}

bool Interactor3D::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView3D* view = ( RenderView3D* )renderview;
	
	m_nMousePosX = event.GetX();
	m_nMousePosY = event.GetY();
	
	view->CancelUpdateMouseRASPosition();
	
	if ( event.MiddleDown() && event.ControlDown() )
	{
		m_bWindowLevel = true;	
	}	
	else
	{				
		return Interactor::ProcessMouseDownEvent( event, renderview );	// pass down the event
	}	
		
	return false;	// do not pass down the event
}

bool Interactor3D::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{	
	RenderView3D* view = ( RenderView3D* )renderview;
	
	m_bWindowLevel = false;
	if ( event.GetX() == m_nMousePosX && event.GetY() == m_nMousePosY )
	{
		if ( event.LeftUp() )
		{
			view->UpdateCursorRASPosition( event.GetX(), event.GetY() );
		}
	}
	else
	{
		m_nMousePosX = event.GetX();
		m_nMousePosY = event.GetY();
	}	
	
	return Interactor::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3D::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView3D* view = ( RenderView3D* )renderview;
	
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	if ( !lcm->HasAnyLayer() )
	{
		return Interactor::ProcessMouseMoveEvent( event, renderview );
	}
	
	int posX = event.GetX();
	int posY = event.GetY();
			
	if ( m_bWindowLevel )
	{
		LayerMRI* layer = ( LayerMRI* )(lcm->GetLayerCollection( "MRI" )->GetActiveLayer());
		if ( layer && layer->IsVisible() )
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
		if ( event.MiddleIsDown() || event.RightIsDown() || event.LeftIsDown() )
		{
			view->CancelUpdateMouseRASPosition();
		}
		else
			view->UpdateMouseRASPosition( posX, posY );

		return Interactor::ProcessMouseMoveEvent( event, view );
	}
	
	UpdateCursor( event, view );
	return false;
}

bool Interactor3D::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{	
	RenderView3D* view = ( RenderView3D* )renderview;
		
	LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
	if ( !lcm->HasAnyLayer() )
	{
		return Interactor::ProcessKeyDownEvent( event, renderview );
	}
		
	int nKeyCode = event.GetKeyCode();
	if  ( nKeyCode == WXK_UP )
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
	else if ( nKeyCode == 'R' || nKeyCode == 'F' )
	{
		// do nothing, just intercept these keycodes
	}
	else
		return Interactor::ProcessKeyDownEvent( event, view );
	
	return false;
}
