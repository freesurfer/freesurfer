/**
 * @file  Interactor2DVoxelEdit.cpp
 * @brief Interactor2DVoxelEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.4.2.1 $
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

#include "Interactor2DVoxelEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>

Interactor2DVoxelEdit::Interactor2DVoxelEdit() : Interactor2D()
{
	m_bEditing = false;
}

Interactor2DVoxelEdit::~Interactor2DVoxelEdit()
{
}

bool Interactor2DVoxelEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;

	if ( event.LeftDown() )
	{
		if ( event.ControlDown() )
			return Interactor2D::ProcessMouseDownEvent( event, renderview );
		
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( "MRI" );
		LayerMRI* mri = ( LayerMRI* )lc->GetActiveLayer();
		if ( (!mri || !mri->IsVisible()) ) //&& ( event.ControlDown() || m_nAction == EM_Polyline ) )
		{
			SendBroadcast( "MRINotVisible", this );
		}
		else if ( !mri->IsEditable() ) //&& ( event.ControlDown() || m_nAction == EM_Polyline ) )
		{
			SendBroadcast( "MRINotEditable", this );
		}
		else
		{
			m_nMousePosX = event.GetX();
			m_nMousePosY = event.GetY();
			
			double ras[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
			if ( m_nAction == EM_Freehand ) //&& ( event.ControlDown() ) )
			{
				mri->SaveForUndo( view->GetViewPlane() );			
				m_bEditing = true;				
				mri->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
			}
			else if ( m_nAction == EM_Fill ) //&& ( event.ControlDown() ) )
			{
				mri->SaveForUndo( view->GetViewPlane() );
				mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
			}
			else if ( m_nAction == EM_Polyline )
			{
				m_bEditing = true;
				double ras2[3];
				view->GetCursor()->GetPosition( ras2 );
				view->GetCursor()->SetPosition( ras );
				view->GetCursor()->SetPosition2( ras );			
				if ( m_dPolylinePoints.size() > 0 )	
				{
					mri->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), !event.ShiftDown() );	
				}			
				else
				{
					mri->SaveForUndo( view->GetViewPlane() );
					m_dPolylinePoints.push_back( ras[0] );
					m_dPolylinePoints.push_back( ras[1] );
					m_dPolylinePoints.push_back( ras[2] );
				}	
				
				view->ReleaseMouse();
				view->CaptureMouse();				
			}
			else
				return Interactor2D::ProcessMouseDownEvent( event, renderview );
		}
		
		return false;
	}
	else if ( m_bEditing )
	{		
		m_bEditing = false;
		if ( m_nAction == EM_Polyline )
		{
			if ( event.MiddleDown() )
			{
				view->GetCursor()->Update();
				view->NeedRedraw();
			}
			else if ( event.RightDown() )
			{
				if ( m_dPolylinePoints.size() > 0 )
				{
					LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
					LayerMRI* mri = ( LayerMRI* )lc->GetActiveLayer();
					
					double ras1[3] = { m_dPolylinePoints[0], m_dPolylinePoints[1], m_dPolylinePoints[2] };
					double ras2[3];
					view->GetCursor()->GetPosition( ras2 );
					view->GetCursor()->SetPosition2( ras2 );
					view->GetCursor()->SetPosition( ras1 );
					mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() );	
				}
			}
		}
		
		m_dPolylinePoints.clear();
		view->ReleaseMouse();
		
		return false;
	}					
	return Interactor2D::ProcessMouseDownEvent( event, renderview );	// pass down the event	
}

bool Interactor2DVoxelEdit::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
//	RenderView2D* view = ( RenderView2D* )renderview;
	
	if ( m_bEditing )
	{
		m_nMousePosX = event.GetX();
		m_nMousePosY = event.GetY();
		
		if ( !event.LeftUp() || m_nAction != EM_Polyline || m_dPolylinePoints.size() == 0 )
			m_bEditing = false;		
		
		return false;	
	}
	else
	{		
		return Interactor2D::ProcessMouseUpEvent( event, renderview );
	}
}

bool Interactor2DVoxelEdit::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
	
	if ( m_bEditing )
	{
		int posX = event.GetX();
		int posY = event.GetY();
		
		if ( m_nAction == EM_Freehand )
		{
			LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
			LayerMRI* mri = ( LayerMRI* )lc->GetActiveLayer();
					
			double ras1[3], ras2[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
			view->MousePositionToRAS( posX, posY, ras2 );
			
			mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() );
		}
		else if ( m_nAction == EM_Polyline )
		{
			double ras[3];
			view->MousePositionToRAS( posX, posY, ras );
			view->GetCursor()->SetPosition2( ras );
			view->GetCursor()->SetPosition( view->GetCursor()->GetPosition(), true );
			view->NeedRedraw();
		}
		
		m_nMousePosX = posX;
		m_nMousePosY = posY;
		
		return false;
	} 
	else
	{	
		return Interactor2D::ProcessMouseMoveEvent( event, renderview );
	}
}

bool Interactor2DVoxelEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
	if ( !m_bEditing )
		return Interactor2D::ProcessKeyDownEvent( event, renderview );
	else
		return false;
}

void Interactor2DVoxelEdit::UpdateCursor( wxWindow* wnd )
{
/*	if ( wnd->FindFocus() == wnd && m_nAction == EM_Freehand )
		wnd->SetCursor( wxCURSOR_PENCIL );	
	else
	wnd->SetCursor( wxNullCursor );*/
}
