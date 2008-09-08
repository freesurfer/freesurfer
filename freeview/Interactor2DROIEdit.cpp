/**
 * @file  Interactor2DROIEdit.cpp
 * @brief Interactor2DROIEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/09/08 16:25:28 $
 *    $Revision: 1.10 $
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

#include "Interactor2DROIEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerROI.h"
#include "LayerMRI.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>
#include <wx/utils.h>

Interactor2DROIEdit::Interactor2DROIEdit() : 
		Interactor2D(),
		m_bEditing( false )
{
}

Interactor2DROIEdit::~Interactor2DROIEdit()
{
}

bool Interactor2DROIEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
		
	if ( event.LeftDown() || ( event.RightDown() && event.LeftIsDown() ) )
	{
		if ( event.ControlDown() && event.ShiftDown() )
			return Interactor2D::ProcessMouseDownEvent( event, renderview );
		
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( "ROI" );
		LayerROI* roi = ( LayerROI* )lc->GetActiveLayer();
		if ( !roi || !roi->IsVisible() )
		{
			SendBroadcast( "ROINotVisible", this );
		}
		else
		{
			m_nMousePosX = event.GetX();
			m_nMousePosY = event.GetY();
			
			double ras[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
			if ( m_nAction == EM_Freehand ) // && ( event.ControlDown() || event.ShiftDown() ) )
			{
				roi->SaveForUndo( view->GetViewPlane() );		
				if ( event.ControlDown() )
				{
					roi->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
				}
				else
				{	
					m_bEditing = true;				
					roi->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
				}
			}
			else if ( m_nAction == EM_Fill ) // && ( event.ControlDown() || event.ShiftDown() ) )
			{
				roi->SaveForUndo( view->GetViewPlane() );
				roi->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
			}
			else if ( m_nAction == EM_Polyline )
			{
				m_bEditing = true;
				double ras2[3];
				view->GetCursor2D()->GetPosition( ras2 );
				view->GetCursor2D()->SetPosition( ras );
				view->GetCursor2D()->SetPosition2( ras );			
				if ( m_dPolylinePoints.size() > 0 )	
				{
					roi->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );	
				}			
				else
				{
					roi->SaveForUndo( view->GetViewPlane() );
					m_dPolylinePoints.push_back( ras[0] );
					m_dPolylinePoints.push_back( ras[1] );
					m_dPolylinePoints.push_back( ras[2] );
				}					
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
				view->GetCursor2D()->Update();
				view->NeedRedraw();
			}
			else if ( event.RightDown() )
			{
				if ( m_dPolylinePoints.size() > 0 )
				{
					LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
					LayerROI* roi = ( LayerROI* )lc->GetActiveLayer();
					
					double ras1[3] = { m_dPolylinePoints[0], m_dPolylinePoints[1], m_dPolylinePoints[2] };
					double ras2[3];
					view->GetCursor2D()->GetPosition( ras2 );
					view->GetCursor2D()->SetPosition2( ras2 );
					view->GetCursor2D()->SetPosition( ras1 );
					roi->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );	
				}
			}
		}
		
		m_dPolylinePoints.clear();
		return false;
	}					
	return Interactor2D::ProcessMouseDownEvent( event, renderview );	// pass down the event	
}

bool Interactor2DROIEdit::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
//	RenderView2D* view = ( RenderView2D* )renderview;
	UpdateCursor( event, renderview );
	
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

bool Interactor2DROIEdit::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;			
	if ( m_bEditing )
	{
		UpdateCursor( event, view );	
		int posX = event.GetX();
		int posY = event.GetY();
		if ( m_nAction == EM_Freehand )
		{
			LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
			LayerROI* roi = ( LayerROI* )lc->GetActiveLayer();
					
			double ras1[3], ras2[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
			view->MousePositionToRAS( posX, posY, ras2 );
			
			roi->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
		}
		else if ( m_nAction == EM_Polyline )
		{
			double ras[3];
			view->MousePositionToRAS( posX, posY, ras );
			view->GetCursor2D()->SetPosition2( ras );
			view->GetCursor2D()->SetPosition( view->GetCursor2D()->GetPosition(), true );
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

bool Interactor2DROIEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
	UpdateCursor( event, renderview );
	
	if ( !m_bEditing )
		return Interactor2D::ProcessKeyDownEvent( event, renderview );
	else
		return false;
}

bool Interactor2DROIEdit::ProcessKeyUpEvent( wxKeyEvent& event, RenderView* renderview )
{
	UpdateCursor( event, renderview );
	
	return Interactor2D::ProcessKeyUpEvent( event, renderview );
}

void Interactor2DROIEdit::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
	if ( wnd->FindFocus() == wnd )
	{
		if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) )
		{ 
			wxMouseEvent* e = ( wxMouseEvent* )&event;
			if ( ( ( e->MiddleDown() || e->RightDown() ) && !m_bEditing ) ||
							  ( e->ControlDown() && e->ShiftDown() ) )
			{
				Interactor2D::UpdateCursor( event, wnd );
				return;
			}
		}
		
		if ( m_nAction == EM_Freehand || m_nAction == EM_Polyline )
		{
			if ( event.IsKindOf( CLASSINFO( wxKeyEvent ) ) )
			{	
				wxKeyEvent* e = ( wxKeyEvent* )&event;
				if ( e->GetEventType() != wxEVT_KEY_UP && ( e->GetKeyCode() == WXK_CONTROL && !e->ShiftDown() ) )
				{
					wnd->SetCursor( CursorFactory::CursorFill );
					return;
				}
			}
			
			if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) && (( wxMouseEvent* )&event)->ControlDown() )
			{
				wnd->SetCursor( CursorFactory::CursorFill );
			}
			else
				wnd->SetCursor( m_nAction == EM_Freehand ? CursorFactory::CursorPencil : CursorFactory::CursorPolyline );
		}
		else if ( m_nAction == EM_Fill )
			wnd->SetCursor( CursorFactory::CursorFill );
		else 
			Interactor2D::UpdateCursor( event, wnd );			
	}	
	else
		Interactor2D::UpdateCursor( event, wnd );
}
