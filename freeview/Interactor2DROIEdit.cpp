/**
 * @file  Interactor2DROIEdit.cpp
 * @brief Interactor2DROIEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/04/18 19:58:18 $
 *    $Revision: 1.3 $
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
#include <vtkRenderer.h>

Interactor2DROIEdit::Interactor2DROIEdit() : Interactor2D()
{
	m_bEditing = false;
}

Interactor2DROIEdit::~Interactor2DROIEdit()
{
}

bool Interactor2DROIEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
	RenderView2D* view = ( RenderView2D* )renderview;
		
	if ( event.LeftDown() )
	{
		LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( "ROI" );
		LayerROI* roi = ( LayerROI* )lc->GetActiveLayer();
		if ( !roi || !roi->IsVisible() )
		{
			SendBroadcast( "ROINotEditable", this );
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
				m_bEditing = true;				
				roi->SetVoxelByRAS( ras, !event.ShiftDown() );
			}
			else if ( m_nAction == EM_Fill ) // && ( event.ControlDown() || event.ShiftDown() ) )
			{
				roi->SaveForUndo( view->GetViewPlane() );
				roi->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
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
					roi->SetVoxelByRAS( ras, ras2, !event.ShiftDown() );	
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
				view->GetCursor()->Update();
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
					view->GetCursor()->GetPosition( ras2 );
					view->GetCursor()->SetPosition2( ras2 );
					view->GetCursor()->SetPosition( ras1 );
					roi->SetVoxelByRAS( ras1, ras2, !event.ShiftDown() );	
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
		int posX = event.GetX();
		int posY = event.GetY();
		
		if ( m_nAction == EM_Freehand )
		{
			LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
			LayerROI* roi = ( LayerROI* )lc->GetActiveLayer();
					
			double ras1[3], ras2[3];
			view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
			view->MousePositionToRAS( posX, posY, ras2 );
			
			roi->SetVoxelByRAS( ras1, ras2, !event.ShiftDown() );
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

bool Interactor2DROIEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
	if ( !m_bEditing )
		return Interactor2D::ProcessKeyDownEvent( event, renderview );
	else
		return false;
}
