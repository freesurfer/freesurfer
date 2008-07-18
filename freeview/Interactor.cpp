/**
 * @file  Interactor.cpp
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/07/18 20:23:24 $
 *    $Revision: 1.4 $
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

#include "Interactor.h"
#include "RenderView.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

Interactor::Interactor() : Broadcaster( "Interactor" )
{
	m_nAction = 0;
}

Interactor::~Interactor()
{
}

int Interactor::GetAction()
{
	return m_nAction;
}

void Interactor::SetAction( int nAction )
{
	m_nAction = nAction;
}

bool Interactor::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view )
{
	if ( event.RightDown() && !event.AltDown() )
	{
		m_nDownPosX = event.GetX();
		m_nDownPosY = event.GetY();
	}
	UpdateCursor( event, view );
	
	return true;
}

bool Interactor::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view )
{
	if ( event.RightUp() && m_nDownPosX == event.GetX() && m_nDownPosY == event.GetY() )
		view->TriggerContextMenu( event.GetPosition() );
	
	UpdateCursor( event, view );
	
	return true;
}

bool Interactor::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
}

bool Interactor::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
}

bool Interactor::ProcessKeyUpEvent( wxKeyEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
}

bool Interactor::ProcessMouseWheelEvent( wxMouseEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
}
	
bool Interactor::ProcessMouseEnterEvent( wxMouseEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
}
	
bool Interactor::ProcessMouseLeaveEvent( wxMouseEvent& event, RenderView* view )
{
	UpdateCursor( event, view );
	
	return true;
} 

void Interactor::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
	if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) && (( wxMouseEvent* )&event)->GetEventType() == wxEVT_MIDDLE_DOWN )
	{
		wnd->SetCursor( CursorFactory::CursorPan );
	}
	else if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) && (( wxMouseEvent* )&event)->GetEventType() == wxEVT_RIGHT_DOWN )
	{
		wnd->SetCursor( CursorFactory::CursorZoom );
	}
	else if ( event.GetEventType() == wxEVT_MOUSEWHEEL )
	{
		wnd->SetCursor( CursorFactory::CursorZoom );
	}
	else
		wnd->SetCursor( wxNullCursor );
}
