/**
 * @file  RenderView.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:15 $
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
 
#include "RenderView.h"
#include <vtkRenderer.h>
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "Interactor.h"
#include "MainWindow.h"

IMPLEMENT_DYNAMIC_CLASS(RenderView, wxVTKRenderWindowInteractor)

BEGIN_EVENT_TABLE(RenderView, wxVTKRenderWindowInteractor)
	EVT_SET_FOCUS	( RenderView::OnFocus )	
	EVT_LEFT_DOWN	( RenderView::OnButtonDown )	
	EVT_MIDDLE_DOWN	( RenderView::OnButtonDown )
	EVT_RIGHT_DOWN	( RenderView::OnButtonDown )
	EVT_LEFT_UP		( RenderView::OnButtonUp )
	EVT_MIDDLE_UP	( RenderView::OnButtonUp )
	EVT_RIGHT_UP	( RenderView::OnButtonUp )
	EVT_MOTION		( RenderView::OnMouseMove )
	EVT_MOUSEWHEEL	( RenderView::OnMouseWheel )
	EVT_ENTER_WINDOW( RenderView::OnMouseEnter )
	EVT_LEAVE_WINDOW( RenderView::OnMouseLeave )
	EVT_KEY_DOWN	( RenderView::OnKeyDown )
	EVT_SIZE        ( RenderView::OnSize )
END_EVENT_TABLE()

RenderView::RenderView() : wxVTKRenderWindowInteractor(), 
			Listener( "RenderView" ), 
			Broadcaster( "RenderView" )
{
	m_renderWindow = this->GetRenderWindow();
	m_renderer = vtkRenderer::New();
	m_renderWindow->AddRenderer( m_renderer );
	m_renderer->SetBackground( 0, 0, 0 );

	m_interactor = new Interactor();
	m_nInteractionMode = 0;
	m_nRedrawCount = 0;
	
	UseCaptureMouseOn();
}

RenderView::RenderView( wxWindow* parent, int id ) : 
		wxVTKRenderWindowInteractor( parent, id, wxDefaultPosition,
									 wxDefaultSize,
									 wxWANTS_CHARS | wxNO_FULL_REPAINT_ON_RESIZE ),
		Listener( "RenderView" ), 
		Broadcaster( "RenderView" )
{
	m_renderWindow = this->GetRenderWindow();
	m_renderer = vtkRenderer::New();
	m_renderWindow->AddRenderer( m_renderer );
//	m_renderer->SetBackground( 0.7, 0.7, 0.9 );

	m_interactor = new Interactor();
	m_nInteractionMode = 0;
	m_nRedrawCount = 0;
	
	UseCaptureMouseOn();
}

RenderView* RenderView::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new RenderView;
}

RenderView::~RenderView()
{
	if (m_renderer)
		m_renderer->Delete();
	
	delete m_interactor;
}

void RenderView::OnFocus( wxFocusEvent& event )
{
	Render();
	event.Skip();
}

void RenderView::OnButtonDown( wxMouseEvent& event )
{
	this->SetFocus();

	if ( m_interactor->ProcessMouseDownEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnButtonDown( event );
}

void RenderView::OnButtonUp( wxMouseEvent& event )
{
	if ( m_interactor->ProcessMouseUpEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnButtonUp( event );	
}

void RenderView::OnMouseMove( wxMouseEvent& event )
{
	if ( m_interactor->ProcessMouseMoveEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnMotion( event );
	
	m_interactor->ProcessPostMouseMoveEvent( event, this );
}

void RenderView::OnMouseWheel( wxMouseEvent& event )
{
	if ( m_interactor->ProcessMouseWheelEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnMouseWheel( event );
	
	m_interactor->ProcessPostMouseWheelEvent( event, this );
}

void RenderView::OnMouseEnter( wxMouseEvent& event )
{
	if ( m_interactor->ProcessMouseEnterEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnEnter( event );
}

void RenderView::OnMouseLeave( wxMouseEvent& event )
{
	if ( m_interactor->ProcessMouseLeaveEvent( event, this ) )
		wxVTKRenderWindowInteractor::OnLeave( event );
}

void RenderView::OnKeyDown( wxKeyEvent& event )
{
	if ( m_interactor->ProcessKeyDownEvent( event, this ) )
		event.Skip();
}

void RenderView::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

void RenderView::SetWorldCoordinateInfo( const double* origin, const double* size )
{
	for ( int i = 0; i < 3; i++ )
	{
		m_dWorldOrigin[i] = origin[i];
		m_dWorldSize[i] = size[i];
	}
	UpdateViewByWorldCoordinate();
}

void RenderView::DoListenToMessage ( std::string const iMsg, void* const iData )
{
	if ( iMsg == "LayerActorUpdated" )
		Render();
	
	else if ( iMsg == "LayerAdded" || iMsg == "LayerMoved" || iMsg == "LayerRemoved" )
		RefreshAllActors();
}

int	RenderView::GetInteractionMode()
{
	return m_nInteractionMode;
}

void RenderView::SetInteractionMode( int nMode )
{
	m_nInteractionMode = nMode;
}
	
int RenderView::GetAction()
{
	return m_interactor->GetAction();
}
	
void RenderView::SetAction( int nAction )
{
	m_interactor->SetAction( nAction );
}

void RenderView::OnSize( wxSizeEvent& event )
{
	wxVTKRenderWindowInteractor::OnSize( event );
	
#ifdef __WXGTK__
	MainWindow::GetMainWindowPointer()->NeedRedraw();
#endif
}

void RenderView::SetBackgroundColor( const wxColour& color )
{
	m_renderer->SetBackground( (double)color.Red()/255.0, (double)color.Green()/255.0, (double)color.Blue()/255.0 );
}

wxColour RenderView::GetBackgroundColor() const
{
	double* rgb = m_renderer->GetBackground();
	return wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) );
}

void RenderView::NeedRedraw()
{
	m_nRedrawCount = 1;
}

void RenderView::OnInternalIdle()
{
	wxWindow::OnInternalIdle();
	
	if ( IsShown() && m_nRedrawCount > 0 )
	{
		Render();
		m_nRedrawCount--;
	}
}
