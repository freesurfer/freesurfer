/**
 * @file  Interactor.cpp
 * @brief Base Interactor class manage mouse and key input in render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
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
{}

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
  {  
    view->TriggerContextMenu( event.GetPosition() );
    return true;
  }
  else
  {
    UpdateCursor( event, view );
    return true;
  }
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
