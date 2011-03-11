/**
 * @file  Interactor3DMeasure.cpp
 * @brief Interactor for measurement in 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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


#include "Interactor3DMeasure.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include "SurfaceRegion.h"
#include "CursorFactory.h"

Interactor3DMeasure::Interactor3DMeasure() :
    Interactor3D(),
    m_bSelectRegion( false )
{}

Interactor3DMeasure::~Interactor3DMeasure()
{}

bool Interactor3DMeasure::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  bool ret = Interactor3D::ProcessMouseDownEvent( event, renderview );

  if ( m_nAction == MM_SurfaceRegion && !Interactor3D::IsInAction() && event.LeftDown() )
  {
    if ( event.CmdDown() && !event.ShiftDown() )
    {
      if ( view->InitializeSelectRegion( event.GetX(), event.GetY() ) )
      {
        m_bSelectRegion = true;
        return false;   // do not pass down the event
      }
    }
    else if ( event.CmdDown() && event.ShiftDown() )
    {
      LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
      if ( mri && mri->GetCurrentSurfaceRegion() )
      {
        mri->GetCurrentSurfaceRegion()->ResetOutline();
        m_bSelectRegion = true;
        return false;
      }
    }
    else if ( !event.CmdDown() && event.ShiftDown() )
    {
      LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetActiveLayer( "MRI" );
      if ( mri && mri->GetCurrentSurfaceRegion() )
      {
        if ( mri->GetCurrentSurfaceRegion()->DeleteCell( view, event.GetX(), event.GetY() ) )
          view->NeedRedraw();
        return false;
      }
    }
  }
   
  return ret;
}

bool Interactor3DMeasure::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bSelectRegion )
  {
     view->CloseSelectRegion();
     m_bSelectRegion = false;
  }

  return Interactor3D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3DMeasure::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bSelectRegion )
  {
    view->AddSelectRegionLoopPoint( event.GetX(), event.GetY() );
    return false;
  }
  else
  {
    return Interactor3D::ProcessMouseMoveEvent( event, view );
  }
}

bool Interactor3DMeasure::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
 
  int nKeyCode = event.GetKeyCode();
  if ( nKeyCode == WXK_DELETE )
  {
    view->DeleteCurrentSelectRegion();
    return false;
  }
  else
    return Interactor3D::ProcessKeyDownEvent( event, view );
}

void Interactor3DMeasure::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) )
  {
    wxMouseEvent* e = ( wxMouseEvent* )&event;
    if ( !e->RightDown() && !e->MiddleDown() )
    {
      if ( e->CmdDown() && !e->ShiftDown() )
      {
        wnd->SetCursor( CursorFactory::CursorContour );
        return;
      }
      else if ( e->CmdDown() && e->ShiftDown() )
      {
        wnd->SetCursor( CursorFactory::CursorAdd );
        return;
      }
      else if ( !e->CmdDown() && e->ShiftDown() )
      {
        wnd->SetCursor( CursorFactory::CursorRemove );
        return;
      }
    }
  }
  else if ( event.IsKindOf( CLASSINFO( wxKeyEvent ) ) )
  {
    wxKeyEvent* e = ( wxKeyEvent* )&event;
    if ( e->GetEventType() != wxEVT_KEY_UP )
    {
      if ( e->GetKeyCode() == WXK_CONTROL && !e->ShiftDown() )
      {
        wnd->SetCursor( CursorFactory::CursorContour );
        return;
      }
      else if ( ( e->GetKeyCode() == WXK_CONTROL && e->ShiftDown() ) ||
                  ( e->GetKeyCode() == WXK_SHIFT && e->CmdDown() ) )
      {
        wnd->SetCursor( CursorFactory::CursorAdd );
        return;
      }
      else if ( e->GetKeyCode() == WXK_SHIFT && !e->CmdDown() )
      {
        wnd->SetCursor( CursorFactory::CursorRemove );
        return;
      }
    }
  }
  
  Interactor3D::UpdateCursor( event, wnd );
}
