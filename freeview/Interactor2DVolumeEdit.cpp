/**
 * @file  Interactor2DVolumeEdit.cpp
 * @brief Interactor2DVolumeEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/01/27 18:27:25 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008-2009,
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

#include "Interactor2DVolumeEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerVolumeBase.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

#define FV_CURSOR_PENCIL wxCursor( wxImage("res/images/cursor_pencil.png") )
#define FV_CURSOR_FILL  wxCursor( wxImage("res/images/cursor_fill.png" ) )

Interactor2DVolumeEdit::Interactor2DVolumeEdit( const char* layerTypeName) :
    Interactor2D(),
    m_bEditing( false )
{
  m_strLayerTypeName = layerTypeName;
}

Interactor2DVolumeEdit::~Interactor2DVolumeEdit()
{}

bool Interactor2DVolumeEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
// UpdateCursor( event, view );

  if ( event.LeftDown() || ( event.RightDown() && event.LeftIsDown() ) )
  {
    if ( event.ControlDown() && event.ShiftDown() )
      return Interactor2D::ProcessMouseDownEvent( event, renderview );

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( m_strLayerTypeName.c_str() );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    if ( (!mri || !mri->IsVisible()) ) //&& ( event.ControlDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( m_strLayerTypeName + "NotVisible", this );
    }
    else if ( !mri->IsEditable() ) //&& ( event.ControlDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( m_strLayerTypeName + "NotEditable", this );
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
        if ( event.ControlDown() )
        {
          mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
        else
        {
          m_bEditing = true;
          mri->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
      }
      else if ( m_nAction == EM_Fill ) //&& ( event.ControlDown() ) )
      {
        mri->SaveForUndo( view->GetViewPlane() );
        mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
      }
      else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
      {
        m_bEditing = true;
        double ras2[3];
        view->GetCursor2D()->ClearInterpolationPoints();
        view->GetCursor2D()->GetPosition( ras2 );
        view->GetCursor2D()->SetPosition( ras );
        view->GetCursor2D()->SetPosition2( ras );
        mri->SaveForUndo( view->GetViewPlane() );
        if ( m_dPolylinePoints.size() > 0 )
        {
          if ( m_nAction == EM_Polyline )
            mri->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
          else
            mri->SetLiveWireByRAS( ras, ras2, view->GetViewPlane() );
        }
        else
        {
          // mri->SaveForUndo( view->GetViewPlane() );
          m_dPolylinePoints.push_back( ras[0] );
          m_dPolylinePoints.push_back( ras[1] );
          m_dPolylinePoints.push_back( ras[2] );
          view->GetCursor2D()->SetPosition( ras );
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
    if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    {
      if ( event.MiddleDown() )
      {
        view->GetCursor2D()->Update();
        view->NeedRedraw();
      }
      else if ( event.RightDown() )
      {
        if ( m_dPolylinePoints.size() > 0 && m_nAction == EM_Polyline )
        {
          LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( m_strLayerTypeName.c_str() );
          LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();

          double ras1[3] = { m_dPolylinePoints[0], m_dPolylinePoints[1], m_dPolylinePoints[2] };
          double ras2[3];
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition2( ras2 );
          view->GetCursor2D()->SetPosition( ras1 );
          mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() );
        }
        else
        {
          // mri->SetLiveWireByRAS( ras1, ras2, view->GetViewPlane() );
          view->GetCursor2D()->Update();
          view->NeedRedraw();
        }
      }
    }

    m_dPolylinePoints.clear();
    view->ReleaseMouse();

    return false;
  }
  return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DVolumeEdit::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
// RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );

  if ( m_bEditing )
  {
    m_nMousePosX = event.GetX();
    m_nMousePosY = event.GetY();

    if ( !event.LeftUp() || (m_nAction != EM_Polyline && m_nAction != EM_Livewire )|| m_dPolylinePoints.size() == 0 )
      m_bEditing = false;

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( m_strLayerTypeName.c_str() );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    mri->SendBroadcast( "LayerEdited", mri );

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2DVolumeEdit::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bEditing )
  {
    UpdateCursor( event, view );
    int posX = event.GetX();
    int posY = event.GetY();

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( m_strLayerTypeName.c_str() );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();

    if ( m_nAction == EM_Freehand )
    {
      double ras1[3], ras2[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
      view->MousePositionToRAS( posX, posY, ras2 );

      mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
    }
    else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    {
      double ras[3];
      view->MousePositionToRAS( posX, posY, ras );
      view->GetCursor2D()->SetPosition2( ras );
      if ( m_nAction == EM_Livewire )
      {
        view->GetCursor2D()->SetInterpolationPoints(
          mri->GetLiveWirePointsByRAS( ras,
                                       view->GetCursor2D()->GetPosition(),
                                       view->GetViewPlane() ) );
      }
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

bool Interactor2DVolumeEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  if ( !m_bEditing )
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  else
    return false;
}

bool Interactor2DVolumeEdit::ProcessKeyUpEvent( wxKeyEvent& event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  return Interactor2D::ProcessKeyUpEvent( event, renderview );
}

void Interactor2DVolumeEdit::UpdateCursor( wxEvent& event, wxWindow* wnd )
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

    if ( m_nAction != EM_Fill )
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
