/**
 * @file  Interactor2DRegionEdit.cpp
 * @brief Interactor2DRegionEdit to manage mouse and key input on render view.
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

#include "Interactor2DRegionEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

Interactor2DRegionEdit::Interactor2DRegionEdit() :
    Interactor2D(),
    m_bEditing( false )
{}

Interactor2DRegionEdit::~Interactor2DRegionEdit()
{}

bool Interactor2DRegionEdit::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
// UpdateCursor( event, view );

  if ( event.LeftDown() )
  {
    // if ( event.CmdDown() )
    //  return Interactor2D::ProcessMouseDownEvent( event, renderview );

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( "MRI" );
    LayerMRI* mri = ( LayerMRI* )lc->GetActiveLayer();
    if ( (!mri || !mri->IsVisible()) ) //&& ( event.CmdDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( "MRINotVisible", this );
    }
    else if ( !mri->IsEditable() ) //&& ( event.CmdDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( "MRINotEditable", this );
    }
    else
    {
      m_nMousePosX = event.GetX();
      m_nMousePosY = event.GetY();

      double ras[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
      if ( m_nAction == EM_Freehand ) //&& ( event.CmdDown() ) )
      {
        mri->SaveForUndo( view->GetViewPlane() );
        if ( event.CmdDown() )
        {
          mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
        }
        else
        {
          m_bEditing = true;
          mri->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
        }
      }
      else if ( m_nAction == EM_Fill ) //&& ( event.CmdDown() ) )
      {
        mri->SaveForUndo( view->GetViewPlane() );
        mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() );
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
          mri->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), !event.ShiftDown() );
        }
        else
        {
          mri->SaveForUndo( view->GetViewPlane() );
          m_dPolylinePoints.push_back( ras[0] );
          m_dPolylinePoints.push_back( ras[1] );
          m_dPolylinePoints.push_back( ras[2] );
        }

        if ( view->GetCapture() == view )
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
        view->GetCursor2D()->Update();
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
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition2( ras2 );
          view->GetCursor2D()->SetPosition( ras1 );
          mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !event.ShiftDown() );
        }
      }
    }

    m_dPolylinePoints.clear();
    if ( view->GetCapture() == view )
      view->ReleaseMouse();

    return false;
  }
  return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DRegionEdit::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
// RenderView2D* view = ( RenderView2D* )renderview;
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

bool Interactor2DRegionEdit::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bEditing )
  {
    UpdateCursor( event, view );
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

bool Interactor2DRegionEdit::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  if ( !m_bEditing )
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  else
    return false;
}

void Interactor2DRegionEdit::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  if ( wnd->FindFocus() == wnd )
  {
    if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) )
         && ( event.GetEventType() == wxEVT_MIDDLE_DOWN || event.GetEventType() == wxEVT_RIGHT_DOWN )
         && !m_bEditing )
    {
      Interactor2D::UpdateCursor( event, wnd );
      return;
    }

    if ( m_nAction == EM_Freehand || m_nAction == EM_Polyline )
    {
      if ( event.IsKindOf( CLASSINFO( wxKeyEvent ) ) && (( wxKeyEvent* )&event)->GetKeyCode() == WXK_CONTROL &&
           (( wxKeyEvent* )&event)->GetEventType() != wxEVT_KEY_UP )
      {
        wnd->SetCursor( CursorFactory::CursorFill );
      }
      else if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) && (( wxMouseEvent* )&event)->CmdDown() )
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
