/**
 * @file  Interactor2DVolumeEdit.cpp
 * @brief Interactor for editing volume in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.19 $
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

#include "Interactor2DVolumeEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "Contour2D.h"
#include "CursorFactory.h"
#include "BrushProperty.h"
#include <vtkRenderer.h>

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
    if ( event.CmdDown() && event.ShiftDown() )
      return Interactor2D::ProcessMouseDownEvent( event, renderview );

    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->GetLayerCollection( m_strLayerTypeName.c_str() );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    if ( (!mri || !mri->IsVisible()) ) //&& ( event.CmdDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( m_strLayerTypeName + "NotVisible", this );
    }
    else if ( !mri->IsEditable() ) //&& ( event.CmdDown() || m_nAction == EM_Polyline ) )
    {
      SendBroadcast( m_strLayerTypeName + "NotEditable", this );
    }
    else if ( m_strLayerTypeName == "MRI" && ((LayerMRI*)mri)->IsTransformed() )
    {
      SendBroadcast( m_strLayerTypeName + "NotEditableForTransformation", this );
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
          mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
        else
        {
          m_bEditing = true;
          mri->SetVoxelByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
      }
      else if ( m_nAction == EM_Fill ) //&& ( event.CmdDown() ) )
      {
        mri->SaveForUndo( view->GetViewPlane() );
        mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
      }
      else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
      {
        mri->SaveForUndo( view->GetViewPlane() );
        if ( event.CmdDown() )
        {
          mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
        else
        {
          m_bEditing = true;
          double ras2[3];
          view->GetCursor2D()->ClearInterpolationPoints();
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition( ras );
          view->GetCursor2D()->SetPosition2( ras );
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
  
          if ( view->GetCapture() == view )
            view->ReleaseMouse();
          view->CaptureMouse();
        }
      }
      else if ( m_nAction == EM_ColorPicker && mri->IsTypeOf( "MRI" ) )
      {
        if ( event.CmdDown() )
        {
          mri->SaveForUndo( view->GetViewPlane() );
          mri->FloodFillByRAS( ras, view->GetViewPlane(), !event.ShiftDown() && !event.RightIsDown() );
        }
        else
        {
          double dValue = ((LayerMRI*)mri)->GetVoxelValue( ras );
          if ( dValue != 0 )
          {
            mri->SetFillValue( (float)dValue );
            mri->SendBroadcast( "LayerActorUpdated", mri );
          }
        }
      }
      else if ( m_nAction == EM_Contour && mri->IsTypeOf( "MRI" ) )
      {
        LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetBrushProperty()->GetReferenceLayer();
        if ( !mri_ref )
        {
          SendBroadcast( m_strLayerTypeName + "ReferenceNotSet", this );
          return false;
        }
          
        Contour2D* c2d = view->GetContour2D();
        if ( event.CmdDown() && event.AltDown() )
        {
          double dValue = mri_ref->GetVoxelValue( ras );
          if ( dValue != 0 )
          {
            m_bEditing = true;
            c2d->SetInput( mri_ref->GetSliceImageData( view->GetViewPlane() ), dValue, ras[view->GetViewPlane()], mri_ref->GetActiveFrame() );
            c2d->SetVisible( true );
            view->NeedRedraw();
          }
          else if ( c2d->IsVisible() )
          {
            m_bEditing = true;
          }
        }
        else if ( event.CmdDown() && !event.AltDown() )
        {
          mri->SaveForUndo( view->GetViewPlane() );
          ((LayerMRI*)mri)->FloodFillByContour2D( ras, c2d );
        }
        else if ( event.ShiftDown() )
        {
          m_bEditing = true;
          c2d->RemoveLine( ras, ras );
          view->NeedRedraw();
        }
        else
        {
          m_bEditing = true;
          c2d->AddLine( ras, ras );
          view->NeedRedraw();
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
    if ( view->GetCapture() == view )
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
    else if ( m_nAction == EM_Contour )
    {
      LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetBrushProperty()->GetReferenceLayer();
      Contour2D* c2d = view->GetContour2D();
      if ( event.ShiftDown() )
      {
        double ras1[3], ras2[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
        view->MousePositionToRAS( posX, posY, ras2 );
        c2d->RemoveLine( ras1, ras2 );
      }
      else if ( event.CmdDown() && event.AltDown() )
      {
        double scale = 0.2;
        if ( mri_ref )
        {
          double dMin = mri_ref->GetProperties()->GetMinValue();
          double dMax = mri_ref->GetProperties()->GetMaxValue();
          scale = ( dMax - dMin ) * 0.0005;
        }
        c2d->SetContourValue( c2d->GetContourValue() + scale * ( posY - m_nMousePosY ) ); 
      }
      else 
      {
        double ras1[3], ras2[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
        view->MousePositionToRAS( posX, posY, ras2 );
        c2d->AddLine( ras1, ras2 );
      }
      
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
  
  RenderView2D* view = ( RenderView2D* )renderview;
  if ( event.GetModifiers() == wxMOD_ALT && event.GetKeyCode() == 'H' )
  {
    Contour2D* c2d = view->GetContour2D();
    c2d->SetVisible( !c2d->IsVisible() );
    view->NeedRedraw();
    return false;
  }
  
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
           ( e->CmdDown() && e->ShiftDown() ) )
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
        if ( e->GetEventType() != wxEVT_KEY_UP && ( e->GetKeyCode() == WXK_CONTROL && !e->ShiftDown() && !e->AltDown() ) )
        {
          wnd->SetCursor( CursorFactory::CursorFill );
          return;
        }
      }

      if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) && (( wxMouseEvent* )&event)->CmdDown()
           && !(( wxMouseEvent* )&event)->ShiftDown() && !(( wxMouseEvent* )&event)->AltDown() )
      {
        wnd->SetCursor( CursorFactory::CursorFill );
      }
      else if ( m_nAction == EM_ColorPicker )
        wnd->SetCursor( CursorFactory::CursorColorPicker );
      else if ( event.IsKindOf( CLASSINFO( wxMouseEvent ) ) )
      {
        switch ( m_nAction )
        {
          case EM_Freehand:
            wnd->SetCursor( CursorFactory::CursorPencil );
            break;
          case EM_Polyline:
            wnd->SetCursor( CursorFactory::CursorPolyline );
            break;
          case EM_Contour:
            wnd->SetCursor( CursorFactory::CursorContour );
            break;
        }
      }
    }
    else 
      wnd->SetCursor( CursorFactory::CursorFill );
  }
  else
    Interactor2D::UpdateCursor( event, wnd );
}
