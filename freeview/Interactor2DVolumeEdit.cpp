/**
 * @file  Interactor2DVolumeEdit.cpp
 * @brief Interactor for editing volume in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/09/04 18:50:45 $
 *    $Revision: 1.27 $
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
 *
 */

#include "Interactor2DVolumeEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "Contour2D.h"
#include "CursorFactory.h"
#include "BrushProperty.h"
#include <vtkRenderer.h>
#include <QDebug>

Interactor2DVolumeEdit::Interactor2DVolumeEdit( const QString& layerTypeName, QObject* parent ) :
  Interactor2D( parent ),
  m_bEditing( false ),
  m_bColorPicking( false )
{
  m_strLayerTypeName = layerTypeName;
}

Interactor2DVolumeEdit::~Interactor2DVolumeEdit()
{}

bool Interactor2DVolumeEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( !view->hasFocus() )
  {
    return Interactor2D::ProcessMouseDownEvent( event, renderview );
  }

  if ( event->button() == Qt::LeftButton ||
       ( event->button() == Qt::RightButton && (event->buttons() & Qt::LeftButton) ) )
  {
    if ( (event->modifiers() & CONTROL_MODIFIER ) && (event->modifiers() & Qt::ShiftModifier) )
    {
      return Interactor2D::ProcessMouseDownEvent( event, renderview );
    }

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    if ( (!mri || !mri->IsVisible()) ) //&& ( event->ControlDown() || m_nAction == EM_Polyline ) )
    {
      emit Error( "Layer Not Visible", mri );
    }
    else if ( !mri->IsEditable() ) //&& ( event->ControlDown() || m_nAction == EM_Polyline ) )
    {
      emit Error( "Layer Not Editable", mri );
    }
    else if ( m_strLayerTypeName == "MRI" && ((LayerMRI*)mri)->IsTransformed() )
    {
      emit Error( "Layer Not Editable For Transformation", mri );
    }
    else
    {
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();
      bool bFill3D = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetFill3D();

      double ras[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
      bool bCondition = !(event->modifiers() & Qt::ShiftModifier) && !(event->buttons() & Qt::RightButton);
      if ( (m_nAction == EM_ColorPicker || m_bColorPicking ) && mri->IsTypeOf( "MRI" ) )
      {
        if ( event->modifiers() & CONTROL_MODIFIER )
        {
          mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
          mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
        }
        else
        {
          double dValue = ((LayerMRI*)mri)->GetVoxelValue( ras );
          if ( dValue != 0 )
          {
            mri->SetFillValue( (float)dValue );
          }
        }
        m_bColorPicking = false;
      }
      else if ( m_nAction == EM_Freehand ) //&& ( event->ControlDown() ) )
      {
        if ( event->modifiers() & CONTROL_MODIFIER )
        {
          mri->SaveForUndo( bFill3D ? -1 : view->GetViewPlane());
          mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
        }
        else
        {
          mri->SaveForUndo( view->GetViewPlane() );
          m_bEditing = true;
          mri->SetVoxelByRAS( ras, view->GetViewPlane(),bCondition );
        }
      }
      else if ( m_nAction == EM_Fill ) //&& ( event->ControlDown() ) )
      {
        mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
        mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
      }
      else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
      {
        if ( event->modifiers() & CONTROL_MODIFIER )
        {
          mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
          mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
        }
        else
        {
          mri->SaveForUndo( view->GetViewPlane() );
          m_bEditing = true;
          double ras2[3];
          view->GetCursor2D()->ClearInterpolationPoints();
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition( ras );
          view->GetCursor2D()->SetPosition2( ras );
          if ( m_dPolylinePoints.size() > 0 )
          {
            if ( m_nAction == EM_Polyline )
            {
              mri->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), bCondition );
            }
            else
            {
              mri->SetLiveWireByRAS( ras, ras2, view->GetViewPlane() );
            }
          }
          else
          {
            // mri->SaveForUndo( view->GetViewPlane() );
            m_dPolylinePoints.push_back( ras[0] );
            m_dPolylinePoints.push_back( ras[1] );
            m_dPolylinePoints.push_back( ras[2] );
            view->GetCursor2D()->SetPosition( ras );
          }

          view->grabMouse();
        }
      }
      else if ( m_nAction == EM_Contour && mri->IsTypeOf( "MRI" ) )
      {
        LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetReferenceLayer();
        if ( !mri_ref )
        {
          emit Error( "Reference Layer Not Set" );
          return false;
        }

        Contour2D* c2d = view->GetContour2D();
        if ( (event->modifiers() & CONTROL_MODIFIER) && (event->modifiers() & Qt::AltModifier) )
        {
          double dValue = mri_ref->GetVoxelValue( ras );
          if ( dValue != 0 )
          {
            m_bEditing = true;
            c2d->SetInput( mri_ref->GetSliceImageData( view->GetViewPlane() ), dValue, ras[view->GetViewPlane()], mri_ref->GetActiveFrame() );
            c2d->SetVisible( true );
            view->RequestRedraw();
          }
          else if ( c2d->IsVisible() )
          {
            m_bEditing = true;
          }
        }
        else if ( (event->modifiers() & CONTROL_MODIFIER) && !(event->modifiers() & Qt::AltModifier) )
        {
          mri->SaveForUndo( view->GetViewPlane() );
          ((LayerMRI*)mri)->FloodFillByContour2D( ras, c2d );
        }
        else if ( event->modifiers() & Qt::ShiftModifier )
        {
          m_bEditing = true;
          c2d->RemoveLine( ras, ras );
          view->RequestRedraw();
        }
        else
        {
          m_bEditing = true;
          c2d->AddLine( ras, ras );
          view->RequestRedraw();
        }
      }
      else
      {
        return Interactor2D::ProcessMouseDownEvent( event, renderview );
      }
    }

    return false;
  }
  else if ( m_bEditing )
  {
    m_bEditing = false;
    if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    {
      if ( event->button() == Qt::MidButton )
      {
        view->GetCursor2D()->Update();
        view->RequestRedraw();
      }
      else if ( event->button() == Qt::RightButton )
      {
        if ( m_dPolylinePoints.size() > 0 && m_nAction == EM_Polyline )
        {
          LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
          LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();

          double ras1[3] = { m_dPolylinePoints[0], m_dPolylinePoints[1], m_dPolylinePoints[2] };
          double ras2[3];
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition2( ras2 );
          view->GetCursor2D()->SetPosition( ras1 );
          mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !(event->modifiers() & Qt::ShiftModifier) );
        }
        else
        {
          // mri->SetLiveWireByRAS( ras1, ras2, view->GetViewPlane() );
          view->GetCursor2D()->Update();
          view->RequestRedraw();
        }
      }
    }

    m_dPolylinePoints.clear();
    view->releaseMouse();

    return false;
  }
  return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DVolumeEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
// RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );

  if ( m_bEditing )
  {
    m_nMousePosX = event->x();
    m_nMousePosY = event->y();

    if ( event->button() != Qt::LeftButton ||
         (m_nAction != EM_Polyline && m_nAction != EM_Livewire ) ||
         m_dPolylinePoints.size() == 0 )
    {
      m_bEditing = false;
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2DVolumeEdit::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bEditing )
  {
    UpdateCursor( event, view );
    int posX = event->x();
    int posY = event->y();

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    if ( m_nAction == EM_Freehand )
    {
      double ras1[3], ras2[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
      view->MousePositionToRAS( posX, posY, ras2 );

      mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(),
                          !(event->modifiers() & Qt::ShiftModifier) && !(event->buttons() & Qt::RightButton) );
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
      view->RequestRedraw();
    }
    else if ( m_nAction == EM_Contour )
    {
      LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetReferenceLayer();
      Contour2D* c2d = view->GetContour2D();
      if ( event->modifiers() & Qt::ShiftModifier )
      {
        double ras1[3], ras2[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
        view->MousePositionToRAS( posX, posY, ras2 );
        c2d->RemoveLine( ras1, ras2 );
      }
      else if ( (event->modifiers() & CONTROL_MODIFIER) && (event->modifiers() & Qt::AltModifier) )
      {
        double scale = 0.2;
        if ( mri_ref )
        {
          double dMin = mri_ref->GetProperty()->GetMinValue();
          double dMax = mri_ref->GetProperty()->GetMaxValue();
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

      view->RequestRedraw();
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

bool Interactor2DVolumeEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  RenderView2D* view = ( RenderView2D* )renderview;
  if ( event->modifiers() & Qt::AltModifier && event->key() == Qt::Key_H )
  {
    Contour2D* c2d = view->GetContour2D();
    c2d->SetVisible( !c2d->IsVisible() );
    view->RequestRedraw();
    return false;
  }
  else if (event->modifiers() & Qt::ShiftModifier && event->key() == Qt::Key_C)
  {
    m_bColorPicking = true;
    return false;
  }
  else if (event->key() == Qt::Key_Escape)
  {
    m_bColorPicking = false;
    return false;
  }

  if ( !m_bEditing )
  {
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  }
  else
  {
    return false;
  }
}

bool Interactor2DVolumeEdit::ProcessKeyUpEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  return Interactor2D::ProcessKeyUpEvent( event, renderview );
}

void Interactor2DVolumeEdit::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( wnd->hasFocus() )
  {
    bool bMouseEvent = false;
    if ( event->type() == QEvent::MouseButtonPress ||
         event->type() == QEvent::MouseButtonRelease ||
         event->type() == QEvent::MouseMove)
    {
      QMouseEvent* e = ( QMouseEvent* )event;
      bMouseEvent = true;
      if ( ( ( e->button() == Qt::MidButton || e->button() == Qt::RightButton ) && !m_bEditing ) ||
           ( ( e->modifiers() & CONTROL_MODIFIER) && (e->modifiers() & Qt::ShiftModifier) ) )
      {
        Interactor2D::UpdateCursor( event, wnd );
        return;
      }
    }

    if ( m_nAction != EM_Fill )
    {
      if ( event->type() == QEvent::KeyPress )
      {
        QKeyEvent* e = ( QKeyEvent* )event;
        if ( e->key() == CONTROL_KEY && !(e->modifiers() & Qt::ShiftModifier) && !(e->modifiers() & Qt::AltModifier) )
        {
          wnd->setCursor( CursorFactory::CursorFill );
          return;
        }
      }

      if ( bMouseEvent && (( QMouseEvent* )event)->modifiers() & CONTROL_MODIFIER
           && !((( QMouseEvent* )event)->modifiers() & Qt::ShiftModifier)
           && !((( QMouseEvent* )event)->modifiers() & Qt::AltModifier) )
      {
        wnd->setCursor( CursorFactory::CursorFill );
      }
      else if ( m_nAction == EM_ColorPicker || m_bColorPicking )
      {
        wnd->setCursor( CursorFactory::CursorColorPicker );
      }
      else //if ( bMouseEvent )
      {
        switch ( m_nAction )
        {
        case EM_Freehand:
          wnd->setCursor( CursorFactory::CursorPencil );
          break;
        case EM_Polyline:
          wnd->setCursor( CursorFactory::CursorPolyline );
          break;
        case EM_Contour:
          wnd->setCursor( CursorFactory::CursorContour );
          break;
        }
      }
    }
    else
    {
      wnd->setCursor( CursorFactory::CursorFill );
    }
  }
  else
  {
    Interactor2D::UpdateCursor( event, wnd );
  }
}
