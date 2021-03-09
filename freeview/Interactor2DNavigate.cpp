/**
 * @brief Interactor for navigating (zoom, pan, etc.) in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#include "Interactor2DNavigate.h"
#include "MainWindow.h"
#include "RenderView2D.h"
#include "LayerLandmarks.h"

Interactor2DNavigate::Interactor2DNavigate(QObject* parent) :
  Interactor2D(parent),
  m_bEditing( false ),
  m_nCurrentLandmark( -1 )
{}

bool Interactor2DNavigate::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  if (m_nCurrentLandmark < 0)
    return Interactor2D::ProcessMouseDownEvent(event, renderview);

  RenderView2D* view = ( RenderView2D* )renderview;

  if ( event->button() == Qt::LeftButton )
  {
    if ( event->modifiers() & CONTROL_MODIFIER && event->modifiers() & Qt::ShiftModifier )
    {
      return Interactor2D::ProcessMouseDownEvent( event, renderview );
    }

    LayerLandmarks* landmarks = (LayerLandmarks*)MainWindow::GetMainWindow()->GetSupplementLayer("Landmarks");
    if ( !landmarks )
    {
      emit Error( "Landmarks non-exist", landmarks );
    }
    else
    {
      m_bEditing = true;
      double pos[3];
      view->MousePositionToRAS( event->x(), event->y(), pos );
      landmarks->SetLandmarkPosition(m_nCurrentLandmark, pos);
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseDownEvent( event, renderview );  // pass down the event
  }
}

bool Interactor2DNavigate::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  m_bEditing = false;

  return Interactor2D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor2DNavigate::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  if ( m_bEditing && m_nCurrentLandmark >= 0 )
  {
    UpdateCursor( event, view );
    LayerLandmarks* landmarks = (LayerLandmarks*)MainWindow::GetMainWindow()->GetSupplementLayer("Landmarks");
    if ( landmarks )
    {
      double pos[3];
      view->MousePositionToRAS( event->x(), event->y(), pos );
      landmarks->SetLandmarkPosition(m_nCurrentLandmark, pos);
    }
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

bool Interactor2DNavigate::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  return Interactor2D::ProcessKeyDownEvent( event, renderview );
}

void Interactor2DNavigate::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( wnd->hasFocus() )
  {
    if ( event->type() == QEvent::MouseButtonPress ||
         event->type() == QEvent::MouseButtonRelease ||
         event->type() == QEvent::MouseMove)
    {
      QMouseEvent* e = ( QMouseEvent* )event;
      if ( e->button() == Qt::MidButton || e->button() == Qt::RightButton ||
           ( ( e->modifiers() & CONTROL_MODIFIER) && (e->modifiers() & Qt::ShiftModifier) ) )
      {
        Interactor2D::UpdateCursor( event, wnd );
      }
      else
      {
        // set own cursor
        wnd->setCursor( QCursor() );
      }
    }
  }
  else
  {
    Interactor2D::UpdateCursor( event, wnd );
  }
}

void Interactor2DNavigate::SetCurrentLandmark(int n)
{
  m_nCurrentLandmark = n;
}
