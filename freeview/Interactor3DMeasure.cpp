/**
 * @brief Interactor for measurement in 3D render view.
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


#include "Interactor3DMeasure.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerPropertyMRI.h"
#include "LayerMRI.h"
#include "SurfaceRegion.h"
#include "CursorFactory.h"
#include <QDebug>

Interactor3DMeasure::Interactor3DMeasure(QObject* parent) :
  Interactor3D(parent),
  m_bSelectRegion( false )
{}

Interactor3DMeasure::~Interactor3DMeasure()
{}

bool Interactor3DMeasure::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  bool ret = Interactor3D::ProcessMouseDownEvent( event, renderview );
#ifdef Q_OS_MAC
  if ( m_nAction == MM_SurfaceRegion && !Interactor3D::IsInAction() &&
       (event->button() == Qt::LeftButton || ( event->button() == Qt::RightButton && (event->buttons() & Qt::LeftButton) ) ) )
#else
  if ( m_nAction == MM_SurfaceRegion && !Interactor3D::IsInAction() &&
       event->button() == Qt::LeftButton )
#endif
  {
    if ( event->modifiers() & CONTROL_MODIFIER && !(event->modifiers() & Qt::ShiftModifier) )
    {
      if ( view->InitializeSelectRegion( event->x(), event->y() ) )
      {
        m_bSelectRegion = true;
        return false;   // do not pass down the event
      }
    }
    else if ( event->modifiers() & CONTROL_MODIFIER && event->modifiers() & Qt::ShiftModifier )
    {
      LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
      if ( mri && mri->GetCurrentSurfaceRegion() )
      {
        mri->GetCurrentSurfaceRegion()->ResetOutline();
        m_bSelectRegion = true;
        return false;
      }
    }
    else if ( !(event->modifiers() & CONTROL_MODIFIER) && event->modifiers() & Qt::ShiftModifier )
    {
      LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
      if ( mri && mri->GetCurrentSurfaceRegion() )
      {
        if ( mri->GetCurrentSurfaceRegion()->DeleteCell( view, event->x(), event->y() ) )
        {
          view->RequestRedraw();
        }
        return false;
      }
    }
  }

  return ret;
}

bool Interactor3DMeasure::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bSelectRegion )
  {
    view->CloseSelectRegion();
    m_bSelectRegion = false;
  }

  return Interactor3D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3DMeasure::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bSelectRegion )
  {
    view->AddSelectRegionLoopPoint( event->x(), event->y() );
    return false;
  }
  else
  {
    return Interactor3D::ProcessMouseMoveEvent( event, view );
  }
}

bool Interactor3DMeasure::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( event->key() == Qt::Key_Delete )
  {
    view->DeleteCurrentSelectRegion();
    return false;
  }
  else
  {
    return Interactor3D::ProcessKeyDownEvent( event, view );
  }
}

void Interactor3DMeasure::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( event->type() == QEvent::MouseButtonPress ||
       event->type() == QEvent::MouseButtonRelease ||
       event->type() == QEvent::MouseMove)
  {
    QMouseEvent* e = ( QMouseEvent* )event;
#ifdef Q_OS_MAC
    if ( (e->button() != Qt::RightButton && e->button() != Qt::MidButton) ||
         ( e->button() == Qt::RightButton && (e->buttons() & Qt::LeftButton) ) )
#else
    if (e->button() != Qt::RightButton && e->button() != Qt::MidButton)
#endif
    {
      if ( e->modifiers() & CONTROL_MODIFIER && !(e->modifiers() & Qt::ShiftModifier) )
      {
        wnd->setCursor( CursorFactory::CursorContour );
        return;
      }
      else if ( e->modifiers() & CONTROL_MODIFIER && e->modifiers() & Qt::ShiftModifier )
      {
        wnd->setCursor( CursorFactory::CursorAdd );
        return;
      }
      else if ( !(e->modifiers() & CONTROL_MODIFIER) && e->modifiers() & Qt::ShiftModifier )
      {
        wnd->setCursor( CursorFactory::CursorRemove );
        return;
      }
    }
  }
  else if ( event->type() == QEvent::KeyPress )
  {
    QKeyEvent* e = ( QKeyEvent* )event;
    if ( e->key() == CONTROL_KEY && !(e->modifiers() & Qt::ShiftModifier) )
    {
      wnd->setCursor( CursorFactory::CursorContour );
      return;
    }
    else if ( ( e->key() == CONTROL_KEY && e->modifiers() & Qt::ShiftModifier ) ||
              ( e->key() == Qt::Key_Shift && e->modifiers() & CONTROL_MODIFIER ) )
    {
      wnd->setCursor( CursorFactory::CursorAdd );
      return;
    }
    else if ( e->key() == Qt::Key_Shift && !(e->modifiers() & CONTROL_MODIFIER) )
    {
      wnd->setCursor( CursorFactory::CursorRemove );
      return;
    }
  }

  Interactor3D::UpdateCursor( event, wnd );
}
