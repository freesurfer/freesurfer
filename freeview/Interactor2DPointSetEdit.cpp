/**
 * @brief Interactor for editing way points in 2D render view.
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

#include "Interactor2DPointSetEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerPropertyMRI.h"
#include "LayerPointSet.h"
#include "LayerMRI.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

Interactor2DPointSetEdit::Interactor2DPointSetEdit(QObject* parent) :
  Interactor2D(parent),
  m_bEditing( false ),
  m_nCurrentIndex( -1 )
{}

Interactor2DPointSetEdit::~Interactor2DPointSetEdit()
{}

bool Interactor2DPointSetEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( event->button() == Qt::LeftButton )
  {
    if ( event->modifiers() & CONTROL_MODIFIER && event->modifiers() & Qt::ShiftModifier )
    {
      view->UpdateCursorRASPosition( event->x(), event->y());
      view->RequestRedraw();
      return false;
      //  return Interactor2D::ProcessMouseDownEvent( event, renderview );
    }

    if ( !(event->modifiers() & Qt::ShiftModifier) )
    {
      view->PickPointSetAtCursor(event->x(), event->y());
    }

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "PointSet" );
    LayerPointSet* wp = ( LayerPointSet* )lc->GetActiveLayer();
    if ( !wp || !wp->IsVisible() )
    {
      emit Error( "PointSetNotVisible", wp );
    }
    else
    {
      wp->SaveForUndo();
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();
      m_bEditing = true;

      double ras[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
      if ( !(event->modifiers() & Qt::ShiftModifier) )
      {
        m_nCurrentIndex = wp->FindPoint( ras );
        if ( m_nCurrentIndex < 0 )
        {
          m_nCurrentIndex = wp->AddPoint( ras, 1, false, property("insert_after").toInt()-1 );
          m_bEditing = false;
        }
      }
      else if (wp->RemovePoint( ras ))
      {
        m_nCurrentIndex = wp->GetNumberOfPoints()-1;
        m_bEditing = false;
      }
    }

    setProperty("insert_after", 0);
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseDownEvent( event, renderview );  // pass down the event
  }
}

bool Interactor2DPointSetEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  // RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );

  m_bEditing = false;
  //  m_nCurrentIndex = -1;

  return Interactor2D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor2DPointSetEdit::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  if ( m_bEditing )
  {
    UpdateCursor( event, view );
    int posX = event->x();
    int posY = event->y();
    m_nMousePosX = posX;
    m_nMousePosY = posY;
    if ( m_nCurrentIndex >= 0 )
    {
      double ras[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
      LayerPointSet* wp = ( LayerPointSet* )MainWindow::GetMainWindow()->GetActiveLayer( "PointSet" );
      wp->UpdatePoint( m_nCurrentIndex, ras );
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

bool Interactor2DPointSetEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  if ( !m_bEditing )
  {
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  }
  else
  {
    return false;
  }
}

void Interactor2DPointSetEdit::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( wnd->hasFocus() )
  {
    if ( event->type() == QEvent::MouseButtonPress ||
         event->type() == QEvent::MouseButtonRelease ||
         event->type() == QEvent::MouseMove)
    {
      QMouseEvent* e = ( QMouseEvent* )event;
      if ( e->button() == Qt::MiddleButton || e->button() == Qt::RightButton ||
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
