/**
 * @brief Interactor for volume cropping in 2D render view.
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

#include "Interactor2DVolumeCrop.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "VolumeCropper.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>

Interactor2DVolumeCrop::Interactor2DVolumeCrop(QObject* parent) :
  Interactor2D(parent),
  m_bSelected( false )
{
}

Interactor2DVolumeCrop::~Interactor2DVolumeCrop()
{}

bool Interactor2DVolumeCrop::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( event->button() == Qt::LeftButton &&
       MainWindow::GetMainWindow()->GetVolumeCropper()
       ->PickActiveBound2D( view, event->x(), event->y() ) )
  {
    UpdateCursor( event, view );
    m_bSelected = true;
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseDownEvent( event, renderview );  // pass down the event
  }
}

bool Interactor2DVolumeCrop::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  if ( m_bSelected )
  {
    MainWindow::GetMainWindow()->GetVolumeCropper()->ReleaseActiveBound();
    MainWindow::GetMainWindow()->GetRenderView( 3 )->RequestRedraw();
    m_bSelected = false;
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2DVolumeCrop::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bSelected )
  {
    UpdateCursor( event, view );

    MainWindow::GetMainWindow()->GetVolumeCropper()
        ->MoveActiveBound2D( view, event->x(), event->y() );

    return false;
  }
  else if ( !(event->buttons() & Qt::LeftButton) && !(event->buttons() & Qt::RightButton) &&
            !(event->buttons() & Qt::MiddleButton) &&
            MainWindow::GetMainWindow()->GetVolumeCropper()
            ->PickActiveBound2D( view, event->x(), event->y() ) )
  {
    UpdateCursor( event, view );
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

void Interactor2DVolumeCrop::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( MainWindow::GetMainWindow()->GetVolumeCropper()->GetActivePlane() >= 0 )
  {
    wnd->setCursor( CursorFactory::CursorGrab );
  }
  else
  {
    Interactor2D::UpdateCursor( event, wnd );
  }
}
