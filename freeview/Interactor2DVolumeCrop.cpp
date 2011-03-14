/**
 * @file  Interactor2DVolumeCrop.cpp
 * @brief Interactor for volume cropping in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:57 $
 *    $Revision: 1.3 $
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
    return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
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
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
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
            !(event->buttons() & Qt::MidButton) &&
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
    wnd->setCursor( CursorFactory::CursorGrab );
  else
    Interactor2D::UpdateCursor( event, wnd );
}
