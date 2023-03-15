/**
 * @brief Base Interactor class manage mouse and key input in render view.
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

#include "Interactor.h"
#include "RenderView.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>
#include <QDebug>

#ifdef Q_OS_MAC
Qt::KeyboardModifier Interactor::CONTROL_MODIFIER = Qt::ControlModifier;
Qt::Key Interactor::CONTROL_KEY = Qt::Key_Meta;
#else
Qt::KeyboardModifier Interactor::CONTROL_MODIFIER = Qt::ControlModifier;
Qt::Key Interactor::CONTROL_KEY = Qt::Key_Control;
#endif

Interactor::Interactor( QObject* parent ) : QObject( parent )
{
  m_nAction = 0;
}

Interactor::~Interactor()
{}

int Interactor::GetAction()
{
  return m_nAction;
}

void Interactor::SetAction( int nAction )
{
  m_nAction = nAction;
}

void Interactor::SetUseCommandControl(bool b)
{
#ifdef Q_OS_MAC
  if (true)
  {
    CONTROL_MODIFIER = Qt::ControlModifier;
    CONTROL_KEY = Qt::Key_Control;
  }
  else
  {
    CONTROL_MODIFIER = Qt::MetaModifier;
    CONTROL_KEY = Qt::Key_Meta;
  }
#else
  Q_UNUSED(b);
#endif
}

bool Interactor::ProcessMouseDownEvent( QMouseEvent* event, RenderView* view )
{
  if ( event->button() == Qt::RightButton && event->modifiers() == Qt::NoModifier) //!(event->modifiers() | Qt::AltModifier) )
  {
    m_nDownPosX = event->x();
    m_nDownPosY = event->y();
  }
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessMouseUpEvent( QMouseEvent* event, RenderView* view )
{
  view->releaseMouse();
  if ( event->button() == Qt::RightButton && m_nDownPosX == event->x() && m_nDownPosY == event->y() )
  {
    view->TriggerContextMenu( event );
    return true;
  }
  else
  {
    UpdateCursor( event, view );
    return true;
  }

  return true;
}

bool Interactor::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessKeyDownEvent( QKeyEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessKeyUpEvent( QKeyEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessMouseWheelEvent( QWheelEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessMouseEnterEvent( QEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

bool Interactor::ProcessMouseLeaveEvent( QEvent* event, RenderView* view )
{
  UpdateCursor( event, view );

  return true;
}

void Interactor::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( (event->type() == QEvent::MouseButtonPress && ((QMouseEvent*)event)->button() == Qt::MiddleButton) ||
       (event->type() == QEvent::MouseMove && ((QMouseEvent*)event)->buttons() & Qt::MiddleButton ))
  {
    wnd->setCursor( CursorFactory::CursorPan );
  }
  else if ((event->type() == QEvent::MouseButtonPress && ((QMouseEvent*)event)->button() == Qt::RightButton) ||
           (event->type() == QEvent::MouseMove && ((QMouseEvent*)event)->buttons() & Qt::RightButton ))
  {
    wnd->setCursor( CursorFactory::CursorZoom );
  }
  else if ( event->type() == QEvent::Wheel )
  {
    wnd->setCursor( CursorFactory::CursorZoom );
  }
  else
  {
    wnd->unsetCursor();
  }
}
