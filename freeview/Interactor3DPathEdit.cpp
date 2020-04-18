#include "Interactor3DPathEdit.h"
#include <QMouseEvent>
#include "LayerSurface.h"
#include "MainWindow.h"
#include "CursorFactory.h"
#include "RenderView3D.h"
#include <QDebug>

Interactor3DPathEdit::Interactor3DPathEdit(QObject* parent) :
  Interactor3D(parent)
{

}

bool Interactor3DPathEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  m_bEditAttempted = false;
  if (event->button() == Qt::LeftButton && event->modifiers() != (Qt::ShiftModifier | Qt::ControlModifier))
  {
    m_bEditAttempted = true;
  }

  return Interactor3D::ProcessMouseDownEvent( event, renderview );
}

bool Interactor3DPathEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  if ( MainWindow::GetMainWindow()->IsEmpty() )
  {
    return Interactor3D::ProcessKeyDownEvent( event, renderview );
  }

  int nKeyCode = event->key();
  if ( nKeyCode == Qt::Key_Delete )
  {
    LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
    if (surf)
      surf->RemoveLastPathPoint();
  }
  else
  {
    return Interactor3D::ProcessKeyDownEvent( event, renderview );
  }

  return false;
}

bool Interactor3DPathEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
  if (m_bEditAttempted && qAbs(event->x() - this->m_nPressedPosX) < 2
      && qAbs(event->y() - this->m_nPressedPosY) < 2)
  {
    LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
    if (surf)
    {
      int nvo = view->PickCurrentSurfaceVertex(event->x(), event->y(), surf);
      if (nvo >= 0)
      {
        if (event->modifiers() == Qt::ShiftModifier)
        {
          surf->RemovePathPoint(nvo);
        }
        else
        {
          int nPath = -1;
          if (event->modifiers() & Qt::ControlModifier)
            nPath = surf->FindPathAt(nvo);
          if (nPath >= 0)
            surf->SetActivePath(nPath);
          else
            surf->AddPathPoint(nvo);
        }
        return true;
      }
    }
  }
  return Interactor3D::ProcessMouseUpEvent(event, renderview);
}
