#include "Interactor3DPointSetEdit.h"
#include <QMouseEvent>
#include "LayerSurface.h"
#include "LayerPointSet.h"
#include "LayerMRI.h"
#include "MainWindow.h"
#include "CursorFactory.h"
#include "RenderView3D.h"
#include "FSSurface.h"
#include <QDebug>

Interactor3DPointSetEdit::Interactor3DPointSetEdit(QObject* parent) :
  Interactor3D(parent)
{

}

bool Interactor3DPointSetEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  m_bEditAttempted = false;
  if (event->button() == Qt::LeftButton && event->modifiers() != Qt::ControlModifier)
  {
    m_bEditAttempted = true;
  }

  return Interactor3D::ProcessMouseDownEvent( event, renderview );
}

bool Interactor3DPointSetEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  if ( MainWindow::GetMainWindow()->IsEmpty() )
  {
    return Interactor3D::ProcessKeyDownEvent( event, renderview );
  }

  return Interactor3D::ProcessKeyDownEvent( event, renderview );
}

bool Interactor3DPointSetEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
  MainWindow* mwnd = MainWindow::GetMainWindow();
  if (m_bEditAttempted && qAbs(event->x() - this->m_nPressedPosX) < 2
      && qAbs(event->y() - this->m_nPressedPosY) < 2)
  {
    LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
    LayerPointSet* wp = ( LayerPointSet* )MainWindow::GetMainWindow()->GetActiveLayer("PointSet");
    if ( !wp || !wp->IsVisible() )
    {
//      emit Error( "PointSetNotVisible", wp );
      return false;
    }

    if (wp)
    {
      if ( !(event->modifiers() & Qt::ShiftModifier) )
      {
        if (!mwnd->GetVisibleLayers("Surf").isEmpty())
        {
          int nvo = -1;
          nvo = view->PickCurrentSurfaceVertex(event->x(), event->y(), surf);
          if (nvo >= 0)
          {
            wp->SaveForUndo();
            double ras[3];
            surf->GetTargetAtVertex(nvo, ras, surf->IsInflated()?FSSurface::SurfaceWhite:-1);
            int nIndex = wp->FindPoint( ras );
            if ( nIndex < 0 )
            {
              nIndex = wp->AddPoint( ras, 1, true );
              return true;
            }
          }
        }
//        else if (!mwnd->GetVisibleLayers("MRI").isEmpty())
//        {
//          double ras[3];
//          if (view->PickVolume(event->x(), event->y(), ras))
//          {
//            int nIndex = wp->FindPoint( ras );
//            if ( nIndex < 0 )
//            {
//              nIndex = wp->AddPoint( ras, 1, true );
//              return true;
//            }
//          }
//        }
      }
      else
      {
        int nPt = view->PickCurrentPointSetPoint(event->x(), event->y());
        if (nPt >= 0)
        {
          wp->SaveForUndo();
          wp->RemovePoint(nPt);
          return true;
        }
      }
    }
  }
  return Interactor3D::ProcessMouseUpEvent(event, renderview);
}
