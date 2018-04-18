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
  if (event->button() == Qt::LeftButton)
  {
    m_bEditAttempted = true;
  }

  return Interactor3D::ProcessMouseDownEvent( event, renderview );
}

//bool Interactor3DPathEdit::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
//{
//  RenderView3D* view = ( RenderView3D* )renderview;
//  if (m_bEditing)
//  {
//    LayerROI* roi = (LayerROI*)MainWindow::GetMainWindow()->GetActiveLayer("ROI");
//    if (roi && roi->GetMappedSurface())
//    {
//      int nvo = view->PickCurrentSurfaceVertex(event->x(), event->y(), roi->GetMappedSurface());
//      if (nvo >= 0)
//      {
//        if (m_nPrevVertex >= 0 && m_nPrevVertex != nvo)
//        {
//          LayerSurface* surf = roi->GetMappedSurface();
//          if (surf)
//          {
//            QList<int> seeds;
//            seeds << m_nPrevVertex << nvo;
//            roi->EditVertex(surf->FindPath(seeds), !(event->modifiers() & Qt::ShiftModifier));
//          }
//        }
//        else
//          roi->EditVertex(nvo, !(event->modifiers() & Qt::ShiftModifier));
//        m_nPrevVertex = nvo;
//      }
//    }
//    renderview->setCursor( CursorFactory::CursorPencil );
//    return false;
//  }
//  bool ret = Interactor3D::ProcessMouseMoveEvent(event, renderview);
//  renderview->setCursor( CursorFactory::CursorPencil );
//  return ret;
//}

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
          int nPath = surf->FindPathAt(nvo);
          if (nPath >= 0)
            surf->SetActivePath(nPath);
          else
            surf->AddPathPoint(nvo);
        }
      }
    }
  }
  return Interactor3D::ProcessMouseUpEvent(event, renderview);
}
