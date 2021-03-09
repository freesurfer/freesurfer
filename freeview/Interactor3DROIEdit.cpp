/**
 * @brief Interactor for navigating in 3D render view.
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


#include "Interactor3DROIEdit.h"
#include "RenderView3D.h"
#include <QMouseEvent>
#include "LayerROI.h"
#include "LayerSurface.h"
#include "MainWindow.h"
#include "CursorFactory.h"

Interactor3DROIEdit::Interactor3DROIEdit(QObject* parent) :
  Interactor3D(parent), m_bEditing(false), m_nPrevVertex(-1)
{

}

bool Interactor3DROIEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
  if (event->button() == Qt::LeftButton && !(event->modifiers() & CONTROL_MODIFIER))
  {
    m_bEditing = true;
    m_nPrevVertex = -1;
    LayerROI* roi = (LayerROI*)MainWindow::GetMainWindow()->GetActiveLayer("ROI");
    if (roi && roi->GetMappedSurface())
    {
      roi->SaveForUndo();
      int nvo = view->PickCurrentSurfaceVertex(event->x(), event->y(), roi->GetMappedSurface());
      if (nvo >= 0)
      {
        m_nPrevVertex = nvo;
        roi->EditVertex(nvo, !(event->modifiers() & Qt::ShiftModifier));
        return false;
      }
    }
  }

  return Interactor3D::ProcessMouseDownEvent( event, renderview );
}

bool Interactor3DROIEdit::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
  if (m_bEditing)
  {
    LayerROI* roi = (LayerROI*)MainWindow::GetMainWindow()->GetActiveLayer("ROI");
    if (roi && roi->GetMappedSurface())
    {
      int nvo = view->PickCurrentSurfaceVertex(event->x(), event->y(), roi->GetMappedSurface());
      if (nvo >= 0)
      {
        if (m_nPrevVertex >= 0 && m_nPrevVertex != nvo)
        {
          LayerSurface* surf = roi->GetMappedSurface();
          if (surf)
          {
            QVector<int> seeds;
            seeds << m_nPrevVertex << nvo;
            roi->EditVertex(surf->FindPath(seeds), !(event->modifiers() & Qt::ShiftModifier));
          }
        }
        else
          roi->EditVertex(nvo, !(event->modifiers() & Qt::ShiftModifier));
        m_nPrevVertex = nvo;
      }
    }
    renderview->setCursor( CursorFactory::CursorPencil );
    return false;
  }
  bool ret = Interactor3D::ProcessMouseMoveEvent(event, renderview);
  renderview->setCursor( CursorFactory::CursorPencil );
  return ret;
}

bool Interactor3DROIEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  m_bEditing = false;
  bool ret = Interactor3D::ProcessMouseUpEvent(event, renderview);
  renderview->setCursor( CursorFactory::CursorPencil );
  return ret;
}
