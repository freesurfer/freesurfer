/**
 * @brief Interactor for volume cropping in 3D render view.
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


#include "Interactor3DVolumeCrop.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerPropertyMRI.h"
#include "LayerMRI.h"
#include "VolumeCropper.h"
#include "CursorFactory.h"

Interactor3DVolumeCrop::Interactor3DVolumeCrop(QObject* parent) :
  Interactor3D(parent),
  m_bCropping( false )
{}

Interactor3DVolumeCrop::~Interactor3DVolumeCrop()
{}

bool Interactor3DVolumeCrop::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  m_nMousePosX = event->x();
  m_nMousePosY = event->y();

  if ( MainWindow::GetMainWindow()->GetVolumeCropper()->IsShown() &&
       view->PickCroppingBound( m_nMousePosX, m_nMousePosY ) )
  {
    UpdateCursor( event, view );
    m_bCropping = true;
    return false;
  }
  else
  {
    return Interactor3D::ProcessMouseDownEvent( event, view );
  }
}

bool Interactor3DVolumeCrop::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bCropping )
  {
    MainWindow::GetMainWindow()->GetVolumeCropper()->ReleaseActiveBound();
    view->RequestRedraw();
    m_bCropping = false;
    return false;
  }
  else
  {
    return Interactor3D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor3DVolumeCrop::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  int posX = event->x();
  int posY = event->y();

  if ( m_bCropping )
  {
    UpdateCursor( event, view );
    view->MoveCroppingBound( posX - m_nMousePosX, posY - m_nMousePosY );
    m_nMousePosX = posX;
    m_nMousePosY = posY;
    return false;
  }
  else
  {
    return Interactor3D::ProcessMouseMoveEvent( event, view );
  }
}

bool Interactor3DVolumeCrop::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  return Interactor3D::ProcessKeyDownEvent( event, view );
}

void Interactor3DVolumeCrop::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( MainWindow::GetMainWindow()->GetVolumeCropper()->GetActivePlane() >= 0 )
  {
    wnd->setCursor( CursorFactory::CursorGrab );
  }
  else
  {
    Interactor3D::UpdateCursor( event, wnd );
  }
}
