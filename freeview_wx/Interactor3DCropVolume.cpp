/**
 * @file  Interactor3DCropVolume.cpp
 * @brief Interactor for volume cropping in 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "Interactor3DCropVolume.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include "VolumeCropper.h"
#include "CursorFactory.h"

Interactor3DCropVolume::Interactor3DCropVolume() :
    Interactor3D(),
    m_bCropping( false )
{}

Interactor3DCropVolume::~Interactor3DCropVolume()
{}

bool Interactor3DCropVolume::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  m_nMousePosX = event.GetX();
  m_nMousePosY = event.GetY();
  
  if ( MainWindow::GetMainWindowPointer()->GetVolumeCropper()->IsShown() && 
            view->PickCroppingBound( m_nMousePosX, m_nMousePosY ) )
  {
    UpdateCursor( event, view );
    m_bCropping = true;
    return false;
  }
  else
    return Interactor3D::ProcessMouseDownEvent( event, view );
}

bool Interactor3DCropVolume::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( m_bCropping )
  {
    MainWindow::GetMainWindowPointer()->GetVolumeCropper()->ReleaseActiveBound();
    view->NeedRedraw();
    m_bCropping = false;
    return false;
  }
  else
    return Interactor3D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3DCropVolume::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  int posX = event.GetX();
  int posY = event.GetY();

  if ( m_bCropping )
  {
    UpdateCursor( event, view );
    view->MoveCroppingBound( posX - m_nMousePosX, posY - m_nMousePosY );
    m_nMousePosX = posX;
    m_nMousePosY = posY;
    return false;
  }
  else
    return Interactor3D::ProcessMouseMoveEvent( event, view );
}

bool Interactor3DCropVolume::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;
 
  return Interactor3D::ProcessKeyDownEvent( event, view );
}

void Interactor3DCropVolume::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  if ( MainWindow::GetMainWindowPointer()->GetVolumeCropper()->GetActivePlane() >= 0 )
    wnd->SetCursor( CursorFactory::CursorGrab );
  else
    Interactor3D::UpdateCursor( event, wnd );
}
