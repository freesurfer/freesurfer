/**
 * @file  Interactor2DCropVolume.cpp
 * @brief Interactor for volume cropping in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
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

#include "Interactor2DCropVolume.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "CursorFactory.h"
#include "VolumeCropper.h"
#include <vtkRenderer.h>

Interactor2DCropVolume::Interactor2DCropVolume() :
    Interactor2D(),
    m_bSelected( false )
{
}

Interactor2DCropVolume::~Interactor2DCropVolume()
{}

bool Interactor2DCropVolume::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
 
  if ( event.LeftDown() &&
       MainWindow::GetMainWindowPointer()->GetVolumeCropper()
       ->PickActiveBound2D( view, event.GetX(), event.GetY() ) )
  {
    UpdateCursor( event, view );
    m_bSelected = true;
    return false;
  }
  else 
    return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DCropVolume::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
  if ( m_bSelected )
  {
    MainWindow::GetMainWindowPointer()->GetVolumeCropper()->ReleaseActiveBound();
    MainWindow::GetMainWindowPointer()->GetRenderView( 3 )->NeedRedraw();
    m_bSelected = false;
    return false;
  }
  else 
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
}

bool Interactor2DCropVolume::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bSelected )
  {
    UpdateCursor( event, view );

    MainWindow::GetMainWindowPointer()->GetVolumeCropper()
        ->MoveActiveBound2D( view, event.GetX(), event.GetY() );
    
    return false;
  }
  else if ( !event.LeftIsDown() && !event.RightIsDown() && !event.MiddleIsDown() &&       
            MainWindow::GetMainWindowPointer()->GetVolumeCropper()
            ->PickActiveBound2D( view, event.GetX(), event.GetY() ) )
  {
    UpdateCursor( event, view );
    return false;
  }
  else 
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

void Interactor2DCropVolume::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  if ( MainWindow::GetMainWindowPointer()->GetVolumeCropper()->GetActivePlane() >= 0 )
    wnd->SetCursor( CursorFactory::CursorGrab );
  else
    Interactor2D::UpdateCursor( event, wnd );
}
