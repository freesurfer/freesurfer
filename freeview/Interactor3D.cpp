/**
 * @file  Interactor3D.cpp
 * @brief Base Interactor class for 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/11/01 21:38:45 $
 *    $Revision: 1.24 $
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

#include "Interactor3D.h"
#include "RenderView3D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include "CursorFactory.h"
#include "VolumeCropper.h"
#include <vtkRenderer.h>

Interactor3D::Interactor3D() :
    Interactor(),
    m_nMousePosX( -1 ),
    m_nMousePosY( -1 ),
    m_bWindowLevel( false ),
    m_bMoveSlice( false )
{}

Interactor3D::~Interactor3D()
{}

bool Interactor3D::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  m_nMousePosX = event.GetX();
  m_nMousePosY = event.GetY();

  view->CancelUpdateMouseRASPosition();

  if ( event.MiddleDown() && event.CmdDown() )
  {
    m_bWindowLevel = true;
  }
  else if ( view->GetHighlightedSlice() >= 0 )
  {
    m_bMoveSlice = true;
  }
  else
  {
    return Interactor::ProcessMouseDownEvent( event, renderview ); // pass down the event
  }

  return false; // do not pass down the event
}

bool Interactor3D::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  if ( event.GetX() == m_nMousePosX && event.GetY() == m_nMousePosY )
  {
    if ( event.LeftUp() )
    {
      view->UpdateCursorRASPosition( event.GetX(), event.GetY() );
      view->UpdateConnectivityDisplay();
    }
  }
  else
  {
    m_nMousePosX = event.GetX();
    m_nMousePosY = event.GetY();
  }
  
  m_bWindowLevel = false;
  m_bMoveSlice = false;

  return Interactor::ProcessMouseUpEvent( event, renderview );
}

bool Interactor3D::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  if ( !lcm->HasAnyLayer() )
  {
    return Interactor::ProcessMouseMoveEvent( event, renderview );
  }

  int posX = event.GetX();
  int posY = event.GetY();

  if ( m_bWindowLevel )
  {
    LayerMRI* layer = ( LayerMRI* )(lcm->GetLayerCollection( "MRI" )->GetActiveLayer());
    if ( layer && layer->IsVisible() )
    {
      double scaleX = 0.002;
      double scaleY = 0.002;
      double w = ( posX - m_nMousePosX ) * scaleX;
      double l = ( posY - m_nMousePosY ) * scaleY;
      double scaleOverall = layer->GetProperties()->GetMaxValue() -
                            layer->GetProperties()->GetMinValue();
      w *= scaleOverall;
      l *= scaleOverall;
      w += layer->GetProperties()->GetWindow();
      l += layer->GetProperties()->GetLevel();
      if ( w < 0 )
        w = 0;
      layer->GetProperties()->SetWindowLevel( w, l );
    }
  }
  else if ( m_bMoveSlice )
  {
    view->MoveSliceToScreenCoord( posX, posY );
  }
  else
  {
    if ( event.MiddleIsDown() || event.RightIsDown() || event.LeftIsDown() )
    {
      view->CancelUpdateMouseRASPosition();
    }
    else
      view->UpdateMouseRASPosition( posX, posY );

    return Interactor::ProcessMouseMoveEvent( event, view );
  }
  
  m_nMousePosX = posX;
  m_nMousePosY = posY;

  UpdateCursor( event, view );
  return false;
}

bool Interactor3D::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  RenderView3D* view = ( RenderView3D* )renderview;

  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  if ( !lcm->HasAnyLayer() )
  {
    return Interactor::ProcessKeyDownEvent( event, renderview );
  }

  int nKeyCode = event.GetKeyCode();
  if  ( nKeyCode == WXK_UP )
  {
    view->MoveUp();
  }
  else if ( nKeyCode == WXK_DOWN )
  {
    view->MoveDown();
  }
  else if ( nKeyCode == WXK_LEFT )
  {
    view->MoveLeft();
  }
  else if ( nKeyCode == WXK_RIGHT )
  {
    view->MoveRight();
  }
  else if ( nKeyCode == 'R' ) //|| nKeyCode == 'F' || nKeyCode == 'S' || nKeyCode == 'W' )
  {
    // do nothing, just intercept these keycodes
  }
  else if ( nKeyCode == WXK_DELETE )
  {
    view->DeleteCurrentSelectRegion();
  }
  else
    return Interactor::ProcessKeyDownEvent( event, view );

  return false;
}

void Interactor3D::UpdateCursor( wxEvent& event, wxWindow* wnd )
{
  RenderView3D* view = ( RenderView3D* )wnd;
  if ( view->GetHighlightedSlice() >= 0 )
  {
    wnd->SetCursor( CursorFactory::CursorGrab );
  }
  else
    Interactor::UpdateCursor( event, wnd );
}

