/**
 * @file  Interactor2D.cpp
 * @brief Base Interactor class to manage mouse and key input on 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.24 $
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

#include "Interactor2D.h"
#include "RenderView2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>

Interactor2D::Interactor2D() : Interactor(),
    m_nMousePosX( -1 ),
    m_nMousePosY( -1 ),
    m_bWindowLevel( false ),
    m_bChangeSlice( false ),
    m_bMovingCursor( false ),
    m_bSelecting( false )
{}

Interactor2D::~Interactor2D()
{}

bool Interactor2D::ProcessMouseDownEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  m_nMousePosX = event.GetX();
  m_nMousePosY = event.GetY();

  view->UpdateAnnotation();

  if ( event.CmdDown() && !event.ShiftDown() )
  {
    if ( event.LeftDown() )
    {
      view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, true );
      return false;
    }
    else if ( event.RightDown() )
    {
      view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, false );
      return false;
    }
  }

  if ( event.LeftDown() )
  {
    m_nDownPosX = m_nMousePosX;
    m_nDownPosY = m_nMousePosY;

    if ( event.ShiftDown() && !event.CmdDown() )
    {
      m_bWindowLevel = true;
    }
    else
    {
      m_bMovingCursor = true;
      view->UpdateCursorRASPosition( m_nMousePosX, m_nMousePosY );
      view->NeedRedraw();
    }
  }
  else if ( event.MiddleDown() && event.ShiftDown() )
  {
      m_bSelecting = true;
      view->StartSelection( m_nMousePosX, m_nMousePosY );
  }
  else if ( event.RightDown() && event.ShiftDown() && !event.CmdDown() )
  {
    m_bWindowLevel = true;
  }
  else
  {
    return Interactor::ProcessMouseDownEvent( event, renderview ); // pass down the event
  }

  return false; // do not pass down the event
}

bool Interactor2D::ProcessMouseUpEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bSelecting )
  {
    view->StopSelection();
    view->NeedRedraw();
  }
  
  m_nMousePosX = event.GetX();
  m_nMousePosY = event.GetY();
  m_bWindowLevel = false;
  m_bChangeSlice = false;
  m_bMovingCursor = false;
  m_bSelecting = false;

  view->UpdateAnnotation();
  view->Update2DOverlay();

  if ( event.LeftUp() )
  {
    return false;
  }
  else
  {
    return Interactor::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2D::ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  if ( !lcm->HasAnyLayer() )
  {
    return Interactor::ProcessMouseMoveEvent( event, renderview );
  }

  int posX = event.GetX();
  int posY = event.GetY();

  if ( m_bChangeSlice )
  {
    double* voxelSize = lcm->GetLayerCollection( "MRI" )->GetWorldVoxelSize();
    int nPlane = view->GetViewPlane();
    double dPixelPer = -0.20;
    double dPosDiff =  ( ( (int)( dPixelPer * ( posY - m_nDownPosY ) ) ) / dPixelPer -
                         ( (int)( dPixelPer * ( m_nMousePosY - m_nDownPosY ) ) ) / dPixelPer )
                       * dPixelPer * voxelSize[nPlane];
    if ( lcm->OffsetSlicePosition( nPlane, dPosDiff ) )
    {
      m_nMousePosX = posX;
      m_nMousePosY = posY;
    }
  }
  else if ( m_bMovingCursor )
  {
    view->UpdateCursorRASPosition( posX, posY );
    view->NeedRedraw();
  }
  else if ( m_bWindowLevel )
  {
    std::vector<Layer*> layers = lcm->GetLayerCollection( "MRI" )->GetLayers();
    LayerMRI* layer = NULL;
    for ( size_t i = 0; i < layers.size(); i++ )
    {
      layer = ( LayerMRI*)layers[i];
      if ( layer->IsVisible() && layer->GetProperties()->GetColorMap() == LayerPropertiesMRI::Grayscale )
        break;
    }
    if ( layer )
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
    m_nMousePosX = posX;
    m_nMousePosY = posY;
  }
  else if ( m_bSelecting )
  {
    view->UpdateSelection( posX, posY );   
    view->NeedRedraw();
  }
  else
  {
    if ( event.MiddleIsDown() || event.RightIsDown() )
    {
      view->UpdateAnnotation();
      view->Update2DOverlay();
      if ( event.RightIsDown() )
        view->SendBroadcast( "Zooming", view );
    }
    else
      view->UpdateMouseRASPosition( posX, posY );

    return Interactor::ProcessMouseMoveEvent( event, renderview );
  }

  return false;
}

void Interactor2D::ProcessPostMouseWheelEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  view->UpdateAnnotation();
  view->Update2DOverlay();
  view->NeedRedraw();
  view->SendBroadcast( "Zooming", view );

  Interactor::ProcessPostMouseWheelEvent( event, renderview );
}

void Interactor2D::ProcessPostMouseMoveEvent( wxMouseEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  if ( event.RightIsDown() )
  {
    view->Update2DOverlay();
    view->NeedRedraw();
  }

  Interactor::ProcessPostMouseMoveEvent( event, renderview );
}

bool Interactor2D::ProcessKeyDownEvent( wxKeyEvent& event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  if ( !lcm->HasAnyLayer() )
  {
    return Interactor::ProcessKeyDownEvent( event, renderview );
  }

  int nKeyCode = event.GetKeyCode();
  if ( nKeyCode == WXK_PAGEUP )
  {
    view->MoveSlice( 1 );
  }
  else if ( nKeyCode == WXK_PAGEDOWN)
  {
    view->MoveSlice( -1 );
  }
  else if ( nKeyCode == WXK_UP )
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
  else if ( nKeyCode == '3' /*|| nKeyCode == 'W' || nKeyCode == 'S'*/ || nKeyCode == 'R' || nKeyCode == 'F' )
  {
    // do nothing, just intercept these vtk default keycodes
  }
  else
    return Interactor::ProcessKeyDownEvent( event, view );

  return false;
}
