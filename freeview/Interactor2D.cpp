/**
 * @brief Base Interactor class to manage mouse and key input on 2D render view.
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

#include "Interactor2D.h"
#include "RenderView2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerPropertyMRI.h"
#include "LayerMRI.h"
#include <vtkRenderer.h>
#include <QDebug>
#include <QTimer>
#include <QApplication>

Interactor2D::Interactor2D( QObject* parent ) : Interactor( parent ),
  m_nMousePosX( -1 ),
  m_nMousePosY( -1 ),
  m_bWindowLevel( false ),
  m_bChangeSlice( false ),
  m_bMovingCursor( false ),
  m_bSelecting( false )
{}

Interactor2D::~Interactor2D()
{}

bool Interactor2D::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  m_nMousePosX = event->x();
  m_nMousePosY = event->y();

  view->UpdateAnnotation();

  view->PickLineProfile(m_nMousePosX, m_nMousePosY);

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if ( ( event->modifiers() & CONTROL_MODIFIER ) &&  !( event->modifiers() & Qt::ShiftModifier ) )
  {
    if ( event->button() == Qt::LeftButton )
    {
#ifndef Q_OS_MAC
      view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, 2.0 );    // zoom in
#endif
      return false;
    }
    else if ( event->button() == Qt::RightButton )
    {
#ifndef Q_OS_MAC
      view->ZoomAtCursor( m_nMousePosX, m_nMousePosY, 0.5 );    // zoom out
#endif
      return false;
    }
  }

  if ( event->button() == Qt::LeftButton )
  {
    m_nDownPosX = m_nMousePosX;
    m_nDownPosY = m_nMousePosY;

    if ( !( event->modifiers() & CONTROL_MODIFIER ) &&  ( event->modifiers() & Qt::ShiftModifier ) &&
         !mainwnd->IsRepositioningSurface())
    {
      // m_bWindowLevel = true;
      return Interactor::ProcessMouseDownEvent( event, renderview );
    }
    else
    {
      view->grabMouse();
      m_bMovingCursor = true;
      view->UpdateCursorRASPosition( m_nMousePosX, m_nMousePosY,
                                     !mainwnd->IsRepositioningSurface() && ( event->modifiers() & CONTROL_MODIFIER ) &&
                                     ( event->modifiers() & Qt::ShiftModifier ) );
      view->RequestRedraw();
      if (mainwnd->IsRepositioningSurface())
      {
        if ( event->modifiers() & CONTROL_MODIFIER &&
             event->modifiers() & Qt::ShiftModifier)
          QTimer::singleShot(0, mainwnd, SIGNAL(SurfaceRepositionIntensityChanged()));
        else if (event->modifiers() & Qt::ShiftModifier)
          QTimer::singleShot(0, mainwnd, SIGNAL(SurfaceRepositionVertexChanged()));
      }
      emit CursorLocationClicked();
    }
  }
  else if ( (event->button() == Qt::MidButton && ( event->modifiers() & Qt::ShiftModifier )) ||
            (event->button() == Qt::RightButton && ( event->modifiers() & CONTROL_MODIFIER ) &&
             ( event->modifiers() & Qt::ShiftModifier )) )
  {
    m_bSelecting = true;
    view->StartSelection( m_nMousePosX, m_nMousePosY );
  }
  else if ( event->button() == Qt::RightButton &&
            !( event->modifiers() & CONTROL_MODIFIER ) &&
            ( event->modifiers() & Qt::ShiftModifier ) )
  {
    m_bWindowLevel = true;
  }
#ifdef Q_OS_MAC
  else if ( event->button() == Qt::RightButton &&
            ( event->modifiers() & CONTROL_MODIFIER ) &&
            ( event->modifiers() & Qt::ShiftModifier ) )
  {
    m_bMovingCursor = true;
    view->UpdateCursorRASPosition( m_nMousePosX, m_nMousePosY );
    view->RequestRedraw();
  }
#endif
  else if (event->button() == Qt::MidButton && (event->modifiers() & CONTROL_MODIFIER ))
  {
    m_bChangeSlice = true;
  }
  else
  {
    return Interactor::ProcessMouseDownEvent( event, renderview ); // pass down the event
  }

  return false; // do not pass down the event
}

bool Interactor2D::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  view->releaseMouse();

  if ( m_bSelecting )
  {
    view->StopSelection();
    view->RequestRedraw();
  }

  m_nMousePosX = event->x();
  m_nMousePosY = event->y();
  m_bWindowLevel = false;
  m_bChangeSlice = false;
  m_bMovingCursor = false;
  m_bSelecting = false;

  view->UpdateAnnotation();
  view->Update2DOverlay();

  if ( event->button() == Qt::LeftButton && !( event->modifiers() & Qt::ShiftModifier ) )
  {
    return false;
  }
  else
  {
    return Interactor::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2D::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if ( mainwnd->IsEmpty() )
  {
    return Interactor::ProcessMouseMoveEvent( event, renderview );
  }

  int posX = event->x();
  int posY = event->y();

  if ( m_bChangeSlice )
  {
    double* voxelSize = mainwnd->GetLayerCollection( "MRI" )->GetWorldVoxelSize();
    int nPlane = view->GetViewPlane();
    double dPixelPer = -0.25;
    double dPosDiff =  ( ( (int)( dPixelPer * ( posY - m_nDownPosY ) ) ) / dPixelPer -
                         ( (int)( dPixelPer * ( m_nMousePosY - m_nDownPosY ) ) ) / dPixelPer )
        * dPixelPer * voxelSize[nPlane];
    if ( mainwnd->OffsetSlicePosition( nPlane, dPosDiff ) )
    {
      m_nMousePosX = posX;
      m_nMousePosY = posY;
    }
  }
  else if ( m_bMovingCursor )
  {
    view->UpdateCursorRASPosition( posX, posY );
    view->RequestRedraw();
  }
  else if ( m_bWindowLevel )
  {
    QList<Layer*> layers = mainwnd->GetLayerCollection( "MRI" )->GetLayers();
    LayerMRI* layer = (LayerMRI*)mainwnd->GetActiveLayer("MRI");
    if (layer && !layer->IsWindowAdjustable())
    {
      layer = NULL;
    }
    for ( int i = 0; i < layers.size(); i++ )
    {
      LayerMRI* mri = ( LayerMRI*)layers[i];
      if (mri->IsWindowAdjustable())
      {
        if (layer == NULL || layer == mri || mri->IsObscuring())
        {
          layer = mri;
          break;
        }
      }
    }
    if ( layer )
    {
      double scaleX = 0.002;
      double scaleY = 0.002;
      double w = ( posX - m_nMousePosX ) * scaleX;
      double l = ( posY - m_nMousePosY ) * scaleY;
      double scaleOverall = layer->GetProperty()->GetMaxValue() -
          layer->GetProperty()->GetMinValue();
      w *= scaleOverall;
      l *= scaleOverall;
      switch ( layer->GetProperty()->GetColorMap() )
      {
      case LayerPropertyMRI::Grayscale:
        w += layer->GetProperty()->GetWindow();
        l += layer->GetProperty()->GetLevel();
        if ( w < 0 )
        {
          w = 0;
        }
        layer->GetProperty()->SetWindowLevel(w, l);
        break;
      case LayerPropertyMRI::Heat:
        w += layer->GetProperty()->GetHeatScaleMaxThreshold() - layer->GetProperty()->GetHeatScaleMinThreshold();
        l += (layer->GetProperty()->GetHeatScaleMaxThreshold() + layer->GetProperty()->GetHeatScaleMinThreshold())/2;
        if ( w < 0 )
        {
          w = 0;
        }
        layer->GetProperty()->SetHeatScale( l-w/2, l, l+w/2 );
        break;
      default:
        w += layer->GetProperty()->GetMaxGenericThreshold() - layer->GetProperty()->GetMinGenericThreshold();
        l += (layer->GetProperty()->GetMaxGenericThreshold() + layer->GetProperty()->GetMinGenericThreshold())/2;
        if ( w < 0 )
        {
          w = 0;
        }
        layer->GetProperty()->SetMinMaxGenericThreshold( l-w/2, l+w/2 );
        break;
      }
    }
    m_nMousePosX = posX;
    m_nMousePosY = posY;
  }
  else if ( m_bSelecting )
  {
    view->UpdateSelection( posX, posY );
    view->RequestRedraw();
  }
  else
  {
    if ( event->buttons() & Qt::MidButton || event->buttons() & Qt::RightButton ||
         ((event->buttons() & Qt::LeftButton) && (event->modifiers() & Qt::ShiftModifier)))
    {
      view->UpdateAnnotation();
      view->Update2DOverlay();
      if ( event->buttons() & Qt::RightButton )
      {
        view->EmitZooming();
      }
    }
    else
    {
      view->UpdateMouseRASPosition( posX, posY );
    }

    return Interactor::ProcessMouseMoveEvent( event, renderview );
  }

  return false;
}

void Interactor2D::ProcessPostMouseWheelEvent( QWheelEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  view->UpdateAnnotation();
  view->Update2DOverlay();
  view->RequestRedraw();
  view->EmitZooming();

  Interactor::ProcessPostMouseWheelEvent( event, renderview );
}

void Interactor2D::ProcessPostMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  if (event->buttons() & Qt::RightButton)
  {
    view->Update2DOverlay();
    view->RequestRedraw();
  }

  Interactor::ProcessPostMouseMoveEvent( event, renderview );
}

bool Interactor2D::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  int nKeyCode = event->key();
  if (nKeyCode == Qt::Key_Escape)
  {
    m_bWindowLevel = false;
    m_bChangeSlice = false;
    m_bMovingCursor = false;
    m_bSelecting = false;
  }

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  if ( mainwnd->IsEmpty() )
  {
    return Interactor::ProcessKeyDownEvent( event, renderview );
  }

  if ( nKeyCode == Qt::Key_Plus )
  {
    LayerMRI* mri = (LayerMRI*)mainwnd->GetActiveLayer("MRI");
    if (mri && mri->GetNumberOfFrames() > 1)
    {
      int n = mri->GetActiveFrame() + 1;
      if (n >= mri->GetNumberOfFrames())
        n = 0;
      mri->SetActiveFrame(n);
    }
  }
  else if ( nKeyCode == Qt::Key_Minus )
  {
    LayerMRI* mri = (LayerMRI*)mainwnd->GetActiveLayer("MRI");
    if (mri && mri->GetNumberOfFrames() > 1)
    {
      int n = mri->GetActiveFrame() - 1;
      if (n < 0)
        n = mri->GetNumberOfFrames()-1;
      mri->SetActiveFrame(n);
    }
  }
  else if ( event->modifiers() & Qt::ShiftModifier )
  {
    if ( nKeyCode == Qt::Key_Up )
    {
      view->Zoom(1.05);
    }
    else if ( nKeyCode == Qt::Key_Down )
    {
      view->Zoom(0.95);
    }
  }
  else if (event->modifiers() & Qt::ControlModifier )
  {
    if ( nKeyCode == Qt::Key_Up )
    {
      view->MoveUp();
    }
    else if ( nKeyCode == Qt::Key_Down )
    {
      view->MoveDown();
    }
    else if ( nKeyCode == Qt::Key_Left )
    {
      view->MoveLeft();
    }
    else if ( nKeyCode == Qt::Key_Right )
    {
      view->MoveRight();
    }
  }
  else if ( nKeyCode == Qt::Key_PageUp || nKeyCode == Qt::Key_Up )
  {
    view->MoveSlice( 1 );
  }
  else if ( nKeyCode == Qt::Key_PageDown || nKeyCode == Qt::Key_Down )
  {
    view->MoveSlice( -1 );
  }
  else if ( nKeyCode == Qt::Key_3 /*|| nKeyCode == 'W' || nKeyCode == 'S'*/ || nKeyCode == Qt::Key_R || nKeyCode == Qt::Key_F )
  {
    // do nothing, just intercept these vtk default keycodes
  }
  else
  {
    return Interactor::ProcessKeyDownEvent( event, view );
  }

  return false;
}
