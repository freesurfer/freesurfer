/**
 * @brief Interactor for measure tool in 2D render view.
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

#include "Interactor2DMeasure.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "Region2DLine.h"
#include "Region2DPolyline.h"
#include "Region2DRectangle.h"
//#include "ToolWindowMeasure.h"
#include "CursorFactory.h"
#include <vtkRenderer.h>
#include <QDebug>

Interactor2DMeasure::Interactor2DMeasure( QObject* parent ) :
  Interactor2D( parent ),
  m_bEditing( false ),
  m_bDrawing( false ),
  m_nPointIndex( -1 ),
  m_region( NULL )
{
}

Interactor2DMeasure::~Interactor2DMeasure()
{}

bool Interactor2DMeasure::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  // UpdateCursor( event, view );

  if ( m_region && !m_bDrawing && !m_bEditing )
  {
    m_region->Highlight( false );
  }

  if ( event->button() == Qt::LeftButton )
  {
    if ( ( event->modifiers() & CONTROL_MODIFIER ) && ( event->modifiers() & Qt::ShiftModifier ) )
    {
      view->UpdateCursorRASPosition( event->x(), event->y());
      view->RequestRedraw();
      return false;
      //  return Interactor2D::ProcessMouseDownEvent( event, renderview );
    }

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();
    if ( mri )
    {
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();

      if ( m_region && m_bDrawing ) // drawing
      {
        ((Region2DPolyline*)m_region)->AddPoint( m_nMousePosX, m_nMousePosY );
      }
      else
      {
        Region2D* reg = view->GetRegion( m_nMousePosX, m_nMousePosY, &m_nPointIndex );
        if ( !reg ) // new region
        {
          if ( m_nAction == MM_Line )
          {
            Region2DLine* reg_line = new Region2DLine( view );
            reg_line->SetLine( m_nMousePosX, m_nMousePosY, m_nMousePosX, m_nMousePosY );
            view->AddRegion( reg_line );
            m_region = reg_line;
          }
          else if ( m_nAction == MM_Spline || m_nAction == MM_Polyline )
          {
            Region2DPolyline* reg_polyline = new Region2DPolyline( view, m_nAction == MM_Spline );
            reg_polyline->AddPoint( m_nMousePosX, m_nMousePosY );
            reg_polyline->AddPoint( m_nMousePosX, m_nMousePosY ); // add second point
            view->AddRegion( reg_polyline );
            m_region = reg_polyline;
          }
          else if ( m_nAction == MM_Rectangle )
          {
            Region2DRectangle* reg_rect = new Region2DRectangle( view );
            reg_rect->SetRect( m_nMousePosX, m_nMousePosY, 1, 1 );
            view->AddRegion( reg_rect );
            m_region = reg_rect;
          }
          m_bDrawing = true;
        }
        else      // editing
        {
          m_region = reg;
          m_bEditing = true;
          m_region->Highlight();
          view->EmitRegionSelected( reg );
          view->RequestRedraw();
        }
      }
      return false;
    }
  }
  else if ( event->button() == Qt::RightButton )
  {
    if ( m_bDrawing && m_region )
    {
      m_bDrawing = false;
      m_bEditing = false;
      if ( m_nAction == MM_Spline || m_nAction == MM_Polyline )
      {
        ((Region2DPolyline*)m_region)->RemoveLastPoint();
      }
      view->RequestRedraw();
      return false;
    }
  }

  return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DMeasure::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );

  if ( m_bDrawing )
  {
    if ( ( m_nAction == MM_Spline || m_nAction == MM_Polyline ) && m_region )
    {
      return false;
    }

    if ( m_nMousePosX != event->x() || m_nMousePosY != event->y() )
    {
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();
      if ( m_region )
      {
        if ( m_nAction == MM_Line )
        {
          ((Region2DLine*)m_region)->SetPoint2( m_nMousePosX, m_nMousePosY );
        }
        else if ( m_nAction == MM_Rectangle )
        {
          ((Region2DRectangle*)m_region)->SetBottomRight( m_nMousePosX, m_nMousePosY );
        }
      }
    }
    else
    {
      if ( m_region && m_nAction != MM_Polyline && m_nAction != MM_Spline )
      {
        view->DeleteRegion( m_region );
        m_region = NULL;
      }
    }

    m_bEditing = false;
    m_bDrawing = false;
    return false;
  }
  else if ( m_bEditing )
  {
    m_bEditing = false;
    m_bDrawing = false;
    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2DMeasure::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  if ( m_bDrawing )
  {
    UpdateCursor( event, view );
    int posX = event->x();
    int posY = event->y();

    //    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( m_strLayerTypeName.c_str() );
    //    LayerVolumeBase* mri = ( LayerVolumeBase* )lc->GetActiveLayer();

    if ( m_region )
    {
      if ( m_nAction == MM_Line )
      {
        ((Region2DLine*)m_region)->SetPoint2( posX, posY );
      }
      else if ( m_nAction == MM_Spline || m_nAction == MM_Polyline )
      {
        ((Region2DPolyline*)m_region)->UpdatePoint( -1, posX, posY );
      }
      else if ( m_nAction == MM_Rectangle )
      {
        ((Region2DRectangle*)m_region)->SetBottomRight( posX, posY );
      }
      view->RequestRedraw();
    }

    return false;
  }
  else if ( m_bEditing )
  {
    UpdateCursor( event, view );
    int offsetX = event->x() - m_nMousePosX;
    int offsetY = event->y() - m_nMousePosY;
    if ( m_region )
    {
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();
      if ( m_nPointIndex >= 0 )
      {
        m_region->UpdatePoint( m_nPointIndex, m_nMousePosX, m_nMousePosY );
      }
      else
      {
        m_region->Offset( offsetX, offsetY );
      }
      view->RequestRedraw();
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

bool Interactor2DMeasure::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  UpdateCursor( event, renderview );
  if ( m_region && (event->key() == Qt::Key_Delete || event->key() == Qt::Key_Backspace) )
  {
    view->DeleteRegion( m_region );
    m_region = NULL;
    return false;
  }
  else
  {
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  }
}

bool Interactor2DMeasure::ProcessKeyUpEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  return Interactor2D::ProcessKeyUpEvent( event, renderview );
}

void Interactor2DMeasure::UpdateCursor( QEvent* event, QWidget* wnd )
{
  RenderView2D* view = ( RenderView2D* )wnd;
  if ( wnd->hasFocus() )
  {
    if ( event->type() == QEvent::MouseButtonPress ||
         event->type() == QEvent::MouseButtonRelease ||
         event->type() == QEvent::MouseMove)
    {
      QMouseEvent* e = ( QMouseEvent* )event;
      if ( ( ( e->button() == Qt::MidButton || e->button() == Qt::RightButton ) && !m_bEditing ) ||
           ( ( e->modifiers() & CONTROL_MODIFIER) && (e->modifiers() & Qt::ShiftModifier) ) )
      {
        Interactor2D::UpdateCursor( event, wnd );
        return;
      }
      else if ( ( !m_bEditing && !m_bDrawing && view->GetRegion( e->x(), e->y() ) ) ||
                ( m_bEditing && m_nPointIndex < 0 ) )
      {
        wnd->setCursor( CursorFactory::CursorGrab );
      }
      else
      {
        if ( m_nAction == MM_Line )
        {
          wnd->setCursor( CursorFactory::CursorMeasureLine );
        }
        else if ( m_nAction == MM_Rectangle )
        {
          wnd->setCursor( CursorFactory::CursorMeasureRectangle );
        }
        else if ( m_nAction == MM_Polyline || m_nAction == MM_Spline )
        {
          wnd->setCursor( CursorFactory::CursorMeasurePolyline );
        }
      }
    }
  }
  else
  {
    Interactor2D::UpdateCursor( event, wnd );
  }
}
