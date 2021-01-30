/**
 * @brief Interactor for editing volume in 2D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#include "Interactor2DVolumeEdit.h"
#include "RenderView2D.h"
#include "Cursor2D.h"
#include "MainWindow.h"
#include "LayerCollection.h"
#include "LayerVolumeBase.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "Contour2D.h"
#include "CursorFactory.h"
#include "BrushProperty.h"
#include <vtkRenderer.h>
#include <QDebug>
#include "LayerROI.h"
#include "LayerPropertyROI.h"
#include <QTimer>

Interactor2DVolumeEdit::Interactor2DVolumeEdit( const QString& layerTypeName, QObject* parent ) :
  Interactor2D( parent ),
  m_bEditing( false ),
  m_bColorPicking( false )
{
  m_strLayerTypeName = layerTypeName;
}

Interactor2DVolumeEdit::~Interactor2DVolumeEdit()
{}

void Interactor2DVolumeEdit::PreprocessMouseEvent(QMouseEvent *event)
{
  if (m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    return;

  bool bRightButtonErase = MainWindow::GetMainWindow()->GetSetting("RightButtonErase").toBool();
  if (bRightButtonErase && event->button() == Qt::RightButton && event->modifiers() == Qt::NoModifier)
  {
    QMouseEvent e(event->type(), event->pos(), Qt::LeftButton, Qt::LeftButton, Qt::ShiftModifier);
    *event = e;
  }
}

#include "LUTDataHolder.h"

bool Interactor2DVolumeEdit::ProcessMouseDownEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  //  if ( !view->hasFocus() )
  //  {
  //    return Interactor2D::ProcessMouseDownEvent( event, renderview );
  //  }

  PreprocessMouseEvent(event);

  if (m_nAction == EM_GeoSeg)
  {
    if ( (event->button() == Qt::LeftButton || event->button() == Qt::RightButton)
         && !(event->modifiers() & CONTROL_MODIFIER))
    {
      BrushProperty* bp = MainWindow::GetMainWindow()->GetBrushProperty();
      LayerMRI* mri = (LayerMRI*)bp->GetReferenceLayer();
      if (!mri)
      {
        emit Error( QString("Must select a reference layer."));
        return false;
      }

      LayerMRI* layer_draw = qobject_cast<LayerMRI*>(MainWindow::GetMainWindow()->FindSupplementLayer("GEOS_DRAW"));
      LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection("Supplement");
      if (!layer_draw)
      {
        layer_draw = new LayerMRI(mri);
        if ( !layer_draw->Create( mri, false, MRI_UCHAR) )
        {
          //  QMessageBox::warning( this, "Error", "Can not create drawing layer." );
          delete layer_draw;
          return false;
        }
        QVariantMap map = bp->GetGeosSettings();
        QMap<int, QColor> colors;
        colors[0] = QColor(0,0,0,0);
        colors[1] = map.contains("ForegroundColor")?map["ForegroundColor"].value<QColor>():QColor(0,255,0);
        colors[2] = map.contains("BackgroundColor")?map["BackgroundColor"].value<QColor>():QColor(255,0,0);
        layer_draw->GetProperty()->SetCustomColors(colors);
        layer_draw->GetProperty()->SetColorMapToLUT();
        layer_draw->SetName( "GEOS_DRAW" );

        LayerMRI* layer_fill = new LayerMRI(mri);
        if ( !layer_fill->Create( mri, false, MRI_UCHAR) )
        {
          //  QMessageBox::warning( this, "Error", "Can not create drawing layer." );
          delete layer_fill;
          return false;
        }
        colors.clear();
        colors[0] = QColor(0,0,0,0);
        colors[1] = map.contains("FillColor")?map["FillColor"].value<QColor>():QColor(255,255,0);
        layer_fill->GetProperty()->SetCustomColors(colors);
        layer_fill->GetProperty()->SetColorMapToLUT();
        layer_fill->SetName( "GEOS_FILL" );
        layer_fill->SetFillValue(1);
        if (map.contains("Opacity"))
          layer_fill->GetProperty()->SetOpacity(map.value("Opacity").toDouble());

        lc->AddLayer(layer_fill);
        lc->AddLayer(layer_draw);
      }

      layer_draw->SetFillValue(event->button() == Qt::RightButton ? 2:1);
      double ras[3];
      m_nMousePosX = event->x();
      m_nMousePosY = event->y();
      view->MousePositionToRAS( event->x(), event->y(), ras );
      layer_draw->SaveForUndo(view->GetViewPlane());
      if (event->modifiers() & Qt::ShiftModifier)
        layer_draw->SetVoxelByRAS(ras, view->GetViewPlane(), false);
      else
        layer_draw->SetVoxelByRAS(ras, view->GetViewPlane());
      m_bEditing = true;
      view->grabMouse();
      return false;
    }
  }

  if ( event->button() == Qt::LeftButton ||
       ( event->button() == Qt::RightButton && (event->buttons() & Qt::LeftButton) ) )
  {
    if ( (event->modifiers() & Qt::ControlModifier ) && (event->modifiers() & Qt::ShiftModifier) )
    {
      if (event->button() == Qt::LeftButton)
      {
        view->UpdateCursorRASPosition( event->x(), event->y());
        view->RequestRedraw();
        return false;
      }
      else
        return Interactor2D::ProcessMouseDownEvent( event, renderview );
    }

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetSelectedLayers(m_strLayerTypeName);
    if (layers.isEmpty())
      layers << lc->GetActiveLayer();
    QList<LayerVolumeBase*> mriLayers;
    foreach (Layer* layer, layers)
      mriLayers << qobject_cast<LayerVolumeBase*>(layer);

    foreach (LayerVolumeBase* mri, mriLayers)
    {
      if ( (!mri || !mri->IsVisible()) ) //&& ( event->ControlDown() || m_nAction == EM_Polyline ) )
      {
        emit Error( QString("Selected layer '%1' is not visible").arg(mri->GetName()), mri );
        break;
      }
      else if ( !mri->IsEditable() ) //&& ( event->ControlDown() || m_nAction == EM_Polyline ) )
      {
        emit Error( QString("Selected layer '%1' is not editable").arg(mri->GetName()), mri );
        break;
      }
      else if ( m_strLayerTypeName == "MRI" && ((LayerMRI*)mri)->IsTransformed() )
      {
        emit Error( QString("Selected layer '%1' is not editable for transformation").arg(mri->GetName()), mri );
        break;
      }
      else
      {
        m_nMousePosX = event->x();
        m_nMousePosY = event->y();
        bool bFill3D = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetFill3D();

        double ras[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras );
        bool bCondition = !(event->modifiers() & Qt::ShiftModifier) && !(event->buttons() & Qt::RightButton);
        if (bCondition && (m_nAction == EM_Freehand || m_nAction == EM_Fill || m_nAction == EM_Polyline))
        {
          if (mri->IsTypeOf("MRI") && ((LayerMRI*)mri)->GetProperty()->GetColorMap() == LayerPropertyMRI::LUT
              && !((LayerMRI*)mri)->GetProperty()->IsValueInColorTable(mri->GetFillValue()))
          {
            emit Error("Brush value is not in the current color table. I don't know what color to use. Drawing cannot continue.");
            return false;
          }
        }

        if ( (m_nAction == EM_ColorPicker || m_bColorPicking ) && mri->IsTypeOf( "MRI" ) )
        {
          if ( event->modifiers() & CONTROL_MODIFIER )
          {
            mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
            mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
          }
          else
          {
            double dValue = ((LayerMRI*)mri)->GetVoxelValue( ras );
            if ( dValue != 0 )
            {
              MainWindow::GetMainWindow()->GetBrushProperty()->SetFillValue( (float)dValue );
            }
          }
          m_bColorPicking = false;
        }
        else if ( m_nAction == EM_Freehand ) //&& ( event->ControlDown() ) )
        {
          if ( event->modifiers() & CONTROL_MODIFIER )
          {
            mri->SaveForUndo( bFill3D ? -1 : view->GetViewPlane());
            mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
          }
          else
          {
            mri->SaveForUndo( view->GetViewPlane() );
            m_bEditing = true;
            mri->SetVoxelByRAS( ras, view->GetViewPlane(), bCondition );
          }
          view->grabMouse();
        }
        else if ( m_nAction == EM_Clone )
        {
          mri->SaveForUndo( view->GetViewPlane() );
          m_bEditing = true;
          mri->CloneVoxelByRAS( ras, view->GetViewPlane() );
        }
        else if ( m_nAction == EM_Shift)
        {
          mri->SaveForUndo(view->GetViewPlane());
          m_bEditing = true;
          mri->PrepareShifting(view->GetViewPlane());
        }
        else if ( m_nAction == EM_Fill )
        {
          mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
          if (event->modifiers() & CONTROL_MODIFIER)
          {
            mri->BorderFillByRAS(ras, view->GetViewPlane(), true);
          }
          else
            mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
        }
        else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
        {
          if ( event->modifiers() & CONTROL_MODIFIER )
          {
            mri->SaveForUndo(bFill3D ? -1 : view->GetViewPlane());
            mri->FloodFillByRAS( ras, view->GetViewPlane(), bCondition, bFill3D );
          }
          else
          {
            mri->SaveForUndo( view->GetViewPlane() );
            m_bEditing = true;
            double ras2[3];
            view->GetCursor2D()->ClearInterpolationPoints();
            view->GetCursor2D()->GetPosition( ras2 );
            view->GetCursor2D()->SetPosition( ras );
            view->GetCursor2D()->SetPosition2( ras );
            if ( m_dPolylinePoints.size() > 0 )
            {
              if ( m_nAction == EM_Polyline )
              {
                mri->SetVoxelByRAS( ras, ras2, view->GetViewPlane(), bCondition );
              }
              else
              {
                mri->SetLiveWireByRAS( ras, ras2, view->GetViewPlane() );
              }
            }
            else
            {
              // mri->SaveForUndo( view->GetViewPlane() );
              m_dPolylinePoints.push_back( ras[0] );
              m_dPolylinePoints.push_back( ras[1] );
              m_dPolylinePoints.push_back( ras[2] );
              view->GetCursor2D()->SetPosition( ras );
            }

            view->grabMouse();
          }
        }
        else if ( m_nAction == EM_Contour && mri->IsTypeOf( "MRI" ) )
        {
          LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetReferenceLayer();
          if ( !mri_ref )
          {
            emit Error( "Reference Layer Not Set" );
            return false;
          }

          Contour2D* c2d = view->GetContour2D();
          if ( (event->modifiers() & CONTROL_MODIFIER) && (event->modifiers() & Qt::AltModifier) )
          {
            double dValue = mri_ref->GetVoxelValue( ras );
            if ( dValue != 0 )
            {
              m_bEditing = true;
              c2d->SetInput( mri_ref->GetSliceImageData( view->GetViewPlane() ), dValue, ras[view->GetViewPlane()], mri_ref->GetActiveFrame() );
              c2d->SetVisible( true );
              view->RequestRedraw();
            }
            else if ( c2d->IsVisible() )
            {
              m_bEditing = true;
            }
          }
          else if ( (event->modifiers() & CONTROL_MODIFIER) && !(event->modifiers() & Qt::AltModifier) )
          {
            mri->SaveForUndo( view->GetViewPlane() );
            ((LayerMRI*)mri)->FloodFillByContour2D( ras, c2d );
          }
          else if ( event->modifiers() & Qt::ShiftModifier )
          {
            m_bEditing = true;
            c2d->RemoveLine( ras, ras );
            view->RequestRedraw();
          }
          else
          {
            m_bEditing = true;
            c2d->AddLine( ras, ras );
            view->RequestRedraw();
          }
        }
        else
        {
          return Interactor2D::ProcessMouseDownEvent( event, renderview );
        }
      }
    }

    return false;
  }
  else if ( m_bEditing )
  {
    m_bEditing = false;
    if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    {
      if ( event->button() == Qt::MidButton )
      {
        view->GetCursor2D()->Update();
        view->RequestRedraw();
      }
      else if ( event->button() == Qt::RightButton )
      {
        if ( m_dPolylinePoints.size() > 0 && m_nAction == EM_Polyline )
        {
          LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
          QList<Layer*> layers = MainWindow::GetMainWindow()->GetSelectedLayers(m_strLayerTypeName);
          if (layers.isEmpty())
            layers << lc->GetActiveLayer();

          double ras1[3] = { m_dPolylinePoints[0], m_dPolylinePoints[1], m_dPolylinePoints[2] };
          double ras2[3];
          view->GetCursor2D()->GetPosition( ras2 );
          view->GetCursor2D()->SetPosition2( ras2 );
          view->GetCursor2D()->SetPosition( ras1 );
          foreach (Layer* layer, layers)
          {
            LayerVolumeBase* mri = qobject_cast<LayerVolumeBase*>(layer);
            if (mri)
              mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), !(event->modifiers() & Qt::ShiftModifier) );
          }
        }
        else
        {
          // mri->SetLiveWireByRAS( ras1, ras2, view->GetViewPlane() );
          view->GetCursor2D()->Update();
          view->RequestRedraw();
        }
      }
    }

    m_dPolylinePoints.clear();
    view->releaseMouse();

    return false;
  }
  return Interactor2D::ProcessMouseDownEvent( event, renderview ); // pass down the event
}

bool Interactor2DVolumeEdit::ProcessMouseUpEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;
  view->releaseMouse();
  PreprocessMouseEvent(event);

  UpdateCursor( event, renderview );

  //  if (m_nAction == EM_GeoSeg)
  //    QTimer::singleShot(0, MainWindow::GetMainWindow(), SIGNAL(SupplementLayerChanged()));

  if ( m_bEditing )
  {
    m_nMousePosX = event->x();
    m_nMousePosY = event->y();

    if ( event->button() != Qt::LeftButton ||
         (m_nAction != EM_Polyline && m_nAction != EM_Livewire ) ||
         m_dPolylinePoints.size() == 0 )
    {
      m_bEditing = false;
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseUpEvent( event, renderview );
  }
}

bool Interactor2DVolumeEdit::ProcessMouseMoveEvent( QMouseEvent* event, RenderView* renderview )
{
  RenderView2D* view = ( RenderView2D* )renderview;

  PreprocessMouseEvent(event);

  if (m_nAction == EM_GeoSeg && m_bEditing)
  {
    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection("Supplement");
    QList<Layer*> layers = lc->GetLayers("MRI");    // Get foreground/background drawing layers
    LayerMRI* layer_draw = NULL;
    foreach (Layer* layer, layers)
    {
      if (layer->GetName() == "GEOS_DRAW")
      {
        layer_draw = (LayerMRI*)layer;
        break;
      }
    }
    if (!layer_draw)
      return false;

    layer_draw->SetFillValue(event->buttons() & Qt::RightButton ? 2:1);
    double ras1[3], ras2[3];
    view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
    view->MousePositionToRAS( event->x(), event->y(), ras2 );
    if (event->modifiers() & Qt::ShiftModifier)
      layer_draw->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(), false);
    else
      layer_draw->SetVoxelByRAS( ras1, ras2, view->GetViewPlane());

    m_nMousePosX = event->x();
    m_nMousePosY = event->y();
    return false;
  }

  if ( m_bEditing )
  {
    UpdateCursor( event, view );
    int posX = event->x();
    int posY = event->y();

    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetSelectedLayers(m_strLayerTypeName);
    if (layers.isEmpty())
      layers << lc->GetActiveLayer();
    QList<LayerVolumeBase*> mriLayers;
    foreach (Layer* layer, layers)
      mriLayers << qobject_cast<LayerVolumeBase*>(layer);

    if ( m_nAction == EM_Freehand )
    {
      double ras1[3], ras2[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
      view->MousePositionToRAS( posX, posY, ras2 );

      foreach (LayerVolumeBase* mri, mriLayers)
        mri->SetVoxelByRAS( ras1, ras2, view->GetViewPlane(),
                            !(event->modifiers() & Qt::ShiftModifier) && !(event->buttons() & Qt::RightButton) );
    }
    else if ( m_nAction == EM_Clone )
    {
      double ras1[3], ras2[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
      view->MousePositionToRAS( posX, posY, ras2 );

      foreach (LayerVolumeBase* mri, mriLayers)
        mri->CloneVoxelByRAS( ras1, ras2, view->GetViewPlane() );
    }
    else if ( m_nAction == EM_Polyline || m_nAction == EM_Livewire )
    {
      double ras[3];
      view->MousePositionToRAS( posX, posY, ras );
      view->GetCursor2D()->SetPosition2( ras );
      if ( m_nAction == EM_Livewire )
      {
        foreach (LayerVolumeBase* mri, mriLayers)
          view->GetCursor2D()->SetInterpolationPoints(
                mri->GetLiveWirePointsByRAS( ras,
                                             view->GetCursor2D()->GetPosition(),
                                             view->GetViewPlane() ) );
      }
      view->GetCursor2D()->SetPosition( view->GetCursor2D()->GetPosition(), true );
      view->RequestRedraw();
    }
    else if (m_nAction == EM_Shift)
    {
      double ras1[3], ras2[3];
      view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
      view->MousePositionToRAS( posX, posY, ras2 );
      for (int i = 0; i < 3; i++)
        ras1[i] = ras2[i] - ras1[i];
      foreach (LayerVolumeBase* mri, mriLayers)
        mri->ShiftVoxelsByRAS(ras1, view->GetViewPlane());
    }
    else if ( m_nAction == EM_Contour )
    {
      LayerMRI* mri_ref = (LayerMRI*)MainWindow::GetMainWindow()->GetBrushProperty()->GetReferenceLayer();
      Contour2D* c2d = view->GetContour2D();
      if ( event->modifiers() & Qt::ShiftModifier )
      {
        double ras1[3], ras2[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
        view->MousePositionToRAS( posX, posY, ras2 );
        c2d->RemoveLine( ras1, ras2 );
      }
      else if ( (event->modifiers() & CONTROL_MODIFIER) && (event->modifiers() & Qt::AltModifier) )
      {
        double scale = 0.2;
        if ( mri_ref )
        {
          double dMin = mri_ref->GetProperty()->GetMinValue();
          double dMax = mri_ref->GetProperty()->GetMaxValue();
          scale = ( dMax - dMin ) * 0.0005;
        }
        c2d->SetContourValue( c2d->GetContourValue() + scale * ( posY - m_nMousePosY ) );
      }
      else
      {
        double ras1[3], ras2[3];
        view->MousePositionToRAS( m_nMousePosX, m_nMousePosY, ras1 );
        view->MousePositionToRAS( posX, posY, ras2 );
        c2d->AddLine( ras1, ras2 );
      }

      view->RequestRedraw();
    }

    if (m_nAction != EM_Shift)
    {
      m_nMousePosX = posX;
      m_nMousePosY = posY;
    }

    return false;
  }
  else
  {
    return Interactor2D::ProcessMouseMoveEvent( event, renderview );
  }
}

bool Interactor2DVolumeEdit::ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  RenderView2D* view = ( RenderView2D* )renderview;
  if ( event->modifiers() & Qt::AltModifier && event->key() == Qt::Key_H )
  {
    Contour2D* c2d = view->GetContour2D();
    c2d->SetVisible( !c2d->IsVisible() );
    view->RequestRedraw();
    return false;
  }
  else if (event->modifiers() & Qt::ShiftModifier && event->key() == Qt::Key_C)
  {
    m_bColorPicking = !m_bColorPicking;
    return false;
  }
  else if (event->key() == Qt::Key_Escape)
  {
    m_bColorPicking = false;
    return false;
  }

  // disable using keyboard to shift voxels
  if (false) // (event->modifiers() & Qt::ShiftModifier) && m_nAction == EM_Shift )
  {
    int nKeyCode = event->key();
    int n[3] = {0, 0, 0};
    int nx = 0, ny = 1;
    int nPlane = view->GetViewPlane();
    int nx_sign = -1;
    switch (nPlane)
    {
    case 0:
      nx = 1;
      ny = 2;
      nx_sign = 1;
      break;
    case 1:
      nx = 0;
      ny = 2;
      break;
    default:
      break;
    }
    if ( nKeyCode == Qt::Key_Up )
    {
      n[ny]+=1;
    }
    else if ( nKeyCode == Qt::Key_Down )
    {
      n[ny]-=1;
    }
    else if ( nKeyCode == Qt::Key_Left )
    {
      n[nx]-=nx_sign;
    }
    else if ( nKeyCode == Qt::Key_Right )
    {
      n[nx]+=nx_sign;
    }
    LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( m_strLayerTypeName );
    LayerVolumeBase* mri = qobject_cast<LayerVolumeBase*>(lc->GetActiveLayer());
    if (mri)
    {
      mri->SaveForUndo(nPlane);
      mri->PrepareShifting(nPlane);
      mri->ShiftVoxels(n, nPlane);
      return false;
    }
  }

  if ( !m_bEditing )
  {
    return Interactor2D::ProcessKeyDownEvent( event, renderview );
  }
  else
  {
    return false;
  }
}

bool Interactor2DVolumeEdit::ProcessKeyUpEvent( QKeyEvent* event, RenderView* renderview )
{
  UpdateCursor( event, renderview );

  return Interactor2D::ProcessKeyUpEvent( event, renderview );
}

void Interactor2DVolumeEdit::UpdateCursor( QEvent* event, QWidget* wnd )
{
  if ( wnd->hasFocus() )
  {
    bool bMouseEvent = false;
    if ( event->type() == QEvent::MouseButtonPress ||
         event->type() == QEvent::MouseButtonRelease ||
         event->type() == QEvent::MouseMove)
    {
      QMouseEvent* e = ( QMouseEvent* )event;
      bMouseEvent = true;
      if ( ( ( e->button() == Qt::MidButton || e->button() == Qt::RightButton ) && !m_bEditing ) ||
           ( ( e->modifiers() & CONTROL_MODIFIER) && (e->modifiers() & Qt::ShiftModifier) ) )
      {
        Interactor2D::UpdateCursor( event, wnd );
        return;
      }
    }

    if ( m_nAction != EM_Fill )
    {
      if ( event->type() == QEvent::KeyPress )
      {
        QKeyEvent* e = ( QKeyEvent* )event;
        if ( e->key() == CONTROL_KEY && !(e->modifiers() & Qt::ShiftModifier) && !(e->modifiers() & Qt::AltModifier) )
        {
          wnd->setCursor( CursorFactory::CursorFill );
          return;
        }
      }

      if ( bMouseEvent && (( QMouseEvent* )event)->modifiers() & CONTROL_MODIFIER
           && !((( QMouseEvent* )event)->modifiers() & Qt::ShiftModifier)
           && !((( QMouseEvent* )event)->modifiers() & Qt::AltModifier) )
      {
        wnd->setCursor( CursorFactory::CursorFill );
      }
      else if ( m_nAction == EM_ColorPicker || m_bColorPicking )
      {
        wnd->setCursor( CursorFactory::CursorColorPicker );
      }
      else //if ( bMouseEvent )
      {
        switch ( m_nAction )
        {
        case EM_Freehand:
          wnd->setCursor( CursorFactory::CursorPencil );
          break;
        case EM_Clone:
          wnd->setCursor( CursorFactory::CursorPencil );
          break;
        case EM_Polyline:
          wnd->setCursor( CursorFactory::CursorPolyline );
          break;
        case EM_Contour:
          wnd->setCursor( CursorFactory::CursorContour );
          break;
        case EM_Shift:
          wnd->setCursor( Qt::OpenHandCursor );
          break;
        }
      }
    }
    else if (m_nAction != EM_GeoSeg)
    {
      wnd->setCursor( CursorFactory::CursorFill );
    }
  }
  else
  {
    Interactor2D::UpdateCursor( event, wnd );
  }
}
