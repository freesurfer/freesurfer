/**
 * @file  RenderView3D.cpp
 * @brief 3D view
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/23 17:35:44 $
 *    $Revision: 1.67 $
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
#include "RenderView3D.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include "SurfaceOverlay.h"
#include "VolumeCropper.h"
#include "Cursor3D.h"
#include "SurfaceRegion.h"
#include "SurfaceROI.h"
#include "SurfaceOverlayProperty.h"
#include <vtkProp.h>
#include <vtkCellPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkActor2D.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkCubeSource.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPlane.h>
#include <vtkCellArray.h>
#include <QtGlobal>
#include <QMenu>
#include <QDebug>
#include "Interactor3DNavigate.h"
#include "Interactor3DMeasure.h"
#include "Interactor3DVolumeCrop.h"
#include "LayerVolumeTrack.h"
#include <vtkScalarBarActor.h>
#include "vtkRGBAColorTransferFunction.h"
#include <vtkAnnotatedCubeActor.h>

#define SLICE_PICKER_PIXEL_TOLERANCE  15

RenderView3D::RenderView3D( QWidget* parent ) : RenderView( parent )
{
  this->GetRenderWindow()->GetInteractor()->SetDesiredUpdateRate(30);
  this->GetRenderWindow()->GetInteractor()->SetStillUpdateRate(0.01);

  m_bShowSlices = true;
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSliceFrames[i] = vtkSmartPointer<vtkActor>::New();
    m_actorSliceFrames[i]->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
    m_actorSliceFrames[i]->GetProperty()->SetRepresentationToWireframe();
    m_actorSliceFrames[i]->GetProperty()->SetDiffuse( 0.0 );
    m_actorSliceFrames[i]->GetProperty()->SetAmbient( 1.0 );
    m_bSliceVisibility[i] = true;

    m_actorSliceBoundingBox[i] = vtkSmartPointer<vtkActor>::New();
    m_cubeSliceBoundingBox[i] = vtkSmartPointer<vtkCubeSource>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_cubeSliceBoundingBox[i]->GetOutputPort() );
    m_actorSliceBoundingBox[i]->SetMapper( mapper );
    mapper->Update();
    m_actorSliceBoundingBox[i]->GetProperty()->SetRepresentationToWireframe();
  }
  HighlightSliceFrame( -1 );

  m_cursor3D = new Cursor3D( this );
  m_interactorNavigate = new Interactor3DNavigate( this );
  m_interactorMeasure = new Interactor3DMeasure( this );
  m_interactorVolumeCrop = new Interactor3DVolumeCrop(this);
  connect(m_cursor3D, SIGNAL(Updated()), this, SLOT(RequestRedraw()));

  m_actorScalarBar->SetNumberOfLabels( 4 );

  SetInteractionMode(IM_Navigate);
}

void RenderView3D::SetInteractionMode( int nMode )
{
  RenderView::SetInteractionMode( nMode );

  switch ( nMode )
  {
  case IM_Measure:
    m_interactor = m_interactorMeasure;
    break;
  case IM_VolumeCrop:
    m_interactor = m_interactorVolumeCrop;
    break;
  default:
    m_interactor = m_interactorNavigate;
    break;
  }
}

void RenderView3D::OnSlicePositionChanged()
{
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  m_cursor3D->SetPosition( lc->GetSlicePosition() );
  UpdateSliceFrames();
  UpdateSurfaceCorrelationData();

  RenderView::OnSlicePositionChanged();
}

void RenderView3D::OnIdle()
{
  if ( m_bToUpdateRASPosition )
  {
    DoUpdateRASPosition( m_nPickCoord[0], m_nPickCoord[1] );
    m_bToUpdateRASPosition = false;
  }
  if ( m_bToUpdateCursorPosition )
  {
    DoUpdateRASPosition( m_nCursorCoord[0], m_nCursorCoord[1], true );
    m_bToUpdateCursorPosition = false;
  }
  if ( m_bToUpdateConnectivity )
  {
    DoUpdateConnectivityDisplay();
    m_bToUpdateConnectivity = false;
  }

  RenderView::OnIdle();
}

void RenderView3D::UpdateSurfaceCorrelationData()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayerCollection( "Surface" )->GetLayers();
  for ( int i = 0; i < layers.size(); i++ )
  {
    ((LayerSurface*)layers[i])->UpdateCorrelationOverlay();
  }
}

void RenderView3D::UpdateSliceFrames()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerCollection* lc = mainwnd->GetLayerCollection( "MRI" );
  if ( lc->IsEmpty() )
  {
    lc = mainwnd->GetLayerCollection( "Surface" );
  }
  double* bounds = m_dBounds;
  double* slicepos = lc->GetSlicePosition();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  points->InsertPoint( 0, slicepos[0], bounds[2], bounds[4] );
  points->InsertPoint( 1, slicepos[0], bounds[2], bounds[5] );
  points->InsertPoint( 2, slicepos[0], bounds[3], bounds[5] );
  points->InsertPoint( 3, slicepos[0], bounds[3], bounds[4] );
  vtkIdType ids[5] = { 0, 1, 2, 3, 0 };
  lines->InsertNextCell( 5, ids );
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[0]->GetMapper())->SetInput( polydata );

  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], slicepos[1], bounds[4] );
  points->InsertPoint( 1, bounds[0], slicepos[1], bounds[5] );
  points->InsertPoint( 2, bounds[1], slicepos[1], bounds[5] );
  points->InsertPoint( 3, bounds[1], slicepos[1], bounds[4] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[1]->GetMapper())->SetInput( polydata );

  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], bounds[2], slicepos[2] );
  points->InsertPoint( 1, bounds[0], bounds[3], slicepos[2] );
  points->InsertPoint( 2, bounds[1], bounds[3], slicepos[2] );
  points->InsertPoint( 3, bounds[1], bounds[2], slicepos[2] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[2]->GetMapper())->SetInput( polydata );

  for ( int i = 0; i < 3; i++ )
  {
    double bounds[6];
    for ( int n = 0; n < 3; n++ )
    {
      bounds[n*2]   = m_dBounds[n*2] - m_dBoundingTolerance;
      bounds[n*2+1] = m_dBounds[n*2+1] + m_dBoundingTolerance;
    }

    bounds[i*2] = slicepos[i] - m_dBoundingTolerance;
    bounds[i*2+1] = slicepos[i] + m_dBoundingTolerance;
    m_cubeSliceBoundingBox[i]->SetBounds( bounds );
    m_actorSliceBoundingBox[i]->GetMapper()->Update();
  }

  RequestRedraw();
}

void RenderView3D::UpdateViewByWorldCoordinate()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double wcenter[3];
  for ( int i = 0; i < 3; i++ )
  {
    wcenter[i] = m_dWorldOrigin[i] + m_dWorldSize[i] / 2;
  }
  cam->SetFocalPoint( wcenter );
  cam->SetPosition( wcenter[0] - ( m_dWorldSize[1] > m_dWorldSize[2] ? m_dWorldSize[1] : m_dWorldSize[2] ) *2.5,
                    wcenter[1],
                    wcenter[2]);
  cam->SetViewUp( 0, 0, 1 );
  m_renderer->ResetCameraClippingRange();
}

void RenderView3D::RefreshAllActors(bool bForScreenShot)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  SettingsScreenshot setting = mainwnd->GetScreenShotSettings();

  m_renderer->RemoveAllViewProps();
  bool b[3] = { m_bShowSlices, m_bShowSlices, m_bShowSlices };
  mainwnd->GetLayerCollection( "MRI" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "ROI" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Surface" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "PointSet" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Track" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Supplement" )->Append3DProps( m_renderer, b );

  if (!mainwnd->IsEmpty())
  {
    if (!bForScreenShot || !setting.HideCursor)
      m_cursor3D->AppendActor( m_renderer );
  }

  // add focus frame
  if (!bForScreenShot)
  {
    m_renderer->AddViewProp( m_actorFocusFrame );
  }
  if ( !mainwnd->GetLayerCollection("MRI")->IsEmpty() ||! mainwnd->GetLayerCollection("Surface")->IsEmpty() )
  {
    m_renderer->AddViewProp( m_actorScalarBar );
    if ( !mainwnd->GetLayerCollection("MRI")->IsEmpty() && m_bShowSlices )
    {
      for ( int i = 0; i < 3; i++ )
      {
        m_renderer->AddViewProp( m_actorSliceFrames[i] );
      }
    }
  }

//    mainwnd->GetConnectivityData()->AppendProps( m_renderer );
  mainwnd->GetVolumeCropper()->Append3DProps( m_renderer );

  m_renderer->ResetCameraClippingRange();
  RenderView::RefreshAllActors(bForScreenShot);
}

void RenderView3D::DoUpdateRASPosition( int posX, int posY, bool bCursor )
{
  LayerCollection* lc_mri = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  LayerCollection* lc_roi = MainWindow::GetMainWindow()->GetLayerCollection( "ROI" );
  LayerCollection* lc_surface = MainWindow::GetMainWindow()->GetLayerCollection( "Surface" );

  this->setToolTip("");
  if ( lc_mri->IsEmpty() && lc_roi->IsEmpty() && lc_surface->IsEmpty() )
  {
    return;
  }

// MousePositionToRAS( posX, posY, pos );
// vtkPointPicker* picker = vtkPointPicker::SafeDownCast( this->GetPicker() );
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetRenderWindow()->GetInteractor()->GetPicker() );
// vtkPropPicker* picker = vtkPropPicker::SafeDownCast( this->GetPicker() );
  if ( picker )
  {
    picker->InitializePickList();

    vtkPropCollection* props = GetRenderer()->GetViewProps();
    if ( props )
    {
      props->InitTraversal();
      vtkProp* prop = props->GetNextProp();
      while ( prop )
      {
        if ( vtkActor::SafeDownCast( prop ) )
        {
          picker->AddPickList( prop );
        }
        prop = props->GetNextProp();
      }
    }
    // add bounding box for slice frame picking
    for ( int i = 0; i < 3; i++ )
    {
      picker->AddPickList( m_actorSliceBoundingBox[i] );
    }

    double pos[3];
    picker->Pick( posX, rect().height() - posY, 0, GetRenderer() );
    picker->GetPickPosition( pos );

    vtkProp* prop = picker->GetViewProp();
    if ( !prop )
    {
      HighlightSliceFrame( -1 );
      return;
    }

    // check slice frame selection first
    bool bFramePicked = false;
    double tolerance = m_dBoundingTolerance * 1.414;
    for ( int i = 0; i < 3; i++ )
    {
      if ( m_actorSliceBoundingBox[i].GetPointer() == prop )
      {
        if ( fabs( pos[0] - m_dBounds[0] ) < tolerance ||
             fabs( pos[0] - m_dBounds[1] ) < tolerance ||
             fabs( pos[1] - m_dBounds[2] ) < tolerance ||
             fabs( pos[1] - m_dBounds[3] ) < tolerance ||
             fabs( pos[2] - m_dBounds[4] ) < tolerance ||
             fabs( pos[2] - m_dBounds[5] ) < tolerance
           )
        {
          // 3D hit test passed, now we check screen distance
          double screen_pt[3];
          double screen_pts[4][3];
          int x, y;
          this->WorldToScreen( pos[0], pos[1], pos[2], x, y );
          screen_pt[0] = x;
          screen_pt[1] = y;
          screen_pt[2] = 0;
          vtkPoints* pts = vtkPolyDataMapper::SafeDownCast( m_actorSliceFrames[i]->GetMapper() )->GetInput()->GetPoints();
          for ( int j = 0; j < 4; j++ )
          {
            double* p = pts->GetPoint( j );
            this->WorldToScreen( p[0], p[1], p[2], x, y );
            screen_pts[j][0] = x;
            screen_pts[j][1] = y;
            screen_pts[j][2] = 0;
          }
          int ids[4][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
          double dMinDist = 1000000000;
          for ( int j = 0; j < 4; j++ )
          {
            double dist = vtkLine::DistanceToLine( screen_pt, screen_pts[ids[j][0]], screen_pts[ids[j][1]]);
            if ( dist < dMinDist )
            {
              dMinDist = dist;
            }
          }
          if ( dMinDist < SLICE_PICKER_PIXEL_TOLERANCE )
          {
            HighlightSliceFrame( i );
            m_dIntersectPoint[0] = pos[0];
            m_dIntersectPoint[1] = pos[1];
            m_dIntersectPoint[2] = pos[2];
            bFramePicked = true;
            break;
          }
        }
      }
    }

    if ( !bFramePicked )
    {
      //  if ( !lc_surface->IsEmpty() && !lc_surface->HasProp( prop ) )
      {
        for ( int i = 0; i < 3; i++ )
        {
          picker->DeletePickList( m_actorSliceBoundingBox[i] );
        }

        picker->Pick( posX, rect().height() - posY, 0, GetRenderer() );
        picker->GetPickPosition( pos );
        prop = picker->GetViewProp();
      }

      if ( lc_mri->HasProp( prop ) || lc_roi->HasProp( prop ) )
      {
        if ( bCursor )
        {
          LayerMRI* mri = (LayerMRI*)lc_mri->HasProp( prop );
          SurfaceRegion* reg = NULL;
          if ( mri )
          {
            reg = mri->SelectSurfaceRegion( pos );
          }
          if ( reg )
          {
            RequestRedraw( true ); // force redraw
            emit SurfaceRegionSelected(reg);
          }
        }
        else
        {
          LayerVolumeTrack* vt = qobject_cast<LayerVolumeTrack*>(lc_mri->HasProp( prop ));
          if (vt)
          {
            QVariantMap info = vt->GetLabelByProp(prop);
            if (!info.isEmpty())
            {
              this->setToolTip(QString("%1 %2").arg(info["label"].toInt()).arg(info["name"].toString()));
              emit VolumeTrackMouseOver(vt, info);
            }
          }
        }
      }
      else if ( Layer* layer = lc_surface->HasProp( prop ) )
      {
        if ( bCursor )
        {
          lc_mri->SetCursorRASPosition( pos );
          MainWindow::GetMainWindow()->SetSlicePosition( pos );
          if (layer)
          {
            lc_surface->SetActiveLayer(layer);
          }
          emit SurfaceVertexClicked();
        }
        else
        {
          lc_mri->SetCurrentRASPosition( pos );
        }
      }

      HighlightSliceFrame( -1 );
    }
  }
}

void RenderView3D::DoUpdateConnectivityDisplay()
{
  /*
  ConnectivityData* conn = MainWindow::GetMainWindowPointer()->GetConnectivityData();
  if ( conn->IsValid() && conn->GetDisplayMode() != ConnectivityData::DM_All )
  {
  conn->BuildConnectivityActors();
  }
  */
}

void RenderView3D::HighlightSliceFrame( int n )
{
  if ( m_nSliceHighlighted == n )
  {
    return;
  }

  double colors[][3] = { { 1, 0.1, 0.1}, { 0.1, 1, 0.1 }, { 0.1, 0.1, 1 } };
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSliceFrames[i]->GetProperty()->SetLineWidth( 2 );
    m_actorSliceFrames[i]->GetProperty()->SetColor( colors[i] );
  }
  if ( n >= 0 && n <= 2 )
  {
    m_actorSliceFrames[n]->GetProperty()->SetLineWidth( 4 );
    m_actorSliceFrames[n]->GetProperty()->SetColor( 1, 1, 1 );
  }
  m_nSliceHighlighted = n;
  if ( m_actorSliceFrames[0]->GetMapper()->GetInput() )
  {
    RequestRedraw();
  }
}

bool RenderView3D::GetShowSliceFrames()
{
  return m_actorSliceFrames[0]->GetVisibility();
}

void RenderView3D::SetShowSliceFrames( bool bShow )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSliceFrames[i]->SetVisibility( bShow?1:0 );
    m_actorSliceBoundingBox[i]->SetPickable( bShow?1:0 );
  }

  RequestRedraw();
}

void RenderView3D::SetShowSlices(bool bShow)
{
  if (bShow != m_bShowSlices)
  {
    m_bShowSlices = bShow;
    RefreshAllActors();
  }
}

void RenderView3D::UpdateCursorRASPosition( int posX, int posY )
{
  m_bToUpdateCursorPosition = true;
  m_nCursorCoord[0] = posX;
  m_nCursorCoord[1] = posY;
}

void RenderView3D::UpdateMouseRASPosition( int posX, int posY )
{
  m_bToUpdateRASPosition = true;
  m_nPickCoord[0] = posX;
  m_nPickCoord[1] = posY;
}

void RenderView3D::MoveSliceToScreenCoord( int x, int y )
{
  if ( m_nSliceHighlighted < 0 )
  {
    return;
  }

  MainWindow* mainwnd = MainWindow::GetMainWindow();
  LayerCollection* lc = mainwnd->GetLayerCollection( "MRI" );
  if ( lc->IsEmpty() )
  {
    lc = mainwnd->GetLayerCollection( "Surface" );
  }
  double* bounds = m_dBounds;
  double slicepos[3];
  lc->GetSlicePosition( slicepos );
  double pt[3];
  pt[0] = m_dIntersectPoint[0];
  pt[1] = m_dIntersectPoint[1];
  pt[2] = m_dIntersectPoint[2];

  double v[3] = { 0, 0, 0 };
  switch ( m_nSliceHighlighted )
  {
  case 0:
    if ( qMin( fabs( pt[1] - bounds[2] ), fabs( pt[1] - bounds[3] ) ) <
         qMin( fabs( pt[2] - bounds[4] ), fabs( pt[2] - bounds[5] ) ) )
    {
      v[1] = 1;
      if ( fabs( pt[1] - bounds[2] ) < fabs( pt[1] - bounds[3] ) )
      {
        pt[1] = bounds[2];
      }
      else
      {
        pt[1] = bounds[3];
      }
    }
    else
    {
      v[2] = 1;
      if ( fabs( pt[2] - bounds[4] ) < fabs( pt[2] - bounds[5] ) )
      {
        pt[2] = bounds[4];
      }
      else
      {
        pt[2] = bounds[5];
      }
    }
    break;
  case 1:
    if ( qMin( fabs( pt[0] - bounds[0] ), fabs( pt[0] - bounds[1] ) ) <
         qMin( fabs( pt[2] - bounds[4] ), fabs( pt[2] - bounds[5] ) ) )
    {
      v[0] = 1;
      if ( fabs( pt[0] - bounds[0] ) < fabs( pt[0] - bounds[1] ) )
      {
        pt[0] = bounds[0];
      }
      else
      {
        pt[0] = bounds[1];
      }
    }
    else
    {
      v[2] = 1;
      if ( fabs( pt[2] - bounds[4] ) < fabs( pt[2] - bounds[5] ) )
      {
        pt[2] = bounds[4];
      }
      else
      {
        pt[2] = bounds[5];
      }
    }
    break;
  case 2:
    if ( qMin( fabs( pt[0] - bounds[0] ), fabs( pt[0] - bounds[1] ) ) <
         qMin( fabs( pt[1] - bounds[2] ), fabs( pt[1] - bounds[3] ) ) )
    {
      v[0] = 1;
      if ( fabs( pt[0] - bounds[0] ) < fabs( pt[0] - bounds[1] ) )
      {
        pt[0] = bounds[0];
      }
      else
      {
        pt[0] = bounds[1];
      }
    }
    else
    {
      v[1] = 1;
      if ( fabs( pt[1] - bounds[2] ) < fabs( pt[1] - bounds[3] ) )
      {
        pt[1] = bounds[2];
      }
      else
      {
        pt[1] = bounds[3];
      }
    }
    break;
  }
  pt[m_nSliceHighlighted] = slicepos[m_nSliceHighlighted];

  double pt1[3], pt2[3];
  this->ScreenToWorld( x, y, -100, pt1[0], pt1[1], pt1[2] );
  this->ScreenToWorld( x, y, 100, pt2[0], pt2[1], pt2[2] );
  double new_pt[3], t = 0;
  vtkPlane::IntersectWithLine( pt1, pt2, v, pt, t, new_pt );
  if ( t > 100000 )
  {
    new_pt[0] = pt[0];
    new_pt[1] = pt[1];
    new_pt[2] = pt[2];
  }
  if ( new_pt[m_nSliceHighlighted] < bounds[m_nSliceHighlighted*2] )
  {
    new_pt[m_nSliceHighlighted] = bounds[m_nSliceHighlighted*2];
  }
  else if ( new_pt[m_nSliceHighlighted] > bounds[m_nSliceHighlighted*2+1] )
  {
    new_pt[m_nSliceHighlighted] = bounds[m_nSliceHighlighted*2+1];
  }
  mainwnd->OffsetSlicePosition( m_nSliceHighlighted, new_pt[m_nSliceHighlighted] - slicepos[m_nSliceHighlighted], false );
  slicepos[m_nSliceHighlighted] = new_pt[m_nSliceHighlighted];
  lc->SetCursorRASPosition( slicepos );
}

vtkProp* RenderView3D::PickProp( int posX, int posY, double* pos_out )
{
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetRenderWindow()->GetInteractor()->GetPicker() );
  if ( !picker )
  {
    return NULL;
  }

  picker->InitializePickList();
  vtkPropCollection* props = GetRenderer()->GetViewProps();
  if ( props )
  {
    props->InitTraversal();
    vtkProp* prop = props->GetNextProp();
    while ( prop )
    {
      if ( vtkActor::SafeDownCast( prop ) )
      {
        picker->AddPickList( prop );
      }
      prop = props->GetNextProp();
    }
  }
  picker->Pick( posX, rect().height() - posY, 0, GetRenderer() );
  if ( pos_out )
  {
    picker->GetPickPosition( pos_out );
  }
  return picker->GetViewProp();
}

bool RenderView3D::InitializeSelectRegion( int posX, int posY )
{
  double pos[3];
  vtkProp* prop = this->PickProp( posX, posY, pos );
  if ( !prop )
  {
    return false;
  }

  LayerCollection* lc_mri = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  LayerMRI* mri = NULL;
  for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri_temp = (LayerMRI*)lc_mri->GetLayer(i);
    if ( mri_temp->HasProp( prop ) && mri_temp->GetProperty()->GetShowAsContour() )
    {
      mri = mri_temp;
      break;
    }
  }

  if ( !mri )
  {
    return false;
  }

  lc_mri->SetActiveLayer( mri );
  SurfaceRegion* reg = mri->CreateNewSurfaceRegion( pos );
  if ( reg )
  {
    emit SurfaceRegionSelected(reg);
  }
  return true;
}

bool RenderView3D::PickSelectRegion( int nId )
{
  LayerCollection* lc_mri = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  LayerMRI* mri = NULL;
  for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri_temp = (LayerMRI*)lc_mri->GetLayer(i);
    if ( mri_temp->GetProperty()->GetShowAsContour() )
    {
      mri = mri_temp;
      break;
    }
  }

  if ( !mri )
  {
    return false;
  }

  SurfaceRegion* reg = mri->SelectSurfaceRegion( nId );
  if ( reg )
  {
    emit SurfaceRegionSelected(reg);
  }
  return true;
}

void RenderView3D::AddSelectRegionLoopPoint( int posX, int posY )
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( mri )
  {
    double pos[3];
    vtkProp* prop = this->PickProp( posX, posY, pos );
    if ( !prop || !mri->HasProp( prop ) )
    {
      return;
    }

    mri->AddSurfaceRegionLoopPoint( pos );
  }
}

void RenderView3D::CloseSelectRegion()
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( mri )
  {
    mri->CloseSurfaceRegion();
  }
}

void RenderView3D::DeleteCurrentSelectRegion()
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindow()->GetActiveLayer( "MRI" );
  if ( mri )
  {
    SurfaceRegion* reg = mri->GetCurrentSurfaceRegion();
    if ( mri->DeleteCurrentSurfaceRegion() )
    {
      emit SurfaceRegionRemoved(reg);
    }
  }
}

SurfaceROI* RenderView3D::InitializeSurfaceROI( int posX, int posY )
{
  double pos[3];
  vtkProp* prop = this->PickProp( posX, posY, pos );
  if ( !prop )
  {
    return NULL;
  }

  LayerCollection* lc_surf = MainWindow::GetMainWindow()->GetLayerCollection( "Surface" );
  LayerSurface* surf = NULL;
  for ( int i = 0; i < lc_surf->GetNumberOfLayers(); i++ )
  {
    LayerSurface* temp = (LayerSurface*)lc_surf->GetLayer(i);
    if ( temp->HasProp( prop ) )
    {
      surf = temp;
      break;
    }
  }

  if ( !surf )
  {
    return false;
  }

  lc_surf->SetActiveLayer( surf );
  SurfaceROI* roi = surf->GetSurfaceROI();
  if ( roi )
  {
    roi->InitializeOutline(pos);
  }
  return roi;
}

void RenderView3D::AddSurfaceROIPoint( int posX, int posY )
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if ( surf )
  {
    double pos[3];
    vtkProp* prop = this->PickProp( posX, posY, pos );
    if ( !prop || !surf->HasProp( prop ) )
    {
      return;
    }

    surf->GetSurfaceROI()->AddPoint( pos );
    RequestRedraw();
  }
}

bool RenderView3D::UpdateBounds()
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  double bounds[6] = { 1000000, -1000000, 1000000, -1000000, 1000000, -1000000 };
  for ( int n = 0; n < 1; n++ )
  {
    LayerCollection* lc = mainwnd->GetLayerCollection( (n == 0 ? "MRI" : "Surface") );
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
    {
      double bd[6];
      lc->GetLayer( i )->GetDisplayBounds( bd );
      for ( int j = 0; j < 3; j++ )
      {
        if ( bounds[j*2] > bd[j*2] )
        {
          bounds[j*2] = bd[j*2];
        }
        if ( bounds[j*2+1] < bd[j*2+1] )
        {
          bounds[j*2+1] = bd[j*2+1];
        }
      }
    }
  }
  for ( int i = 0; i < 6; i++ )
  {
    m_dBounds[i] = bounds[i];
  }

  double dMaxLength = 0;
  for ( int i = 0; i < 3; i++ )
  {
    if ( dMaxLength < ( bounds[i*2+1]-bounds[i*2] ) )
    {
      dMaxLength = bounds[i*2+1]-bounds[i*2];
    }
  }

  m_dBoundingTolerance = dMaxLength * 0.02;
  UpdateSliceFrames();

  return true;
}


bool RenderView3D::PickCroppingBound( int nX, int nY )
{
  vtkProp* prop = PickProp( nX, nY );
  if ( prop && MainWindow::GetMainWindow()->GetVolumeCropper()->PickActiveBound( prop ) )
  {
    RequestRedraw( true );
    return true;
  }
  else
  {
    return false;
  }
}

void RenderView3D::MoveCroppingBound( int nX, int nY )
{
  MainWindow::GetMainWindow()->GetVolumeCropper()
  ->MoveActiveBound( this, nX, nY );
}

void RenderView3D::SnapToNearestAxis()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3], v_up[3];

  cam->OrthogonalizeViewUp();

  cam->GetDirectionOfProjection(v);

  cam->GetViewUp(v_up);

  double wcenter[3];
  for ( int i = 0; i < 3; i++ )
  {
    wcenter[i] = m_dWorldOrigin[i] + m_dWorldSize[i] / 2;
  }
  cam->SetFocalPoint( wcenter );

  if ( fabs(v[0]) > fabs(v[1]) && fabs(v[0]) > fabs(v[2]) )
  {
    v[0] = ( v[0] > 0 ? 1 : -1 );
    v[1] = v[2] = 0;
  }
  else if ( fabs(v[1]) > fabs(v[2]) )
  {
    v[1] = ( v[1] > 0 ? 1 : -1 );
    v[0] = v[2] = 0;
  }
  else
  {
    v[2] = ( v[2] > 0 ? 1 : -1 );
    v[0] = v[1] = 0;
  }

  if ( fabs(v_up[0]) > fabs(v_up[1]) && fabs(v_up[0]) > fabs(v_up[2]) )
  {
    v_up[0] = ( v_up[0] > 0 ? 1 : -1 );
    v_up[1] = v_up[2] = 0;
  }
  else if ( fabs(v_up[1]) > fabs(v_up[2]) )
  {
    v_up[1] = ( v_up[1] > 0 ? 1 : -1 );
    v_up[0] = v_up[2] = 0;
  }
  else
  {
    v_up[2] = ( v_up[2] > 0 ? 1 : -1 );
    v_up[0] = v_up[1] = 0;
  }

  double pos[3];
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = wcenter[i] - ( m_dWorldSize[i] * v[i] * 2.5 );
  }
  cam->SetPosition( pos );
  cam->SetViewUp( v_up );
  cam->OrthogonalizeViewUp();
  m_renderer->ResetCameraClippingRange();

  RequestRedraw();
}

void RenderView3D::UpdateScalarBar()
{
  /*
  LayerSurface* surf = (LayerSurface*) MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  if ( surf && surf->GetActiveOverlay() )
  {
    m_actorScalarBar->SetLookupTable( surf->GetActiveOverlay()->GetProperty()->GetLookupTable() );
  }
  else
  */
  {
    RenderView::UpdateScalarBar();
  }
}

void RenderView3D::TriggerContextMenu( QMouseEvent* event )
{
  QMenu* menu = new QMenu(this);
  menu->addAction(MainWindow::GetMainWindow()->ui->actionShowSlices);
  menu->addAction(MainWindow::GetMainWindow()->ui->actionShowSliceFrames);
  bool bShowBar = this->GetShowScalarBar();
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainwnd->GetLayers("Surface");
  layers << mainwnd->GetLayers("MRI");
  layers << mainwnd->GetLayers("PointSet");
  if (!layers.isEmpty())
  {
    QMenu* menu2 = menu->addMenu("Show Color Bar");
    QActionGroup* ag = new QActionGroup(this);
    ag->setExclusive(false);
    foreach (Layer* layer, layers)
    {
      QAction* act = new QAction(layer->GetName(), this);
      act->setCheckable(true);
      act->setChecked(bShowBar && layer == m_layerScalarBar);
      act->setData(QVariant::fromValue((QObject*)layer));
      menu2->addAction(act);
      ag->addAction(act);
    }
    connect(ag, SIGNAL(triggered(QAction*)), this, SLOT(SetScalarBarLayer(QAction*)));
  }
  LayerMRI* mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (mri && mri->GetProperty()->GetShowAsContour())
  {
    menu->addSeparator();
    QAction* act = new QAction("Save IsoSurface As...", this);
    menu->addAction(act);
    connect(act, SIGNAL(triggered()), mainwnd, SLOT(OnSaveIsoSurface()));
  }
  menu->exec(event->globalPos());
}
