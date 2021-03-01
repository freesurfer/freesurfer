/**
 * @brief 3D view
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
#include "Interactor3DROIEdit.h"
#include "LayerVolumeTrack.h"
#include <vtkScalarBarActor.h>
#include "vtkRGBAColorTransferFunction.h"
#include <vtkAnnotatedCubeActor.h>
#include <vtkCubeAxesActor.h>
#include <vtkTextProperty.h>
#include <QFileInfo>
#include "MyUtils.h"
#include "FSSurface.h"
#include "Interactor3DPathEdit.h"
#include <QElapsedTimer>
#include "vtkInteractorStyleMyTrackballCamera.h"
#include <vtkCubeAxesActor.h>

#define SLICE_PICKER_PIXEL_TOLERANCE  15

RenderView3D::RenderView3D( QWidget* parent ) : RenderView( parent )
{
  m_interactorStyle = vtkSmartPointer<vtkInteractorStyleMyTrackballCamera>::New();
  this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(m_interactorStyle);
  this->GetRenderWindow()->GetInteractor()->SetDesiredUpdateRate(30);
  this->GetRenderWindow()->GetInteractor()->SetStillUpdateRate(0.01);

  GetRenderWindow()->SetAlphaBitPlanes(0);
  GetRenderWindow()->SetMultiSamples(0);
  GetRenderer()->SetUseDepthPeeling(true);
  GetRenderer()->SetMaximumNumberOfPeels(4);
  GetRenderer()->SetOcclusionRatio(0);

  m_bShowSliceFrames = true;
  m_bShowAxes = false;
  m_bShowCursor = true;
  m_bFocalPointAtCursor = false;
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
  m_interactorROIEdit = new Interactor3DROIEdit(this);
  m_interactorPathEdit = new Interactor3DPathEdit(this);
  connect(m_cursor3D, SIGNAL(Updated()), this, SLOT(RequestRedraw()));

  m_cursorInflatedSurf = new Cursor3D(this);
  m_cursorInflatedSurf->Hide();
  connect(m_cursorInflatedSurf, SIGNAL(Updated()), this, SLOT(RequestRedraw()));

  m_actorScalarBar->SetNumberOfLabels( 4 );

  m_actorAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
  m_actorAxesActor->XAxisLabelVisibilityOn();
  m_actorAxesActor->YAxisLabelVisibilityOn();
  m_actorAxesActor->ZAxisLabelVisibilityOn();
  m_actorAxesActor->SetXTitle("");
  m_actorAxesActor->SetYTitle("");
  m_actorAxesActor->SetZTitle("");
  m_actorAxesActor->SetFlyModeToClosestTriad();
  m_actorAxesActor->SetCamera(m_renderer->GetActiveCamera());

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
  case IM_ROIEdit:
    m_interactor = m_interactorROIEdit;
    break;
  case IM_SurfacePath:
  case IM_SurfaceCut:
    m_interactor = m_interactorPathEdit;
    break;
  default:
    m_interactor = m_interactorNavigate;
    break;
  }
}

void RenderView3D::OnSlicePositionChanged(bool bCenter)
{
  Q_UNUSED(bCenter);
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  //  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
  //  m_cursorInflatedSurf->Show(m_cursor3D->IsShown() && surf && surf->IsInflated());
  m_cursor3D->SetPosition( lc->GetSlicePosition() );
  UpdateSliceFrames();
  UpdateSurfaceCorrelationData();
  if (m_bFocalPointAtCursor)
    m_interactorStyle->SetRotateByPoint(true, lc->GetSlicePosition());

  RenderView::OnSlicePositionChanged();
}

void RenderView3D::OnIdle()
{
  if ( m_bToUpdateRASPosition )
  {
    if (QApplication::mouseButtons() == Qt::NoButton)
      DoUpdateRASPosition( m_nPickCoord[0], m_nPickCoord[1], false, m_bSlicePickOnly);
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
#if VTK_MAJOR_VERSION > 5
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[0]->GetMapper())->SetInputData( polydata );
#else
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[0]->GetMapper())->SetInput( polydata );
#endif

  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], slicepos[1], bounds[4] );
  points->InsertPoint( 1, bounds[0], slicepos[1], bounds[5] );
  points->InsertPoint( 2, bounds[1], slicepos[1], bounds[5] );
  points->InsertPoint( 3, bounds[1], slicepos[1], bounds[4] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
#if VTK_MAJOR_VERSION > 5
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[1]->GetMapper())->SetInputData( polydata );
#else
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[1]->GetMapper())->SetInput( polydata );
#endif

  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], bounds[2], slicepos[2] );
  points->InsertPoint( 1, bounds[0], bounds[3], slicepos[2] );
  points->InsertPoint( 2, bounds[1], bounds[3], slicepos[2] );
  points->InsertPoint( 3, bounds[1], bounds[2], slicepos[2] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
#if VTK_MAJOR_VERSION > 5
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[2]->GetMapper())->SetInputData( polydata );
#else
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[2]->GetMapper())->SetInput( polydata );
#endif

  for ( int i = 0; i < 3; i++ )
  {
    double bds[6];
    for ( int n = 0; n < 3; n++ )
    {
      bds[n*2]   = m_dBounds[n*2] - m_dBoundingTolerance;
      bds[n*2+1] = m_dBounds[n*2+1] + m_dBoundingTolerance;
    }

    bds[i*2] = slicepos[i] - m_dBoundingTolerance;
    bds[i*2+1] = slicepos[i] + m_dBoundingTolerance;
    m_cubeSliceBoundingBox[i]->SetBounds( bds );
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
  cam->SetPosition( wcenter[0] - qMax(m_dWorldSize[1], m_dWorldSize[2]) *2.5,
      wcenter[1],
      wcenter[2]);
  cam->SetViewUp( 0, 0, 1 );
  m_renderer->GetActiveCamera()->SetViewAngle(30);
  m_renderer->ResetCameraClippingRange();
}

void RenderView3D::ResetViewAnterior()
{
  Reset();
  Azimuth(-90);
}

void RenderView3D::ResetViewPosterior()
{
  Reset();
  Azimuth(90);
}

void RenderView3D::ResetViewInferior()
{
  ResetViewAnterior();
  Elevation(-90);
}

void RenderView3D::ResetViewSuperior()
{
  ResetViewPosterior();
  Elevation(90);
}

void RenderView3D::ResetViewLeft()
{
  Reset();
}

void RenderView3D::ResetViewRight()
{
  Reset();
  Azimuth(180);
}

void RenderView3D::ResetViewLateral()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("Surface");
  foreach (Layer* layer, layers)
  {
    if (layer->IsVisible())
    {
      if (((LayerSurface*)layer)->GetHemisphere() == 0)
        ResetViewLeft();
      else
        ResetViewRight();
      break;
    }
  }
}

void RenderView3D::ResetViewMedial()
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("Surface");
  foreach (Layer* layer, layers)
  {
    if (layer->IsVisible())
    {
      if (((LayerSurface*)layer)->GetHemisphere() == 0)
        ResetViewRight();
      else
        ResetViewLeft();
      break;
    }
  }
}

void RenderView3D::RefreshAllActors(bool bForScreenShot)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  SettingsScreenshot setting = mainwnd->GetScreenShotSettings();

  m_renderer->RemoveAllViewProps();
  bool* b = m_bSliceVisibility;
  mainwnd->GetLayerCollection( "MRI" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "ROI" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Surface" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "ODF" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "PointSet" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Tract" )->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "CMAT")->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "FCD")->Append3DProps( m_renderer, b );
  mainwnd->GetLayerCollection( "Supplement" )->Append3DProps( m_renderer, b );

  if (!mainwnd->IsEmpty())
  {
    if (!bForScreenShot || !setting.HideCursor)
    {
      m_cursor3D->AppendActor( m_renderer );
      m_cursorInflatedSurf->AppendActor(m_renderer);
    }
  }

  // add focus frame
  if (!bForScreenShot)
  {
    m_renderer->AddViewProp( m_actorFocusFrame );
  }
  if ( !mainwnd->GetLayerCollection("MRI")->IsEmpty() ||! mainwnd->GetLayerCollection("Surface")->IsEmpty() )
  {
    m_renderer->AddViewProp( m_actorScalarBar );
    if ( !mainwnd->GetLayerCollection("MRI")->IsEmpty() )
    {
      for ( int i = 0; i < 3; i++ )
      {
        m_renderer->AddViewProp( m_actorSliceFrames[i] );
      }
    }
  }

  if (m_bShowAxes)
    m_renderer->AddViewProp(m_actorAxesActor);

  //    mainwnd->GetConnectivityData()->AppendProps( m_renderer );
  mainwnd->GetVolumeCropper()->Append3DProps( m_renderer );


  m_renderer->ResetCameraClippingRange();
  RenderView::RefreshAllActors(bForScreenShot);
}

void RenderView3D::DoUpdateRASPosition( int posX, int posY, bool bCursor, bool bSlicePickOnly )
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

    if (!bSlicePickOnly)
    {
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
    }
    // add bounding box for slice frame picking
    for ( int i = 0; i < 3; i++ )
    {
      picker->AddPickList( m_actorSliceBoundingBox[i] );
    }

    double pos[3];
    posY = this->rect().height() - posY;
  #if VTK_MAJOR_VERSION > 7
    if (devicePixelRatio() > 1)
    {
        posX *= devicePixelRatio();
        posY *= devicePixelRatio();
    }
  #endif
    picker->Pick( posX, posY, 0, GetRenderer() );
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

    if ( !bFramePicked)
      HighlightSliceFrame( -1 );

    if ( !bFramePicked && !bSlicePickOnly )
    {
      for ( int i = 0; i < 3; i++ )
      {
        picker->DeletePickList( m_actorSliceBoundingBox[i] );
      }

      picker->Pick( posX, posY, 0, GetRenderer() );
      picker->GetPickPosition( pos );
      prop = picker->GetViewProp();

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
        if (layer)
        {
          if (bCursor)
            lc_surface->SetActiveLayer(layer);
          LayerSurface* surf = (LayerSurface*)layer;
          QVariantMap settings = MainWindow::GetMainWindow()->GetGeneralSettings();

          MapInflatedCoords(surf, pos, pos,
                            ((!bCursor) || GetInteractionMode() == IM_ROIEdit)?false:settings["AutoReorientView"].toBool(),
                            bCursor);
        }

        if(bCursor)
        {
          lc_mri->SetCursorRASPosition( pos );
          MainWindow::GetMainWindow()->SetSlicePosition( pos );
          if (layer)
            emit SurfaceVertexClicked((LayerSurface*)layer);
        }
        else
          lc_mri->SetCurrentRASPosition( pos );
      }
    }
  }
}

bool RenderView3D::MapInflatedCoords(LayerSurface *surf, double *pos_in, double *pos_out, bool bAutoOrient, bool bCursor)
{
  int nVertex = surf->GetVertexIndexAtTarget(pos_in, NULL);
  if (bCursor)
    surf->SetCurrentVertex(nVertex);
  else
    surf->SetMouseVertex(nVertex);
  if (QFileInfo(surf->GetFileName()).fileName().contains("inflated") || surf->GetActiveSurface() == FSSurface::SurfaceInflated)
  {
    if (m_cursor3D->IsShown())
      m_cursorInflatedSurf->Show();
    if (bCursor)
      m_cursorInflatedSurf->SetPosition(pos_in);
    FSSurface* fsurf = surf->GetSourceSurface();
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("Surface");
    bool bFoundMappingSurface = false;
    foreach (Layer* s, layers)
    {
      LayerSurface* f = (LayerSurface*)s;
      if (f != surf && QFileInfo(f->GetFileName()).fileName().contains(surf->GetMappingSurfaceName()))
      {
        if (f->GetHemisphere() == surf->GetHemisphere() &&
            QFileInfo(f->GetFileName()).absolutePath() == QFileInfo(surf->GetFileName()).absolutePath())
        {
          f->SetCurrentVertex(nVertex);
          f->GetTargetAtVertex(nVertex, pos_out);
          if (bAutoOrient)
          {
            double v[3];
            surf->GetSmoothedVertexNormal(nVertex, v);
            this->AlignViewToNormal(v);
          }
          bFoundMappingSurface = true;
          return true;
        }
      }
    }
    if (!bFoundMappingSurface && fsurf->IsSurfaceLoaded( FSSurface::SurfaceWhite ))
    {
      fsurf->GetSurfaceRASAtVertex(nVertex, pos_out, FSSurface::SurfaceWhite);
      surf->GetTargetAtSurfaceRAS(pos_out, pos_out);
      if (bAutoOrient)
      {
        double v[3];
        surf->GetSmoothedVertexNormal(nVertex, v);
        this->AlignViewToNormal(v);
      }
      return true;
    }
  }
  else
    m_cursorInflatedSurf->Hide();
  return false;
}

void RenderView3D::MapToInflatedCoords(double *pos_in)
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("Surface");
  LayerSurface* inflated = NULL;
  foreach (Layer* s, layers)
  {
    LayerSurface* f = (LayerSurface*)s;
    if (f->IsInflated())
      inflated = f;
  }

  if (!inflated)
  {
    m_cursorInflatedSurf->Hide();
    RequestRedraw();
    return;
  }

  foreach (Layer* s, layers)
  {
    LayerSurface* f = (LayerSurface*)s;
    if (f != inflated && QFileInfo(f->GetFileName()).fileName().contains(inflated->GetMappingSurfaceName()))
    {
      if (f->GetHemisphere() == inflated->GetHemisphere() &&
          QFileInfo(f->GetFileName()).absolutePath() == QFileInfo(inflated->GetFileName()).absolutePath())
      {
        int nvo = f->GetVertexIndexAtTarget(pos_in, NULL);
        if (nvo >= 0)
        {
          double pos[3];
          inflated->GetTargetAtVertex(nvo, pos);
          inflated->SetCurrentVertex(nvo);
          if (m_cursor3D->IsShown())
            m_cursorInflatedSurf->Show();
          m_cursorInflatedSurf->SetPosition(pos);
          RequestRedraw();
        }
        return;
      }
    }
  }
}

int RenderView3D::PickCurrentSurfaceVertex(int posX, int posY, LayerSurface* curSurf)
{
  LayerCollection* lc_surface = MainWindow::GetMainWindow()->GetLayerCollection( "Surface" );

  this->setToolTip("");
  if ( lc_surface->IsEmpty() )
  {
    return -1;
  }

  posY = this->rect().height() - posY;
  int nvo = -1;
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetRenderWindow()->GetInteractor()->GetPicker() );
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

    double pos[3];
  #if VTK_MAJOR_VERSION > 7
    if (devicePixelRatio() > 1)
    {
        posX *= devicePixelRatio();
        posY *= devicePixelRatio();
    }
  #endif
    picker->Pick( posX, posY, 0, GetRenderer() );
    picker->GetPickPosition( pos );

    vtkProp* prop = picker->GetViewProp();
    Layer* layer = lc_surface->HasProp( prop );
    if (layer)
    {
      LayerSurface* surf = (LayerSurface*)layer;
      if (curSurf == NULL || curSurf == surf)
        nvo = surf->GetVertexIndexAtTarget(pos, NULL);
      if (curSurf == NULL)
      {
        surf->SetCurrentVertex(nvo);
        emit SurfaceVertexClicked(surf);
      }
    }
  }
  return nvo;
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
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    m_actorSliceFrames[i]->GetProperty()->SetLineWidth( 2*ratio );
    m_actorSliceFrames[i]->GetProperty()->SetColor( colors[i] );
  }
  if ( n >= 0 && n <= 2 )
  {
    m_actorSliceFrames[n]->GetProperty()->SetLineWidth( 4*ratio );
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
  return m_bShowSliceFrames;
}

void RenderView3D::SetShowSliceFrames( bool bShow )
{
  for ( int i = 0; i < 3; i++ )
  {
    int nFlag = (m_bSliceVisibility[i] && bShow ? 1: 0);
    m_actorSliceFrames[i]->SetVisibility( nFlag );
    m_actorSliceBoundingBox[i]->SetVisibility( nFlag );
    m_actorSliceBoundingBox[i]->SetPickable( nFlag );
  }
  m_bShowSliceFrames = bShow;

  RequestRedraw();
}

void RenderView3D::SetShowAllSlices(bool bShow)
{
  for (int i = 0; i < 3; i++)
    m_bSliceVisibility[i] = bShow;
  SetShowSliceFrames(m_bShowSliceFrames);
  RefreshAllActors();
}

void RenderView3D::ShowSlice(int nPlane, bool bshow)
{
  m_bSliceVisibility[nPlane] = bshow;
  SetShowSliceFrames(m_bShowSliceFrames);
  RefreshAllActors();
}

void RenderView3D::UpdateCursorRASPosition( int posX, int posY )
{
  m_bToUpdateCursorPosition = true;
  m_nCursorCoord[0] = posX;
  m_nCursorCoord[1] = posY;
}

void RenderView3D::UpdateMouseRASPosition( int posX, int posY, bool bSlicePickOnly )
{
  m_bToUpdateRASPosition = true;
  m_bSlicePickOnly = bSlicePickOnly;
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
  posY = this->rect().height() - posY;
#if VTK_MAJOR_VERSION > 7
  if (devicePixelRatio() > 1)
  {
      posX *= devicePixelRatio();
      posY *= devicePixelRatio();
  }
#endif
  picker->Pick( posX, posY, 0, GetRenderer() );
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
    return NULL;
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
  LayerCollection* lc = mainwnd->GetLayerCollection("MRI");
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
  if (dMaxLength > 0)
    m_cursor3D->RebuildActor(dMaxLength/256);

  return true;
}

void RenderView3D::UpdateAxesActor()
{
  LayerSurface* surf = (LayerSurface*)MainWindow::GetMainWindow()->GetActiveLayer("Surface");
  if (surf && surf->IsVisible())
  {
    vtkActor* prop = surf->GetMainActor();
    double bounds[6];
    prop->GetBounds(bounds);
    m_actorAxesActor->SetXAxisRange(0, bounds[1]-bounds[0]);
    m_actorAxesActor->SetYAxisRange(0, bounds[3]-bounds[2]);
    m_actorAxesActor->SetZAxisRange(0, bounds[5]-bounds[4]);
    m_actorAxesActor->SetBounds(bounds);
    m_actorAxesActor->VisibilityOn();
  }
  else
  {
    m_actorAxesActor->VisibilityOff();
  }
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

void RenderView3D::Azimuth(double degrees)
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->OrthogonalizeViewUp();
  cam->Azimuth(degrees);
  m_renderer->ResetCameraClippingRange();
  RequestRedraw();
}

void RenderView3D::Elevation(double degrees)
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->OrthogonalizeViewUp();
  cam->Elevation(degrees);
  cam->OrthogonalizeViewUp();
  m_renderer->ResetCameraClippingRange();
  RequestRedraw();
}

void RenderView3D::Rotate90()
{
  Azimuth(90);
}

void RenderView3D::Rotate180()
{
  Azimuth(180);
}

void RenderView3D::UpdateScalarBar()
{
  //    LayerSurface* surf = (LayerSurface*) MainWindow::GetMainWindow()->GetActiveLayer( "Surface" );
  //    if ( surf && surf->GetActiveOverlay() )
  //    {
  //        m_actorScalarBar->SetLookupTable( surf->GetActiveOverlay()->GetProperty()->GetLookupTable() );
  //    }
  //    else
  {
    RenderView::UpdateScalarBar();
  }
}

void RenderView3D::TriggerContextMenu( QMouseEvent* event )
{
  QMenu* menu = new QMenu(this);
  bool bShowBar = this->GetShowScalarBar();
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  QList<Layer*> layers = mainwnd->GetLayers("Surface");
  layers << mainwnd->GetLayers("MRI");
  layers << mainwnd->GetLayers("PointSet");

  QMenu* submenu = menu->addMenu("Reset View");
  //  QAction* act = new QAction("Right", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewRight()));
  //  submenu->addAction(act);
  //  act = new QAction("Left", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewLeft()));
  //  submenu->addAction(act);
  //  act = new QAction("Anterior", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewAnterior()));
  //  submenu->addAction(act);
  //  act = new QAction("Posterior", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewPosterior()));
  //  submenu->addAction(act);
  //  act = new QAction("Superior", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewSuperior()));
  //  submenu->addAction(act);
  //  act = new QAction("Inferior", this);
  //  connect(act, SIGNAL(triggered()), this, SLOT(ResetViewInferior()));
  //  submenu->addAction(act);
  submenu->addAction(mainwnd->ui->actionResetViewLeft);
  submenu->addAction(mainwnd->ui->actionResetViewRight);
  submenu->addAction(mainwnd->ui->actionResetViewAnterior);
  submenu->addAction(mainwnd->ui->actionResetViewPosterior);
  submenu->addAction(mainwnd->ui->actionResetViewSuperior);
  submenu->addAction(mainwnd->ui->actionResetViewInferior);
  menu->addSeparator();

  QAction* act = new QAction("Rotate Around Cursor", this);
  act->setCheckable(true);
  act->setChecked(GetFocalPointAtCursor());
  connect(act, SIGNAL(toggled(bool)), SLOT(SetFocalPointAtCursor(bool)));
  menu->addAction(act);
  menu->addSeparator();

  act = new QAction("Show All Slices", this);
  act->setData(3);
  menu->addAction(act);
  connect(act, SIGNAL(triggered()), this, SLOT(OnShowSlice()));
  act = new QAction("Hide All Slices", this);
  act->setData(-1);
  menu->addAction(act);
  connect(act, SIGNAL(triggered()), this, SLOT(OnShowSlice()));
  menu->addSeparator();

  QStringList slice_names;
  slice_names << "Show Sagittal Slice" << "Show Coronal Slice" << "Show Axial Slice";
  for (int i = 0; i < 3; i++)
  {
    act = new QAction(slice_names[i], this);
    act->setData(i);
    act->setCheckable(true);
    act->setChecked(m_bSliceVisibility[i]);
    menu->addAction(act);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(OnShowSlice(bool)));
  }
  menu->addSeparator();

  if (!mainwnd->GetLayers("MRI").isEmpty())
  {
    menu->addAction(MainWindow::GetMainWindow()->ui->actionShowSliceFrames);
  }

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

    QAction* act = new QAction("Show 3D Scale", this);
    act->setCheckable(true);
    act->setChecked(GetShowAxes());
    connect(act, SIGNAL(toggled(bool)), SLOT(SetShowAxes(bool)));
    menu->addAction(act);
  }
  LayerMRI* mri = qobject_cast<LayerMRI*>(mainwnd->GetActiveLayer("MRI"));
  if (mri && mri->GetProperty()->GetShowAsContour())
  {
    menu->addSeparator();
    QAction* act = new QAction("Save IsoSurface As...", this);
    menu->addAction(act);
    connect(act, SIGNAL(triggered()), mainwnd, SLOT(OnSaveIsoSurface()));
  }

  LayerSurface* surf = (LayerSurface*)mainwnd->GetActiveLayer("Surface");
  if ( (surf && surf->IsContralateralPossible() ) || !mainwnd->GetLayers("FCD").isEmpty())
  {
    menu->addSeparator();
    QAction* act = new QAction("Go To Contralateral Point", this);
    menu->addAction(act);
    connect(act, SIGNAL(triggered()), mainwnd, SLOT(GoToContralateralPoint()));
  }

  if (!mainwnd->IsEmpty() && mainwnd->GetMainView() == this)
  {
      menu->addSeparator();
      QAction* action = new QAction("Copy", this);
      connect(action, SIGNAL(triggered(bool)), mainwnd, SLOT(OnCopyView()));
      menu->addAction(action);
      menu->addAction(mainwnd->ui->actionSaveScreenshot);
  }
  menu->exec(event->globalPos());
}

void RenderView3D::OnShowSlice(bool bShow)
{
  QAction* act = qobject_cast<QAction*>(sender());
  if (act)
  {
    int n = act->data().toInt();
    if (n < 0)
    {
      for (int i = 0; i < 3; i++)
        m_bSliceVisibility[i] = false;
    }
    else if (n > 2)
    {
      for (int i = 0; i < 3; i++)
        m_bSliceVisibility[i] = true;
    }
    else
      m_bSliceVisibility[n] = bShow;
    SetShowSliceFrames(m_bShowSliceFrames);
    RefreshAllActors();
  }
}

void RenderView3D::HideSlices()
{
  for (int i = 0; i < 3; i++)
    m_bSliceVisibility[i] = false;
  SetShowSliceFrames(false);
  RefreshAllActors();
}

void RenderView3D::SetCamera(const QVariantMap &info)
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  if (cam)
  {
    if (!info.contains("Position") || !info.contains("FocalPoint") ||
        !info.contains("ViewUp") || !info.contains("ViewAngle"))
      return;
    double pos[3], focal_pt[3], view_up[3];
    QVariantMap map = info["Position"].toMap();
    pos[0] = map["x"].toDouble();
    pos[1] = map["y"].toDouble();
    pos[2] = map["z"].toDouble();
    cam->SetPosition(pos);
    map = info["FocalPoint"].toMap();
    focal_pt[0] = map["x"].toDouble();
    focal_pt[1] = map["y"].toDouble();
    focal_pt[2] = map["z"].toDouble();
    cam->SetFocalPoint(focal_pt);
    map = info["ViewUp"].toMap();
    view_up[0] = map["x"].toDouble();
    view_up[1] = map["y"].toDouble();
    view_up[2] = map["z"].toDouble();
    cam->SetViewUp(view_up);
    cam->SetViewAngle(info["ViewAngle"].toDouble());
    //    map = info["ClippingRange"].toMap();
    //    clip_range[0] = map["near"].toDouble();
    //    clip_range[1] = map["far"].toDouble();
    //    cam->SetClippingRange(clip_range);
    m_renderer->ResetCameraClippingRange();
    Render();
    if (info.contains("ViewSize"))
    {
      map = info["ViewSize"].toMap();
      MainWindow::GetMainWindow()->SetViewSize(map["width"].toInt(), map["height"].toInt());
    }
  }
}

QVariantMap RenderView3D::GetCamera()
{
  QVariantMap info;
  vtkCamera* cam = m_renderer->GetActiveCamera();
  if (cam)
  {
    double pos[3], focal_pt[3], view_up[3];
    cam->GetPosition(pos);
    cam->GetFocalPoint(focal_pt);
    cam->GetViewUp(view_up);
    QVariantMap map;
    map["x"] = pos[0];
    map["y"] = pos[1];
    map["z"] = pos[2];
    info["Position"] = map;
    map.clear();
    map["x"] = focal_pt[0];
    map["y"] = focal_pt[1];
    map["z"] = focal_pt[2];
    info["FocalPoint"] = map;
    map.clear();
    map["x"] = view_up[0];
    map["y"] = view_up[1];
    map["z"] = view_up[2];
    info["ViewUp"] = map;
    info["ViewAngle"] = cam->GetViewAngle();
    //    map.clear();
    //    cam->GetClippingRange(clip_range);
    //    map["near"] = clip_range[0];
    //    map["far"] = clip_range[1];
    //    info["ClippingRange"] = map;
    map.clear();
    map["width"] = width();
    map["height"] = height();
    info["ViewSize"] = map;
  }
  return info;
}

void RenderView3D::ZoomAtCursor(int x, int y, double factor)
{
  double pt[3];
  this->ScreenToWorld( x, y, 0, pt[0], pt[1], pt[2] );
  vtkCamera* cam = m_renderer->GetActiveCamera();
  if (cam)
  {
    double pos[3], vproj[3];
    cam->GetPosition(pos);
    MyUtils::GetVector(pos, pt, vproj);
    double dist = cam->GetDistance();
    for (int i = 0; i < 3; i++)
    {
      pt[i] = pos[i] + vproj[i]*dist;
    }
    PanToWorld(pt);
    if (factor >= 1 || cam->GetViewAngle() < 90)
      cam->Zoom(factor);
    RequestRedraw();
  }
}

void RenderView3D::ShowCursor(bool bshow)
{
  m_cursor3D->Show(bshow);
  m_cursorInflatedSurf->Show(bshow);
  m_bShowCursor = bshow;
}

void RenderView3D::OnLayerVisibilityChanged()
{
  LayerCollection* lc_surface = MainWindow::GetMainWindow()->GetLayerCollection( "Surface" );
  bool bShowCursor = m_cursor3D->IsShown();
  bool bShowInflated = m_cursorInflatedSurf->IsShown();
  if (lc_surface->IsEmpty())
    m_cursorInflatedSurf->Hide();
  else
  {
    QList<Layer*> list = lc_surface->GetLayers();
    m_cursor3D->Hide();
    m_cursorInflatedSurf->Hide();
    foreach (Layer* layer, list)
    {
      LayerSurface* surf = (LayerSurface*)layer;
      if (!surf->IsInflated() && surf->IsVisible() && surf->GetVisibleIn3D() && m_bShowCursor)
        m_cursor3D->Show();
      else if (surf->IsInflated() && surf->IsVisible() && surf->GetVisibleIn3D() && m_bShowCursor)
        m_cursorInflatedSurf->Show();
    }
  }
  if (bShowCursor != m_cursor3D->IsShown() || bShowInflated != m_cursorInflatedSurf->IsShown())
    RequestRedraw();
}

void RenderView3D::SetFocalPointAtCursor(bool b)
{
  m_bFocalPointAtCursor = b;
  m_interactorStyle->SetRotateByPoint(b, b?m_cursor3D->GetPosition():NULL);
}

void RenderView3D::SetShowAxes(bool b)
{
  m_bShowAxes = b;
  RefreshAllActors();
}
