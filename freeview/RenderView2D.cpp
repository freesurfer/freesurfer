/**
 * @file  RenderView2D.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/05/06 19:43:02 $
 *    $Revision: 1.46 $
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
#include "RenderView2D.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "Contour2D.h"
#include "VolumeCropper.h"
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include "Annotation2D.h"
#include "Interactor2DNavigate.h"
#include "Interactor2DMeasure.h"
#include "Interactor2DVoxelEdit.h"
#include "Interactor2DROIEdit.h"
#include "Interactor2DPointSetEdit.h"
#include "Interactor2DVolumeCrop.h"
#include <vtkActor2D.h>
#include <vtkScalarBarActor.h>
#include "Region2DRectangle.h"
#include "Cursor2D.h"
#include "MyUtils.h"
#include <QDebug>

RenderView2D::RenderView2D( QWidget* parent ) : RenderView( parent )
{
  m_renderer->GetActiveCamera()->ParallelProjectionOn();
  m_contour2D = new Contour2D( this );
  m_cursor2D = new Cursor2D( this );
  m_annotation2D = new Annotation2D( this );
  m_selection2D = new Region2DRectangle( this );
  m_selection2D->SetEnableStats( false );
  connect(m_cursor2D, SIGNAL(Updated()), this, SLOT(RequestRedraw()));

  m_interactorNavigate = new Interactor2DNavigate( this );
  m_interactorMeasure = new Interactor2DMeasure( this );
  m_interactorVoxelEdit = new Interactor2DVoxelEdit( this );
  m_interactorROIEdit = new Interactor2DROIEdit( this );
  m_interactorPointSetEdit = new Interactor2DPointSetEdit( this );
  m_interactorVolumeCrop = new Interactor2DVolumeCrop( this );
  SetInteractionMode( IM_Navigate );
}

void RenderView2D::SetViewPlane( int nPlane )
{
  m_nViewPlane = nPlane;
  m_contour2D->SetPlane( nPlane );
}

void RenderView2D::SetInteractionMode( int nMode )
{
  RenderView::SetInteractionMode( nMode );

  switch ( nMode )
  {
  case IM_Measure:
    m_interactor = m_interactorMeasure;
    break;
  case IM_VoxelEdit:
    m_interactor = m_interactorVoxelEdit;
    break;
  case IM_ROIEdit:
    m_interactor = m_interactorROIEdit;
    break;
  case IM_PointSetEdit:
    m_interactor = m_interactorPointSetEdit;
    break;
  case IM_VolumeCrop:
    m_interactor = m_interactorVolumeCrop;
    break;
  default:
    m_interactor = m_interactorNavigate;
    break;
  }
}

void RenderView2D::RefreshAllActors(bool bForScreenShot)
{
  MainWindow* mainwnd = MainWindow::GetMainWindow();
  SettingsScreenshot setting = mainwnd->GetScreenShotSettings();
  LayerCollection* lc = mainwnd->GetLayerCollection( "MRI" );
  m_renderer->RemoveAllViewProps();
  lc->Append2DProps( m_renderer, m_nViewPlane );

  LayerMRI* mri = NULL;
  if ( !lc->IsEmpty() )
  {
    mri = (LayerMRI*)lc->GetLayer(0);
    mri->Remove2DProps( m_renderer, m_nViewPlane );

    m_renderer->AddViewProp( m_contour2D->GetActor() );

    mri->Append2DProps( m_renderer, m_nViewPlane );

    // add annotation and cursor
    if (!bForScreenShot || !setting.HideCursor)
    {
      m_cursor2D->AppendActor( m_renderer );
    }
    if (!bForScreenShot || !setting.HideCoords)
    {
      m_annotation2D->AppendAnnotations( m_renderer );
    }
    m_selection2D->AppendProp( m_renderer );

    // add scalar bar
    m_renderer->AddViewProp( m_actorScalarBar );
  }

  mainwnd->GetLayerCollection( "ROI" )->Append2DProps( m_renderer, m_nViewPlane );
  mainwnd->GetLayerCollection( "Surface" )->Append2DProps( m_renderer, m_nViewPlane );
  mainwnd->GetLayerCollection( "PointSet" )->Append2DProps( m_renderer, m_nViewPlane );

  mainwnd->GetLayerCollection( "Supplement" )->Append2DProps( m_renderer, m_nViewPlane );

  mainwnd->GetVolumeCropper()->Append2DProps( m_renderer, m_nViewPlane );

  // add regions
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    m_regions[i]->AppendProp( m_renderer );
  }

  // add focus frame
  if (!bForScreenShot)
  {
    m_renderer->AddViewProp( m_actorFocusFrame );
  }

  if (!lc->IsEmpty())
  {
    double* orig = lc->GetWorldOrigin();
    double* size = lc->GetWorldSize();
    m_renderer->ResetCameraClippingRange(orig[0], orig[0]+size[0],
                                         orig[1], orig[1]+size[1],
                                         orig[2], orig[2]+size[2]);
  }
  RenderView::RefreshAllActors(bForScreenShot);
}

void RenderView2D::UpdateViewByWorldCoordinate()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double wcenter[3];
  for ( int i = 0; i < 3; i++ )
  {
    wcenter[i] = m_dWorldOrigin[i] + m_dWorldSize[i] / 2;
  }
  cam->SetFocalPoint( wcenter );
  switch ( m_nViewPlane )
  {
  case 0:
    cam->SetPosition( wcenter[0] + m_dWorldSize[0], wcenter[1], wcenter[2] );
    cam->SetViewUp( 0, 0, 1 );
    break;
  case 1:
    cam->SetPosition( wcenter[0], wcenter[1] + m_dWorldSize[1], wcenter[2] );
    cam->SetViewUp( 0, 0, 1 );
    break;
  case 2:
    cam->SetPosition( wcenter[0], wcenter[1], wcenter[2] - m_dWorldSize[2] );
    break;
  }
//  m_renderer->ResetCameraClippingRange();
  cam->SetParallelScale( qMax( qMax(m_dWorldSize[0], m_dWorldSize[1]), m_dWorldSize[2]) );
}

void RenderView2D::UpdateAnnotation()
{
  m_annotation2D->Update( m_renderer, m_nViewPlane );
}

void RenderView2D::UpdateMouseRASPosition( int posX, int posY )
{
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  if ( !lc )
  {
    return;
  }

  double pos[3];
  MousePositionToRAS( posX, posY, pos );

  lc->SetCurrentRASPosition( pos );
}

void RenderView2D::UpdateCursorRASPosition( int posX, int posY )
{
  LayerCollection* lc = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );
  if ( !lc )
  {
    return;
  }

  double pos[3];
  MousePositionToRAS( posX, posY, pos );

  lc->SetCursorRASPosition( pos );
  MainWindow::GetMainWindow()->SetSlicePosition( pos );
}

void RenderView2D::Update2DOverlay()
{
  m_cursor2D->Update();
  m_selection2D->Update();
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    m_regions[i]->Update();
  }

  double slicePos[3];
  MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
  m_contour2D->UpdateSliceLocation( slicePos[m_nViewPlane] );
}

void RenderView2D::OnSlicePositionChanged()
{
  double slicePos[3];
  MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    m_regions[i]->UpdateSlicePosition( m_nViewPlane, slicePos[m_nViewPlane] );
  }
  m_cursor2D->SetPosition( slicePos );
  Update2DOverlay();
  UpdateAnnotation();

  RenderView::OnSlicePositionChanged();
}

void RenderView2D::MousePositionToRAS( int posX, int posY, double* pos )
{
  pos[0] = posX;
  pos[1] = rect().height() - posY;
  pos[2] = 0;
  m_renderer->ViewportToNormalizedViewport( pos[0], pos[1] );
  m_renderer->NormalizedViewportToView( pos[0], pos[1], pos[2] );
  m_renderer->ViewToWorld( pos[0], pos[1], pos[2] );

  double slicePos[3];
  MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
  pos[m_nViewPlane] = slicePos[m_nViewPlane];
}

LayerMRI* RenderView2D::GetFirstNonLabelVolume( )
{
  QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" )->GetLayers();
  for ( int i = 0; i < layers.size(); i++ )
  {
    LayerMRI* layer = ( LayerMRI*)layers[i];
    if ( layer->IsVisible() && layer->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT )
    {
      return layer;
    }
  }
  return NULL;
}

void RenderView2D::StartSelection( int nX, int nY )
{
  m_selection2D->SetRect( nX, nY, 1, 1 );
  m_selection2D->Show( true );
}

void RenderView2D::UpdateSelection( int nX, int nY )
{
  m_selection2D->SetBottomRight( nX, nY );
}

void RenderView2D::StopSelection()
{
  m_selection2D->Show( false );

  LayerMRI* layer = GetFirstNonLabelVolume();
  if ( layer )
  {
    double range[2], m_dPt0[3], m_dPt2[3];
    m_selection2D->GetWorldPoint( 0, m_dPt0 );
    m_selection2D->GetWorldPoint( 2, m_dPt2 );
    if ( layer->GetVoxelValueRange( m_dPt0, m_dPt2, m_nViewPlane, range ) )
    {
      switch ( layer->GetProperty()->GetColorMap() )
      {
      case LayerPropertyMRI::Grayscale:
        layer->GetProperty()->SetMinMaxGrayscaleWindow( range[0], range[1] );
        break;
      case LayerPropertyMRI::Heat:
        layer->GetProperty()->SetHeatScale( range[0], (range[1]-range[0])/2, range[1] );
        break;
      default:
        layer->GetProperty()->SetMinMaxGenericThreshold( range[0], range[1] );
        break;
      }
    }
  }
}

Region2D* RenderView2D::GetRegion( int nX, int nY, int* index_out )
{
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i]->Contains( nX, nY, index_out ) )
    {
      return m_regions[i];
    }
  }
  return NULL;
}

void RenderView2D::AddRegion( Region2D* region )
{
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i] == region )
    {
      return;
    }
  }

  m_regions.push_back( region );
  emit RegionSelected( region );
  RefreshAllActors();
}

void RenderView2D::DeleteRegion( Region2D* region )
{
  for ( int i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i] == region )
    {
      m_regions.erase( m_regions.begin() + i );
      emit RegionRemoved( region );
      delete region;
      RefreshAllActors();
      return;
    }
  }
}

void RenderView2D::MoveSlice( int nStep )
{
  MainWindow* mainWnd = MainWindow::GetMainWindow();
  LayerCollection* lc_mri = mainWnd->GetLayerCollection( "MRI" );
  double* voxelSize = lc_mri->GetWorldVoxelSize();
  int nPlane = GetViewPlane();
  mainWnd->OffsetSlicePosition( nPlane, voxelSize[nPlane]*nStep );
  lc_mri->SetCursorRASPosition( lc_mri->GetSlicePosition() );
  UpdateAnnotation();
}

void RenderView2D::SyncZoomTo( RenderView2D* view )
{
  m_renderer->GetActiveCamera()->SetParallelScale( view->m_renderer->GetActiveCamera()->GetParallelScale() );
// PanToWorld( GetCursor2D()->GetPosition() );
  EnsureCursor2DVisible();
  Update2DOverlay();
  UpdateAnnotation();
  RequestRedraw();
}

bool RenderView2D::EnsureCursor2DVisible()
{
  double* pos = GetCursor2D()->GetPosition();
  double x = pos[0], y = pos[1], z = pos[2];
  m_renderer->WorldToView( x, y, z );
  m_renderer->ViewToNormalizedViewport( x, y, z );
  if ( x < 0 || x > 1 || y < 0 || y > 1 )
  {
    PanToWorld( pos );
    return true;
  }
  else
  {
    return false;
  }
}

void RenderView2D::PanToWorld( double* pos )
{
  double focalPt[3], camPos[3], vproj[3];
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->GetFocalPoint( focalPt );
  cam->GetPosition( camPos );
  cam->GetDirectionOfProjection( vproj );
  double camDist = cam->GetDistance();

  double dist = MyUtils::GetDistance( pos, focalPt );
  double v[3];
  MyUtils::GetVector( pos, focalPt, v );
  dist *= MyUtils::Dot( vproj, v );

  for ( int i = 0; i < 3; i++ )
  {
    focalPt[i] = pos[i] + vproj[i] * dist;
    camPos[i] = focalPt[i] - vproj[i] * camDist;
  }

  cam->SetFocalPoint( focalPt );
  cam->SetPosition( camPos );
}

void RenderView2D::ZoomAtCursor( int nX, int nY, double factor )
{
  double pos[3];
  pos[0] = nX;
  pos[1] = rect().height() - nY;
  pos[2] = 0;
  m_renderer->ViewportToNormalizedViewport( pos[0], pos[1] );
  m_renderer->NormalizedViewportToView( pos[0], pos[1], pos[2] );
  m_renderer->ViewToWorld( pos[0], pos[1], pos[2] );

  // first move the click point to the center of viewport
  PanToWorld( pos );

  // then zoom
  m_renderer->GetActiveCamera()->Zoom( factor );
  Update2DOverlay();
  UpdateAnnotation();
  RequestRedraw();
}

void RenderView2D::resizeEvent(QResizeEvent *event)
{
  RenderView::resizeEvent(event);

  Update2DOverlay();
  UpdateAnnotation();
}

void RenderView2D::ShowCoordinateAnnotation( bool bShow )
{
  m_annotation2D->Show( bShow );
}

bool RenderView2D::GetShowCoordinateAnnotation()
{
  return m_annotation2D->IsVisible();
}

bool RenderView2D::SetSliceNumber( int nNum )
{
  LayerCollection* lc_mri = MainWindow::GetMainWindow()->GetLayerCollection( "MRI" );

  LayerMRI* mri = (LayerMRI*)lc_mri->GetActiveLayer();
  if ( !mri )
  {
    return false;
  }

  vtkImageData* imagedata = mri->GetImageData();
  int nPlane = GetViewPlane();
  int* dim = imagedata->GetDimensions();
  double* voxelsize = imagedata->GetSpacing();
  double* orig = imagedata->GetOrigin();
  if ( nNum < 0 || nNum >= dim[nPlane] )
  {
    return false;
  }

  double pos[3];
  lc_mri->GetSlicePosition( pos );
  pos[nPlane] = orig[nPlane] + nNum * voxelsize[nPlane];
  MainWindow::GetMainWindow()->SetSlicePosition( pos );
  lc_mri->SetCursorRASPosition( pos );
  return true;
}
