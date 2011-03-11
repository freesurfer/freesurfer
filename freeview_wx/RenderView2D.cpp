/**
 * @file  RenderView2D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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

#include <wx/xrc/xmlres.h>
#include "RenderView2D.h"
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkActor2D.h>
#include <vtkScalarBarActor.h>
#include <vtkPropCollection.h>
#include <vtkCoordinate.h>
#include <vtkImageActor.h>
#include "LayerCollection.h"
#include "LayerCollectionManager.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include "Annotation2D.h"
#include "Cursor2D.h"
#include "Interactor2DNavigate.h"
#include "Interactor2DROIEdit.h"
#include "Interactor2DVoxelEdit.h"
#include "Interactor2DWayPointsEdit.h"
#include "Interactor2DMeasure.h"
#include "Interactor2DCropVolume.h"
#include "MyUtils.h"
#include "Region2DRectangle.h"
#include "ToolWindowMeasure.h"
#include "Contour2D.h"
#include "VolumeCropper.h"

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))


IMPLEMENT_DYNAMIC_CLASS(RenderView2D, RenderView)

BEGIN_EVENT_TABLE(RenderView2D, RenderView)
EVT_SIZE        ( RenderView2D::OnSize )
END_EVENT_TABLE()

RenderView2D::RenderView2D( int nPlane ) : RenderView(), m_nViewPlane( nPlane )
{
  Initialize2D();
}

RenderView2D::RenderView2D( int nPlane, wxWindow* parent, int id ) : RenderView( parent, id ), m_nViewPlane( nPlane )
{
  Initialize2D();
}

void RenderView2D::Initialize2D()
{
  m_renderer->GetActiveCamera()->ParallelProjectionOn();
  m_annotation2D = new Annotation2D;
  m_cursor2D = new Cursor2D( this );
  m_selection2D = new Region2DRectangle( this );
  m_selection2D->SetEnableStats( false );
  m_contour2D = new Contour2D( this );

  m_interactorNavigate  = new Interactor2DNavigate();
  m_interactorMeasure  = new Interactor2DMeasure();
  m_interactorVoxelEdit  = new Interactor2DVoxelEdit();
  m_interactorROIEdit  = new Interactor2DROIEdit();
  m_interactorWayPointsEdit = new Interactor2DWayPointsEdit();
  m_interactorCropVolume = new Interactor2DCropVolume();

  m_interactorNavigate->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorMeasure->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorVoxelEdit->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorROIEdit->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorWayPointsEdit->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorCropVolume->AddListener( MainWindow::GetMainWindowPointer() );

  SetInteractionMode( IM_Navigate );
}

RenderView2D* RenderView2D::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new RenderView2D();
}

RenderView2D::~RenderView2D()
{
  m_interactor = NULL;
  delete m_annotation2D;
  delete m_cursor2D;
  delete m_contour2D;
  
  delete m_interactorNavigate;
  delete m_interactorMeasure;
  delete m_interactorVoxelEdit;
  delete m_interactorROIEdit;
  delete m_interactorWayPointsEdit;
  delete m_interactorCropVolume;
  
  for ( size_t i = 0; i < m_regions.size(); i++ )
  {
    delete m_regions[i];
  }
  m_regions.clear();
}

void RenderView2D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void RenderView2D::SetInteractionMode( int nMode )
{
  RenderView::SetInteractionMode( nMode );

  switch ( nMode )
  {
  case IM_Measure:
    m_interactor = m_interactorMeasure;
    break;
  case IM_ROIEdit:
    m_interactor = m_interactorROIEdit;
    break;
  case IM_VoxelEdit:
    m_interactor = m_interactorVoxelEdit;
    break;
  case IM_WayPointsEdit:
    m_interactor = m_interactorWayPointsEdit;
    break;
  case IM_VolumeCrop:    
    m_interactor = m_interactorCropVolume;
    break;     
  default:
    m_interactor = m_interactorNavigate;
    break;
  }

// m_interactor->AddListener( MainWindow::GetMainWindowPointer() );
}

void RenderView2D::RefreshAllActors()
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  
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
    m_cursor2D->AppendActor( m_renderer );
    m_annotation2D->AppendAnnotations( m_renderer );
    m_selection2D->AppendProp( m_renderer );

    // add scalar bar
    m_renderer->AddViewProp( m_actorScalarBar );
  }
  
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" )->Append2DProps( m_renderer, m_nViewPlane );
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->Append2DProps( m_renderer, m_nViewPlane );
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "WayPoints" )->Append2DProps( m_renderer, m_nViewPlane );
  
  MainWindow::GetMainWindowPointer()->GetVolumeCropper()->Append2DProps( m_renderer, m_nViewPlane );
  
  // add regions
  for ( size_t i = 0; i < m_regions.size(); i++ )
    m_regions[i]->AppendProp( m_renderer );
  
  // add focus frame
  m_renderer->AddViewProp( m_actorFocusFrame );

  NeedRedraw();
  //Render();
}

int RenderView2D::GetViewPlane()
{
  return m_nViewPlane;
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
  cam->SetParallelScale( max( max(m_dWorldSize[0], m_dWorldSize[1]), m_dWorldSize[2]) / 2 );
//  m_renderer->ResetCameraClippingRange();
}

void RenderView2D::UpdateAnnotation()
{
  m_annotation2D->Update( m_renderer, m_nViewPlane );
}

void RenderView2D::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerActorUpdated" ||
       iMsg == "LayerAdded" ||
       iMsg == "LayerMoved" ||
       iMsg == "LayerRemoved" )
  {
    if ( iMsg == "LayerRemoved" )
    {
      if ( MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->IsEmpty() && m_regions.size() > 0 )
      {
        for ( size_t i = 0; i < m_regions.size(); i++ )
        {
          delete m_regions[i];
        }
        m_regions.clear();
        RefreshAllActors();
      }
    }
    UpdateAnnotation();
  }
  else if ( iMsg == "CursorRASPositionChanged" )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
    m_cursor2D->SetPosition( lc->GetCursorRASPosition() );
    // if ( EnsureCursor2DVisible() )
    Update2DOverlay();
  }
  else if ( iMsg == "Zooming" )
  {
    Settings2D s = MainWindow::GetMainWindowPointer()->Get2DSettings();
    if ( s.SyncZoomFactor )
    {
      RenderView2D* view = (RenderView2D*)iData;
      if ( view )
      {
        SyncZoomTo( view );
      }
    }
  }
  else if ( iMsg == "SlicePositionChanged" )
  {
    double slicePos[3];
    MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
    for ( size_t i = 0; i < m_regions.size(); i++ )
    {
      m_regions[i]->UpdateSlicePosition( m_nViewPlane, slicePos[m_nViewPlane] );
    }
  }

  RenderView::DoListenToMessage( iMsg, iData, sender );
}

void RenderView2D::SyncZoomTo( RenderView2D* view )
{
  m_renderer->GetActiveCamera()->SetParallelScale( view->m_renderer->GetActiveCamera()->GetParallelScale() );
// PanToWorld( GetCursor2D()->GetPosition() );
  EnsureCursor2DVisible();
  Update2DOverlay();
  UpdateAnnotation();
  NeedRedraw();
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

void RenderView2D::TriggerContextMenu( const wxPoint& pos )
{
  /*
  if ( m_interactor != m_interactorVoxelEdit )
    return;
  
  double ras[3];
  MousePositionToRAS( pos.x, pos.y, ras );
  
  wxMenu* menu = new wxMenu;
  menu->Append( XRCID("ID_EDIT_COPY"),  _T("Copy") );
  menu->Append( XRCID("ID_EDIT_PASTE"), _T("Paste") );
  menu->Append( XRCID("ID_FILE_EXIT"), _T("Exit") );
  this->PopupMenu( menu );
  delete menu;
  */
}


void RenderView2D::UpdateMouseRASPosition( int posX, int posY )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  if ( !lc )
    return;

  double pos[3];
  MousePositionToRAS( posX, posY, pos );

  lc->SetCurrentRASPosition( pos );
}

void RenderView2D::UpdateCursorRASPosition( int posX, int posY, bool bConnectPrevious )
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  if ( !lc )
    return;

  double pos[3];
  MousePositionToRAS( posX, posY, pos );

  lc->SetCursorRASPosition( pos );
  MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->SetSlicePosition( pos );

// m_cursor2D->SetPosition( pos, bConnectPrevious );
}

void RenderView2D::MousePositionToRAS( int posX, int posY, double* pos )
{
  wxSize sz = GetClientSize();
  pos[0] = posX;
  pos[1] = sz.GetHeight() - posY;
  pos[2] = 0;
// m_renderer->DisplayToNormalizedDisplay( pos[0], pos[1] );
// m_renderer->NormalizedDisplayToViewport( pos[0], pos[1] );
  m_renderer->ViewportToNormalizedViewport( pos[0], pos[1] );
  m_renderer->NormalizedViewportToView( pos[0], pos[1], pos[2] );
  m_renderer->ViewToWorld( pos[0], pos[1], pos[2] );

  double slicePos[3];
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
  pos[m_nViewPlane] = slicePos[m_nViewPlane];
}

void RenderView2D::ZoomAtCursor( int nX, int nY, bool ZoomIn, double factor )
{
  wxSize sz = GetClientSize();
  double pos[3];
  pos[0] = nX;
  pos[1] = sz.GetHeight() - nY;
  pos[2] = 0;
  m_renderer->ViewportToNormalizedViewport( pos[0], pos[1] );
  m_renderer->NormalizedViewportToView( pos[0], pos[1], pos[2] );
  m_renderer->ViewToWorld( pos[0], pos[1], pos[2] );

  // first move the click point to the center of viewport
  PanToWorld( pos );

  // then zoom
  m_renderer->GetActiveCamera()->Zoom( ZoomIn ? factor : 1.0/factor );
  Update2DOverlay();
  UpdateAnnotation();
  NeedRedraw();
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

void RenderView2D::Update2DOverlay()
{
  m_cursor2D->Update();
  m_selection2D->Update();
  for ( size_t i = 0; i < m_regions.size(); i++ )
    m_regions[i]->Update();
  
  double slicePos[3];
  MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetSlicePosition( slicePos );
  m_contour2D->UpdateSliceLocation( slicePos[m_nViewPlane] );
}

void RenderView2D::MoveLeft()
{
  RenderView::MoveLeft();

  UpdateAnnotation();
  Update2DOverlay();
}

void RenderView2D::MoveRight()
{
  RenderView::MoveRight();

  UpdateAnnotation();
  Update2DOverlay();
}

void RenderView2D::MoveUp()
{
  RenderView::MoveUp();

  UpdateAnnotation();
  Update2DOverlay();
}

void RenderView2D::MoveDown()
{
  RenderView::MoveDown();

  UpdateAnnotation();
  Update2DOverlay();
}

void RenderView2D::OnSize( wxSizeEvent& event )
{
  RenderView::OnSize( event );
  
  Update2DOverlay();
  UpdateAnnotation();
}

void RenderView2D::PreScreenshot()
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();

  m_renderer->RemoveAllViewProps();
  lcm->Append2DProps( m_renderer, m_nViewPlane );

  // add coordinate annotation
  SettingsScreenshot s = MainWindow::GetMainWindowPointer()->GetScreenshotSettings();
  if ( !s.HideCoords )
    m_annotation2D->AppendAnnotations( m_renderer );
  if ( !s.HideCursor )
    m_cursor2D->AppendActor( m_renderer );

  // add scalar bar
  m_renderer->AddViewProp( m_actorScalarBar );
}


void RenderView2D::PostScreenshot()
{
  RefreshAllActors();
}

void RenderView2D::ShowCoordinateAnnotation( bool bShow )
{
  m_annotation2D->Show( bShow );
}

bool RenderView2D::GetShowCoordinateAnnotation()
{
  return m_annotation2D->IsVisible();
}

void RenderView2D::MoveSlice( int nStep )
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  LayerCollection* lc_mri = lcm->GetLayerCollection( "MRI" );

  double* voxelSize = lc_mri->GetWorldVoxelSize();
  int nPlane = GetViewPlane();
  lcm->OffsetSlicePosition( nPlane, voxelSize[nPlane]*nStep );
  lc_mri->SetCursorRASPosition( lc_mri->GetSlicePosition() );
}

bool RenderView2D::SetSliceNumber( int nNum )
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  LayerCollection* lc_mri = lcm->GetLayerCollection( "MRI" );
  
  LayerMRI* mri = (LayerMRI*)lc_mri->GetActiveLayer();
  if ( !mri )
    return false;
  
  vtkImageData* imagedata = mri->GetImageData();
  int nPlane = GetViewPlane();
  int* dim = imagedata->GetDimensions();
  double* voxelsize = imagedata->GetSpacing();
  double* orig = imagedata->GetOrigin();
  if ( nNum < 0 || nNum >= dim[nPlane] )
    return false;
  
  double pos[3];
  lc_mri->GetSlicePosition( pos );
  pos[nPlane] = orig[nPlane] + nNum * voxelsize[nPlane];
  lcm->SetSlicePosition( pos );
  lc_mri->SetCursorRASPosition( pos );
  return true;
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
      switch ( layer->GetProperties()->GetColorMap() )
      {
        case LayerPropertiesMRI::Grayscale:   
          layer->GetProperties()->SetMinMaxGrayscaleWindow( range[0], range[1] );
          break;
        case LayerPropertiesMRI::Heat:
          layer->GetProperties()->SetHeatScale( range[0], (range[1]-range[0])/2, range[1] );
          break;
        default:
          layer->GetProperties()->SetMinMaxGenericThreshold( range[0], range[1] );
          break;
      }
    }
  }
}

LayerMRI* RenderView2D::GetFirstNonLabelVolume( )
{
  std::vector<Layer*> layers = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetLayers();
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    LayerMRI* layer = ( LayerMRI*)layers[i];
    if ( layer->IsVisible() && layer->GetProperties()->GetColorMap() != LayerPropertiesMRI::LUT )
      return layer;
  }
  return NULL;
}

Region2D* RenderView2D::GetRegion( int nX, int nY, int* index_out )
{
  for ( size_t i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i]->Contains( nX, nY, index_out ) )
      return m_regions[i];
  }
  return NULL;
}

void RenderView2D::AddRegion( Region2D* region )
{
  for ( size_t i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i] == region )
      return;
  }
  
  m_regions.push_back( region );
  this->SendBroadcast( "RegionSelected", region );
  RefreshAllActors();
}

void RenderView2D::DeleteRegion( Region2D* region )
{
  for ( size_t i = 0; i < m_regions.size(); i++ )
  {
    if ( m_regions[i] == region )
    {
      m_regions.erase( m_regions.begin() + i );
      this->SendBroadcast( "RegionRemoved", region );
      delete region;
      RefreshAllActors();
      return;
    }
  }
}
