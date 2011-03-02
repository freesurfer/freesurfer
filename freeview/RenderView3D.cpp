/**
 * @file  RenderView3D.cpp
 * @brief View for rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 22:00:37 $
 *    $Revision: 1.50 $
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
#include "ConnectivityData.h"
#include "LayerCollection.h"
#include "LayerMRI.h"
#include "LayerPropertiesMRI.h"
#include <vtkRenderer.h>
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkActor2D.h"
#include "vtkCellPicker.h"
#include "vtkPointPicker.h"
#include "vtkPropPicker.h"
#include "vtkProp3DCollection.h"
#include "vtkScalarBarActor.h"
#include "vtkPlane.h"
#include "vtkMath.h"
#include "vtkTubeFilter.h"
#include "vtkCubeSource.h"
#include "vtkLine.h"
#include "Interactor3DNavigate.h"
#include "Interactor3DMeasure.h"
#include "Interactor3DCropVolume.h"
#include "LayerSurface.h"
#include "SurfaceOverlayProperties.h"
#include "SurfaceOverlay.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkBoundingBox.h"
#include "SurfaceRegion.h"
#include "Cursor3D.h"
#include "VolumeCropper.h"

#define SLICE_PICKER_PIXEL_TOLERANCE  15

IMPLEMENT_DYNAMIC_CLASS(RenderView3D, RenderView)

BEGIN_EVENT_TABLE(RenderView3D, RenderView)

END_EVENT_TABLE()

RenderView3D::RenderView3D() : RenderView()
{
  InitializeRenderView3D();
}

RenderView3D::RenderView3D( wxWindow* parent, int id ) : RenderView( parent, id )
{
  InitializeRenderView3D();
}

void RenderView3D::InitializeRenderView3D()
{
  this->SetDesiredUpdateRate( 5000 );

  if ( m_interactor )
    delete m_interactor;

  m_interactor = NULL;
  m_interactorNavigate = new Interactor3DNavigate();
  m_interactorMeasure = new Interactor3DMeasure();  
  m_interactorCropVolume = new Interactor3DCropVolume();
  m_interactorNavigate->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorMeasure->AddListener( MainWindow::GetMainWindowPointer() );
  m_interactorCropVolume->AddListener( MainWindow::GetMainWindowPointer() );

  m_bToUpdateRASPosition = false;
  m_bToUpdateCursorPosition = false;
  m_bToUpdateConnectivity = false;
  m_dBoundingTolerance = 4;

  vtkCellPicker* picker = vtkCellPicker::New();
// vtkPointPicker* picker = vtkPointPicker::New();
// vtkPropPicker* picker = vtkPropPicker::New();
  picker->SetTolerance( 0.0005 );
  picker->PickFromListOn();
  this->SetPicker( picker );
  picker->Delete();

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
  
  m_actorScalarBar->SetNumberOfLabels( 4 );
  
  SetInteractionMode( IM_Navigate );
}

RenderView3D* RenderView3D::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new RenderView3D;
}

RenderView3D::~RenderView3D()
{
  m_interactor = NULL;
  
  delete m_interactorNavigate;
  delete m_interactorMeasure;
  delete m_interactorCropVolume;
  
  delete m_cursor3D;
}

void RenderView3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
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
      m_interactor = m_interactorCropVolume;
      break;
    default:
      m_interactor = m_interactorNavigate;
      break;
  }
}

void RenderView3D::RefreshAllActors()
{
  if ( m_bDisabled)
    return;
  
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();

  m_renderer->RemoveAllViewProps();
  bool b[3] = { true, true, true };
  lcm->Append3DProps( m_renderer, b );

  m_cursor3D->AppendActor( m_renderer );

  // add focus frame
  m_renderer->AddViewProp( m_actorFocusFrame );

  if ( lcm->HasLayer( "MRI" ) || lcm->HasLayer( "Surface" ) )
  {
    m_renderer->AddViewProp( m_actorScalarBar ); 
    if ( lcm->HasLayer( "MRI" ) )
    {
      for ( int i = 0; i < 3; i++ )
        m_renderer->AddViewProp( m_actorSliceFrames[i] );  
    }
  }
  
  MainWindow::GetMainWindowPointer()->GetConnectivityData()->AppendProps( m_renderer );
  MainWindow::GetMainWindowPointer()->GetVolumeCropper()->Append3DProps( m_renderer );
  
  m_renderer->ResetCameraClippingRange();

  NeedRedraw();
  // Render();
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

// snap the camera to the nearest axis
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
  m_renderer->ResetCameraClippingRange();
  
  NeedRedraw();
}


void RenderView3D::UpdateMouseRASPosition( int posX, int posY )
{
  m_bToUpdateRASPosition = true;
  m_nPickCoord[0] = posX;
  m_nPickCoord[1] = posY;
}

void RenderView3D::CancelUpdateMouseRASPosition()
{
  m_bToUpdateRASPosition = false;
}

void RenderView3D::DoUpdateRASPosition( int posX, int posY, bool bCursor )
{
  LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerCollection* lc_roi = MainWindow::GetMainWindowPointer()->GetLayerCollection( "ROI" );
  LayerCollection* lc_surface = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" );

  if ( lc_mri->IsEmpty() && lc_roi->IsEmpty() && lc_surface->IsEmpty() )
    return;
  
// MousePositionToRAS( posX, posY, pos );
// vtkPointPicker* picker = vtkPointPicker::SafeDownCast( this->GetPicker() );
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetPicker() );
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
          picker->AddPickList( prop );
        prop = props->GetNextProp();
      }
    }
    // add bounding box for slice frame picking
    for ( int i = 0; i < 3; i++ )
    {
      picker->AddPickList( m_actorSliceBoundingBox[i] );
    }
    
    double pos[3];
    picker->Pick( posX, GetClientSize().GetHeight() - posY, 0, GetRenderer() );
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
              dMinDist = dist;
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
        
        picker->Pick( posX, GetClientSize().GetHeight() - posY, 0, GetRenderer() );
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
            reg = mri->SelectSurfaceRegion( pos );
          if ( reg )
          {
            NeedRedraw( true ); // force redraw
            this->SendBroadcast( "SurfaceRegionSelected", reg );
          }
        }
      }
      else if ( lc_surface->HasProp( prop ) )
      {  
        if ( bCursor )
        {
          lc_mri->SetCursorRASPosition( pos );
          MainWindow::GetMainWindowPointer()->GetLayerCollectionManager()->SetSlicePosition( pos );
          this->SendBroadcast( "SurfaceVertexClicked", this );
        }
        else
          lc_mri->SetCurrentRASPosition( pos );
      }
        
      HighlightSliceFrame( -1 );
    }
  }
}

void RenderView3D::UpdateSurfaceCorrelationData()
{
  std::vector<Layer*> layers = MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetLayers();
  for ( size_t i = 0; i < layers.size(); i++ )
  {
    ((LayerSurface*)layers[i])->UpdateCorrelationOverlay();
  }
}

vtkProp* RenderView3D::PickProp( int posX, int posY, double* pos_out )
{
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetPicker() );
  if ( !picker )
    return NULL;
  
  picker->InitializePickList();
  vtkPropCollection* props = GetRenderer()->GetViewProps();
  if ( props )
  {
    props->InitTraversal();
    vtkProp* prop = props->GetNextProp();
    while ( prop )
    {
      if ( vtkActor::SafeDownCast( prop ) )
        picker->AddPickList( prop );
      prop = props->GetNextProp();
    }
  }
  picker->Pick( posX, GetClientSize().GetHeight() - posY, 0, GetRenderer() );
  if ( pos_out )
    picker->GetPickPosition( pos_out );
  return picker->GetViewProp();
}

int RenderView3D::PickCell( vtkProp* prop, int posX, int posY, double* pos_out )
{
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( this->GetPicker() );
  if ( !picker )
    return -1; 
  
  picker->InitializePickList();
  picker->AddPickList( prop );
  picker->Pick( posX, GetClientSize().GetHeight() - posY, 0, GetRenderer() );
  if ( pos_out )
    picker->GetPickPosition( pos_out );
  return picker->GetCellId();
}

bool RenderView3D::InitializeSelectRegion( int posX, int posY )
{
  double pos[3];
  vtkProp* prop = this->PickProp( posX, posY, pos );  
  if ( !prop )
    return false;
  
  LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerMRI* mri = NULL;
  for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri_temp = (LayerMRI*)lc_mri->GetLayer(i);
    if ( mri_temp->HasProp( prop ) && mri_temp->GetProperties()->GetShowAsContour() )
    {
      mri = mri_temp;
      break;
    }
  }
  
  if ( !mri )
    return false;
  
  lc_mri->SetActiveLayer( mri );
  SurfaceRegion* reg = mri->CreateNewSurfaceRegion( pos );
  if ( reg )
    this->SendBroadcast( "SurfaceRegionSelected", reg );
  return true;
}
  
bool RenderView3D::PickSelectRegion( int nId )
{
  LayerCollection* lc_mri = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerMRI* mri = NULL;
  for ( int i = 0; i < lc_mri->GetNumberOfLayers(); i++ )
  {
    LayerMRI* mri_temp = (LayerMRI*)lc_mri->GetLayer(i);
    if ( mri_temp->GetProperties()->GetShowAsContour() )
    {
      mri = mri_temp;
      break;
    }
  }
  
  if ( !mri )
    return false;
  
  SurfaceRegion* reg = mri->SelectSurfaceRegion( nId );
  if ( reg )
    this->SendBroadcast( "SurfaceRegionSelected", reg );
  return true;
}
  
void RenderView3D::AddSelectRegionLoopPoint( int posX, int posY )
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    double pos[3];
    vtkProp* prop = this->PickProp( posX, posY, pos );  
    if ( !prop || !mri->HasProp( prop ) )
      return;
    
    mri->AddSurfaceRegionLoopPoint( pos );
  }
}

void RenderView3D::CloseSelectRegion()
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    mri->CloseSurfaceRegion();
  }
}

void RenderView3D::DeleteCurrentSelectRegion()
{
  LayerMRI* mri = (LayerMRI*)MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" )->GetActiveLayer();
  if ( mri )
  {
    SurfaceRegion* reg = mri->GetCurrentSurfaceRegion();
    if ( mri->DeleteCurrentSurfaceRegion() )
      this->SendBroadcast( "SurfaceRegionRemoved", reg );
  }
}

void RenderView3D::DoUpdateConnectivityDisplay()
{
  ConnectivityData* conn = MainWindow::GetMainWindowPointer()->GetConnectivityData();
  if ( conn->IsValid() && conn->GetDisplayMode() != ConnectivityData::DM_All )
  {
    conn->BuildConnectivityActors();
//    NeedRedraw();
  }
}

void RenderView3D::UpdateCursorRASPosition( int posX, int posY )
{
  m_bToUpdateCursorPosition = true;
  m_nCursorCoord[0] = posX;
  m_nCursorCoord[1] = posY;
}


void RenderView3D::UpdateConnectivityDisplay()
{
  m_bToUpdateConnectivity = true;
}


void RenderView3D::OnInternalIdle()
{
  RenderView::OnInternalIdle();

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
}

void RenderView3D::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "CursorRASPositionChanged" )
  {
    LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
    m_cursor3D->SetPosition( lc->GetCursorRASPosition() );
    UpdateSliceFrames();
  }
  else if ( iMsg == "ConnectivityActorUpdated" )
  {
    m_renderer->ResetCameraClippingRange();
    NeedRedraw();
  }
  else if ( iMsg == "SlicePositionChanged" )
  {  
    UpdateSliceFrames();
    UpdateSurfaceCorrelationData();
  }
  else if ( iMsg == "LayerAdded" || iMsg == "LayerRemoved" || iMsg == "LayerTransformed" )
  {
    UpdateBounds();
    UpdateSliceFrames();
  }
  else if ( iMsg == "SurfaceRegionAdded" || iMsg == "SurfaceRegionRemoved" )
  {
    RefreshAllActors();
  }
  else if ( iMsg == "SurfaceRegionUpdated" || iMsg == "SurfaceRegionColorChanged" )
  {
    NeedRedraw( true );
  }

  RenderView::DoListenToMessage( iMsg, iData, sender );
}

void RenderView3D::ShowVolumeSlice( int nPlane, bool bShow )
{
  m_bSliceVisibility[nPlane] = bShow;
  RefreshAllActors();
}

void RenderView3D::PreScreenshot()
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();

  m_renderer->RemoveAllViewProps();
  lcm->Append3DProps( m_renderer );

  // add coordinate annotation
  SettingsScreenshot s = MainWindow::GetMainWindowPointer()->GetScreenshotSettings();
  if ( !s.HideCursor )
    m_cursor3D->AppendActor( m_renderer );

  MainWindow::GetMainWindowPointer()->GetConnectivityData()->AppendProps( m_renderer );
  
  // add scalar bar
  m_renderer->AddViewProp( m_actorScalarBar );
}

void RenderView3D::PostScreenshot()
{
  RefreshAllActors();
}

void RenderView3D::UpdateScalarBar()
{
  LayerSurface* surf = (LayerSurface*) MainWindow::GetMainWindowPointer()->GetLayerCollection( "Surface" )->GetActiveLayer();
  if ( surf && surf->GetActiveOverlay() )
  {
    m_actorScalarBar->SetLookupTable( surf->GetActiveOverlay()->GetProperties()->GetLookupTable() );
  }
  else
    RenderView::UpdateScalarBar();
}

void RenderView3D::Azimuth( double angle )
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->Azimuth( angle );
  NeedRedraw();
}

bool RenderView3D::UpdateBounds()
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  double bounds[6] = { 1000000, -1000000, 1000000, -1000000, 1000000, -1000000 };
  for ( int n = 0; n < 1; n++ )
  {
    LayerCollection* lc = lcm->GetLayerCollection( (n == 0 ? "MRI" : "Surface") );
    for ( int i = 0; i < lc->GetNumberOfLayers(); i++ )
    {
      double bd[6];
      lc->GetLayer( i )->GetDisplayBounds( bd );
      for ( int j = 0; j < 3; j++ )
      {
        if ( bounds[j*2] > bd[j*2] )
          bounds[j*2] = bd[j*2];
        if ( bounds[j*2+1] < bd[j*2+1] )
          bounds[j*2+1] = bd[j*2+1];
      }
    } 
  }
  for ( int i = 0; i < 6; i++ )
    m_dBounds[i] = bounds[i];
  
  double dMaxLength = 0;
  for ( int i = 0; i < 3; i++ )
  {
    if ( dMaxLength < ( bounds[i*2+1]-bounds[i*2] ) )
      dMaxLength = bounds[i*2+1]-bounds[i*2];
  }
  
  m_dBoundingTolerance = dMaxLength * 0.02;
  UpdateSliceFrames();
  
  return true;
}

void RenderView3D::UpdateSliceFrames()
{
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  LayerCollection* lc = lcm->GetLayerCollection( "MRI" );
  if ( lc->IsEmpty() )
    lc = lcm->GetLayerCollection( "Surface" );
  double* bounds = m_dBounds;
  double* slicepos = lc->GetSlicePosition();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  points->InsertPoint( 0, slicepos[0], bounds[2], bounds[4] );
  points->InsertPoint( 1, slicepos[0], bounds[2], bounds[5] );
  points->InsertPoint( 2, slicepos[0], bounds[3], bounds[5] );
  points->InsertPoint( 3, slicepos[0], bounds[3], bounds[4] );
  vtkIdType ids[5] = { 0, 1, 2, 3, 0 };
  lines->InsertNextCell( 5, ids );
  polys->InsertNextCell( 4, ids );
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
//  polydata->SetPolys( polys );
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[0]->GetMapper())->SetInput( polydata );
  
  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], slicepos[1], bounds[4] );
  points->InsertPoint( 1, bounds[0], slicepos[1], bounds[5] );
  points->InsertPoint( 2, bounds[1], slicepos[1], bounds[5] );
  points->InsertPoint( 3, bounds[1], slicepos[1], bounds[4] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
//  polydata->SetPolys( polys );
  vtkPolyDataMapper::SafeDownCast(m_actorSliceFrames[1]->GetMapper())->SetInput( polydata );
  
  points = vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint( 0, bounds[0], bounds[2], slicepos[2] );
  points->InsertPoint( 1, bounds[0], bounds[3], slicepos[2] );
  points->InsertPoint( 2, bounds[1], bounds[3], slicepos[2] );
  points->InsertPoint( 3, bounds[1], bounds[2], slicepos[2] );
  polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( points );
  polydata->SetLines( lines );
//  polydata->SetPolys( polys );
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
  
  NeedRedraw();
}

void RenderView3D::HighlightSliceFrame( int n )
{
  if ( m_nSliceHighlighted == n )
    return;
  
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
    NeedRedraw();
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
  
  NeedRedraw();
}

#define min(a, b) (a < b ? a : b)
// move slice from (x1, y1) to (x2, y2) in screen coordinate
void RenderView3D::MoveSliceToScreenCoord( int x, int y )
{
  if ( m_nSliceHighlighted < 0 )
    return;
  
  LayerCollectionManager* lcm = MainWindow::GetMainWindowPointer()->GetLayerCollectionManager();
  LayerCollection* lc = lcm->GetLayerCollection( "MRI" );
  if ( lc->IsEmpty() )
    lc = lcm->GetLayerCollection( "Surface" );
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
      if ( min( fabs( pt[1] - bounds[2] ), fabs( pt[1] - bounds[3] ) ) <
           min( fabs( pt[2] - bounds[4] ), fabs( pt[2] - bounds[5] ) ) )
      {
        v[1] = 1;
        if ( fabs( pt[1] - bounds[2] ) < fabs( pt[1] - bounds[3] ) )
          pt[1] = bounds[2];
        else
          pt[1] = bounds[3];
      }
      else
      {
        v[2] = 1;
        if ( fabs( pt[2] - bounds[4] ) < fabs( pt[2] - bounds[5] ) )
          pt[2] = bounds[4];
        else
          pt[2] = bounds[5];
      }
      break;
    case 1:
      if ( min( fabs( pt[0] - bounds[0] ), fabs( pt[0] - bounds[1] ) ) <
           min( fabs( pt[2] - bounds[4] ), fabs( pt[2] - bounds[5] ) ) )
      {
        v[0] = 1;
        if ( fabs( pt[0] - bounds[0] ) < fabs( pt[0] - bounds[1] ) )
          pt[0] = bounds[0];
        else
          pt[0] = bounds[1];
      }
      else
      {
        v[2] = 1;
        if ( fabs( pt[2] - bounds[4] ) < fabs( pt[2] - bounds[5] ) )
          pt[2] = bounds[4];
        else
          pt[2] = bounds[5];
      }
      break;
    case 2:
      if ( min( fabs( pt[0] - bounds[0] ), fabs( pt[0] - bounds[1] ) ) <
           min( fabs( pt[1] - bounds[2] ), fabs( pt[1] - bounds[3] ) ) )
      {
        v[0] = 1;
        if ( fabs( pt[0] - bounds[0] ) < fabs( pt[0] - bounds[1] ) )
          pt[0] = bounds[0];
        else
          pt[0] = bounds[1];
      }
      else
      {
        v[1] = 1;
        if ( fabs( pt[1] - bounds[2] ) < fabs( pt[1] - bounds[3] ) )
          pt[1] = bounds[2];
        else
          pt[1] = bounds[3];
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
    new_pt[m_nSliceHighlighted] = bounds[m_nSliceHighlighted*2];
  else if ( new_pt[m_nSliceHighlighted] > bounds[m_nSliceHighlighted*2+1] )
    new_pt[m_nSliceHighlighted] = bounds[m_nSliceHighlighted*2+1];
  lcm->OffsetSlicePosition( m_nSliceHighlighted, new_pt[m_nSliceHighlighted] - slicepos[m_nSliceHighlighted], false );
  slicepos[m_nSliceHighlighted] = new_pt[m_nSliceHighlighted];
  lc->SetCursorRASPosition( slicepos );
}

bool RenderView3D::PickCroppingBound( int nX, int nY )
{
  vtkProp* prop = PickProp( nX, nY );
  if ( prop && MainWindow::GetMainWindowPointer()->GetVolumeCropper()->PickActiveBound( prop ) )
  {
    NeedRedraw( true );
    return true;
  }
  else
    return false;
}

void RenderView3D::MoveCroppingBound( int nX, int nY )
{
  MainWindow::GetMainWindowPointer()->GetVolumeCropper()
      ->MoveActiveBound( this, nX, nY );
}

