/**
 * @file  RenderView.cpp
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

#include "RenderView.h"
#include "MainWindow.h"
#include <vtkRenderer.h>
#include <vtkActor2D.h>
#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "Interactor.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerCollection.h"
#include "LayerPropertiesMRI.h"
#include <vtkCoordinate.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkProperty2D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCamera.h>
#include <vtkMath.h>
#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>
#include <vtkLookupTable.h>
#include <vtkLightKit.h>
#include "MyUtils.h"


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define SCALE_FACTOR  200


IMPLEMENT_DYNAMIC_CLASS(RenderView, wxVTKRenderWindowInteractor)

BEGIN_EVENT_TABLE(RenderView, wxVTKRenderWindowInteractor)
  EVT_SET_FOCUS   ( RenderView::OnSetFocus )
  EVT_KILL_FOCUS  ( RenderView::OnKillFocus )
  EVT_LEFT_DOWN   ( RenderView::OnButtonDown )
  EVT_MIDDLE_DOWN ( RenderView::OnButtonDown )
  EVT_RIGHT_DOWN  ( RenderView::OnButtonDown )
  EVT_LEFT_UP     ( RenderView::OnButtonUp )
  EVT_MIDDLE_UP   ( RenderView::OnButtonUp )
  EVT_RIGHT_UP    ( RenderView::OnButtonUp )
  EVT_MOTION      ( RenderView::OnMouseMove )
  EVT_MOUSEWHEEL  ( RenderView::OnMouseWheel )
  EVT_ENTER_WINDOW( RenderView::OnMouseEnter )
  EVT_LEAVE_WINDOW( RenderView::OnMouseLeave )
  EVT_KEY_DOWN    ( RenderView::OnKeyDown )
  EVT_KEY_UP      ( RenderView::OnKeyUp )
  EVT_SIZE        ( RenderView::OnSize )
END_EVENT_TABLE()

RenderView::RenderView() : wxVTKRenderWindowInteractor(),
    Listener( "RenderView" ),
    Broadcaster( "RenderView" )
{
  InitializeRenderView();
}

RenderView::RenderView( wxWindow* parent, int id ) :
    wxVTKRenderWindowInteractor( parent, id, wxDefaultPosition,
                                 wxDefaultSize,
                                 wxWANTS_CHARS | wxNO_FULL_REPAINT_ON_RESIZE ),
    Listener( "RenderView" ),
    Broadcaster( "RenderView" )
{
  InitializeRenderView();
}

void RenderView::InitializeRenderView()
{
  vtkInteractorStyleTrackballCamera* istyle = vtkInteractorStyleTrackballCamera::New();
  this->SetInteractorStyle(istyle);
  istyle->Delete();

  m_renderWindow = this->GetRenderWindow();
  m_renderer = vtkRenderer::New();
  m_renderWindow->AddRenderer( m_renderer );
// m_renderer->SetBackground( 0.7, 0.7, 0.9 );
// m_renderWindow->SetDesiredUpdateRate( 5000 );

  m_interactor = NULL;
  m_nInteractionMode = 0;
  m_nRedrawCount = 0;
  m_bDisabled = false;

  // initialize focus frame
  m_actorFocusFrame = vtkActor2D::New();
  m_actorFocusFrame->VisibilityOff();
  vtkPoints* Pts = vtkPoints::New();
  Pts->InsertNextPoint( 0, 0, 0 );
  Pts->InsertNextPoint( 0, 1, 0 );
  Pts->InsertNextPoint( 1, 1, 0 );
  Pts->InsertNextPoint( 1, 0, 0 );
  vtkCellArray* Lines = vtkCellArray::New();
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 0 );
  Lines->InsertCellPoint( 1 );
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 1 );
  Lines->InsertCellPoint( 2 );
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 2 );
  Lines->InsertCellPoint( 3 );
  Lines->InsertNextCell( 2 );
  Lines->InsertCellPoint( 3 );
  Lines->InsertCellPoint( 0 );

  vtkPolyData* Grid = vtkPolyData::New();
  Grid->SetPoints(Pts);
  Grid->SetLines(Lines);
  Pts->Delete();
  Lines->Delete();

  vtkCoordinate* normCoords = vtkCoordinate::New();
  normCoords->SetCoordinateSystemToNormalizedDisplay();

  vtkPolyDataMapper2D* pMapper = vtkPolyDataMapper2D::New();
  pMapper->SetInput(Grid);
  pMapper->SetTransformCoordinate(normCoords);
  Grid->Delete();
  normCoords->Delete();

  m_actorFocusFrame->SetMapper(pMapper);
  m_actorFocusFrame->GetProperty()->SetColor( 0.9, 0.9, 0);
  m_actorFocusFrame->GetProperty()->SetLineWidth( 5 );
  pMapper->Delete();

  // scalar bar actor
  m_actorScalarBar = vtkScalarBarActor::New();
  m_actorScalarBar->SetOrientationToVertical();
  m_actorScalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  m_actorScalarBar->GetPositionCoordinate()->SetValue(0.9, 0.7);

// qDebug() << m_actorScalarBar->GetLabelTextProperty()->GetFontSize();

// m_actorScalarBar->GetTitleTextProperty()->SetFontSize(5);
  m_actorScalarBar->SetHeight(0.3);
  m_actorScalarBar->SetWidth(0.09);
  m_actorScalarBar->VisibilityOff();
  vtkLookupTable* lut = vtkLookupTable::New();
  m_actorScalarBar->SetLookupTable(lut);
  lut->Delete();

  vtkSmartPointer<vtkScalarBarWidget> barWidget = vtkSmartPointer<vtkScalarBarWidget>::New();
  barWidget->SetScalarBarActor( m_actorScalarBar );
  barWidget->RepositionableOn();
  barWidget->SetInteractor( this );

  UseCaptureMouseOn();
  
  // light kit
  vtkSmartPointer<vtkLightKit> lights = vtkSmartPointer<vtkLightKit>::New();
  lights->AddLightsToRenderer( m_renderer );
//  lights->SetKeyLightIntensity( lights->GetKeyLightIntensity() * 1.1 );
//  lights->SetKeyLightAngle( 35, 10 );
//  lights->SetKeyLightWarmth( 0.5 );
}

RenderView* RenderView::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new RenderView;
}

RenderView::~RenderView()
{
  if (m_renderer)
    m_renderer->Delete();

  if ( m_interactor )
    delete m_interactor;

  m_actorFocusFrame->Delete();
  m_actorScalarBar->Delete();
}

void RenderView::OnSetFocus( wxFocusEvent& event )
{
  m_actorFocusFrame->VisibilityOn();
  NeedRedraw();
  event.Skip();
}

void RenderView::OnKillFocus( wxFocusEvent& event )
{
  if ( m_actorFocusFrame->GetVisibility() )
  {
    m_actorFocusFrame->VisibilityOff();
    NeedRedraw();
  }
  event.Skip();
}

void RenderView::SetFocusFrameColor( double r, double g, double b )
{
  m_actorFocusFrame->GetProperty()->SetColor( r, g, b);
}

void RenderView::OnButtonDown( wxMouseEvent& event )
{
  if ( FindFocus() != this )
  {
    this->SetFocus();
    // if left button down, do not pass along mouse event further down
    if ( MainWindow::GetMainWindowPointer()->GetPreviousActiveView() != this && event.LeftDown() )
      return;
  }

  if ( m_interactor->ProcessMouseDownEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnButtonDown( event );
}

void RenderView::OnButtonUp( wxMouseEvent& event )
{
  if ( m_interactor->ProcessMouseUpEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnButtonUp( event );
}

void RenderView::OnMouseMove( wxMouseEvent& event )
{
// if ( FindFocus() != this )
//  this->SetFocus();

  if ( m_interactor->ProcessMouseMoveEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnMotion( event );

  m_interactor->ProcessPostMouseMoveEvent( event, this );
}

void RenderView::OnMouseWheel( wxMouseEvent& event )
{
  if ( m_interactor->ProcessMouseWheelEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnMouseWheel( event );

  m_interactor->ProcessPostMouseWheelEvent( event, this );
}

void RenderView::OnMouseEnter( wxMouseEvent& event )
{
  if ( FindFocus() != this )
  {
    this->SetFocus();
  }

  if ( m_interactor->ProcessMouseEnterEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnEnter( event );
}

void RenderView::OnMouseLeave( wxMouseEvent& event )
{
  if ( m_interactor->ProcessMouseLeaveEvent( event, this ) )
    wxVTKRenderWindowInteractor::OnLeave( event );
}

void RenderView::OnKeyDown( wxKeyEvent& event )
{
  if ( m_interactor->ProcessKeyDownEvent( event, this ) )
    event.Skip();
}

void RenderView::OnKeyUp( wxKeyEvent& event )
{
  if ( m_interactor->ProcessKeyUpEvent( event, this ) )
    event.Skip();
}

void RenderView::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void RenderView::SetWorldCoordinateInfo( const double* origin, const double* size )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldOrigin[i] = origin[i];
    m_dWorldSize[i] = size[i];
  }
  UpdateViewByWorldCoordinate();
}

void RenderView::ResetView()
{
  UpdateViewByWorldCoordinate();
  NeedRedraw();
}

void RenderView::DoListenToMessage ( std::string const iMsg, void* iData, void* sender )
{
  if ( iMsg == "LayerActorUpdated" )
  {
    UpdateScalarBar();
    NeedRedraw();
  }
  else if ( iMsg == "LayerAdded" || iMsg == "LayerMoved" || iMsg == "LayerRemoved" || iMsg == "LayerTransformed" ||
            iMsg == "LayerContourShown" || iMsg == "DisplayModeChanged" )
  {
    UpdateScalarBar();
    RefreshAllActors();
  }
}

int RenderView::GetInteractionMode()
{
  return m_nInteractionMode;
}

void RenderView::SetInteractionMode( int nMode )
{
  m_nInteractionMode = nMode;
}

int RenderView::GetAction()
{
  return m_interactor->GetAction();
}

void RenderView::SetAction( int nAction )
{
  m_interactor->SetAction( nAction );
}

void RenderView::OnSize( wxSizeEvent& event )
{
  wxVTKRenderWindowInteractor::OnSize( event );

#ifdef __WXGTK__
  MainWindow::GetMainWindowPointer()->NeedRedraw();
#endif
}

void RenderView::SetBackgroundColor( const wxColour& color )
{
  m_renderer->SetBackground( (double)color.Red()/255.0, (double)color.Green()/255.0, (double)color.Blue()/255.0 );
}

wxColour RenderView::GetBackgroundColor() const
{
  double* rgb = m_renderer->GetBackground();
  return wxColour( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) );
}

void RenderView::NeedRedraw( bool bForce )
{
  if ( bForce )
  {
    m_nRedrawCount = 0;
    Render();
  }
  else
    m_nRedrawCount = 1;
}

void RenderView::OnInternalIdle()
{
  wxWindow::OnInternalIdle();

  if ( IsShown() && m_nRedrawCount > 0 )
  {
    wxCursor cursor = GetCursor();
    SetCursor( wxCURSOR_WAIT ); 
    Render();
    SetCursor( cursor );
    m_nRedrawCount--;
  }
}

bool RenderView::SaveScreenshot( const wxString& fn, int nMagnification, bool bAntiAliasing )
{
  PreScreenshot();
  bool ret = MyUtils::VTKScreenCapture( GetRenderWindow(), m_renderer, fn.char_str(), bAntiAliasing, nMagnification );
  PostScreenshot();

  return ret;
}


void RenderView::MoveLeft()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double viewup[3], proj[3], v[3];
  cam->GetViewUp( viewup );
  cam->GetDirectionOfProjection( proj );
  vtkMath::Cross( viewup, proj, v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = min( min( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] -= v[i] * scale;
    cam_pos[i] -= v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  NeedRedraw();
}

void RenderView::MoveRight()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double viewup[3], proj[3], v[3];
  cam->GetViewUp( viewup );
  cam->GetDirectionOfProjection( proj );
  vtkMath::Cross( viewup, proj, v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = min( min( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] += v[i] * scale;
    cam_pos[i] += v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  NeedRedraw();
}

void RenderView::MoveUp()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3];
  cam->GetViewUp( v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = min( min( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] -= v[i] * scale;
    cam_pos[i] -= v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  NeedRedraw();
}

void RenderView::MoveDown()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3];
  cam->GetViewUp( v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = min( min( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] += v[i] * scale;
    cam_pos[i] += v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  NeedRedraw();
}

void RenderView::Zoom( double dFactor )
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->Zoom( dFactor );
  
  Render();
}

void RenderView::ShowScalarBar( bool bShow )
{
  m_actorScalarBar->SetVisibility( bShow?1:0 );
}

bool RenderView::GetShowScalarBar()
{
  return m_actorScalarBar->GetVisibility();
}

void RenderView::UpdateScalarBar()
{
  LayerCollection* lc = MainWindow::GetMainWindowPointer()->GetLayerCollection( "MRI" );
  LayerMRI* mri = (LayerMRI*)lc->GetActiveLayer();
  if ( mri )
    m_actorScalarBar->SetLookupTable( mri->GetProperties()->GetActiveLookupTable() );
}

void RenderView::ViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z )
{
  MyUtils::ViewportToWorld( m_renderer, x, y, world_x, world_y, world_z );
}

void RenderView::NormalizedViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z )
{
  MyUtils::NormalizedViewportToWorld( m_renderer, x, y, world_x, world_y, world_z );
}

void RenderView::WorldToViewport( double world_x, double world_y, double world_z, double& x, double& y, double& z )
{
  MyUtils::WorldToViewport( m_renderer, world_x, world_y, world_z, x, y, z );
}

void RenderView::WorldToScreen( double world_x, double world_y, double world_z, int& x, int& y )
{
  double dx, dy, dz;
  MyUtils::WorldToViewport( m_renderer, world_x, world_y, world_z, dx, dy, dz );
  x = (int)dx;
  y = (int)( GetClientSize().GetHeight()-dy );
}

void RenderView::ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z )
{
  MyUtils::ViewportToWorld( m_renderer, x, GetClientSize().GetHeight()-y, z, world_x, world_y, world_z );
}

void RenderView::GetWorldBound( double* bound )
{
  for ( int i = 0; i < 3; i++ )
  {
    bound[i*2] = m_dWorldOrigin[i];
    bound[i*2+1] = m_dWorldOrigin[i] + m_dWorldSize[i];
  }
}
