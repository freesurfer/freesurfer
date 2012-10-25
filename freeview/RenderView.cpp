/**
 * @file  RenderView.cpp
 * @brief View class for rendering 2D and 3D actors
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/25 00:36:17 $
 *    $Revision: 1.47 $
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
#include "Interactor.h"
#include "MainWindow.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerPointSet.h"
#include "LayerPropertyPointSet.h"
#include "LayerSurface.h"
#include "SurfaceOverlay.h"
#include "SurfaceOverlayProperty.h"
#include <QTimer>
#include <QApplication>
#include "MyVTKUtils.h"
#include <QDebug>
#include "vtkActor2D.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkProperty2D.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkMath.h"
#include "vtkScalarBarActor.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include <QPainter>
#include <QAction>
#include <vtkCellPicker.h>
#include <vtkRenderWindow.h>

#define SCALE_FACTOR  200

RenderView::RenderView( QWidget* parent ) : GenericRenderView( parent),
  m_bNeedRedraw( false ),
  m_nInteractionMode( IM_Navigate )
{
  m_interactor = new Interactor( this );

  // initialize focus frame
  m_actorFocusFrame = vtkSmartPointer<vtkActor2D>::New();
  m_actorFocusFrame->VisibilityOff();
  vtkSmartPointer<vtkPoints> Pts = vtkSmartPointer<vtkPoints>::New();
  Pts->InsertNextPoint( 0, 0, 0 );
  Pts->InsertNextPoint( 0, 1, 0 );
  Pts->InsertNextPoint( 1, 1, 0 );
  Pts->InsertNextPoint( 1, 0, 0 );
  vtkSmartPointer<vtkCellArray> Lines = vtkSmartPointer<vtkCellArray>::New();
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

  vtkCellPicker* picker = vtkCellPicker::New();
  picker->SetTolerance( 0.0005 );
  picker->PickFromListOn();
  this->GetRenderWindow()->GetInteractor()->SetPicker( picker );
  picker->Delete();

  vtkSmartPointer<vtkPolyData> Grid = vtkSmartPointer<vtkPolyData>::New();
  Grid->SetPoints(Pts);
  Grid->SetLines(Lines);

  vtkSmartPointer<vtkCoordinate> normCoords = vtkSmartPointer<vtkCoordinate>::New();
  normCoords->SetCoordinateSystemToNormalizedDisplay();

  vtkSmartPointer<vtkPolyDataMapper2D> pMapper = vtkSmartPointer<vtkPolyDataMapper2D>::New();
  pMapper->SetInput(Grid);
  pMapper->SetTransformCoordinate(normCoords);

  m_actorFocusFrame->SetMapper(pMapper);
  m_actorFocusFrame->GetProperty()->SetColor( 0.9, 0.9, 0);
  m_actorFocusFrame->GetProperty()->SetLineWidth( 5 );

  m_actorScalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  m_actorScalarBar->SetOrientationToVertical();
  m_actorScalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  m_actorScalarBar->GetPositionCoordinate()->SetValue(0.9, 0.7);
  m_actorScalarBar->SetHeight(0.3);
  m_actorScalarBar->SetWidth(0.09);
  m_actorScalarBar->VisibilityOff();
  m_actorScalarBar->SetLookupTable(vtkSmartPointer<vtkLookupTable>::New());

  this->setMouseTracking( true );
  connect( this, SIGNAL(BackgroundColorChanged(QColor)), this, SLOT(RequestRedraw()));

  QTimer* timer = new QTimer( this );
  connect( timer, SIGNAL(timeout()), this, SLOT(OnIdle()) );
  connect( this, SIGNAL(ActorsUpdated()), this, SLOT(RequestRedraw()), Qt::QueuedConnection );
  timer->start( 50 );
}

int RenderView::GetInteractionMode()
{
  return m_nInteractionMode;
}

void RenderView::SetInteractionMode( int nMode )
{
  m_nInteractionMode = nMode;
}

void RenderView::RequestRedraw(bool bForce)
{
  if ( bForce )
  {
    Render();
  }
  else
  {
    m_bNeedRedraw = true;
  }
}

void RenderView::OnIdle()
{
//   if ( qApp->hasPendingEvents() )
//       return;

  if ( m_bNeedRedraw )
  {
    Render();
    m_bNeedRedraw = false;
  }
}

void RenderView::SetFocusFrameColor( double r, double g, double b )
{
  m_actorFocusFrame->GetProperty()->SetColor( r, g, b);
}

void RenderView::paintEvent(QPaintEvent *event)
{
  GenericRenderView::paintEvent( event );
}

void RenderView::focusInEvent( QFocusEvent* event )
{
  m_actorFocusFrame->VisibilityOn();
  RequestRedraw();
  GenericRenderView::focusInEvent( event );
}

void RenderView::focusOutEvent( QFocusEvent* event )
{
  if ( m_actorFocusFrame->GetVisibility() )
  {
    m_actorFocusFrame->VisibilityOff();
    RequestRedraw();
  }
  GenericRenderView::focusOutEvent( event );
}

void RenderView::mousePressEvent( QMouseEvent* event )
{
  if ( !hasFocus() )
  {
    this->setFocus();
    // if left button down, do not pass along mouse event further down
    //  if ( MainWindow::GetMainWindowPointer()->GetPreviousActiveView() != this && event.LeftDown() )
    //    return;
  }

  if ( m_interactor->ProcessMouseDownEvent( event, this ) )
  {
    GenericRenderView::mousePressEvent( event );
  }
  else
  {
    event->ignore();
  }
}

void RenderView::mouseReleaseEvent( QMouseEvent* event )
{
  if ( m_interactor->ProcessMouseUpEvent( event, this ) )
  {
    GenericRenderView::mouseReleaseEvent( event );
  }
}

void RenderView::mouseMoveEvent( QMouseEvent* event )
{
// if ( FindFocus() != this )
//  this->SetFocus();

  if ( m_interactor->ProcessMouseMoveEvent( event, this ) )
  {
    GenericRenderView::mouseMoveEvent( event );
  }
  else
  {
    event->ignore();
  }

  m_interactor->ProcessPostMouseMoveEvent( event, this );
}

void RenderView::wheelEvent( QWheelEvent* event )
{
  if ( m_interactor->ProcessMouseWheelEvent( event, this ) )
  {
    GenericRenderView::wheelEvent( event );
  }

  m_interactor->ProcessPostMouseWheelEvent( event, this );
}

void RenderView::enterEvent( QEvent* event )
{
  if ( !hasFocus() )
  {
    this->setFocus();
  }

  if ( m_interactor->ProcessMouseEnterEvent( event, this ) )
  {
    GenericRenderView::enterEvent( event );
  }
}

void RenderView::leaveEvent( QEvent* event )
{
  if ( m_interactor->ProcessMouseLeaveEvent( event, this ) )
  {
    GenericRenderView::leaveEvent( event );
  }
}

void RenderView::keyPressEvent( QKeyEvent* event )
{
  if ( m_interactor->ProcessKeyDownEvent( event, this ) )
  {
    GenericRenderView::keyPressEvent( event );
  }
}

void RenderView::keyReleaseEvent( QKeyEvent* event )
{
  if ( m_interactor->ProcessKeyUpEvent( event, this ) )
  {
    GenericRenderView::keyReleaseEvent( event );
  }
}

void RenderView::SetWorldCoordinateInfo( const double* origin, const double* size, bool bResetView )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_dWorldOrigin[i] = origin[i];
    m_dWorldSize[i] = size[i];
  }
  if (bResetView)
    UpdateViewByWorldCoordinate();
}

void RenderView::ViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z )
{
  MyVTKUtils::ViewportToWorld( m_renderer, x, y, world_x, world_y, world_z );
}

void RenderView::NormalizedViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z )
{
  MyVTKUtils::NormalizedViewportToWorld( m_renderer, x, y, world_x, world_y, world_z );
}

void RenderView::WorldToViewport( double world_x, double world_y, double world_z, double& x, double& y, double& z )
{
  MyVTKUtils::WorldToViewport( m_renderer, world_x, world_y, world_z, x, y, z );
}

void RenderView::WorldToScreen( double world_x, double world_y, double world_z, int& x, int& y )
{
  double dx, dy, dz;
  MyVTKUtils::WorldToViewport( m_renderer, world_x, world_y, world_z, dx, dy, dz );
  x = (int)dx;
  y = (int)( this->rect().height()-dy );
}

void RenderView::ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z )
{
  MyVTKUtils::ViewportToWorld( m_renderer, x, rect().height()-y, z, world_x, world_y, world_z );
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
  double scale = qMax( qMax( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] -= v[i] * scale;
    cam_pos[i] -= v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  RequestRedraw();
  emit ViewChanged();
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
  double scale = qMax( qMax( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] += v[i] * scale;
    cam_pos[i] += v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  RequestRedraw();
  emit ViewChanged();
}

void RenderView::MoveUp()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3];
  cam->GetViewUp( v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = qMax( qMax( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] -= v[i] * scale;
    cam_pos[i] -= v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  RequestRedraw();
  emit ViewChanged();
}

void RenderView::MoveDown()
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3];
  cam->GetViewUp( v );
  double focal_pt[3], cam_pos[3];
  cam->GetFocalPoint( focal_pt );
  cam->GetPosition( cam_pos );
  double scale = qMax( qMax( m_dWorldSize[0], m_dWorldSize[1] ), m_dWorldSize[2] ) / SCALE_FACTOR;
  for ( int i = 0; i < 3; i++ )
  {
    focal_pt[i] += v[i] * scale;
    cam_pos[i] += v[i] * scale;
  }
  cam->SetFocalPoint( focal_pt );
  cam->SetPosition( cam_pos );

  RequestRedraw();
  emit ViewChanged();
}

void RenderView::Zoom( double dFactor )
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  cam->Zoom( dFactor );

  emit ViewChanged();
  Render();
}

void RenderView::CenterAtWorldPosition(double *pos)
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double v[3], cam_pos[3];
  cam->GetDirectionOfProjection( v );
  double dist = cam->GetDistance();
  for ( int i = 0; i < 3; i++ )
  {
    cam_pos[i] = pos[i] - v[i] * dist;
  }
  cam->SetFocalPoint( pos );
  cam->SetPosition( cam_pos );

  RequestRedraw();
  emit ViewChanged();
}

int RenderView::GetAction()
{
  return m_interactor->GetAction();
}

void RenderView::SetAction( int nAction )
{
  m_interactor->SetAction( nAction );
}

void RenderView::Reset()
{
  UpdateViewByWorldCoordinate();
  RequestRedraw();
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
  if (m_layerScalarBar.isNull())
  {
    QList<Layer*> layers = MainWindow::GetMainWindow()->GetLayers("MRI");
    foreach (Layer* layer, layers)
    {
      LayerMRI* mri = qobject_cast<LayerMRI*>(layer);
      if (mri && mri->GetProperty()->GetColorMap() != LayerPropertyMRI::LUT)
      {
        m_actorScalarBar->SetLookupTable( mri->GetProperty()->GetActiveLookupTable() );
        m_layerScalarBar = mri;
        break;
      }
    }
  }
  else
  {
    if (m_layerScalarBar->IsTypeOf("MRI"))
    {
      LayerMRI* mri = qobject_cast<LayerMRI*>(m_layerScalarBar);
      m_actorScalarBar->SetLookupTable( mri->GetProperty()->GetActiveLookupTable() );
    }
    else if (m_layerScalarBar->IsTypeOf("PointSet"))
    {
      LayerPointSet* ps = qobject_cast<LayerPointSet*>(m_layerScalarBar);
      if (ps->GetProperty()->GetColorMap() == LayerPropertyPointSet::HeatScale)
      {
        m_actorScalarBar->SetLookupTable(ps->GetProperty()->GetHeatScaleLUT());
      }
    }
    else if (m_layerScalarBar->IsTypeOf("Surface"))
    {
      LayerSurface* surf = qobject_cast<LayerSurface*>(m_layerScalarBar);
      if (surf->GetActiveOverlay())
        m_actorScalarBar->SetLookupTable( surf->GetActiveOverlay()->GetProperty()->GetLookupTable() );
    }
  }
}

void RenderView::SetScalarBarLayer(Layer *layer)
{
  m_layerScalarBar = layer;
  UpdateScalarBar();
  if (!GetShowScalarBar())
    ShowScalarBar(true);
}

void RenderView::SetScalarBarLayer(QAction *act)
{
  Layer* layer = qobject_cast<Layer*>(act->data().value<QObject*>());
  if (layer)
  {
    if (act->isChecked())
      SetScalarBarLayer(layer);
    else
      ShowScalarBar(false);
  }
}

bool RenderView::SaveScreenShot(const QString& filename, bool bAntiAliasing, int nMag)
{
  blockSignals(true);
  RefreshAllActors(true);
  blockSignals(false);
  bool ret = SaveImage(filename, bAntiAliasing, nMag);
  RefreshAllActors(false);
  return ret;
}

int RenderView::PickCell( vtkProp* prop, int posX, int posY, double* pos_out )
{
  vtkCellPicker* picker = vtkCellPicker::SafeDownCast( GetRenderWindow()->GetInteractor()->GetPicker() );
  if ( !picker )
  {
    return -1;
  }

  picker->InitializePickList();
  picker->AddPickList( prop );
  picker->Pick( posX, this->rect().height() - posY, 0, GetRenderer() );
  if ( pos_out )
  {
    picker->GetPickPosition( pos_out );
  }
  return picker->GetCellId();
}
