/**
 * @brief View class for rendering 2D and 3D actors
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
#include "vtkRenderWindowInteractor.h"
#include <QPainter>
#include <QAction>
#include <vtkCellPicker.h>
#include <vtkRenderWindow.h>
#include "MyUtils.h"
#include <QDir>
#include <QFileInfo>

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
#if VTK_MAJOR_VERSION > 5
  pMapper->SetInputData(Grid);
#else
  pMapper->SetInput(Grid);
#endif
  pMapper->SetTransformCoordinate(normCoords);

  m_actorFocusFrame->SetMapper(pMapper);
  m_actorFocusFrame->GetProperty()->SetColor( 0.9, 0.9, 0);
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = devicePixelRatio();
#endif
  m_actorFocusFrame->GetProperty()->SetLineWidth( 5*ratio );

  m_actorScalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
  m_actorScalarBar->SetOrientationToVertical();
  m_actorScalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
  m_actorScalarBar->GetPositionCoordinate()->SetValue(0.9, 0.7);
  m_actorScalarBar->SetHeight(0.3);
  m_actorScalarBar->SetWidth(0.09);
  m_actorScalarBar->VisibilityOff();
  m_actorScalarBar->SetLookupTable(vtkSmartPointer<vtkLookupTable>::New());
  m_actorScalarBar->SetUseOpacity(1);

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

  if ( m_bNeedRedraw && isVisible())
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

  emit MouseIn();

  if ( m_interactor->ProcessMouseEnterEvent( event, this ) )
  {
    GenericRenderView::enterEvent( event );
  }
}

void RenderView::leaveEvent( QEvent* event )
{
  emit MouseOut();

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
#if VTK_MAJOR_VERSION > 7
  if (devicePixelRatio() > 1)
  {
      dx /= devicePixelRatio();
      dy /= devicePixelRatio();
  }
#endif
  x = (int)dx;
  y = (int)( this->rect().height()-dy );
}

void RenderView::ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z )
{
  y = rect().height()-y;
#if VTK_MAJOR_VERSION > 7
  if (devicePixelRatio() > 1)
  {
      x *= devicePixelRatio();
      y *= devicePixelRatio();
  }
#endif
  MyVTKUtils::ViewportToWorld( m_renderer, x, y, z, world_x, world_y, world_z );
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

void RenderView::PanToWorld( double* pos )
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

  ResetCameraClippingRange();
  RequestRedraw();
  emit ViewChanged();
}

void RenderView::AlignViewToNormal(double *v)
{
  vtkCamera* cam = m_renderer->GetActiveCamera();
  double f_pos[3], dist;
  cam->GetFocalPoint(f_pos);
  dist = cam->GetDistance();
  for (int i = 0; i < 3; i++)
    f_pos[i] += dist*v[i];
  cam->SetPosition(f_pos);
  cam->OrthogonalizeViewUp();
  ResetCameraClippingRange();
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
  emit RequestRedraw();
}

void RenderView::SetScalarBarLayer(Layer *layer)
{
  m_layerScalarBar = layer;
  UpdateScalarBar();
}

void RenderView::SetScalarBarLayer(QAction *act)
{
  Layer* layer = qobject_cast<Layer*>(act->data().value<QObject*>());
  if (layer)
  {
    if (act->isChecked())
    {
      SetScalarBarLayer(layer);
      if (!GetShowScalarBar())
        ShowScalarBar(true);
    }
    else
      ShowScalarBar(false);
  }
}

bool RenderView::SaveScreenShot(const QString& filename, bool bAntiAliasing, int nMag, bool bAutoTrim)
{
  blockSignals(true);
  RefreshAllActors(true);
  blockSignals(false);
  QString fn = filename;
//  if (bAutoTrim)
//  {
//    fn = QFileInfo(QDir::temp(), QString::number(qrand()) + "." + QFileInfo(filename).suffix()).absoluteFilePath();
//  }
  bool ret = SaveImage(fn, bAntiAliasing, nMag);
  if (bAutoTrim)
  {
//    system(QString("convert -trim %1 %2").arg(fn).arg(filename).toLatin1().data());
    TrimImageFiles(QStringList(fn));
  }
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
  return picker->GetCellId();
}

void RenderView::mouseDoubleClickEvent(QMouseEvent *e)
{
  emit DoubleClicked();
  e->accept();
}

void RenderView::TrimImageFiles(const QStringList &files)
{
  if (files.isEmpty())
      return;

  int x0 = 1e6, y0 = 1e6, x1 = 0, y1 = 0;
  foreach (QString fn, files)
  {
    QImage image(fn);
    QRgb* rgb = (QRgb*)image.constBits();
    QRgb bg_val = rgb[0];
    for (int i = 0; i < image.height(); i++)
    {
        for (int j = 0; j < image.width(); j++)
        {
            if (i >= y0 || rgb[i*image.width()+j] != bg_val)
            {
                if (i < y0)
                    y0 = i;
                break;
            }
        }
    }

    for (int i = image.height()-1; i >= 0; i--)
    {
        for (int j = 0; j < image.width(); j++)
        {
            if (i <= y1 || rgb[i*image.width()+j] != bg_val)
            {
                if (i > y1)
                    y1 = i;
                break;
            }
        }
    }

    for (int i = 0; i < image.width(); i++)
    {
        for (int j = 0; j < image.height(); j++)
        {
            if (i >= x0 || rgb[j*image.width()+i] != bg_val)
            {
                if (i < x0)
                    x0 = i;
                break;
            }
        }
    }

    for (int i = image.width()-1; i >= 0; i--)
    {
        for (int j = 0; j < image.height(); j++)
        {
            if (i <= x1 || rgb[j*image.width()+i] != bg_val)
            {
                if (i > x1)
                    x1 = i;
                break;
            }
        }
    }
  }

  QRect rc(x0, y0, x1-x0+1, y1-y0+1);
  foreach (QString fn, files)
  {
    QImage image = QImage(fn).copy(rc);
    image.save(fn);
  }
}
