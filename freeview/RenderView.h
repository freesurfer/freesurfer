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
#ifndef RENDERVIEW_H
#define RENDERVIEW_H

#include "GenericRenderView.h"
#include "vtkSmartPointer.h"
#include <QPointer>

class Interactor;
class QEvent;
class QMouseEvent;
class QKeyEvent;
class QWheelEvent;
class QFocusEvent;
class vtkActor2D;
class vtkScalarBarActor;
class vtkProp;
class Layer;

class RenderView : public GenericRenderView
{
  Q_OBJECT
public:
  RenderView( QWidget* parent = NULL );

  enum InteractionMode { IM_Navigate = 0, IM_Measure, IM_VoxelEdit, IM_ReconEdit,
                         IM_ROIEdit, IM_PointSetEdit, IM_VolumeCrop, IM_SurfaceCut, IM_SurfacePath };

  void SetWorldCoordinateInfo( const double* origin, const double* size, bool bResetView = true );
  virtual void UpdateViewByWorldCoordinate() {}

  int PickCell( vtkProp* prop, int posX, int posY, double* pos_out = NULL );

  int GetInteractionMode();
  virtual void SetInteractionMode( int nMode );

  int GetAction();

  void SetFocusFrameColor( double r, double g, double b );

  void ViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
  void NormalizedViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
  void WorldToViewport( double world_x, double world_y, double world_z, double& x, double& y, double& z );
  void WorldToScreen( double world_x, double world_y, double world_z, int& x, int& y );
  void ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z );

  virtual void focusInEvent     ( QFocusEvent* event);
  virtual void focusOutEvent    ( QFocusEvent* event );
  virtual void mousePressEvent  ( QMouseEvent* event );
  virtual void mouseReleaseEvent  ( QMouseEvent* event );
  virtual void mouseMoveEvent   ( QMouseEvent* event );
  virtual void wheelEvent       ( QWheelEvent* event );
  virtual void enterEvent       ( QEvent* event );
  virtual void leaveEvent       ( QEvent* event );
  virtual void keyPressEvent    ( QKeyEvent* event );
  virtual void keyReleaseEvent  ( QKeyEvent* event );
  virtual void mouseDoubleClickEvent( QMouseEvent* e);

  bool GetShowScalarBar();

  bool SaveScreenShot(const QString& filename, bool bAntiAliasing, int nMag = 1, bool bAutoTrim = false);

  void TrimImageFiles(const QStringList& files);

  virtual void TriggerContextMenu( QMouseEvent* event ) { Q_UNUSED(event); }

signals:
  void ViewChanged();
  void MouseIn();
  void MouseOut();
  void DoubleClicked();

public slots:
  void RequestRedraw( bool bForce = false );
  void MoveUp();
  void MoveDown();
  void MoveLeft();
  void MoveRight();
  void Zoom( double factor );
  void PanToWorld(double* pos);
  void Reset();
  void SetAction( int nAction );
  void ShowScalarBar( bool bShow );
  void SetScalarBarLayer( Layer* layer );
  void SetScalarBarLayer( QAction* act );
  void CenterAtWorldPosition( double* pos );
  void AlignViewToNormal(double* v);
  virtual void UpdateScalarBar();
  void SetParallelProjection(bool bParallel);

protected:
  virtual void paintEvent(QPaintEvent *event);

protected slots:
  virtual void OnIdle();
  virtual void OnSlicePositionChanged(bool bCenterView = false)
  {
    Q_UNUSED(bCenterView);
    RequestRedraw();
  }

protected:
  bool    m_bNeedRedraw;
  double  m_dWorldOrigin[3];
  double  m_dWorldSize[3];

  Interactor*   m_interactor;
  int           m_nInteractionMode;

  vtkSmartPointer<vtkActor2D>   m_actorFocusFrame;
  vtkSmartPointer<vtkScalarBarActor>  m_actorScalarBar;
  QPointer<Layer>        m_layerScalarBar;
};

#endif // RENDERVIEW_H
