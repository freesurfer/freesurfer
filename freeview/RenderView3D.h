/**
 * @brief 3D view
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
#ifndef RENDERVIEW3D_H
#define RENDERVIEW3D_H

#include "RenderView.h"
#include <vtkSmartPointer.h>
#include <QVariantMap>
#include <QThread>

class vtkActor;
class vtkProp;
class vtkCubeSource;
class vtkCubeAxesActor;
class Cursor3D;
class SurfaceRegion;
class Interactor3DNavigate;
class Interactor3DMeasure;
class Interactor3DROIEdit;
class Interactor3DPointSetEdit;
class Interactor3DVolumeCrop;
class vtkAnnotatedCubeActor;
class Layer;
class LayerSurface;
class LayerPointSet;
class SurfaceROI;
class Interactor3DPathEdit;
class RenderView3D;
class vtkInteractorStyleMyTrackballCamera;
class Region3D;
class SurfacePath;

class RenderView3D : public RenderView
{
  friend class PropPickingThread;

  Q_OBJECT
public:
  RenderView3D( QWidget* parent );

  void SetInteractionMode(int nMode);

  void UpdateViewByWorldCoordinate();

  bool GetShowSliceFrames();

  void CancelUpdateMouseRASPosition()
  {
    m_bToUpdateRASPosition = false;
  }

  inline int GetHighlightedSlice()
  {
    return m_nSliceHighlighted;
  }

  void UpdateConnectivityDisplay()
  {
    m_bToUpdateConnectivity = true;
  }

  void MoveSliceToScreenCoord( int x, int y );

  void UpdateCursorRASPosition( int posX, int posY );
  void UpdateMouseRASPosition( int posX, int posY, bool bSlicePickOnly = false );
  vtkProp* InitializeSelectRegion( int posX, int poboolsY, int nDrawMode );

  void AddSelectRegionLoopPoint( int posX, int posY, vtkProp* prop_in, int nAction );

  void CloseSelectRegion(int nDrawMode);

  bool PickSelectRegion( int nId );

  bool PickCroppingBound( int nX, int nY );
  void MoveCroppingBound( int nX, int nY );

  Cursor3D* GetCursor3D()
  {
    return m_cursor3D;
  }

  Cursor3D* GetInflatedSurfCursor()
  {
    return m_cursorInflatedSurf;
  }

  void TriggerContextMenu( QMouseEvent* event );

  bool GetShowSlices()
  {
    return m_bSliceVisibility[0] || m_bSliceVisibility[1] || m_bSliceVisibility[2];
  }

  SurfaceROI* InitializeSurfaceROI( int posX, int posY );

  void AddSurfaceROIPoint( int posX, int posY );

  int PickCurrentSurfaceVertex(int posX, int posY, LayerSurface* curSurf = NULL);

  int PickCurrentPointSetPoint(int posX, int posY, LayerPointSet* curPointSet = NULL);

  void ShowSlice(int nPlane, bool bshow);

  QVariantMap GetCameraInfo();

  void SetCamera(const QVariantMap& cam);

  void ZoomAtCursor(int x, int y, double factor);

  bool MapInflatedCoords(LayerSurface* surf, double* pos_in, double* pos_out, bool AutoOrient, bool bCursor = true);

  void MapToInflatedCoords(double* pos_in);

  bool GetFocalPointAtCursor()
  {
    return m_bFocalPointAtCursor;
  }

  bool GetShowAxes()
  {
    return m_bShowAxes;
  }

  int GetAxesFlyMode();

signals:
  void SurfaceVertexClicked(LayerSurface* surf);
  void SurfaceRegionSelected(SurfaceRegion*);
  void SurfaceRegionRemoved(SurfaceRegion*);
  void VolumeTrackMouseOver(Layer* layer, const QVariantMap& info);
  void Region3DSelected(Region3D*);
  void Region3DRemoved(Region3D*);

public slots:
  void RefreshAllActors(bool bForScreenShot = false);
  void SetShowSliceFrames( bool bShow );
  void UpdateSliceFrames();
  bool UpdateBounds();
  void SnapToNearestAxis();
  void Rotate90();
  void Rotate180();
  void UpdateSurfaceCorrelationData();
  void SetShowAllSlices(bool bShow);
  void OnShowSlice(bool bShow = true);
  void HideSlices();
  void ResetViewLeft();
  void ResetViewRight();
  void ResetViewSuperior();
  void ResetViewInferior();
  void ResetViewAnterior();
  void ResetViewPosterior();
  void ResetViewLateral();
  void ResetViewMedial();
  void ShowCursor(bool bshow);
  void OnLayerVisibilityChanged();
  void Azimuth(double degrees);
  void Elevation(double degrees);
  void UpdateScalarBar();
  void SetFocalPointAtCursor(bool b);
  void UpdateAxesActor();
  void SetShowAxes(bool b);  
  void DeleteCurrentSelectRegion();
  void SetAxesFlyMode(int n);
  void DeleteCurrent3DRegion();
  void DeleteAll3DRegions();
  void SavePathAsControlPoints();
  void SaveMarksAsControlPoints();

protected:
  void DoUpdateRASPosition( int posX, int posY, bool bCursor = false, bool bSlicePickOnly = false );
  void DoUpdateConnectivityDisplay();
  void SavePathAsControlPoints(SurfacePath* sp);

  void HighlightSliceFrame( int n );

  virtual void OnSlicePositionChanged(bool bCenterView = false);
  virtual void OnIdle();

  vtkProp* PickProp( int posX, int posY, double* pos_out = NULL, vtkPropCollection* props_in = NULL );

private:
  int  m_nPickCoord[2];
  int  m_nCursorCoord[2];
  bool m_bToUpdateRASPosition;
  bool m_bToUpdateCursorPosition;
  bool m_bToUpdateConnectivity;
  bool m_bSlicePickOnly;

  Cursor3D* m_cursor3D;
  Cursor3D* m_cursorInflatedSurf;
  bool m_bSliceVisibility[3];
  vtkSmartPointer<vtkActor> m_actorSliceFrames[3];
  vtkSmartPointer<vtkActor> m_actorSliceBoundingBox[3];
  vtkSmartPointer<vtkCubeSource>  m_cubeSliceBoundingBox[3];
  vtkSmartPointer<vtkAnnotatedCubeActor> m_actorAnnotatedCube;
  vtkSmartPointer<vtkCubeAxesActor>          m_actorAxesActor;
  vtkActor*     m_actorForAxes;

  double  m_dBounds[6];
  double  m_dBoundingTolerance;
  int     m_nSliceHighlighted;

  bool    m_bShowSliceFrames;
  bool    m_bShowAxes;
  bool    m_bShowCursor;
  bool    m_bFocalPointAtCursor;

  double  m_dIntersectPoint[3];
  Interactor3DNavigate*   m_interactorNavigate;
  Interactor3DMeasure*    m_interactorMeasure;
  Interactor3DVolumeCrop* m_interactorVolumeCrop;
  Interactor3DROIEdit*    m_interactorROIEdit;
  Interactor3DPathEdit*   m_interactorPathEdit;
  Interactor3DPointSetEdit* m_interactorPointSetEdit;

  vtkSmartPointer<vtkInteractorStyleMyTrackballCamera>  m_interactorStyle;
};

#endif // RENDERVIEW3D_H
