/**
 * @file  RenderView3D.h
 * @brief 3D view
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/23 17:35:44 $
 *    $Revision: 1.37 $
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
#ifndef RENDERVIEW3D_H
#define RENDERVIEW3D_H

#include "RenderView.h"
#include <vtkSmartPointer.h>
#include <QVariantMap>

class vtkActor;
class vtkProp;
class vtkCubeSource;
class Cursor3D;
class SurfaceRegion;
class Interactor3DNavigate;
class Interactor3DMeasure;
class Interactor3DVolumeCrop;
class vtkAnnotatedCubeActor;
class Layer;
class SurfaceROI;

class RenderView3D : public RenderView
{
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
  void UpdateMouseRASPosition( int posX, int posY );
  bool InitializeSelectRegion( int posX, int posY );

  void AddSelectRegionLoopPoint( int posX, int posY );

  void CloseSelectRegion();

  void DeleteCurrentSelectRegion();

  bool PickSelectRegion( int nId );

  bool PickCroppingBound( int nX, int nY );
  void MoveCroppingBound( int nX, int nY );

  Cursor3D* GetCursor3D()
  {
    return m_cursor3D;
  }

  void UpdateScalarBar();

  void TriggerContextMenu( QMouseEvent* event );

  bool GetShowSlices()
  {
    return m_bShowSlices;
  }

  SurfaceROI* InitializeSurfaceROI( int posX, int posY );

  void AddSurfaceROIPoint( int posX, int posY );

signals:
  void SurfaceVertexClicked();
  void SurfaceRegionSelected(SurfaceRegion*);
  void SurfaceRegionRemoved(SurfaceRegion*);
  void VolumeTrackMouseOver(Layer* layer, const QVariantMap& info);

public slots:
  void RefreshAllActors(bool bForScreenShot = false);
  void SetShowSliceFrames( bool bShow );
  void UpdateSliceFrames();
  bool UpdateBounds();
  void SnapToNearestAxis();
  void UpdateSurfaceCorrelationData();
  void SetShowSlices(bool bShow = true);

protected:
  void DoUpdateRASPosition( int posX, int posY, bool bCursor = false );
  void DoUpdateConnectivityDisplay();

  void HighlightSliceFrame( int n );

  virtual void OnSlicePositionChanged();
  virtual void OnIdle();

  vtkProp* PickProp( int posX, int posY, double* pos_out = NULL );

private:
  int  m_nPickCoord[2];
  int  m_nCursorCoord[2];
  bool m_bToUpdateRASPosition;
  bool m_bToUpdateCursorPosition;
  bool m_bToUpdateConnectivity;

  Cursor3D* m_cursor3D;
  bool m_bSliceVisibility[3];
  vtkSmartPointer<vtkActor> m_actorSliceFrames[3];
  vtkSmartPointer<vtkActor> m_actorSliceBoundingBox[3];
  vtkSmartPointer<vtkCubeSource>  m_cubeSliceBoundingBox[3];
  vtkSmartPointer<vtkAnnotatedCubeActor> m_actorAnnotatedCube;

  double  m_dBounds[6];
  double  m_dBoundingTolerance;
  int     m_nSliceHighlighted;

  bool    m_bShowSlices;

  double  m_dIntersectPoint[3];
  Interactor3DNavigate*   m_interactorNavigate;
  Interactor3DMeasure*    m_interactorMeasure;
  Interactor3DVolumeCrop* m_interactorVolumeCrop;
};

#endif // RENDERVIEW3D_H
