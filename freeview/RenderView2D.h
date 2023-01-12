/**
 * @brief 2D slice view
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
#ifndef RENDERVIEW2D_H
#define RENDERVIEW2D_H

#include "RenderView.h"
#include <QList>

class Region2D;
class Region2DRectangle;
class Contour2D;
class Cursor2D;
class Annotation2D;
class Interactor2DNavigate;
class Interactor2DMeasure;
class Interactor2DVoxelEdit;
class Interactor2DROIEdit;
class Interactor2DPointSetEdit;
class Interactor2DVolumeCrop;
class LayerMRI;
class LayerLineProfile;
class LayerPointSet;

class RenderView2D : public RenderView
{
  //    friend class Interactor2D;

  Q_OBJECT
public:
  RenderView2D( QWidget* parent );

  void SetViewPlane( int nPlane );

  int GetViewPlane()
  {
    return m_nViewPlane;
  }

  void UpdateViewByWorldCoordinate();

  Contour2D* GetContour2D()
  {
    return m_contour2D;
  }

  Cursor2D* GetCursor2D()
  {
    return m_cursor2D;
  }

  Annotation2D* GetAnnotation2D()
  {
    return m_annotation2D;
  }

  void UpdateMouseRASPosition( int posX, int posY );
  void UpdateCursorRASPosition( int posX, int posY, bool bSnapToVertex = false );
  void MoveSlice( int nStep );

  void SetInteractionMode( int nMode );

  void MousePositionToRAS( int posX, int posY, double* pos );
  LayerMRI* GetFirstNonLabelVolume();

  void StartSelection( int nX, int nY );
  void UpdateSelection( int nX, int nY );

  Region2D* GetRegion( int nX, int nY, int* index_out = NULL );
  void AddRegion( Region2D* region );
  void DeleteRegion( Region2D* region );

  void EmitZooming()
  {
    emit Zooming(this);
  }

  void EmitRegionSelected(Region2D* reg)
  {
    emit RegionSelected(reg);
  }

  bool GetShowCoordinateAnnotation();

  void ZoomAtCursor( int nX, int nY, double factor);

  bool SetSliceNumber( int nNum );

  Interactor2DNavigate* GetInteractorNavigate()
  {
    return m_interactorNavigate;
  }

  void TriggerContextMenu( QMouseEvent* event );

  bool PickLineProfile(int x, int y);

  bool GetAutoScaleText()
  {
    return m_bAutoScaleText;
  }

  int GetTextSize()
  {
    return m_nTextSize;
  }

  LayerPointSet* PickPointSetAtCursor(int nX, int nY);

  QList<Region2D*> GetRegions()
  {
    return m_regions;
  }

  bool GetNeurologicalView()
  {
    return m_bNeurologicalView;
  }

public slots:
  void RefreshAllActors(bool bForScreenShot = false);
  void StopSelection();
  void UpdateAnnotation();
  void Update2DOverlay();
  void ShowCoordinateAnnotation( bool bShow );
  void CenterAtCursor();
  void SetAutoScaleText(bool b);
  void SetTextSize(int nsize);
  void OnMovePointToLocalMaximum(bool bUseLast = false);
  void OnMovePointToLocalMaximumDefault();
  void OnMoveAllPointsToLocalMaximum();
  void SetNeurologicalView(bool b);
  void OnInsertPointAfter();

signals:
  void RegionSelected( Region2D* );
  void RegionRemoved( Region2D* );
  void Zooming(RenderView2D* view);
  void LineProfileIdPicked(LayerLineProfile* lp, int nId);
  void CursorLocationClicked();
  void PointSetPicked(LayerPointSet* wp, int nIndex);

protected slots:
  virtual void OnSlicePositionChanged(bool bCenterView = false);
  void SyncZoomTo(RenderView2D* view);
  void OnDuplicateRegion();
  void OnInteractorError(const QString& msg);
  void OnCopyVoxelValue();
  void OnCopyLabelVolume();
  void OnCopyRegionValue();

protected:
  virtual void resizeEvent(QResizeEvent *event);
  bool EnsureCursor2DVisible();

private:
  int m_nViewPlane;
  double    m_dPreSlicePosition;
  Cursor2D*       m_cursor2D;
  Contour2D*      m_contour2D;
  Annotation2D*   m_annotation2D;
  Region2DRectangle*    m_selection2D;
  QList<Region2D*>      m_regions;

  Interactor2DNavigate*   m_interactorNavigate;
  Interactor2DMeasure*    m_interactorMeasure;
  Interactor2DVoxelEdit*  m_interactorVoxelEdit;
  Interactor2DROIEdit*    m_interactorROIEdit;
  Interactor2DPointSetEdit*   m_interactorPointSetEdit;
  Interactor2DVolumeCrop* m_interactorVolumeCrop;

  bool      m_bNeurologicalView;
  bool      m_bAutoScaleText;
  int       m_nTextSize;
};

#endif // RENDERVIEW2D_H
