/**
 * @file  RenderView2D.h
 * @brief 2D slice view
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/23 17:35:44 $
 *    $Revision: 1.32 $
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

  void UpdateMouseRASPosition( int posX, int posY );
  void UpdateCursorRASPosition( int posX, int posY );
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

public slots:
  void RefreshAllActors(bool bForScreenShot = false);
  void StopSelection();
  void UpdateAnnotation();
  void Update2DOverlay();
  void ShowCoordinateAnnotation( bool bShow );

signals:
  void RegionSelected( Region2D* );
  void RegionRemoved( Region2D* );
  void Zooming(RenderView2D* view);
  void LineProfileIdPicked(LayerLineProfile* lp, int nId);

protected slots:
  virtual void OnSlicePositionChanged();
  void SyncZoomTo(RenderView2D* view);
  void OnDuplicateRegion();
  void OnInteractorError(const QString& msg);

protected:
  virtual void resizeEvent(QResizeEvent *event);
  bool EnsureCursor2DVisible();
  void PanToWorld( double* pos );

private:
  int m_nViewPlane;
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
};

#endif // RENDERVIEW2D_H
