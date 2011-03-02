/**
 * @file  RenderView2D.h
 * @brief View class for 2D image rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.21 $
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

#ifndef RenderView2D_h
#define RenderView2D_h

#include "RenderView.h"
#include <vector>

class Annotation2D;
class Cursor2D;
class Interactor2DNavigate;
class Interactor2DMeasure;
class Interactor2DVoxelEdit;
class Interactor2DROIEdit;
class Interactor2DWayPointsEdit;
class Interactor2DCropVolume;
class Region2DRectangle;
class Region2D;
class Contour2D;
class LayerMRI;

class VTK_RENDERING_EXPORT RenderView2D : public RenderView
{
  friend class Interactor2D;

  DECLARE_DYNAMIC_CLASS(RenderView2D)

public:
  RenderView2D( int nPlane = 0 );
  RenderView2D( int nPlane, wxWindow *parent, int id = wxID_ANY );
  virtual ~RenderView2D();

  void OnSize( wxSizeEvent& event );

  static RenderView2D * New();
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void RefreshAllActors();

  virtual void TriggerContextMenu( const wxPoint& pos );

  int GetViewPlane();

  void UpdateViewByWorldCoordinate();

  void UpdateMouseRASPosition( int nX, int nY );
  void UpdateCursorRASPosition( int nX, int nY, bool bConnectPrevious = false );

  void SetInteractionMode( int nMode );

  void MousePositionToRAS( int nX, int nY, double* ras );

  Cursor2D* GetCursor2D()
  {
    return m_cursor2D;
  }

  void Update2DOverlay();
  bool EnsureCursor2DVisible();

  void MoveLeft();
  void MoveRight();
  void MoveUp();
  void MoveDown();
  void MoveSlice( int nStep );
  void ZoomAtCursor( int nX, int nY, bool ZoomIn, double factor = 2 );
  void PanToWorld( double* pos );
  void SyncZoomTo( RenderView2D* view );

  bool SetSliceNumber( int nSliceNumber );
  
  void PreScreenshot();
  void PostScreenshot();

  void ShowCoordinateAnnotation( bool bShow );
  bool GetShowCoordinateAnnotation();
  
  void StartSelection( int nX, int nY );
  void UpdateSelection( int nX, int nY );
  void StopSelection();
  
  Region2D* GetRegion( int nX, int nY, int* index_out = NULL );
  void AddRegion( Region2D* region );
  void DeleteRegion( Region2D* region );
  
  LayerMRI* GetFirstNonLabelVolume();
  
  Contour2D* GetContour2D()
  {
    return m_contour2D;
  }
  
protected:
  void Initialize2D();
  void UpdateAnnotation();
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  int  m_nViewPlane;

  int  m_nCursorPosX;
  int  m_nCursorPosY;

  Annotation2D*   m_annotation2D;
  Cursor2D*       m_cursor2D;
  Contour2D*      m_contour2D;
  Region2DRectangle*     m_selection2D;
  
  Interactor2DNavigate*       m_interactorNavigate;
  Interactor2DMeasure*        m_interactorMeasure;
  Interactor2DVoxelEdit*      m_interactorVoxelEdit;
  Interactor2DROIEdit*        m_interactorROIEdit;
  Interactor2DWayPointsEdit*  m_interactorWayPointsEdit;
  Interactor2DCropVolume*     m_interactorCropVolume;

  std::vector<Region2D*>      m_regions;
  
  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


