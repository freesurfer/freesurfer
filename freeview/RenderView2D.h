/**
 * @file  RenderView2D.h
 * @brief View class for 2D image rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/07 22:05:04 $
 *    $Revision: 1.12 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef RenderView2D_h
#define RenderView2D_h

#include "RenderView.h"

class Annotation2D;
class Cursor2D;
class Interactor2DNavigate;
class Interactor2DVoxelEdit;
class Interactor2DROIEdit;
class Interactor2DWayPointsEdit;

class VTK_RENDERING_EXPORT RenderView2D : public RenderView
{
  friend class Interactor2D;

  DECLARE_DYNAMIC_CLASS(RenderView2D)

public:
  RenderView2D();
  RenderView2D(wxWindow *parent, int id);
  virtual ~RenderView2D();

  void OnSize( wxSizeEvent& event );

  enum InteractionMode { IM_Navigate = 0, IM_VoxelEdit, IM_ROIEdit, IM_WayPointsEdit };

  static RenderView2D * New();
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void RefreshAllActors();

  virtual void TriggerContextMenu( const wxPoint& pos );

  void SetViewPlane( int nPlane );
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

  void UpdateCursor2D();
  bool EnsureCursor2DVisible();

  void MoveLeft();
  void MoveRight();
  void MoveUp();
  void MoveDown();
  void MoveSlice( int nStep );
  void ZoomAtCursor( int nX, int nY, bool ZoomIn, double factor = 2 );
  void PanToWorld( double* pos );
  void SyncZoomTo( RenderView2D* view );

  void PreScreenshot();
  void PostScreenshot();

  void ShowCoordinateAnnotation( bool bShow );
  bool GetShowCoordinateAnnotation();

protected:
  void Initialize2D();
  void UpdateAnnotation();
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  int  m_nViewPlane;

  int  m_nCursorPosX;
  int  m_nCursorPosY;

  Annotation2D*  m_annotation2D;
  Cursor2D*   m_cursor2D;
  Interactor2DNavigate* m_interactorNavigate;
  Interactor2DVoxelEdit* m_interactorVoxelEdit;
  Interactor2DROIEdit* m_interactorROIEdit;
  Interactor2DWayPointsEdit* m_interactorWayPointsEdit;

  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


