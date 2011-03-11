/**
 * @file  RenderView.h
 * @brief View class for rendering.
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

#ifndef RenderView_h
#define RenderView_h

#include "wxVTKRenderWindowInteractor.h"
#include "Listener.h"
#include "Broadcaster.h"

class vtkRenderer;
class vtkRenderWindow;
class vtkActor2D;
class LayerCollection;
class Interactor;
class vtkScalarBarActor;

class VTK_RENDERING_EXPORT RenderView : public wxVTKRenderWindowInteractor, public Listener, public Broadcaster
{
  DECLARE_DYNAMIC_CLASS( RenderView )

public:
  RenderView();
  RenderView( wxWindow *parent, int id );
  virtual ~RenderView();

  enum InteractionMode { IM_Navigate = 0, IM_Measure, IM_VoxelEdit, IM_ROIEdit, IM_WayPointsEdit, IM_VolumeCrop };
  
  static RenderView * New();
  void PrintSelf( ostream& os, vtkIndent indent );

  void OnSetFocus  ( wxFocusEvent& event);
  void OnKillFocus ( wxFocusEvent& event );
  void OnButtonDown ( wxMouseEvent& event );
  void OnButtonUp  ( wxMouseEvent& event );
  void OnMouseMove ( wxMouseEvent& event );
  void OnMouseWheel ( wxMouseEvent& event );
  void OnMouseEnter ( wxMouseEvent& event );
  void OnMouseLeave ( wxMouseEvent& event );
  void OnKeyDown  ( wxKeyEvent& event );
  void OnKeyUp  ( wxKeyEvent& event );
  void OnSize   ( wxSizeEvent& event );

  void SetWorldCoordinateInfo( const double* origin, const double* size );
  
  void GetWorldBound( double* bound );
  
  virtual void UpdateViewByWorldCoordinate()
  {}

  virtual void RefreshAllActors()
  {}

  virtual void TriggerContextMenu( const wxPoint& pos )
  {}

  virtual void DoListenToMessage( std::string const iMessage, void* iData, void* sender );

  void ViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
  void NormalizedViewportToWorld( double x, double y, double& world_x, double& world_y, double& world_z );
  void WorldToViewport( double world_x, double world_y, double world_z, double& x, double& y, double& z );
  void WorldToScreen( double world_x, double world_y, double world_z, int& x, int& y ); 
  void ScreenToWorld( int x, int y, int z, double& world_x, double& world_y, double& world_z );

  int GetInteractionMode();
  virtual void SetInteractionMode( int nMode );

  int GetAction();
  void SetAction( int nAction );

  wxColour GetBackgroundColor() const;
  void SetBackgroundColor( const wxColour& color );

  virtual void MoveUp();
  virtual void MoveDown();
  virtual void MoveLeft();
  virtual void MoveRight();
  
  void Zoom( double dFactor );

  void NeedRedraw( bool bForce = false );

  void ResetView();

  vtkRenderer* GetRenderer()
  {
    return m_renderer;
  }

  bool SaveScreenshot( const wxString& fn, int nMagnification = 1, bool bAntiAliasing = false );

  virtual void PreScreenshot()
  {}
  virtual void PostScreenshot()
  {}

  void ShowScalarBar( bool bShow );
  bool GetShowScalarBar();

  void SetFocusFrameColor( double r, double g, double b );
  
  virtual void UpdateScalarBar();
  
  void SetRenderDisabled(bool bDisable)
  {
    m_bDisabled = true;
    RefreshAllActors();
  }
  
  bool GetRenderDisabled()
  {
    return m_bDisabled;
  }

protected:
  void InitializeRenderView();
  virtual void OnInternalIdle();

protected:
  vtkRenderer*        m_renderer;
  vtkRenderWindow*    m_renderWindow;
  vtkActor2D*         m_actorFocusFrame;
  vtkScalarBarActor*  m_actorScalarBar;

  double        m_dWorldOrigin[3];
  double        m_dWorldSize[3];

  Interactor*   m_interactor;
  int           m_nInteractionMode;

  int           m_nRedrawCount;
  
  bool          m_bDisabled;

  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


