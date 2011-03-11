/**
 * @file  RenderView3D.h
 * @brief View class for 2D image rendering.
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

#ifndef RenderView3D_h
#define RenderView3D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"

class Cursor3D;
class vtkActor;
class vtkCubeSource;
class vtkProp;
class Interactor3DNavigate;
class Interactor3DMeasure;
class Interactor3DCropVolume;

class VTK_RENDERING_EXPORT RenderView3D : public RenderView
{
  DECLARE_DYNAMIC_CLASS(RenderView3D)

public:
  RenderView3D();
  RenderView3D(wxWindow *parent, int id = wxID_ANY);
  virtual ~RenderView3D();

  static RenderView3D * New();
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void RefreshAllActors();

  void UpdateMouseRASPosition( int posX, int posY );
  void CancelUpdateMouseRASPosition();

  void UpdateCursorRASPosition( int posX, int posY );
  void UpdateConnectivityDisplay();

  void UpdateViewByWorldCoordinate();

  Cursor3D* GetCursor3D()
  {
    return m_cursor3D;
  }

  void ShowVolumeSlice( int nPlane, bool bShow = true );

  bool GetShowVolumeSlice( int nPlane );

  void UpdateScalarBar();
  
  void Azimuth( double angle );
  
  void SnapToNearestAxis();
  
  inline int GetHighlightedSlice()
  {
    return m_nSliceHighlighted;
  }
  
  void MoveSliceToScreenCoord( int x, int y );
  
  bool GetShowSliceFrames();
  
  void SetShowSliceFrames( bool bShow );
  
  bool InitializeSelectRegion( int posX, int posY );
  
  void AddSelectRegionLoopPoint( int posX, int posY );
  
  void CloseSelectRegion();
  
  void DeleteCurrentSelectRegion();
  
  bool PickSelectRegion( int nId );
  
  void SetInteractionMode( int nMode );
  
  bool SaveAllSurfaceRegions( wxString& fn );
  
  bool PickCroppingBound( int nX, int nY );
  
  void MoveCroppingBound( int nX, int nY );
  
  int PickCell( vtkProp* prop, int posX, int posY, double* pos_out = NULL );
  
  void UpdateSurfaceCorrelationData();
  
protected:
  void OnInternalIdle();
  void DoUpdateRASPosition( int posX, int posY, bool bCursor = false );
  void DoUpdateConnectivityDisplay();
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  void UpdateSliceFrames();
  void HighlightSliceFrame( int n );
  bool UpdateBounds();
  
  void PreScreenshot();
  void PostScreenshot();
  

private:
  void InitializeRenderView3D();
  vtkProp* PickProp( int posX, int posY, double* pos_out = NULL );

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
  
  double  m_dBounds[6];
  double  m_dBoundingTolerance;
  int     m_nSliceHighlighted;

  double  m_dIntersectPoint[3];
  
  Interactor3DNavigate*       m_interactorNavigate;
  Interactor3DMeasure*        m_interactorMeasure;
  Interactor3DCropVolume*     m_interactorCropVolume;
  
  // any class wishing to process wxWindows events must use this macro
  DECLARE_EVENT_TABLE()

};

#endif


