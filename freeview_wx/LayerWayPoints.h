/**
 * @file  LayerWayPoints.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
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

#ifndef LayerWayPoints_h
#define LayerWayPoints_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include "FSWayPoints.h"
#include <string>
#include <vector>


class LayerMRI;
class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;
class LayerPropertiesWayPoints;
class wxWindow;
class wxCommandEvent;

class LayerWayPoints : public LayerEditable
{
public:
  LayerWayPoints( LayerMRI* mri, int type = 0 );
  virtual ~LayerWayPoints();

  bool LoadFromFile( std::string filename );
  bool Save();

  bool HasUndo();
  bool HasRedo();

  void Undo();
  void Redo();

  void SaveForUndo();

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  bool HasProp( vtkProp* prop );

  int FindPoint( double* ras, double tolerance = -1 );

  int AddPoint( double* ras, double value = 1 );

  bool RemovePoint( double* ras, double tolerance = -1 );

  bool RemovePoint( int nIndex );

  void UpdatePoint( double* ras, int nIndex, bool bRebuildActors = true );

  void UpdateLabelData();

  LayerMRI* GetReferenceVolume()
  {
    return m_layerRef;
  }

  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

  inline LayerPropertiesWayPoints* GetProperties()
  {
    return (LayerPropertiesWayPoints*)mProperties;
  }

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected:
  void RebuildActors( bool bRebuild3D = true );
  void UpdateColorMap();
  void UpdateOpacity();
  void UpdateScalars();

  virtual void OnSlicePositionChanged( int nPlane );

  std::vector<WayPoint> m_points;
  vtkActor*    m_actorBalls;
  vtkActor*    m_actorSpline;
  vtkActor*    m_actorSlice[3];

  vtkSmartPointer<vtkPolyDataMapper>  m_mapper;

  LayerMRI*    m_layerRef;
  std::vector< std::vector<WayPoint> > m_bufferUndo;
  std::vector< std::vector<WayPoint> > m_bufferRedo;

  FSWayPoints*   m_waypointsSource;
};

#endif


