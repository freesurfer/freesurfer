/**
 * @file  LayerPointSet.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:49 $
 *    $Revision: 1.1 $
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

#ifndef LayerPointSet_h
#define LayerPointSet_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include "FSPointSet.h"

class LayerMRI;
class vtkActor;
class vtkPolyDataMapper;
class vtkPolyData;
class LayerPropertyPointSet;
class wxWindow;
class wxCommandEvent;

class LayerPointSet : public LayerEditable
{
    Q_OBJECT
public:
  LayerPointSet( LayerMRI* mri, int type = 0, QObject* parent = NULL );
  virtual ~LayerPointSet();

  bool LoadFromFile( const QString& filename );
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

  inline LayerPropertyPointSet* GetProperty()
  {
    return (LayerPropertyPointSet*)mProperty;
  }

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected slots:
  void UpdateColorMap();
  void UpdateOpacity();
  void UpdateScalars();
  void UpdateSnapToVoxelCenter();
  void RebuildActors( bool bRebuild3D = true );

protected:
  virtual void OnSlicePositionChanged( int nPlane );

  PointSet m_points;
  vtkActor*    m_actorBalls;
  vtkActor*    m_actorSpline;
  vtkActor*    m_actorSlice[3];

  vtkSmartPointer<vtkPolyDataMapper>  m_mapper;

  LayerMRI*    m_layerRef;
  QList< PointSet > m_bufferUndo;
  QList< PointSet > m_bufferRedo;

  FSPointSet*   m_pointSetSource;
};

#endif


