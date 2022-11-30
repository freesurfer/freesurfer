/**
 * @brief Layer data object for MRI volume.
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
class vtkPoints;
class LayerSurface;

class LayerPointSet : public LayerEditable
{
  Q_OBJECT
public:
  LayerPointSet( LayerMRI* mri, int type = 0, QObject* parent = NULL );
  LayerPointSet( LayerSurface* surf, int type = 0, QObject* parent = NULL );
  virtual ~LayerPointSet();

  bool LoadFromFile( const QString& filename );
  bool LoadFromJsonFile(const QString& filename);
  bool LoadFromString( const QString& content);
  bool Save();
  bool SaveAsJson(const QString& filename);

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

  int AddPoint( double* ras, double value = 1, bool bNotToVoxelCenter = false );

  bool RemovePoint( double* ras, double tolerance = -1 );

  bool RemovePoint( int nIndex );

  void UpdatePoint( int nIndex, double* ras, bool bRebuildActors = true );

  void UpdatePoint( int nIndex, const QString& key, const QVariant& value);

  void UpdateLabelData();

  LayerMRI* GetReferenceVolume()
  {
    return m_layerRef;
  }

  LayerSurface* GetReferenceSurface()
  {
    return m_surfaceRef;
  }

  inline LayerPropertyPointSet* GetProperty()
  {
    return (LayerPropertyPointSet*)mProperty;
  }

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

  std::vector<double> GetPoints();

  void GetPoint(int nIndex, double* pt_out);

  ControlPoint GetPoint(int nIndex);

  FSPointSet* GetPointSetData()
  {
    return m_pointSetSource;
  }

  int GetNumberOfPoints();

  bool GetCentroidPosition(double *pos);

  bool IsEnhanced();

  double GetEndPointDistance();

  vtkPoints* GetSplinedPoints();

  void GetNormalAtPoint(int nIndex, double* vnorm, int nPlane);

  QVariant GetEnhancedData(const QString& key);

  void SetEnhancedData(const QString& key, const QVariant& val);

protected slots:
  void UpdateColorMap();
  void UpdateOpacity();
  void UpdateScalars(vtkPolyData* polydata);
  void UpdateSnapToVoxelCenter();
  void UpdateSplineVisibility();
  void RebuildActors( bool bRebuild3D = true );

signals:
  void PointAdded(int nIndex);
  void PointRemoved(int nIndex);
  void PointSelected(int nIndex);

protected:
  virtual void OnSlicePositionChanged( int nPlane );
  void Initialize(int nType);

  PointSet m_points;
  vtkActor*    m_actorBalls;
  vtkActor*    m_actorSpline;
  vtkActor*    m_actorSlice[3];
  vtkActor*    m_actorSplineSlice[3];

  vtkSmartPointer<vtkPolyDataMapper>  m_mapper;

  LayerMRI*    m_layerRef;
  LayerSurface*   m_surfaceRef;
  QList< PointSet > m_bufferUndo;
  QList< PointSet > m_bufferRedo;

  FSPointSet*   m_pointSetSource;
  QVariantMap   m_mapEnhancedData;

  vtkSmartPointer<vtkPoints>  m_splinedPoints;
};

#endif


