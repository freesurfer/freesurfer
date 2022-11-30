/**
 * @brief Collection of layers of the same type.
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

#ifndef LayerCollection_h
#define LayerCollection_h

#include <QObject>
#include <QString>
#include <QList>
#include <QMap>

class Layer;
class vtkRenderer;
class vtkProp;

class LayerCollection : public QObject
{
  Q_OBJECT
public:
  LayerCollection( const QString& type, QObject* parent = NULL );
  virtual ~LayerCollection();

  int GetNumberOfLayers();
  Layer* GetLayer( int i );
  int GetLayerIndex( Layer* layer );

  bool AddLayer( Layer* layer, bool initializeCoordinate = false );
  bool RemoveLayer( Layer* layer, bool deleteObject = true );
  bool RemoveLayers( QList<Layer*> layers);
  bool MoveLayerUp( Layer* layer );
  bool MoveLayerDown( Layer* layer );
  bool MoveToTop( Layer* layer );
  bool CycleLayer( bool bMoveUp = true, bool bChangeActiveLayer = false );
  void ReorderLayers( const QList<Layer*>& layers);
  void UpdateLayerOrder(const QList<int>& layer_ids);

  void Append2DProps( vtkRenderer* renderer, int nImagePlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool Contains( Layer* layer );
  bool IsEmpty();

  void SetActiveLayer( Layer* layer );
  Layer* GetActiveLayer();

  Layer* GetFirstVisibleLayer();

  Layer* GetLayer(const QString& type);
  QList<Layer*> GetLayers(const QString& type);

  Layer* GetLayerByName( const QString& name );

  Layer* GetLayerById( int nId );

  double* GetSlicePosition();
  void GetSlicePosition( double* slicePos );
  bool SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid = true );
  bool OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid = true );
  bool SetSlicePosition( double* slicePos );

  double* GetCurrentRASPosition();
  void GetCurrentRASPosition( double* pos );
  void SetCurrentRASPosition( double* pos );

  double* GetCursorRASPosition();
  void GetCursorRASPosition( double* pos );
  void SetCursorRASPosition( double* pos );

  void GetCurrentRASIndex( int* nIdx );

  double* GetWorldOrigin();
  void GetWorldOrigin( double* dWorldOrigin_out );
  void SetWorldOrigin( double* dWorldOrigin );

  double* GetWorldSize();
  void GetWorldSize( double* dWorldSize_out );
  void SetWorldSize( double* dWorldSize );

  double* GetWorldVoxelSize();
  void GetWorldVoxelSize( double* dVoxelSize_out );
  void SetWorldVoxelSize( double* dVoxelSize );

  void GetWorldCenter( double* pos );

  QList<Layer*> GetLayers();

  Layer* HasProp( vtkProp* prop );

  QString GetType();

signals:
  void ActiveLayerChanged( Layer* );
  void LayerAdded   ( Layer* );
  void LayerRemoved ( Layer* );
  void LayerCycled  ( Layer* );
  void LayerMoved   ( Layer* );
  void LayersReordered();
  void LayerActorUpdated();
  void LayerActorChanged();
  void LayerPropertyChanged();
  void LayerVisibilityChanged();
  void LayerShowInfoChanged();
  void LayerModified();
  void LayerNameChanged();
  void LayerTransformed();
  void MouseRASPositionChanged();
  void CursorRASPositionChanged();

public slots:
  void LockCurrent( bool bLock );
  void MoveLayerUp();
  void MoveLayerDown();
  void SetMouseRASPosition(double x, double y, double z)
  {
    double ras[3] = {x, y, z};
    SetCurrentRASPosition(ras);
  }
  void ClearLayerIndices();
  void Clear();

protected:
  QList<Layer*> m_layers;

  double      m_dSlicePosition[3];
  double      m_dWorldOrigin[3];
  double      m_dWorldSize[3];
  double      m_dWorldVoxelSize[3];

  double      m_dCurrentRASPosition[3];
  int         m_nCurrentRASIndex[3];

  double      m_dCursorRASPosition[3];

  Layer*      m_layerActive;
  QString     m_strType;
};

#endif


