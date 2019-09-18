/**
 * @file  LayerTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/05 20:03:33 $
 *    $Revision: 1.6 $
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
#ifndef LAYERTRACK_H
#define LAYERTRACK_H

#include "Layer.h"
#include <QColor>
#include <QVariantMap>

class FSTrack;
class LayerMRI;
class vtkActor;
class vtkPoints;
class vtkCellArray;
class LayerPropertyTrack;
class vtkUnsignedCharArray;

class LayerTrack : public Layer
{
  Q_OBJECT
public:
  LayerTrack(LayerMRI* ref, QObject* parent = NULL, bool bCluster = false);
  ~LayerTrack();

  bool LoadTrackFromFiles();

  void Append2DProps(vtkRenderer *renderer, int nPlane);

  void Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility);

  bool HasProp(vtkProp *prop);

  bool IsVisible();

  inline LayerPropertyTrack* GetProperty()
  {
    return (LayerPropertyTrack*)mProperty;
  }

  virtual void SetVisible( bool bVisible = true );

  void SetFileName(const QString& filename);

  void SetFileNames(const QStringList& filenames);

  void SetClusterData(const QVariantMap& data);

  bool IsCluster();

  QVariantMap GetClusterData()
  {
    return m_mapCluster;
  }

  bool HasEmbeddedColor();

signals:
  void Progress(int n);

public slots:
  void RebuildActors();
  void UpdateColor(bool emitSignal = true);
  void LoadTrackFromFiles(const QStringList& filenames)
  {
    SetFileNames(filenames);
    LoadTrackFromFiles();
  }
  void UpdateOpacity(double val);

protected:
  virtual void OnSlicePositionChanged(int nPlane);

  vtkActor* ConstructActor(vtkPoints* points, vtkCellArray* lines, vtkUnsignedCharArray* scalars);
  void VectorToColor(float* pt1, float* pt2, float* c_out, int nMappingType);

  FSTrack*    m_trackData;
  LayerMRI*   m_layerMRIRef;
  QList<vtkActor*>  m_actors;
  QStringList m_listFilenames;
  QVariantMap m_mapCluster;
};

#endif // LAYERTRACK_H
