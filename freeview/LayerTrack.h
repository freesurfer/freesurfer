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
  LayerTrack(LayerMRI* ref, QObject* parent = NULL);
  ~LayerTrack();

  bool LoadTrackFromFile();

  void Append2DProps(vtkRenderer *renderer, int nPlane);

  void Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility);

  bool HasProp(vtkProp *prop);

  bool IsVisible();

  inline LayerPropertyTrack* GetProperty()
  {
    return (LayerPropertyTrack*)mProperty;
  }

  virtual void SetVisible( bool bVisible = true );

signals:
  void Progress(int n);

public slots:
  void RebuildActors();
  void UpdateColor(bool emitSignal = true);

protected:
  virtual void OnSlicePositionChanged(int nPlane);

  vtkActor* ConstructActor(vtkPoints* points, vtkCellArray* lines, vtkUnsignedCharArray* scalars);
  void VectorToColor(float* pt1, float* pt2, float* c_out, int nMappingType);

  FSTrack*    m_trackData;
  LayerMRI*   m_layerMRIRef;
  QList<vtkActor*>  m_actors;
};

#endif // LAYERTRACK_H
