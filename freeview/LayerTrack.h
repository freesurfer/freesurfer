/**
 * @file  LayerTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/11/14 16:30:24 $
 *    $Revision: 1.5 $
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

class FSTrack;
class LayerMRI;
class vtkActor;
class LayerPropertyTrack;

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

signals:
  void Progress(int n);

public slots:
  void RebuildActors();

protected:
  virtual void OnSlicePositionChanged(int nPlane);

  FSTrack*    m_trackData;
  LayerMRI*   m_layerMRIRef;
  QList<vtkActor*>  m_actors;
};

#endif // LAYERTRACK_H
