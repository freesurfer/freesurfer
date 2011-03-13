/**
 * @file  LayerTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.2 $
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
class TrackGroup;
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

protected:
  virtual void OnSlicePositionChanged(int nPlane);

  FSTrack*    m_trackData;
  LayerMRI*   m_layerMRIRef;
  QList<TrackGroup*>  m_trackGroups;
};

#endif // LAYERTRACK_H
