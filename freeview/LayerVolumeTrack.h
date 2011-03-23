/**
 * @file  LayerVolumeTrack.h
 * @brief Layer class for tracks saved in a multi-frame volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/23 21:36:50 $
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
 *
 */

#ifndef LAYERVOLUMETRACK_H
#define LAYERVOLUMETRACK_H

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <QList>

class vtkActor;

class LayerVolumeTrack : public LayerMRI
{
  Q_OBJECT
public:
  LayerVolumeTrack( LayerMRI* ref, QObject* parent = NULL );
  virtual ~LayerVolumeTrack();

  bool LoadFromFile();

  void SetVisible(bool bVisible);

  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );

  virtual COLOR_TABLE* GetEmbeddedColorTable()
  {
    return m_ctabStripped;
  }

  virtual void UpdateOpacity();

protected:
  void UpdateColorMap();
  void UpdateData();
  void RebuildActors();

  QList< vtkSmartPointer<vtkActor> >  m_actors;
  COLOR_TABLE* m_ctabStripped;
};

#endif // LAYERVOLUMETRACK_H
