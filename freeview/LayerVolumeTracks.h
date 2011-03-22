/**
 * @file  LayerVolumeTracks.h
 * @brief Layer class for tracks saved in a multi-frame volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/22 21:21:26 $
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
 *
 */

#ifndef LAYERVOLUMETRACKS_H
#define LAYERVOLUMETRACKS_H

#include "LayerMRI.h"
#include <QList>

class vtkActor;

class LayerVolumeTracks : public LayerMRI
{
  Q_OBJECT
public:
  LayerVolumeTracks( LayerMRI* ref, QObject* parent = NULL );
  virtual ~LayerVolumeTracks();

  bool LoadFromFile();

  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );

protected:
  void UpdateColorMap();
  void UpdateData();

  QList<vtkActor*>  m_actors;
};

#endif // LAYERVOLUMETRACKS_H
