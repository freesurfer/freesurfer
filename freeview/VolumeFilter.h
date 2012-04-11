/**
 * @file  VolumeFilter.h
 * @brief Base VolumeFilter class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/04/11 19:46:21 $
 *    $Revision: 1.7.2.2 $
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

#ifndef VolumeFilter_h
#define VolumeFilter_h

#include <QObject>
#include "CommonDataStruct.h"

extern "C"
{
#include "mri.h"
}

class LayerMRI;

class VolumeFilter : public QObject
{
  Q_OBJECT
public:
  VolumeFilter( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );
  virtual ~VolumeFilter();

  bool Update();

  bool ReadyToUpdate();

  void SetInputOutputVolumes( LayerMRI* input, LayerMRI* output );

  MRI* CreateMRIFromVolume( LayerMRI* layer );

  void MapMRIToVolume( MRI* mri, LayerMRI* layer );

  int GetKernelSize()
  {
    return m_nKernelSize;
  }

  void SetKernelSize( int nKernelSize )
  {
    m_nKernelSize = nKernelSize;
  }

  virtual QString GetName() = 0;

protected slots:
  void OnLayerObjectDeleted();

protected:
  virtual bool Execute() = 0;

  int         m_nKernelSize;
  LayerMRI*   m_volumeInput;
  LayerMRI*   m_volumeOutput;
};

#endif


