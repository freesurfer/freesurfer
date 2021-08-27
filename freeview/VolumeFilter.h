/**
 * @brief Base VolumeFilter class.
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

#ifndef VolumeFilter_h
#define VolumeFilter_h

#include <QObject>
#include "CommonDataStruct.h"



#include "mri.h"


class LayerMRI;
class QTimer;

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

  LayerMRI* GetVolumeInput()
  {
    return m_volumeInput;
  }

  void SetResetWindowLevel()
  {
    m_bResetWindowLevel = true;
  }

  virtual QString GetName() = 0;

signals:
  void Progress(int n);

protected slots:
  void OnLayerObjectDeleted();
  void OnTimeout();
  void TriggerFakeProgress(int interval);

protected:
  virtual bool Execute() = 0;

  int         m_nKernelSize;
  LayerMRI*   m_volumeInput;
  LayerMRI*   m_volumeOutput;
  int         m_nTimerCount;
  QTimer*     m_timerProgress;
  bool        m_bResetWindowLevel;
};

#endif


