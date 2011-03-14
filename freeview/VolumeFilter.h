/**
 * @file  VolumeFilter.h
 * @brief Base VolumeFilter class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:59 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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


