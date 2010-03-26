/**
 * @file  VolumeFilter.h
 * @brief Base VolumeFilter class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/26 19:04:05 $
 *    $Revision: 1.2 $
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

#include "Listener.h"
#include "Broadcaster.h"
#include "CommonDataStruct.h"
#include <string>
#include <vector>

extern "C"
{
#include "mri.h"
}

class LayerMRI;

class VolumeFilter : public Listener, public Broadcaster
{
public:
  VolumeFilter( LayerMRI* input = 0, LayerMRI* output = 0 );
  virtual ~VolumeFilter();

  bool Update();
  
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
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
  
  virtual std::string GetName() = 0;
  
protected:
  virtual bool Execute() = 0;

  int         m_nKernelSize;
  LayerMRI*   m_volumeInput;
  LayerMRI*   m_volumeOutput;
};

#endif


