/**
 * @file  VolumeFilter.h
 * @brief Base VolumeFilter class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.3 $
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


