/**
 * @file  VolumeFilter.h
 * @brief Base VolumeFilter class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/20 19:34:09 $
 *    $Revision: 1.1 $
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

class LayerMRI;

class VolumeFilter : public Listener, public Broadcaster
{
public:
  VolumeFilter( LayerMRI* input, LayerMRI* output );
  virtual ~VolumeFilter();

  virtual void Update() = 0;

  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
  bool ReadyToUpdate();
  
  void SetVolumes( LayerMRI* input, LayerMRI* output );
  
protected:
  LayerMRI*   m_MRIInput;
  LayerMRI*   m_MRIOutput;
};

#endif


