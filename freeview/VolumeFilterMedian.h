/**
 * @file  VolumeFilterMedian.h
 * @brief Base VolumeFilterMedian class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/26 19:04:05 $
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

#ifndef VolumeFilterMedian_h
#define VolumeFilterMedian_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterMedian : public VolumeFilter
{
public:
  VolumeFilterMedian( LayerMRI* input = 0, LayerMRI* output = 0 );
  virtual ~VolumeFilterMedian();

  std::string GetName()
  {
    return "Median";
  }
  
protected:
  bool Execute();
};

#endif


