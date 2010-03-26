/**
 * @file  VolumeFilterMean.h
 * @brief Base VolumeFilterMean class. 
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

#ifndef VolumeFilterMean_h
#define VolumeFilterMean_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterMean : public VolumeFilter
{
public:
  VolumeFilterMean( LayerMRI* input = 0, LayerMRI* output = 0 );
  virtual ~VolumeFilterMean();
  
  std::string GetName()
  {
    return "Mean";
  }
  
protected:
  bool Execute();
};

#endif


