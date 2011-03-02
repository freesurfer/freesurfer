/**
 * @file  VolumeFilterMean.h
 * @brief Base VolumeFilterMean class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
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


