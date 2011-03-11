/**
 * @file  VolumeFilterMedian.h
 * @brief Base VolumeFilterMedian class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
 *    $Revision: 1.1 $
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


