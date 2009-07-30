/**
 * @file  VolumeFilterGradient.h
 * @brief Base VolumeFilterGradient class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/07/30 16:56:08 $
 *    $Revision: 1.1.2.2 $
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

#ifndef VolumeFilterGradient_h
#define VolumeFilterGradient_h

#include "VolumeFilter.h"
#include <string>
#include <vector>

class LayerMRI;

class VolumeFilterGradient : public VolumeFilter
{
public:
  VolumeFilterGradient( LayerMRI* input, LayerMRI* output );
  virtual ~VolumeFilterGradient();

  void SetSmoothing( bool bSmooth )
  {
    m_bSmoothing = bSmooth;
  }
  
  void SetStandardDeviation( double sd )
  {
    m_dSD = sd;
  }
  
  void Update();

private:
  bool    m_bSmoothing;
  double  m_dSD;    
};

#endif


