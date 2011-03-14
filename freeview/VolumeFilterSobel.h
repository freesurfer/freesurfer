/**
 * @file  VolumeFilterSobel.h
 * @brief Base VolumeFilterSobel class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:59 $
 *    $Revision: 1.3 $
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

#ifndef VolumeFilterSobel_h
#define VolumeFilterSobel_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterSobel : public VolumeFilter
{
public:
  VolumeFilterSobel( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  void SetSmoothing( bool bSmooth )
  {
    m_bSmoothing = bSmooth;
  }
  
  void SetStandardDeviation( double sd )
  {
    m_dSD = sd;
  }

  QString GetName()
  {
    return "Gradient";
  }
  
protected:
  bool Execute();
  
  bool    m_bSmoothing;
  double  m_dSD;    
};

#endif


