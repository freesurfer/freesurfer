/**
 * @file  VolumeFilterGradient.h
 * @brief Base VolumeFilterGradient class. 
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

#ifndef VolumeFilterGradient_h
#define VolumeFilterGradient_h

#include "VolumeFilter.h"

class LayerMRI;

class VolumeFilterGradient : public VolumeFilter
{
public:
  VolumeFilterGradient( LayerMRI* input = 0, LayerMRI* output = 0, QObject* parent = 0 );

  void SetSmoothing( bool bSmooth )
  {
    m_bSmoothing = bSmooth;
  }
  
  bool GetSmoothing()
  {
      return m_bSmoothing;
  }

  void SetStandardDeviation( double sd )
  {
    m_dSD = sd;
  }

  double GetStandardDeviation()
  {
      return m_dSD;
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


