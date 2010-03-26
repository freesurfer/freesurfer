/**
 * @file  VolumeFilterConvolve.cpp
 * @brief Base VolumeFilterConvolve class. 
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

#include "VolumeFilterConvolve.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>

VolumeFilterConvolve::VolumeFilterConvolve( LayerMRI* input, LayerMRI* output ) : 
    VolumeFilter( input, output ),
    m_dSigma( 1.0 )
{
}

VolumeFilterConvolve::~VolumeFilterConvolve()
{
}

bool VolumeFilterConvolve::Execute()
{
  MRI* mri_src = CreateMRIFromVolume( m_volumeInput );
  MRI* mri_g = MRIgaussian1d( m_dSigma, m_nKernelSize );
  if ( !mri_src || !mri_g ) 
    return false;

  MRI* mri_dest = MRIconvolveGaussian( mri_src, NULL, mri_g );
  if ( !mri_dest )
    return false;
  
  MapMRIToVolume( mri_dest, m_volumeOutput );
  MRIfree( &mri_src );
  MRIfree( &mri_g );
  MRIfree( &mri_dest );
  
  return true;
}

