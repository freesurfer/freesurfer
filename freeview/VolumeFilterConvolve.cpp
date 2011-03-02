/**
 * @file  VolumeFilterConvolve.cpp
 * @brief Base VolumeFilterConvolve class. 
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

