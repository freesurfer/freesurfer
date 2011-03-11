/**
 * @file  VolumeFilterMean.cpp
 * @brief Base VolumeFilterMean class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
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

#include "VolumeFilterMean.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>

VolumeFilterMean::VolumeFilterMean( LayerMRI* input, LayerMRI* output ) : 
    VolumeFilter( input, output )
{
}

VolumeFilterMean::~VolumeFilterMean()
{
}

bool VolumeFilterMean::Execute()
{
  MRI* mri_src = CreateMRIFromVolume( m_volumeInput );
  if ( !mri_src ) 
    return false;
  
  // clone the src first because MRImean does not keep the src data type if let it do it!
  MRI* mri_dest = MRIclone( mri_src, NULL ) ;
  if ( !mri_dest )
    return false;
  
  MRImean( mri_src, mri_dest, m_nKernelSize );
  MapMRIToVolume( mri_dest, m_volumeOutput );
  MRIfree( &mri_src );
  MRIfree( &mri_dest );
  
  return true;
}

