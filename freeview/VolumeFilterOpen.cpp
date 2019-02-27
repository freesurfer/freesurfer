/**
 * @file  VolumeFilterOpen.cpp
 * @brief Base VolumeFilterOpen class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2014/03/01 04:50:50 $
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
 *
 */

#include "VolumeFilterOpen.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include "ProgressCallback.h"



#include "utils.h"


VolumeFilterOpen::VolumeFilterOpen( LayerMRI* input, LayerMRI* output, QObject* parent ) :
  VolumeFilter( input, output, parent )
{
}

bool VolumeFilterOpen::Execute()
{
  ::SetProgressCallback(ProgressCallback, 0, 50);
  MRI* mri_src = CreateMRIFromVolume( m_volumeInput );
  if ( !mri_src )
  {
    return false;
  }

  ::SetProgressCallback(ProgressCallback, 50, 60);
  MRI* mri_dest = MRIopen( mri_src, NULL );
  if ( !mri_dest )
  {
    return false;
  }

  ::SetProgressCallback(ProgressCallback, 60, 100);
  MapMRIToVolume( mri_dest, m_volumeOutput );
  MRIfree( &mri_src );
  MRIfree( &mri_dest );

  return true;
}

