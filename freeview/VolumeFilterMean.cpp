/**
 * @file  VolumeFilterMean.cpp
 * @brief Base VolumeFilterMean class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:54 $
 *    $Revision: 1.5 $
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

#include "VolumeFilterMean.h"
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>

VolumeFilterMean::VolumeFilterMean( LayerMRI* input, LayerMRI* output, QObject* parent ) :
    VolumeFilter( input, output, parent )
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

