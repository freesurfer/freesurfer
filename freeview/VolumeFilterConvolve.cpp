/**
 * @brief Base VolumeFilterConvolve class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "VolumeFilterConvolve.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>
#include "ProgressCallback.h"
#include <vtkImageGaussianSmooth.h>
#include <QDebug>
#include "utils.h"


VolumeFilterConvolve::VolumeFilterConvolve( LayerMRI* input, LayerMRI* output, QObject* parent ) :
  VolumeFilter( input, output, parent ),
  m_dSigma( 1.0 )
{
}

bool VolumeFilterConvolve::Execute()
{
  if (true) // m_volumeInput->GetDataType() == MRI_RGB)
  {
    TriggerFakeProgress(100);
    vtkSmartPointer<vtkImageGaussianSmooth> filter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    filter->SetStandardDeviation(m_dSigma);
    filter->SetRadiusFactor(m_nKernelSize);
  #if VTK_MAJOR_VERSION > 5
    filter->SetInputData( m_volumeInput->GetImageData() );
  #else
    filter->SetInput( m_volumeInput->GetImageData() );
  #endif
    filter->Update();
    m_volumeOutput->GetImageData()->DeepCopy( filter->GetOutput() );
  }
  else
  {
    ::SetProgressCallback(ProgressCallback, 0, 50);
    MRI* mri_src = CreateMRIFromVolume( m_volumeInput );
    MRI* mri_g = MRIgaussian1d( m_dSigma, m_nKernelSize );
    if ( !mri_src || !mri_g )
    {
      return false;
    }

    ::SetProgressCallback(ProgressCallback, 50, 60);
    MRI* mri_dest = MRIconvolveGaussian( mri_src, NULL, mri_g );
    if ( !mri_dest )
    {
      return false;
    }

    ::SetProgressCallback(ProgressCallback, 60, 100);
    MapMRIToVolume( mri_dest, m_volumeOutput );
    MRIfree( &mri_src );
    MRIfree( &mri_g );
    MRIfree( &mri_dest );
  }

  return true;
}

