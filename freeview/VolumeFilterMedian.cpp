/**
 * @brief Base VolumeFilterMedian class.
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

#include "VolumeFilterMedian.h"
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>

VolumeFilterMedian::VolumeFilterMedian( LayerMRI* input, LayerMRI* output, QObject* parent ) :
  VolumeFilter( input, output, parent )
{
}

bool VolumeFilterMedian::Execute()
{
  TriggerFakeProgress(100);
  vtkSmartPointer<vtkImageMedian3D> filter = vtkSmartPointer<vtkImageMedian3D>::New();
  filter->SetKernelSize( m_nKernelSize, m_nKernelSize, m_nKernelSize );
#if VTK_MAJOR_VERSION > 5
  filter->SetInputData( m_volumeInput->GetImageData() );
#else
  filter->SetInput( m_volumeInput->GetImageData() );
#endif
  filter->Update();
  m_volumeOutput->GetImageData()->DeepCopy( filter->GetOutput() );
  return true;
}

