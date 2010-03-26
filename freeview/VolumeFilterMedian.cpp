/**
 * @file  VolumeFilterMedian.cpp
 * @brief Base VolumeFilterMedian class. 
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

#include "VolumeFilterMedian.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageMedian3D.h>

VolumeFilterMedian::VolumeFilterMedian( LayerMRI* input, LayerMRI* output ) : 
    VolumeFilter( input, output )
{
}

VolumeFilterMedian::~VolumeFilterMedian()
{
}

bool VolumeFilterMedian::Execute()
{
  vtkSmartPointer<vtkImageMedian3D> filter = vtkSmartPointer<vtkImageMedian3D>::New();
  filter->SetKernelSize( m_nKernelSize, m_nKernelSize, m_nKernelSize );
  filter->SetInput( m_volumeInput->GetImageData() );  
  filter->Update();
  m_volumeOutput->GetImageData()->DeepCopy( filter->GetOutput() );
  return true;
}

