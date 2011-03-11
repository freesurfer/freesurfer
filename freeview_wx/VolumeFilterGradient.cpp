/**
 * @file  VolumeFilterGradient.cpp
 * @brief Base VolumeFilterGradient class. 
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

#include "VolumeFilterGradient.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkImageShiftScale.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageAnisotropicDiffusion3D.h>
#include <vtkImageGaussianSmooth.h>

VolumeFilterGradient::VolumeFilterGradient( LayerMRI* input, LayerMRI* output ) : 
    VolumeFilter( input, output ),
    m_bSmoothing( false ),
    m_dSD( 1.0 )
{
}

VolumeFilterGradient::~VolumeFilterGradient()
{
}

bool VolumeFilterGradient::Execute()
{
  
  /*
  vtkSmartPointer<vtkImageClip> clip = vtkSmartPointer<vtkImageClip>::New();
  clip->SetInput( image_in );
  int ext[6];
  image_in->GetExtent( ext );
  ext[m_nPlane*2] = ext[m_nPlane*2 + 1] = m_nSlice;
  clip->SetOutputWholeExtent( ext );
  clip->ClipDataOn();
  clip->ReleaseDataFlagOff();
  
  vtkSmartPointer<vtkImageAnisotropicDiffusion3D> smooth = vtkSmartPointer<vtkImageAnisotropicDiffusion3D>::New();
  smooth->SetInput( m_MRIInput->GetImageData() );
  smooth->SetDiffusionFactor( 0.75 );
  smooth->SetDiffusionThreshold( 50.0 );
  smooth->SetNumberOfIterations( 5 );
*/
  
  vtkSmartPointer<vtkImageGradientMagnitude> grad = vtkSmartPointer<vtkImageGradientMagnitude>::New();
  grad->SetDimensionality( 3 );
  grad->HandleBoundariesOn();
  
  if ( m_bSmoothing )
  {
    vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetInput( m_volumeInput->GetImageData() );
    smooth->SetStandardDeviations( m_dSD, m_dSD, m_dSD );
    smooth->SetRadiusFactors( 1, 1, 1 );
    grad->SetInputConnection( smooth->GetOutputPort() );
  }
  else
    grad->SetInput( m_volumeInput->GetImageData() );
  
  grad->Update();
  m_volumeOutput->GetImageData()->DeepCopy( grad->GetOutput() );
  
  return true;
}

