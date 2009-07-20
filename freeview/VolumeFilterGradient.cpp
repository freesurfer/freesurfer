/**
 * @file  VolumeFilterGradient.cpp
 * @brief Base VolumeFilterGradient class. 
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/07/20 19:34:09 $
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

void VolumeFilterGradient::Update()
{
  if ( !ReadyToUpdate() )
  {
    cerr << "Volume has been removed. Please start over." << endl;
    return;
  }
  
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
    smooth->SetInput( m_MRIInput->GetImageData() );
    smooth->SetStandardDeviations( m_dSD, m_dSD, m_dSD );
    smooth->SetRadiusFactors( 1, 1, 1 );
    grad->SetInput( smooth->GetOutput() );
  }
  else
    grad->SetInput( m_MRIInput->GetImageData() );
  
  grad->Update();
  m_MRIOutput->GetImageData()->DeepCopy( grad->GetOutput() );
  
  m_MRIOutput->SendBroadcast( "LayerActorUpdated", this );
}

