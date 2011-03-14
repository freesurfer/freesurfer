/**
 * @file  VolumeFilterSobel.cpp
 * @brief Base VolumeFilterSobel class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:59 $
 *    $Revision: 1.3 $
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

#include "VolumeFilterSobel.h"
#include <math.h>
#include "LayerMRI.h"
#include <vtkImageData.h>
#include <vtkImageSobel3D.h>
#include <vtkImageMagnitude.h>
#include <vtkImageCast.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkImageShiftScale.h>

VolumeFilterSobel::VolumeFilterSobel( LayerMRI* input, LayerMRI* output, QObject* parent ) :
    VolumeFilter( input, output, parent ),
    m_bSmoothing( false ),
    m_dSD( 1.0 )
{
}

bool VolumeFilterSobel::Execute()
{ 
  vtkSmartPointer<vtkImageSobel3D> filter = vtkSmartPointer<vtkImageSobel3D>::New();
  filter->SetInput( m_volumeInput->GetImageData() );
  vtkSmartPointer<vtkImageMagnitude> mag = vtkSmartPointer<vtkImageMagnitude>::New();
  mag->SetInputConnection(filter->GetOutputPort());
  mag->Update();
  double* orig_range = m_volumeInput->GetImageData()->GetPointData()->GetScalars()->GetRange();
  vtkImageData* img = mag->GetOutput();
  double* range = img->GetPointData()->GetScalars()->GetRange();
  double scale = orig_range[1]/range[1];
  if (scale < 0)
      scale = -scale;
  vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
  scaler->SetInput(img);
  scaler->SetShift(0);
  scaler->SetScale(scale);
  scaler->SetOutputScalarType(m_volumeInput->GetImageData()->GetScalarType());
  scaler->Update();
  m_volumeOutput->GetImageData()->DeepCopy( scaler->GetOutput() );
  
  return true;
}

