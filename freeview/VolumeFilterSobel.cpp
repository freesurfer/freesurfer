/**
 * @brief Base VolumeFilterSobel class.
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
  TriggerFakeProgress(50);
  vtkSmartPointer<vtkImageSobel3D> filter = vtkSmartPointer<vtkImageSobel3D>::New();
#if VTK_MAJOR_VERSION > 5
  filter->SetInputData( m_volumeInput->GetImageData() );
#else
  filter->SetInput( m_volumeInput->GetImageData() );
#endif
  vtkSmartPointer<vtkImageMagnitude> mag = vtkSmartPointer<vtkImageMagnitude>::New();
  mag->SetInputConnection(filter->GetOutputPort());
  mag->Update();
  double* orig_range = m_volumeInput->GetImageData()->GetPointData()->GetScalars()->GetRange();
  vtkImageData* img = mag->GetOutput();
  double* range = img->GetPointData()->GetScalars()->GetRange();
  double scale = orig_range[1]/range[1];
  if (scale < 0)
  {
    scale = -scale;
  }
  vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
#if VTK_MAJOR_VERSION > 5
  scaler->SetInputData(img);
#else
  scaler->SetInput(img);
#endif
  scaler->SetShift(0);
  scaler->SetScale(scale);
  scaler->SetOutputScalarType(m_volumeInput->GetImageData()->GetScalarType());
  scaler->Update();
  m_volumeOutput->GetImageData()->DeepCopy( scaler->GetOutput() );

  return true;
}

