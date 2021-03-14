/**
 * @brief Base VolumeFilterThreshold class.
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

#include "VolumeFilterThreshold.h"
#include "LayerMRI.h"
#include <math.h>
#include <vtkDataArray.h>
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkImageThreshold.h>
#include <vtkPointData.h>

VolumeFilterThreshold::VolumeFilterThreshold(LayerMRI *input, LayerMRI *output,
                                             QObject *parent)
    : VolumeFilter(input, output, parent), m_bReplaceIn(false),
      m_bReplaceOut(false) {}

bool VolumeFilterThreshold::Execute() {
  TriggerFakeProgress(50);
  vtkSmartPointer<vtkImageThreshold> filter =
      vtkSmartPointer<vtkImageThreshold>::New();
  filter->SetInput(m_volumeInput->GetImageData());
  filter->SetReplaceIn(m_bReplaceIn ? 1 : 0);
  filter->SetReplaceOut(m_bReplaceOut ? 1 : 0);
  filter->SetInValue(m_dInValue);
  filter->SetOutValue(m_dOutValue);
  filter->ThresholdBetween(m_dThreshold[0], m_dThreshold[1]);
  filter->SetOutputScalarType(m_volumeInput->GetImageData()->GetScalarType());
  filter->Update();
  m_volumeOutput->GetImageData()->DeepCopy(filter->GetOutput());

  return true;
}
