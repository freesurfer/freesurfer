/**
 * @file  ThreadBuildContour.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/02/05 20:51:41 $
 *    $Revision: 1.7 $
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
#include "ThreadBuildContour.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerCollection.h"
#include "MainWindow.h"
#include "MyVTKUtils.h"
#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageExtractComponents.h"
#include <QDebug>

ThreadBuildContour::ThreadBuildContour(QObject *parent) :
  QThread(parent),
  m_mri( NULL )
{
}

void ThreadBuildContour::BuildContour( LayerMRI* mri, int nSegValue, int nThreadID )
{
  m_mri = mri;
  m_nSegValue = nSegValue;
  m_nThreadID = nThreadID;
  start();
}

void ThreadBuildContour::run()
{
  if (!m_mri)
  {
    return;
  }
  double dTh1 = m_mri->GetProperty()->GetContourMinThreshold();
  double dTh2 = m_mri->GetProperty()->GetContourMaxThreshold();
  bool bExtractAllRegions = m_mri->GetProperty()->GetContourExtractAllRegions();
  bool bLabelContour = m_mri->GetProperty()->GetShowAsLabelContour();
  bool bUpsampleContour = m_mri->GetProperty()->GetContourUpsample();
  if (bLabelContour)
  {
    m_mri->GetProperty()->GetLabelContourRange(&dTh1, &dTh2);
    bExtractAllRegions = true;
  }
  int nSmoothFactor = m_mri->GetProperty()->GetContourSmoothIterations();
  if ( m_nSegValue >= 0 )
  {
    dTh1 = m_nSegValue - 0.5;
    dTh2 = m_nSegValue + 0.5;
  }

  vtkActor* actor = vtkActor::New();
  actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
  //  int ext[6] = { 370, 520, 0, 140, 0, 280 };
  if (m_mri->GetNumberOfFrames() == 1)
  {
    if (bLabelContour)
      MyVTKUtils::BuildLabelContourActor( m_mri->GetImageData(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
    else
      MyVTKUtils::BuildContourActor( m_mri->GetImageData(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
  }
  else
  {
    vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract->SetComponents(m_mri->GetActiveFrame());
    extract->SetInput(m_mri->GetImageData());
    extract->Update();
    if (bLabelContour)
      MyVTKUtils::BuildLabelContourActor( extract->GetOutput(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
    else
      MyVTKUtils::BuildContourActor( extract->GetOutput(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
  }
  m_mri->m_actorContourTemp = actor;
  actor->Delete();

  emit Finished(m_nThreadID);
}
