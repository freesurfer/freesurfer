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
#include <QMap>

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
  bool bLabelContour = m_mri->GetProperty()->GetShowAsLabelContour();
  bool bExtractAllRegions = (bLabelContour || m_mri->GetProperty()->GetContourExtractAllRegions());
  bool bUpsampleContour = m_mri->GetProperty()->GetContourUpsample();
  QList<int> labelList;
  int nSmoothFactor = m_mri->GetProperty()->GetContourSmoothIterations();
  if ( m_nSegValue >= 0 )
  {
    dTh1 = m_nSegValue - 0.5;
    dTh2 = m_nSegValue + 0.5;
  }

  vtkSmartPointer<vtkImageData> imagedata = m_mri->GetImageData();
  if (m_mri->GetNumberOfFrames() > 1)
  {
    vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract->SetComponents(m_mri->GetActiveFrame());
#if VTK_MAJOR_VERSION > 5
    extract->SetInputData(m_mri->GetImageData());
#else
    extract->SetInput(m_mri->GetImageData());
#endif
    extract->Update();
    imagedata = extract->GetOutput();
  }
  QMap<int, vtkActor*> map = m_mri->m_labelActors;
  labelList = m_mri->GetAvailableLabels();
  if (bLabelContour && !labelList.isEmpty())
  {
    foreach (int i, labelList)
    {
      if (!map.contains(i))
      {
        vtkActor* actor = vtkActor::New();
#if VTK_MAJOR_VERSION > 5
        actor->ForceOpaqueOn();
#endif
        actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
        actor->GetMapper()->ScalarVisibilityOn();
        MyVTKUtils::BuildLabelContourActor(imagedata, i, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour,
                                           m_mri->GetProperty()->GetShowVoxelizedContour(), m_mri->GetProperty()->GetContourDilateFirst());
        map[i] = actor;
      }
    }
  }
  else
  {
    vtkActor* actor = vtkActor::New();
    actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
    MyVTKUtils::BuildContourActor( imagedata, dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour, m_mri->GetProperty()->GetContourDilateFirst());
    m_mri->m_actorContourTemp = actor;
    actor->Delete();
  }
  m_mri->m_labelActorsTemp = map;

  emit Finished(m_nThreadID);
}
