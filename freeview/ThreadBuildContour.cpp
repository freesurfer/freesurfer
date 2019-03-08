/**
 * @file  ThreadBuildContour.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2016/06/24 17:11:04 $
 *    $Revision: 1.8 $
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
  bool bExtractAllRegions = m_mri->GetProperty()->GetContourExtractAllRegions();
  bool bLabelContour = m_mri->GetProperty()->GetShowAsLabelContour();
  bool bUpsampleContour = m_mri->GetProperty()->GetContourUpsample();
  QList<int> labelList;
  if (bLabelContour)
  {
    QStringList list = m_mri->GetProperty()->GetLabelContourRange().split(",", QString::SkipEmptyParts);
    foreach (QString strg, list)
    {
      if (strg.contains("-"))
      {
        QStringList sublist = strg.split("-", QString::SkipEmptyParts);
        if (sublist.size() > 1)
        {
          int n1 = qMax(1, sublist[0].toInt());
          int n2 = sublist[1].toInt();
          for (int i = n1; i <= n2; i++)
          {
            if (!labelList.contains(i))
              labelList << i;
          }
        }
      }
      else
      {
        int n = -1;
        n = strg.toInt();
        if (n > 0 && !labelList.contains(n))
          labelList << n;
      }
    }

    bExtractAllRegions = true;
  }
  int nSmoothFactor = m_mri->GetProperty()->GetContourSmoothIterations();
  if ( m_nSegValue >= 0 )
  {
    dTh1 = m_nSegValue - 0.5;
    dTh2 = m_nSegValue + 0.5;
  }

  vtkImageData* imagedata = m_mri->GetImageData();
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
        actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
        actor->GetMapper()->ScalarVisibilityOn();
        MyVTKUtils::BuildLabelContourActor(imagedata, i, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour, false );
        map[i] = actor;
      }
    }
  }
  else
  {
    vtkActor* actor = vtkActor::New();
    actor->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
    MyVTKUtils::BuildContourActor( imagedata, dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
    m_mri->m_actorContourTemp = actor;
    actor->Delete();
  }
  m_mri->m_labelActorsTemp = map;

  emit Finished(m_nThreadID);
}


void ThreadBuildContour::run_old()
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
  QList<int> labelList;
  if (bLabelContour)
  {
    QStringList list = m_mri->GetProperty()->GetLabelContourRange().split(",", QString::SkipEmptyParts);
    foreach (QString strg, list)
    {
      if (strg.contains("-"))
      {
        QStringList sublist = strg.split("-", QString::SkipEmptyParts);
        if (sublist.size() > 1)
        {
          int n1 = qMax(1, sublist[0].toInt());
          int n2 = sublist[1].toInt();
          for (int i = n1; i <= n2; i++)
          {
            if (!labelList.contains(i))
              labelList << i;
          }
        }
      }
      else
      {
        int n = -1;
        n = strg.toInt();
        if (n > 0 && !labelList.contains(n))
          labelList << n;
      }
    }

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
    {
      if (!labelList.isEmpty())
        MyVTKUtils::BuildLabelContourActor( m_mri->GetImageData(), labelList, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
    }
    else
      MyVTKUtils::BuildContourActor( m_mri->GetImageData(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
  }
  else
  {
    vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract->SetComponents(m_mri->GetActiveFrame());
#if VTK_MAJOR_VERSION > 5
    extract->SetInputData(m_mri->GetImageData());
#else
    extract->SetInput(m_mri->GetImageData());
#endif
    extract->Update();
    if (bLabelContour)
    {
      if (!labelList.isEmpty())
        MyVTKUtils::BuildLabelContourActor( extract->GetOutput(), labelList, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
    }
    else
      MyVTKUtils::BuildContourActor( extract->GetOutput(), dTh1, dTh2, actor, nSmoothFactor, NULL, bExtractAllRegions, bUpsampleContour );
  }
  m_mri->m_actorContourTemp = actor;
  actor->Delete();

  emit Finished(m_nThreadID);
}
