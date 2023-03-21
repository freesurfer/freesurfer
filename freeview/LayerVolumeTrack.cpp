/**
 * @brief Layer class for tracks saved in a multi-frame volume.
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

#include "LayerVolumeTrack.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "FSVolume.h"
#include "vtkImageExtractComponents.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "LayerPropertyMRI.h"
#include "vtkProperty.h"
#include "MyVTKUtils.h"
#include <QTimer>
#include <QDebug>


#include "cma.h"


LayerVolumeTrack::LayerVolumeTrack( LayerMRI* ref, QObject* parent ) :
  LayerMRI( ref, parent ),
  m_ctabStripped(NULL)
{
  m_strTypeNames.push_back( "VolumeTrack" );
  SetEditable(false);
  connect(GetProperty(), SIGNAL(ContourSmoothIterationChanged(int)), this, SLOT(RebuildActors()));
}

LayerVolumeTrack::~LayerVolumeTrack()
{
  if (m_ctabStripped)
    CTABfree(&m_ctabStripped);
}

bool LayerVolumeTrack::LoadFromFile()
{
  if ( !LayerMRI::LoadVolumeFromFile() )
  {
    return false;
  }

  UpdateData();

  GetProperty()->SetColorMap(LayerPropertyMRI::Heat);
  return true;
}

void LayerVolumeTrack::UpdateData()
{
  if (!m_volumeSource->GetEmbeddedColorTable())
  {
    cerr << "Did not find color table in track volume.\n";
    m_ctabStripped = CTABdeepCopy(this->GetProperty()->GetLUTCTAB());
  }
  else
    m_ctabStripped = CTABdeepCopy(m_volumeSource->GetEmbeddedColorTable());

  if (m_ctabStripped)
  {
    MRI* mri = m_volumeSource->GetMRI();
    QList<int> list;
    for (int i = 0; i < mri->nframes; i++)
      list << mri->frames[i].label;

    int nTotalCount;
    int nValid = 0;
    CTABgetNumberOfTotalEntries( m_ctabStripped, &nTotalCount );
    for ( int i = 0; i < nTotalCount; i++ )
    {
      CTABisEntryValid( m_ctabStripped, i, &nValid );
      if (nValid)
      {
        if (list.contains(i))
        {
          // update name
          QString name = mri->frames[list.indexOf(i)].name;
          if (name.isEmpty())
            name = cma_label_to_name(i);
          strcpy(m_ctabStripped->entries[i]->name, qPrintable(name));
        }
        else
        {
          // remove this entry
          free(m_ctabStripped->entries[i]);
          m_ctabStripped->entries[i] = 0;
        }
      }
    }
  }

  RebuildActors();
}

void LayerVolumeTrack::RebuildActors()
{
  m_actors.clear();
  MRI* mri = m_volumeSource->GetMRI();
  COLOR_TABLE* ct = m_ctabStripped;
  int nSmoothIterations = GetProperty()->GetContourSmoothIterations();
  for (int i = 0; i < mri->nframes; i++)
  {
    vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract->SetComponents(i);
#if VTK_MAJOR_VERSION > 5
    extract->SetInputData(m_imageData);
#else
    extract->SetInput(m_imageData);
#endif
    extract->Update();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    actor->SetMapper(mapper);
    MyVTKUtils::BuildContourActor(extract->GetOutput(), mri->frames[i].thresh, 1e8, actor, nSmoothIterations, NULL, false, false);
    mapper->ScalarVisibilityOff();
    if (ct)
    {
      int nr, ng, nb;
      CTABrgbAtIndexi( ct, mri->frames[i].label, &nr, &ng, &nb );
      actor->GetProperty()->SetColor(nr/255.0, ng/255.0, nb/255.0);
    }
    m_actors << actor;
  }
  if (m_bVisiblities.isEmpty())
  {
    for (int i = 0; i < m_actors.size(); i++)
      m_bVisiblities << true;
  }
  else
  {
    for (int i = 0; i < m_actors.size(); i++)
      m_actors[i]->SetVisibility(m_bVisiblities[i]?1:0);
  }
  emit ActorChanged();
}

void LayerVolumeTrack::RestoreColors()
{
  MRI* mri = m_volumeSource->GetMRI();
  COLOR_TABLE* ct = m_ctabStripped;
  if (!ct)
    return;
  for (int i = 0; i < mri->nframes; i++)
  {
    int nr, ng, nb;
    CTABrgbAtIndexi( ct, mri->frames[i].label, &nr, &ng, &nb );
    m_actors[i]->GetProperty()->SetColor(nr/255.0, ng/255.0, nb/255.0);
  }
  emit ActorUpdated();
}

void LayerVolumeTrack::SetThreshold(int nLabel, double th)
{
  MRI* mri = m_volumeSource->GetMRI();
  for (int i = 0; i < mri->nframes; i++)
  {
    if (nLabel == mri->frames[i].label)
    {
      mri->frames[i].thresh = th;
      UpdateFrameActor(i);
      GetProperty()->SetHeatScaleMinThreshold(th);
      emit ActorUpdated();
      return;
    }
  }
}

int LayerVolumeTrack::GetFrameLabel(int nFrame)
{
  MRI* mri = m_volumeSource->GetMRI();
  return mri->frames[nFrame].label;
}

void LayerVolumeTrack::UpdateFrameActor(int n)
{
  vtkActor* actor = m_actors[n];
  vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
  extract->SetComponents(n);
#if VTK_MAJOR_VERSION > 5
  extract->SetInputData(m_imageData);
#else
  extract->SetInput(m_imageData);
#endif
  extract->Update();
  MRI* mri = m_volumeSource->GetMRI();
  MyVTKUtils::BuildContourActor(extract->GetOutput(), mri->frames[n].thresh, 1e8, actor);
  actor->GetMapper()->ScalarVisibilityOff();
}

double LayerVolumeTrack::GetThreshold(int nLabel)
{
  MRI* mri = m_volumeSource->GetMRI();
  for (int i = 0; i < mri->nframes; i++)
  {
    if (nLabel == mri->frames[i].label)
      return mri->frames[i].thresh;
  }
  return 0;
}

void LayerVolumeTrack::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  Q_UNUSED(bPlaneVisibility);
  for (int i = 0; i < m_actors.size(); i++)
    renderer->AddViewProp(m_actors[i]);
}

void LayerVolumeTrack::UpdateColorMap()
{
  LayerMRI::UpdateColorMap();
}

void LayerVolumeTrack::SetVisible(bool bVisible)
{
  for (int i = 0; i < m_actors.size(); i++)
  {
    m_actors[i]->SetVisibility(bVisible?m_bVisiblities[i]:false);
  }
  LayerMRI::SetVisible(bVisible);
}

void LayerVolumeTrack::UpdateOpacity()
{
  for (int i = 0; i < m_actors.size(); i++)
    m_actors[i]->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  LayerMRI::UpdateOpacity();
}

bool LayerVolumeTrack::HasProp(vtkProp *prop)
{
  for (int i = 0; i < m_actors.size(); i++)
  {
    if (m_actors[i].GetPointer() == prop)
      return true;
  }
  return false;
}

QVariantMap LayerVolumeTrack::GetLabelByProp(vtkProp* prop)
{
  QVariantMap map;
  for (int i = 0; i < m_actors.size(); i++)
  {
    if (m_actors[i].GetPointer() == prop)
    {
      MRI* mri = m_volumeSource->GetMRI();
      map["label"] = mri->frames[i].label;
      map["name"] = QString(mri->frames[i].name);
    }
  }
  return map;
}

void LayerVolumeTrack::Highlight(int nLabel)
{
  MRI* mri = m_volumeSource->GetMRI();
  for (int i = 0; i < mri->nframes; i++)
  {
    if (nLabel == mri->frames[i].label)
    {
      this->SetActiveFrame(i);
      m_actors[i]->GetProperty()->SetColor(1.5, 1.5, 1.5);  // let it over-flow
      emit ActorUpdated();
      QTimer::singleShot(300, this, SLOT(RestoreColors()));
    }
  }
}

void LayerVolumeTrack::ShowAllLabels(bool bShow)
{
  for (int i = 0; i < m_actors.size(); i++)
  {
    m_actors[i]->SetVisibility(bShow?1:0);
    m_bVisiblities[i] = bShow;
  }
  emit ActorUpdated();
}

void LayerVolumeTrack::SetLabelVisible(int nLabel, bool bVisible)
{
  MRI* mri = m_volumeSource->GetMRI();
  for (int i = 0; i < mri->nframes; i++)
  {
    if (nLabel == mri->frames[i].label)
    {
      m_actors[i]->SetVisibility(bVisible?1:0);
      m_bVisiblities[i] = bVisible;
      emit ActorUpdated();
    }
  }
}

void LayerVolumeTrack::SetFrameVisible(int nFrame, bool bVisible)
{
    if (m_actors.size() > nFrame)
    {
      m_actors[nFrame]->SetVisibility(bVisible?1:0);
      m_bVisiblities[nFrame] = bVisible;
      emit ActorUpdated();
    }
}

QList<int> LayerVolumeTrack::GetVisibleLabels()
{
  QList<int> list;
  MRI* mri = m_volumeSource->GetMRI();
  for (int i = 0; i < mri->nframes; i++)
  {
    if (m_bVisiblities[i])
    {
      list << mri->frames[i].label;
    }
  }
  return list;
}
