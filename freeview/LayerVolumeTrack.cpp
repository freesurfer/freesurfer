/**
 * @file  LayerVolumeTrack.cpp
 * @brief Layer class for tracks saved in a multi-frame volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/23 21:36:50 $
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
#include <QDebug>

LayerVolumeTrack::LayerVolumeTrack( LayerMRI* ref, QObject* parent ) :
    LayerMRI( ref, parent ),
    m_ctabStripped(NULL)
{
  m_strTypeNames.push_back( "VolumeTrack" );
  SetEditable(false);
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
    cerr << "Did not find color table in track volume.\n";
  else
  {
    MRI* mri = m_volumeSource->GetMRI();
    QList<int> list;
    for (int i = 0; i < mri->nframes; i++)
      list << mri->frames[i].label;
    m_ctabStripped = CTABdeepCopy(m_volumeSource->GetEmbeddedColorTable());
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
          strcpy(m_ctabStripped->entries[i]->name, mri->frames[list.indexOf(i)].name);
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
  for (int i = 0; i < mri->nframes; i++)
  {
    vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract->SetComponents(i);
    extract->SetInput(m_imageData);
    extract->Update();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    actor->SetMapper(mapper);
    MyVTKUtils::BuildContourActor(extract->GetOutput(), mri->frames[i].thresh, 10000000000, actor);
    mapper->ScalarVisibilityOff();
    if (ct)
    {
      int nr, ng, nb;
      CTABrgbAtIndexi( ct, mri->frames[i].label, &nr, &ng, &nb );
      actor->GetProperty()->SetColor(nr/255.0, ng/255.0, nb/255.0);
    }
    m_actors << actor;
  }
}

void LayerVolumeTrack::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
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
    m_actors[i]->SetVisibility(bVisible);
  }
  LayerMRI::SetVisible(bVisible);
}

void LayerVolumeTrack::UpdateOpacity()
{
  for (int i = 0; i < m_actors.size(); i++)
    m_actors[i]->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  LayerMRI::UpdateOpacity();
}
