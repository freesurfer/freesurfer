/**
 * @file  LayerTrack.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/11/14 16:30:24 $
 *    $Revision: 1.5 $
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
#include "LayerTrack.h"
#include "FSTrack.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "LayerPropertyTrack.h"
#include <QFileInfo>
#include <QDebug>

LayerTrack::LayerTrack(LayerMRI* ref, QObject* parent) : Layer(parent),
  m_trackData(0),
  m_layerMRIRef(ref)
{
  this->m_strTypeNames << "Track";

  mProperty = new LayerPropertyTrack( this );
}

LayerTrack::~LayerTrack()
{
  if (m_trackData)
  {
    delete m_trackData;
  }
}

bool LayerTrack::LoadTrackFromFile()
{
  if (this->m_sFilename.isEmpty())
  {
    return false;
  }
  FSVolume* refVol = NULL;
  if (m_layerMRIRef)
  {
    refVol = m_layerMRIRef->GetSourceVolume();
  }
  m_trackData = new FSTrack(refVol);
  connect(m_trackData, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)));
  if (!m_trackData->LoadFromFile(m_sFilename))
  {
    delete m_trackData;
    m_trackData = 0;
    qDebug() << "Failed to load from file " << m_sFilename << ".";
    return false;
  }
  SetName(QFileInfo(m_sFilename).completeBaseName());

  RebuildActors();

  return true;
}

void LayerTrack::Append2DProps(vtkRenderer *renderer, int nPlane)
{

}

void LayerTrack::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{

}

bool LayerTrack::HasProp(vtkProp *prop)
{
  return true;
}

bool LayerTrack::IsVisible()
{
  return true;
}

void LayerTrack::OnSlicePositionChanged(int nPlane)
{

}

void LayerTrack::RebuildActors()
{

}
