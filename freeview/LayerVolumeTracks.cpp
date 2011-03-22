/**
 * @file  LayerVolumeTracks.cpp
 * @brief Layer class for tracks saved in a multi-frame volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/22 23:38:45 $
 *    $Revision: 1.1.2.3 $
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

#include "LayerVolumeTracks.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "FSVolume.h"
#include <QDebug>

LayerVolumeTracks::LayerVolumeTracks( LayerMRI* ref, QObject* parent ) :
    LayerMRI( ref, parent )
{
  m_strTypeNames.push_back( "VolumeTracks" );
}

LayerVolumeTracks::~LayerVolumeTracks()
{

}

bool LayerVolumeTracks::LoadFromFile()
{
  if ( !LayerMRI::LoadVolumeFromFile() )
  {
    return false;
  }

  UpdateData();

  return true;
}

void LayerVolumeTracks::UpdateData()
{
  if (!m_volumeSource->GetEmbeddedColorTable())
    cerr << "Did not find color table in track volume.\n";
}

void LayerVolumeTracks::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{
  for (int i = 0; i < m_actors.size(); i++)
    renderer->AddViewProp(m_actors[i]);
}

void LayerVolumeTracks::UpdateColorMap()
{
  LayerMRI::UpdateColorMap();
}
