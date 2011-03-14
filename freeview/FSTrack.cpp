/**
 * @file  FSTrack.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.4 $
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
#include "FSTrack.h"
#include "FSVolume.h"
#include "Track.h"

FSTrack::FSTrack(FSVolume* ref, QObject *parent) :
  TrackData(parent),
  m_volumeRef(ref)
{
}

bool FSTrack::LoadFromFile(const QString &filename)
{
  if (!TrackData::LoadFromFile(filename))
  {
    return false;
  }

  if (m_volumeRef)
  {
    for (int i = 0; i < m_tracks.size(); i++)
    {
      for (int j = 0; j < m_tracks[i].nNum; j++)
      {
        m_volumeRef->RASToTarget(m_tracks[i].fPts + j*3, m_tracks[i].fPts + j*3);
      }
    }
  }
  return true;
}
