/**
 * @file  FSTrack.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/11/14 16:30:23 $
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
#include "FSTrack.h"
#include "FSVolume.h"
#include "Track.h"
#include <QDebug>
#include <vtkMatrix4x4.h>

FSTrack::FSTrack(FSVolume* ref, QObject *parent) :
  TrackData(parent),
  m_volumeRef(ref)
{
}

bool FSTrack::LoadFromFile(const QString &filename, const QString &ref_fn)
{
  if (!TrackData::LoadFromFile(filename))
  {
    return false;
  }

  if (!ref_fn.isEmpty())
  {
    MRI* mri_ref = ::MRIreadHeader(qPrintable(ref_fn), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_ref)
    {
      qDebug() << QString("Could not read reference volume %1.").arg(ref_fn);
      return false;
    }
    MATRIX* m = MRIgetVoxelToRasXform( mri_ref );
    for ( int i = 0; i < 16; i++ )
    {
      m_dVoxToRas[i/4][i%4] = *MATRIX_RELT(m, (i/4)+1, (i%4)+1);
    }
    ::MRIfree(&mri_ref);
    ::MatrixFree(&m);
  }

  float pt[4] = { 0, 0, 0, 1 };
  double mat[16];
  for (int i = 0; i < 16; i++)
    mat[i] = m_dVoxToRas[i/4][i%4];
  for (int i = 0; i < m_tracks.size(); i++)
  {
    for (int j = 0; j < m_tracks[i].nNum; j++)
    {
      pt[0] = m_tracks[i].fPts[j*3];
      pt[1] = m_tracks[i].fPts[j*3+1];
      pt[2] = m_tracks[i].fPts[j*3+2];
      vtkMatrix4x4::MultiplyPoint(mat, pt, m_tracks[i].fPts + j*3);
      if (m_volumeRef)
        m_volumeRef->RASToTarget(m_tracks[i].fPts + j*3, m_tracks[i].fPts + j*3);
    }
  }
  return true;
}
