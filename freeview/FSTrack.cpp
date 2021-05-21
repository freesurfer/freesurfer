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
  QStringList list;
  list << filename;
  return LoadFromFiles(list, ref_fn);
}

bool FSTrack::LoadFromFiles(const QStringList &filenames, const QString &ref_fn)
{
  if (!TrackData::LoadFromFiles(filenames))
  {
    return false;
  }

  if (!ref_fn.isEmpty())
  {
    MRI* mri_ref = ::MRIreadHeader(qPrintable(ref_fn), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_ref)
    {
      cout << "Could not read reference volume " << qPrintable(ref_fn) << endl;
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
      pt[0] = m_tracks[i].fPts[j*3]/m_dVoxelSize[0]-0.5;
      pt[1] = m_tracks[i].fPts[j*3+1]/m_dVoxelSize[1]-0.5;
      pt[2] = m_tracks[i].fPts[j*3+2]/m_dVoxelSize[2]-0.5;
      vtkMatrix4x4::MultiplyPoint(mat, pt, pt);
      if (m_volumeRef)
        m_volumeRef->RASToTarget(pt, pt);
      m_tracks[i].fPts[j*3]   = pt[0];
      m_tracks[i].fPts[j*3+1] = pt[1];
      m_tracks[i].fPts[j*3+2] = pt[2];
    }
  }

  m_dRASBounds[0] = m_dRASBounds[2] = m_dRASBounds[4] = 1e10;
  m_dRASBounds[1] = m_dRASBounds[3] = m_dRASBounds[5] = -1e10;
  for (int i = 0; i <=1; i++)
  {
    for (int j = 0; j <= 1; j++)
    {
      for (int k = 0; k <= 1; k++)
      {
        pt[0] = i*m_nDim[0];
        pt[1] = j*m_nDim[1];
        pt[2] = k*m_nDim[2];
        vtkMatrix4x4::MultiplyPoint(mat, pt, pt);
        if (pt[0] < m_dRASBounds[0])
          m_dRASBounds[0] = pt[0];
        else if (pt[0] > m_dRASBounds[1])
          m_dRASBounds[1] = pt[0];

        if (pt[1] < m_dRASBounds[2])
          m_dRASBounds[2] = pt[1];
        else if (pt[1] > m_dRASBounds[3])
          m_dRASBounds[3] = pt[1];

        if (pt[2] < m_dRASBounds[4])
          m_dRASBounds[4] = pt[2];
        else if (pt[2] > m_dRASBounds[5])
          m_dRASBounds[5] = pt[2];
      }
    }
  }

  return true;
}

void FSTrack::GetRASBounds(double bounds[])
{
  for (int i = 0; i < 6; i++)
    bounds[i] = m_dRASBounds[i];
}
