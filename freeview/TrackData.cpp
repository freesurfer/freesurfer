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
#include "TrackData.h"
#include "track_io/TrackIO.h"
#include <iostream>
#include <stdlib.h>
#include <QDebug>

TrackData::TrackData(QObject *parent) :
  QObject(parent), m_bHasEmbeddedColor(false)
{
}

TrackData::~TrackData()
{
  Clear();
}

void TrackData::Clear()
{
  for (int i = 0; i < m_tracks.size(); i++)
  {
    m_tracks[i].Delete();
  }
  m_tracks.clear();
}

bool TrackData::LoadFromFiles(const QStringList &filenames)
{
  Clear();
  CTrackReader reader;
  TRACK_HEADER header;
  Track track;
  for (int n = 0; n < filenames.size(); n++)
  {
    if (!reader.Open(filenames[n].toLatin1().constData(), &header))
    {
      return false;
    }

    if (n == 0)
    {
      for (int i = 0; i < 3; i++)
      {
        m_nDim[i] = header.dim[i];
        m_dVoxelSize[i] = header.voxel_size[i];
      }
      m_nNumberOfScalars = header.n_scalars;
      m_nNumberOfProperties = header.n_properties;
      for (int i = 0; i < m_nNumberOfScalars; i++)
      {
        m_scalarNames << header.scalar_name[i];
      }
      for (int i = 0; i < m_nNumberOfProperties; i++)
      {
        m_propertyNames << header.property_name[i];
      }

      if (header.vox_to_ras[3][3] == 0)
      {
        m_bValidVoxToRas = false;
        for (int i = 0; i < 4; i++)
        {
          memset(m_dVoxToRas[i], 0, sizeof(double)*4);
        }
        for (int i = 0; i < 3; i++)
        {
          m_dVoxToRas[i][i] = 1.0/m_dVoxelSize[i];
        }
        m_dVoxToRas[3][3] = 1;
      }
      else
      {
        m_bValidVoxToRas = true;
        for (int i = 0; i < 4; i++)
        {
          for (int j = 0; j < 4; j++)
          {
            m_dVoxToRas[i][j] = header.vox_to_ras[i][j];
          }
        }
      }

      for (int i = 0; i < m_nNumberOfScalars; i++)
      {
        m_rangeScalar << qMakePair(1e12, -1e12);
      }
      for (int i = 0; i < m_nNumberOfProperties; i++)
      {
        m_rangeProperty << qMakePair(1e12, -1e12);
      }
      m_nNumberOfPoints = 0;
      m_nNumberOfSegs = 0;
      m_bHasEmbeddedColor = (header.reserved[0] == 'C');
    }

    int nScalars = m_nNumberOfScalars;
    int nProperties = m_nNumberOfProperties;
    short dim[3] = {(short)m_nDim[0], (short)m_nDim[1], (short)m_nDim[2]};
    while (reader.GetNextPointCount(&track.nNum))
    {
      track.fPts = new float[track.nNum*3];
      if (!track.fPts)
      {
        std::cerr << "Can not allocate memory." << std::endl;
        return false;
      }
      float* f = NULL;
      if (nScalars > 0)
      {
        f = new float[track.nNum*nScalars];
        if (!f)
        {
          std::cerr << "Can not allocate memory." << std::endl;
          return false;
        }
      }
      track.fProperty = NULL;
      if (nProperties > 0)
      {
        track.fProperty = new float[nProperties];
        if (!track.fProperty)
        {
          std::cerr << "Can not allocate memory." << std::endl;
          return false;
        }
      }
      reader.GetNextTrackData(track.nNum, track.fPts, f, track.fProperty);
      for (int i = 0; i < nScalars; i++)
      {
        float* p = new float[track.nNum];
        if (!p)
        {
          std::cerr << "Can not allocate memory." << std::endl;
          return false;
        }
        track.fScalars.push_back(p);
        for (int j = 0; j < track.nNum; j++)
        {
          p[j] = f[j*nScalars+i];
          if (m_rangeScalar[i].first > p[j])
          {
            m_rangeScalar[i].first = p[j];
          }
          else if (m_rangeScalar[i].second < p[j])
          {
            m_rangeScalar[i].second = p[j];
          }
        }
      }
      delete[] f;
      for (int i = 0; i < nProperties; i++)
      {
        if (m_rangeProperty[i].first > track.fProperty[i])
        {
          m_rangeProperty[i].first = track.fProperty[i];
        }
        else if (m_rangeProperty[i].second < track.fProperty[i])
        {
          m_rangeProperty[i].second = track.fProperty[i];
        }
      }

      if (m_bHasEmbeddedColor)
        memcpy(track.charColor, header.reserved+1, 3);

      track.Update(dim, m_dVoxelSize);
      m_tracks.push_back(track);
      m_nNumberOfPoints += track.nNum;
      m_nNumberOfSegs += (track.nNum-1);

      if ( m_tracks.size()%100 == 1 )
      {
        emit Progress(reader.GetProgress());
      }
      // must call reset if track instance is to be repeatedly used
      track.Reset();
    }
  }
  m_nNumberOfTracks = m_tracks.size();

  return true;
}
