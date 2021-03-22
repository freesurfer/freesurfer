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
#ifndef TRACKDATA_H
#define TRACKDATA_H

#include <QObject>
#include <QStringList>
#include "Track.h"
#include <QPair>
#include <QColor>

class TrackData : public QObject
{
  friend class LayerTrack;

  Q_OBJECT
public:
  TrackData(QObject *parent = 0);
  ~TrackData();

  bool LoadFromFiles(const QStringList& filenames);

  int GetNumberOfTracks()
  {
    return m_nNumberOfTracks;
  }

  bool HasEmbeddedColor()
  {
    return m_bHasEmbeddedColor;
  }

signals:
  void Progress(int n);

public slots:
  void Clear();

protected:
  int     m_nDim[3];
  float   m_dVoxelSize[3];
  int     m_nNumberOfScalars;
  QStringList m_scalarNames;
  int     m_nNumberOfProperties;
  QStringList m_propertyNames;
  double  m_dVoxToRas[4][4];

  int     m_nNumberOfTracks;
  int     m_nNumberOfPoints;
  int     m_nNumberOfSegs;

  bool    m_bValidVoxToRas;
  bool    m_bHasEmbeddedColor;

  QList<Track>    m_tracks;
  QList< QPair<double, double> > m_rangeScalar;
  QList< QPair<double, double> > m_rangeProperty;
};

#endif // TRACKDATA_H
