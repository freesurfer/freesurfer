/**
 * @file  TrackData.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/12/05 20:03:33 $
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
#ifndef TRACKDATA_H
#define TRACKDATA_H

#include <QObject>
#include <QStringList>
#include "Track.h"
#include <QPair>

class TrackData : public QObject
{
  friend class LayerTrack;

  Q_OBJECT
public:
  TrackData(QObject *parent = 0);
  ~TrackData();

  bool LoadFromFile(const QString& filename);

  int GetNumberOfTracks()
  {
    return m_nNumberOfTracks;
  }

signals:
  void Progress(int n);

public slots:

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
  QString m_sFileName;

  QList<Track>    m_tracks;
  QList< QPair<double, double> > m_rangeScalar;
  QList< QPair<double, double> > m_rangeProperty;
};

#endif // TRACKDATA_H
