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
#ifndef TRACKGROUP_H
#define TRACKGROUP_H

#include <QObject>

class TrackGroup : public QObject
{
  Q_OBJECT
public:
  explicit TrackGroup(QObject *parent = 0);

signals:

public slots:

};

#endif // TRACKGROUP_H
