/**
 * @file  FSTrack.h
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
#ifndef FSTRACK_H
#define FSTRACK_H

#include "TrackData.h"

class FSVolume;

class FSTrack : public TrackData
{
  Q_OBJECT
public:
  FSTrack(FSVolume* ref = 0, QObject *parent = 0);
  bool LoadFromFile(const QString &filename, const QString& ref_fn = QString());

signals:

public slots:

protected:

  FSVolume* m_volumeRef;
};

#endif // FSTRACK_H
