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
#ifndef ThreadBuildContour_H
#define ThreadBuildContour_H

#include <QThread>
#include <QMutex>
#include <QVariantMap>

class LayerMRI;

class ThreadBuildContour : public QThread
{
  Q_OBJECT
public:
  explicit ThreadBuildContour(QObject *parent = 0);

  void BuildContour( LayerMRI* mri, int nSegValue, int nThreadID );

signals:
  void Finished( int nThreadID );

public slots:

protected:
  void run();

  LayerMRI*   m_mri;
  int         m_nSegValue;
  int         m_nThreadID;
};

#endif // ThreadBuildContour_H
