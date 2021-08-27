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
#ifndef ThreadIOWorker_H
#define ThreadIOWorker_H

#include <QThread>
#include <QMutex>
#include <QVariantMap>

class Layer;
class LayerFCD;

class ThreadIOWorker : public QThread
{
  Q_OBJECT
public:
  explicit ThreadIOWorker(QObject *parent = 0);

  enum JobType { JT_LoadVolume = 0, JT_SaveVolume, JT_LoadSurface, JT_SaveSurface, JT_LoadSurfaceOverlay, JT_LoadTrack,
                 JT_LoadConnectome, JT_LoadFCD, JT_LoadODF, JT_TransformVolume };

  void LoadVolume( Layer* layer, const QVariantMap& args = QVariantMap() );
  void SaveVolume( Layer* layer, const QVariantMap& args = QVariantMap() );
  void LoadSurface( Layer* layer, const QVariantMap& args = QVariantMap() );
  void SaveSurface( Layer* layer, const QVariantMap& args = QVariantMap() );
  void LoadSurfaceOverlay( Layer* layer, const QVariantMap& args = QVariantMap() );
  void LoadTrack( Layer* layer, const QVariantMap& args = QVariantMap() );
  void LoadConnectomeMatrix( Layer* layer, const QVariantMap& args = QVariantMap() );
  void LoadFCD(Layer* layer, const QVariantMap& args = QVariantMap());
  void LoadODF(Layer* layer, const QVariantMap& args = QVariantMap());
  void TransformVolume( Layer* layer, const QVariantMap& args = QVariantMap() );

signals:
  void Progress( int n );
  void Error( Layer*, int jobtype );
  void Finished( Layer*, int jobtype );
  void FCDLoadFinished(LayerFCD*);

public slots:

protected:
  void run();

private:
  QMutex    mutex;
  // to abort
  int         m_nJobType;
  Layer*      m_layer;
  QVariantMap m_args;
};

#endif // ThreadIOWorker_H

