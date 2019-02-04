#ifndef GEOSWORKER_H
#define GEOSWORKER_H

#include <QObject>
#include <QThread>

class LayerMRI;

class GeoSWorker : public QObject
{
  Q_OBJECT
public:
  explicit GeoSWorker(QObject *parent = nullptr);
  ~GeoSWorker();

signals:
  void ComputeTriggered();
  void Finished(bool bSuccess);

public slots:
  void Compute(LayerMRI* mri, LayerMRI* seg, LayerMRI* seeds, int max_distance = -1);

private slots:
  void DoCompute();

private:
  LayerMRI* m_seeds;
  LayerMRI* m_mri;
  LayerMRI* m_seg;
  int     m_nMaxDistance;

  QThread   m_thread;
};

#endif // GEOSWORKER_H
