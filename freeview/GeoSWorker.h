#ifndef GEOSWORKER_H
#define GEOSWORKER_H

#include <QObject>
#include <QThread>

class LayerMRI;
class GeodesicMatting;

class GeoSWorker : public QObject
{
  Q_OBJECT
public:
  explicit GeoSWorker(QObject *parent = nullptr);
  ~GeoSWorker();

  QString GetErrorMessage();

signals:
  void ComputeTriggered();
  void ApplyTriggered();
  void ApplyFinished();
  void ComputeFinished(double time_in_secs);  // -1 means fail
  void Progress(double val);

public slots:
  void Compute(LayerMRI* mri, LayerMRI* seg, LayerMRI* seeds, int max_distance = -1, double smoothing = 0, LayerMRI* mask = NULL, double fill_val = -1, int max_foreground_dist = 0);
  void Apply(LayerMRI* seg, LayerMRI* filled);
  void Abort();

private slots:
  void DoCompute();
  void DoApply();

private:
  LayerMRI* m_seeds;
  LayerMRI* m_mri;
  LayerMRI* m_seg;
  LayerMRI* m_filled;
  LayerMRI* m_mask;
  int     m_nMaxDistance;
  int     m_nMaxForegroundDistance;
  double    m_dSmoothing;
  double    m_dFillValue;
  GeodesicMatting*  m_geos;

  QThread   m_thread;
};

#endif // GEOSWORKER_H
