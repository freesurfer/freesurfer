#ifndef SCRIBBLEPROMPTWORKER_H
#define SCRIBBLEPROMPTWORKER_H

#include <QObject>

class TorchScriptModule;
class LayerMRI;

class ScribblePromptWorker : public QObject
{
  Q_OBJECT
public:
  explicit ScribblePromptWorker(QObject *parent = nullptr);
  ~ScribblePromptWorker();

  void Initialize(const QString& fn)
  {
    emit InitializationTriggered(fn);
  }

signals:
  void InitializationTriggered(const QString& fn);
  void ComputeTriggered();
  void ApplyTriggered();
  void ComputeFinished(double elapsed_time);
  void ApplyFinished();

public slots:
  void Compute(LayerMRI *mri, LayerMRI* seg, LayerMRI* seeds, int nPlane, int nSlice, double fill_val);
  void Apply(LayerMRI *seg, LayerMRI *filled);

private slots:
  void DoInitialization(const QString& fn);
  void DoCompute();
  void DoApply();

private:
  LayerMRI* m_seeds;
  LayerMRI* m_mri;
  LayerMRI* m_seg;
  LayerMRI* m_filled;
  int     m_nInputPlane;
  int     m_nInputSlice;
  double    m_dFillValue;
  TorchScriptModule* m_module;
};

#endif // SCRIBBLEPROMPTWORKER_H
