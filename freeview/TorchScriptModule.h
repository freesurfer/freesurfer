#ifndef TORCHSCRIPTMODULE_H
#define TORCHSCRIPTMODULE_H

#include <QObject>

class TorchScriptModule : public QObject
{
  Q_OBJECT
public:
  explicit TorchScriptModule(QObject *parent = nullptr);
  ~TorchScriptModule();
  void Load(const QString& fn);
  void Run(QVector<float*> in_ptr, float* output);

signals:
  void Finished();

private:
  void* m_module;
};

#endif // TORCHSCRIPTMODULE_H
