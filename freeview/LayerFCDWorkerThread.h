#ifndef LayerFCDWorkerThread_H
#define LayerFCDWorkerThread_H

#include <QThread>
#include <QList>
#include <QMutex>

#ifndef IntList
typedef QList<int> IntList;
#endif

class LayerFCD;

class LayerFCDWorkerThread : public QThread
{
  Q_OBJECT
public:
  explicit LayerFCDWorkerThread(LayerFCD *fcd);

signals:
  void Progress(int);

public slots:
  void Abort();

protected:
  void run();

  bool m_bAbort;
  QMutex mutex;
};

#endif // LayerFCDWorkerThread_H
