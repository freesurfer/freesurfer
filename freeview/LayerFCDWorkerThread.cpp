#include "LayerFCDWorkerThread.h"
#include <cstddef>
#include "LayerFCD.h"
#include <QMutexLocker>
#include <QApplication>

LayerFCDWorkerThread::LayerFCDWorkerThread(LayerFCD *fcd) :
  QThread(fcd), m_bAbort(false)
{
}

void LayerFCDWorkerThread::Abort()
{
  QMutexLocker locker(&mutex);
  m_bAbort = true;
}

void LayerFCDWorkerThread::run()
{
  connect(qApp, SIGNAL(GlobalProgress(int)), this, SIGNAL(Progress(int)), Qt::UniqueConnection);
  LayerFCD* fcd = qobject_cast<LayerFCD*>(parent());
  fcd->DoCompute();
  disconnect(qApp, SIGNAL(GlobalProgress(int)), this, SIGNAL(Progress(int)));
}
