#include "VolumeFilterWorkerThread.h"
#include "VolumeFilter.h"

VolumeFilterWorkerThread::VolumeFilterWorkerThread(QObject *parent) :
    QThread(parent), m_filter(NULL)
{
  connect(this, SIGNAL(finished()), this, SLOT(OnFinished()));
}

void VolumeFilterWorkerThread::ExecuteFilter(VolumeFilter *filter)
{
  m_filter = filter;
  connect(filter, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)), Qt::UniqueConnection);
  start(QThread::HighPriority);
}

void VolumeFilterWorkerThread::run()
{
  if (m_filter)
    m_filter->Update();
}

void VolumeFilterWorkerThread::OnFinished()
{
  if (m_filter)
  {
    emit Progress(100);
    emit Finished(m_filter);
    disconnect(m_filter, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)));
  }
}
