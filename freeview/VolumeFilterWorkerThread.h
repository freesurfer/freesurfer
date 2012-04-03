#ifndef VOLUMEFILTERWORKERTHREAD_H
#define VOLUMEFILTERWORKERTHREAD_H

#include <QThread>

class VolumeFilter;

class VolumeFilterWorkerThread : public QThread
{
    Q_OBJECT
public:
    explicit VolumeFilterWorkerThread(QObject *parent = 0);

signals:
    void Progress(int n);
    void Finished(VolumeFilter* filter);

public slots:
    void ExecuteFilter(VolumeFilter* filter);

protected slots:
    void OnFinished();

protected:
    virtual void run();

private:
    VolumeFilter* m_filter;
};

#endif // VOLUMEFILTERWORKERTHREAD_H
