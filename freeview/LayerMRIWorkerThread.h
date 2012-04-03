#ifndef LAYERMRIWORKERTHREAD_H
#define LAYERMRIWORKERTHREAD_H

#include <QThread>
#include <QList>
#include <QMutex>

#ifndef IntList
typedef QList<int> IntList;
#endif

class LayerMRI;

class LayerMRIWorkerThread : public QThread
{
    Q_OBJECT
public:
    explicit LayerMRIWorkerThread(LayerMRI *mri);

signals:
    void AvailableLabels(const IntList& list);

public slots:
    void Abort();

protected:
    void run();

    bool m_bAbort;
    QMutex mutex;
};

#endif // LAYERMRIWORKERTHREAD_H
