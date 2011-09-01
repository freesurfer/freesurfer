#ifndef LAYERMRIWORKERTHREAD_H
#define LAYERMRIWORKERTHREAD_H

#include <QThread>
#include <QList>

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

protected:
    void run();
};

#endif // LAYERMRIWORKERTHREAD_H
