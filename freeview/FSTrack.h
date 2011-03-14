#ifndef FSTRACK_H
#define FSTRACK_H

#include "TrackData.h"

class FSVolume;

class FSTrack : public TrackData
{
    Q_OBJECT
public:
    FSTrack(FSVolume* ref = 0, QObject *parent = 0);
    bool LoadFromFile(const QString &filename);

signals:

public slots:

protected:

    FSVolume* m_volumeRef;
};

#endif // FSTRACK_H
