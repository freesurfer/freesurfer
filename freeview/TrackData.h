#ifndef TRACKDATA_H
#define TRACKDATA_H

#include <QObject>
#include <QStringList>
#include "Track.h"
#include <QPair>

class TrackData : public QObject
{
    Q_OBJECT
public:
    TrackData(QObject *parent = 0);
    ~TrackData();

    bool LoadFromFile(const QString& filename);

    int GetNumberOfTracks()
    {
        return m_nNumberOfTracks;
    }

signals:
    void Progress(int n);

public slots:

protected:
    int     m_nDim[3];
    float   m_dVoxelSize[3];
    int     m_nNumberOfScalars;
    QStringList m_scalarNames;
    int     m_nNumberOfProperties;
    QStringList m_propertyNames;
    double  m_dVoxToRas[4][4];

    int     m_nNumberOfTracks;
    int     m_nNumberOfPoints;
    int     m_nNumberOfSegs;

    bool    m_bValidVoxToRas;
    QString m_sFileName;

    QList<Track>    m_tracks;
    QList< QPair<double, double> > m_rangeScalar;
    QList< QPair<double, double> > m_rangeProperty;
};

#endif // TRACKDATA_H
