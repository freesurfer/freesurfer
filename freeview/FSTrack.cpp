#include "FSTrack.h"
#include "FSVolume.h"
#include "Track.h"

FSTrack::FSTrack(FSVolume* ref, QObject *parent) :
    TrackData(parent),
    m_volumeRef(ref)
{
}

bool FSTrack::LoadFromFile(const QString &filename)
{
    if (!TrackData::LoadFromFile(filename))
        return false;

    if (m_volumeRef)
    {
        for (int i = 0; i < m_tracks.size(); i++)
        {
            for (int j = 0; j < m_tracks[i].nNum; j++)
                m_volumeRef->RASToTarget(m_tracks[i].fPts + j*3, m_tracks[i].fPts + j*3);
        }
    }
    return true;
}
