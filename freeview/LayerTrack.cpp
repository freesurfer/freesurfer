#include "LayerTrack.h"
#include "FSTrack.h"
#include "LayerMRI.h"
#include "FSVolume.h"
#include "TrackGroup.h"
#include "LayerPropertyTrack.h"
#include <QFileInfo>
#include <QDebug>

LayerTrack::LayerTrack(LayerMRI* ref, QObject* parent) : Layer(parent),
    m_trackData(0),
    m_layerMRIRef(ref)
{
    this->m_strTypeNames << "Track";

    mProperty = new LayerPropertyTrack( this );
}

LayerTrack::~LayerTrack()
{
    if (m_trackData)
        delete m_trackData;
}

bool LayerTrack::LoadTrackFromFile()
{
    if (this->m_sFilename.isEmpty())
        return false;
    FSVolume* refVol = NULL;
    if (m_layerMRIRef)
        refVol = m_layerMRIRef->GetSourceVolume();
    m_trackData = new FSTrack(refVol);
    connect(m_trackData, SIGNAL(Progress(int)), this, SIGNAL(Progress(int)));
    if (!m_trackData->LoadFromFile(m_sFilename))
    {
        delete m_trackData;
        m_trackData = 0;
        qDebug() << "Failed to load from file " << m_sFilename << ".";
        return false;
    }
    SetName(QFileInfo(m_sFilename).completeBaseName());

    return true;
}

void LayerTrack::Append2DProps(vtkRenderer *renderer, int nPlane)
{

}

void LayerTrack::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{

}

bool LayerTrack::HasProp(vtkProp *prop)
{
    return true;
}

bool LayerTrack::IsVisible()
{
    return true;
}

void LayerTrack::OnSlicePositionChanged(int nPlane)
{

}
