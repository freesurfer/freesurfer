#ifndef LAYERTRACK_H
#define LAYERTRACK_H

#include "Layer.h"

class FSTrack;
class LayerMRI;
class TrackGroup;
class LayerPropertyTrack;

class LayerTrack : public Layer
{
    Q_OBJECT
public:
    LayerTrack(LayerMRI* ref, QObject* parent = NULL);
    ~LayerTrack();

    bool LoadTrackFromFile();

    void Append2DProps(vtkRenderer *renderer, int nPlane);

    void Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility);

    bool HasProp(vtkProp *prop);

    bool IsVisible();

    inline LayerPropertyTrack* GetProperty()
    {
      return (LayerPropertyTrack*)mProperty;
    }
signals:
    void Progress(int n);

protected:
    virtual void OnSlicePositionChanged(int nPlane);

    FSTrack*    m_trackData;
    LayerMRI*   m_layerMRIRef;
    QList<TrackGroup*>  m_trackGroups;
};

#endif // LAYERTRACK_H
