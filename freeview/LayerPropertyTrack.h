#ifndef LAYERPROPERTYTRACK_H
#define LAYERPROPERTYTRACK_H

#include "LayerProperty.h"

class LayerPropertyTrack : public LayerProperty
{
    Q_OBJECT
public:
    LayerPropertyTrack(QObject* parent = 0);
};

#endif // LAYERPROPERTYTRACK_H
