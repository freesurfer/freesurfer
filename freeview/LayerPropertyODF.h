#ifndef LAYERPROPERTYODF_H
#define LAYERPROPERTYODF_H

#include "LayerPropertyMRI.h"

class LayerPropertyODF : public LayerPropertyMRI
{
  Q_OBJECT
public:
  explicit LayerPropertyODF(QObject *parent = 0);

protected:
  int m_nSkip;
};

#endif // LAYERPROPERTYODF_H
