#ifndef LAYERPROPERTYCONNECTOMEMATRIX_H
#define LAYERPROPERTYCONNECTOMEMATRIX_H

#include "LayerProperty.h"

class LayerPropertyConnectomeMatrix : public LayerProperty
{
  Q_OBJECT
public:
  explicit LayerPropertyConnectomeMatrix(QObject *parent = 0);

  double GetSplineRadius()
  {
    return m_dSplineRadius;
  }

  int GetNumberOfSides()
  {
    return m_nNumberOfSides;
  }

signals:

public slots:


private:
  double  m_dSplineRadius;
  int     m_nNumberOfSides;
};

#endif // LAYERPROPERTYCONNECTOMEMATRIX_H
