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

  double GetFromLabelOpacity()
  {
    return m_dFromLabelOpacity;
  }

  double GetToLabelOpacity()
  {
    return m_dToLabelOpacity;
  }

signals:
  void OpacityChanged();

public slots:
  void SetFromLabelOpacity(double dval);
  void SetToLabelOpacity(double dval);

private:
  double  m_dSplineRadius;
  int     m_nNumberOfSides;
  double  m_dFromLabelOpacity;
  double  m_dToLabelOpacity;
};

#endif // LAYERPROPERTYCONNECTOMEMATRIX_H
