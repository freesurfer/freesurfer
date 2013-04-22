#ifndef LAYERPROPERTYCONNECTOMEMATRIX_H
#define LAYERPROPERTYCONNECTOMEMATRIX_H

#include "LayerProperty.h"
#include <QColor>

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

  QColor GetSplineColor()
  {
    return m_colorSpline;
  }

signals:
  void OpacityChanged();
  void SplineRadiusChanged();
  void SplineColorChanged();

public slots:
  void SetFromLabelOpacity(double dval);
  void SetToLabelOpacity(double dval);
  void SetSplineRadius(double val);
  void SetSplineColor(const QColor& c);

private:
  double  m_dSplineRadius;
  int     m_nNumberOfSides;
  double  m_dFromLabelOpacity;
  double  m_dToLabelOpacity;
  QColor  m_colorSpline;
};

#endif // LAYERPROPERTYCONNECTOMEMATRIX_H
