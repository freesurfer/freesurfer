#include "LayerPropertyConnectomeMatrix.h"

LayerPropertyConnectomeMatrix::LayerPropertyConnectomeMatrix(QObject *parent) :
  LayerProperty(parent),
  m_dSplineRadius(0.5),
  m_nNumberOfSides(5),
  m_dFromLabelOpacity(1.0),
  m_dToLabelOpacity(0.0),
  m_colorSpline(Qt::yellow)
{
  connect(this, SIGNAL(OpacityChanged()), this, SIGNAL(PropertyChanged()));
}

void LayerPropertyConnectomeMatrix::SetFromLabelOpacity(double dval)
{
  if (m_dFromLabelOpacity != dval)
  {
    m_dFromLabelOpacity = dval;
    emit OpacityChanged();
  }
}

void LayerPropertyConnectomeMatrix::SetToLabelOpacity(double dval)
{
  if (m_dToLabelOpacity != dval)
  {
    m_dToLabelOpacity = dval;
    emit OpacityChanged();
  }
}


void LayerPropertyConnectomeMatrix::SetSplineRadius(double val)
{
  if (m_dSplineRadius != val)
  {
    m_dSplineRadius = val;
    emit SplineRadiusChanged();
  }
}

void LayerPropertyConnectomeMatrix::SetSplineColor(const QColor &c)
{
  m_colorSpline = c;
  emit SplineColorChanged();
}
