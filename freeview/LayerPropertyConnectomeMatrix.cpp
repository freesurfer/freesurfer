#include "LayerPropertyConnectomeMatrix.h"

LayerPropertyConnectomeMatrix::LayerPropertyConnectomeMatrix(QObject *parent) :
  LayerProperty(parent),
  m_dSplineRadius(0.5),
  m_nNumberOfSides(5),
  m_dFromLabelOpacity(1.0),
  m_dToLabelOpacity(1.0)
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
