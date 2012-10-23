#include "LayerPropertyLineProfile.h"

LayerPropertyLineProfile::LayerPropertyLineProfile(QObject *parent) :
    LayerProperty(parent),
    m_color(QColor(255, 50, 50)),
    m_dOpacity(0.7),
    m_dRadius(0.25)
{
}

void LayerPropertyLineProfile::SetColor(const QColor &color)
{
  if (m_color != color)
  {
    m_color = color;
    emit ColorChanged(color);
  }
}

void LayerPropertyLineProfile::SetOpacity(double val)
{
  if (m_dOpacity != val)
  {
    m_dOpacity = val;
    emit OpacityChanged(val);
  }
}

void LayerPropertyLineProfile::SetRadius(double val)
{
  if (m_dRadius != val)
  {
    m_dRadius = val;
    emit RadiusChanged(val);
  }
}
