/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "LayerPropertyTrack.h"

LayerPropertyTrack::LayerPropertyTrack(QObject* parent) :
  LayerProperty(parent),
  m_nColorCode(Directional),
  m_nDirectionScheme(EndPoints),
  m_nDirectionMapping(RGB),
  m_color(Qt::yellow),
  m_nRenderRep(Line),
  m_dTubeRadius(0.2),
  m_nNumberOfSides(5),
  m_dOpacity(1)
{
  connect(this, SIGNAL(ColorCodeChanged(int)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(DirectionSchemeChanged(int)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(DirectionMappingChanged(int)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(SolidColorChanged(QColor)), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(RenderRepChanged()), this, SIGNAL(PropertyChanged()));
  connect(this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()));
}

void LayerPropertyTrack::SetColorCode(int nCode)
{
  if (m_nColorCode != nCode)
  {
    m_nColorCode = nCode;
    emit ColorCodeChanged(nCode);
  }
}

void LayerPropertyTrack::SetDirectionScheme(int nVal)
{
  if (m_nDirectionScheme != nVal)
  {
    m_nDirectionScheme = nVal;
    emit DirectionSchemeChanged(nVal);
  }
}

void LayerPropertyTrack::SetDirectionMapping(int nVal)
{
  if (m_nDirectionMapping != nVal)
  {
    m_nDirectionMapping = nVal;
    emit DirectionMappingChanged(nVal);
  }
}

void LayerPropertyTrack::SetSolidColor(const QColor &c)
{
  if (m_color != c)
  {
    m_color = c;
    emit SolidColorChanged(c);
  }
}

void LayerPropertyTrack::SetRenderRep(int nVal)
{
  m_nRenderRep = nVal;
  emit RenderRepChanged();
}

void LayerPropertyTrack::SetTubeRadius(double dVal)
{
  if (dVal != m_dTubeRadius)
  {
    m_dTubeRadius = dVal;
    SetRenderRep(m_nRenderRep);
  }
}

void LayerPropertyTrack::SetNumberOfSides(int nVal)
{
  if (nVal != m_nNumberOfSides)
  {
    m_nNumberOfSides = nVal;
    SetRenderRep(m_nRenderRep);
  }
}

void LayerPropertyTrack::SetOpacity(double val)
{
  if (val != m_dOpacity)
  {
    m_dOpacity = val;
    emit OpacityChanged(val);
  }
}
