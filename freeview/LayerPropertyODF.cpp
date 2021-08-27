#include "LayerPropertyODF.h"
#include "vtkColorTransferFunction.h"
#include <QDebug>

LayerPropertyODF::LayerPropertyODF(QObject *parent) :
  LayerPropertyMRI(parent)
{
  m_dOdfScale = 1;
  m_nOdfSkip = 0;
  m_nOdfInversion = 0;
  m_nOdfColorCode = Directional;
  m_odfLut = vtkSmartPointer<vtkColorTransferFunction>::New();
  m_odfLut->ClampingOn();

  m_spectrum[0] = Qt::blue;
  m_spectrum[0.33] = Qt::cyan;
  m_spectrum[0.5] = Qt::green;
  m_spectrum[0.66] = Qt::yellow;
  m_spectrum[1] = Qt::red;
}

void LayerPropertyODF::SetOdfInversion(int n)
{
  m_nOdfInversion = n;
  emit OdfPropertyChanged();
}

void LayerPropertyODF::SetOdfSkip(int n)
{
  if (n >= 0)
  {
    m_nOdfSkip = n;
    emit OdfPropertyChanged();
  }
}

void LayerPropertyODF::SetOdfScale(double val)
{
  m_dOdfScale = val;
  emit OdfPropertyChanged();
}

void LayerPropertyODF::SetOdfColorCode(int n)
{
  m_nOdfColorCode = n;
  emit ColorCodeChanged();
}

void LayerPropertyODF::GetMagnitudeThreshold(double *th)
{
  th[0] = m_dMagnitudeThreshold[0];
  th[1] = m_dMagnitudeThreshold[1];
}

void LayerPropertyODF::SetMagnitudeThreshold(double *th)
{
  m_dMagnitudeThreshold[0] = th[0];
  m_dMagnitudeThreshold[1] = th[1];
  UpdateOdfLut();
  emit ColorCodeChanged();
}

vtkColorTransferFunction* LayerPropertyODF::GetOdfLut()
{
  return m_odfLut;
}

void LayerPropertyODF::UpdateOdfLut()
{
  double* th = m_dMagnitudeThreshold;
  if (th[0] >= th[1])
    return;

  m_odfLut->RemoveAllPoints();
  QList<double> keys = m_spectrum.keys();
  foreach (double val, keys)
  {
    QColor c = m_spectrum[val];
    m_odfLut->AddRGBPoint(th[0]+(th[1]-th[0])*val,
        c.redF(), c.greenF(), c.blueF());
  }
  m_odfLut->Build();
}
