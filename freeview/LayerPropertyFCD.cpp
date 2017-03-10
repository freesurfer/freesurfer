#include "LayerPropertyFCD.h"
#include "vtkRGBAColorTransferFunction.h"

LayerPropertyFCD::LayerPropertyFCD(QObject *parent) :
  LayerProperty(parent),
  m_dThreshold(2),
  m_dSigma(10),
  m_nMinArea(10),
  mOpacity(0.8)
{
  mLUTTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  // mLUTTable->ClampingOff();
  this->SetColorMapChanged();

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
}

vtkRGBAColorTransferFunction* LayerPropertyFCD::GetLookupTable () const
{
  return mLUTTable;
}

void LayerPropertyFCD::SetLookupTable(vtkRGBAColorTransferFunction *lut)
{
  mLUTTable = lut;
}

void LayerPropertyFCD::SetColorMapChanged()
{
  // Notify the layers that use the color map stuff.
  emit ColorMapChanged();
}

void LayerPropertyFCD::SetOpacity( double opacity )
{
  if ( mOpacity != opacity )
  {
    mOpacity = opacity;
    emit OpacityChanged( opacity );
  }
}

void LayerPropertyFCD::SetThicknessThreshold(double dThreshold)
{
  if (m_dThreshold != dThreshold)
  {
    m_dThreshold = dThreshold;
    emit ThicknessThresholdChanged(dThreshold);
  }
}

void LayerPropertyFCD::SetSigma(double dSigma)
{
  if (m_dSigma != dSigma)
  {
    m_dSigma = dSigma;
    emit SigmaChanged(dSigma);
  }
}

void LayerPropertyFCD::SetMinArea(int val)
{
  if (m_nMinArea != val)
  {
    m_nMinArea = val;
    emit MinAreaChanged(val);
  }
}
