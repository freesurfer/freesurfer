#ifndef LAYERPROPERTYFCD_H
#define LAYERPROPERTYFCD_H

#include "LayerProperty.h"
#include <vtkSmartPointer.h>

class vtkRGBAColorTransferFunction;

class LayerPropertyFCD : public LayerProperty
{
  Q_OBJECT
public:
  explicit LayerPropertyFCD(QObject *parent = 0);

  vtkRGBAColorTransferFunction* GetLookupTable() const;

  double GetOpacity()
  {
    return mOpacity;
  }

  double GetThicknessThreshold()
  {
    return m_dThreshold;
  }

  double GetSigma()
  {
    return m_dSigma;
  }

  void SetLookupTable(vtkRGBAColorTransferFunction* lut);

signals:
  void ThicknessThresholdChanged(double);
  void SigmaChanged(double);
  void OpacityChanged( double opacity );
  void ColorMapChanged();

public slots:
  void SetThicknessThreshold(double dThreshold);
  void SetSigma(double dSigma);
  void SetOpacity(double val);

private:
  void SetColorMapChanged();

  double m_dThreshold;
  double m_dSigma;
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  double mOpacity;
};

#endif // LAYERPROPERTYFCD_H
