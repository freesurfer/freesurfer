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

  double GetMinArea()
  {
    return m_dMinArea;
  }

  void SetLookupTable(vtkRGBAColorTransferFunction* lut);

signals:
  void ThicknessThresholdChanged(double);
  void SigmaChanged(double);
  void OpacityChanged( double opacity );
  void MinAreaChanged(double);
  void ColorMapChanged();

public slots:
  void SetThicknessThreshold(double dThreshold);
  void SetSigma(double dSigma);
  void SetOpacity(double val);
  void SetMinArea(double val);

private:
  void SetColorMapChanged();

  double m_dThreshold;
  double m_dSigma;
  double m_dMinArea;
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  double mOpacity;
};

#endif // LAYERPROPERTYFCD_H
