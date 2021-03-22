/**
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
 * Reimplemented by: Ruopeng Wang
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
 *
 */

#ifndef LayerPropertyROI_h
#define LayerPropertyROI_h

#include "vtkSmartPointer.h"
#include "LayerProperty.h"
#include <QColor>



#include "colortab.h"


class vtkRGBAColorTransferFunction;

class LayerPropertyROI : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyROI ( QObject* parent = NULL );
  ~LayerPropertyROI ();

  double GetOpacity() const;

  double* GetColor()
  {
    return mRGB;
  }
  void SetColor( double r, double g, double b );

  double GetThreshold() const
  {
    return m_dThreshold;
  }

  int GetColorCode()
  {
    return m_nColorCode;
  }

  double GetHeatscaleMin()
  {
    return m_dHeatscaleMin;
  }

  double GetHeatscaleMax()
  {
    return m_dHeatscaleMax;
  }

  void GetValueRange(double* range)
  {
    range[0] = m_dValueRange[0];
    range[1] = m_dValueRange[1];
  }

  vtkRGBAColorTransferFunction* GetLookupTable() const;

  enum ColorCode { SolidColor = 0, Heatscale };

public slots:
  void SetOpacity( double opacity );
  void SetColor( const QColor& c )
  {
    SetColor( c.redF(), c.greenF(), c.blueF() );
  }
  void SetThreshold(double th);

  void SetHeatscaleMin(double val);
  void SetHeatscaleMax(double val);
  void SetHeatscaleValues(double dMin, double dMax);

  void SetValueRange(double range[2])
  {
    m_dValueRange[0] = range[0];
    m_dValueRange[1] = range[1];
    SetColorMapChanged();
  }

  void SetColorCode(int nCode);

  void UpdateLUTTable();

signals:
  void OpacityChanged( double opacity );
  void ColorMapChanged();
  void ThresholdChanged(double th);

private:
  void SetColorMapChanged();

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  // ---------------------------------------------------------------------

  double mOpacity;
  double mRGB[3];
  double m_dThreshold;
  int    m_nColorCode;

  double m_dHeatscaleMin;
  double m_dHeatscaleMax;
  double m_dValueRange[2];
};

#endif
