/**
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
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
 *
 */

#ifndef LayerPropertyPointSet_h
#define LayerPropertyPointSet_h

#include "vtkSmartPointer.h"
#include "LayerProperty.h"
#include <QColor>
#include <vector>
#include <QVariantMap>



#include "colortab.h"



struct ScalarValues
{
  QString strName;
  int  nNum;
  double* dValue;
  double dMin;
  double dMax;
};

class vtkRGBAColorTransferFunction;
class LayerMRI;

class LayerPropertyPointSet : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyPointSet ( QObject* parent = NULL );
  ~LayerPropertyPointSet ();

  enum ColorMapType { SolidColor = 0, HeatScale };

  enum ScalarType { ScalarStat = 0, ScalarLayer, ScalarSet };

  enum PointSetType { WayPoint = 0, ControlPoint, Enhanced };

  double GetOpacity() const;

  double* GetColor()
  {
    return mRGB;
  }

  void SetColor( double r, double g, double b );

  double* GetSplineColor()
  {
    return mRGBSpline;
  }

  void SetSplineColor( double r, double g, double b );

  double GetRadius()
  {
    return m_dRadius;
  }

  double GetSplineRadius()
  {
    return m_dSplineRadius;
  }

  int GetColorMap()
  {
    return m_nColorMap;
  }

  double GetHeatScaleMin()
  {
    return m_dHeatScaleMin;
  }

  double GetHeatScaleMid ()
  {
    return m_dHeatScaleMid;
  }

  double GetHeatScaleMax ()
  {
    return m_dHeatScaleMax;
  }

  double GetHeatScaleOffset ()
  {
    return m_dHeatScaleOffset;
  }

  void SetScalarLayer( LayerMRI* layer );
  LayerMRI* GetScalarLayer()
  {
    return m_layerScalar;
  }


  vtkRGBAColorTransferFunction* GetHeatScaleLUT();

  double GetScalarMinValue();

  double GetScalarMaxValue();

  std::vector< ScalarValues > GetScalarSets()
  {
    return m_scalarSets;
  }
  int GetNumberOfScalarSets()
  {
    return m_scalarSets.size();
  }

  bool LoadScalarsFromFile( const QString& filename );

  int GetScalarType()
  {
    return m_nScalarType;
  }

  int GetScalarSet()
  {
    return m_nScalarSet;
  }
  ScalarValues GetActiveScalarSet()
  {
    return m_scalarSets[m_nScalarSet];
  }

  bool GetSnapToVoxelCenter()
  {
    return m_bSnapToVoxelCenter;
  }

  bool GetShowSpline()
  {
    return m_bShowSpline;
  }

  int  GetType()
  {
    return m_nType;
  }

  void SetStatRange(double dMin, double dMax);

  void SetScalarToStat();

signals:
  void SnapToVoxelCenterChanged( bool bSnag );
  void SplineVisibilityChanged( bool bSpline );
  void ScalarChanged();
  void ColorMapChanged();
  void ColorChanged();
  void OpacityChanged( double );
  void RadiusChanged( double );
  void SplineRadiusChanged( double );

public slots:
  void SetOpacity( double opacity );
  void SetRadius( double r );
  void SetSplineRadius( double r );
  void SetSnapToVoxelCenter( bool bSnap );
  void SetShowSpline( bool bSpline );
  void SetType( int nType );
  void SetScalarSet( int n );
  void SetHeatScaleOffset ( double iValue );
  void SetHeatScaleMid ( double iValue );
  void SetHeatScaleMin( double iValue );
  void SetHeatScaleMax ( double iValue );
  void SetColorMap( int nColorMap );
  void SetColor( const QColor& c )
  {
    SetColor( c.redF(), c.greenF(), c.blueF() );
  }
  void SetSplineColor( const QColor& c )
  {
    SetSplineColor( c.redF(), c.greenF(), c.blueF() );
  }

private:
  void SetColorMapChanged ();
  void UpdateScalarValues();

  void LoadSettings();
  void SaveSettings();

  bool RestoreHeatscaleSettings(const QString& name);
  void SaveHeatscaleSettings();

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  // ---------------------------------------------------------------------

  double  m_dRadius;
  double  m_dSplineRadius;

  double  mOpacity;
  double  mRGB[3];
  double  mRGBSpline[3];

  int     m_nColorMap;
  vtkSmartPointer<vtkRGBAColorTransferFunction> m_lutHeatScale;
  double  m_dHeatScaleMin;
  double  m_dHeatScaleMid;
  double  m_dHeatScaleMax;
  double  m_dHeatScaleOffset;

  double m_dStatMin;
  double m_dStatMax;

  int     m_nType;
  bool    m_bSnapToVoxelCenter;
  bool    m_bShowSpline;

  int         m_nScalarType;
  LayerMRI*   m_layerScalar;
  int         m_nScalarSet;
  std::vector< ScalarValues >  m_scalarSets;
  QVariantMap m_mapHeatscaleSettings;
};

#endif
