/**
 * @file  LayerPropertyPointSet.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:50 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2007-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef LayerPropertyPointSet_h
#define LayerPropertyPointSet_h

#include "vtkSmartPointer.h"
#include "LayerProperty.h"
#include <QColor>
#include <vector>

extern "C"
{
#include "colortab.h"
}


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

  enum ScalarType { ScalarLayer = 0, ScalarSet };
  
  enum PointSetType { WayPoint = 0, ControlPoint };

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

signals:
  void SnapToVoxelCenterChanged( bool bSnag );
  void SplineVisibilityChanged( bool bSpline );
  void ScalarLayerChanged( LayerMRI* );
  void ScalarSetChanged();
  void ColorMapChanged();
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
  
  int     m_nType;
  bool    m_bSnapToVoxelCenter;
  bool    m_bShowSpline;

  int         m_nScalarType;
  LayerMRI*   m_layerScalar;
  int         m_nScalarSet;
  std::vector< ScalarValues >  m_scalarSets;
};

#endif
