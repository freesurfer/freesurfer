/**
 * @file  LayerPropertiesWayPoints.h
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
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef LayerPropertiesWayPoints_h
#define LayerPropertiesWayPoints_h

#include "vtkSmartPointer.h"
#include "LayerProperties.h"

extern "C"
{
#include "colortab.h"
}


struct ScalarValues
{
  std::string strName;
  int  nNum;
  double* dValue;
  double dMin;
  double dMax;
};

class vtkRGBAColorTransferFunction;
class LayerMRI;

class LayerPropertiesWayPoints : public LayerProperties
{
public:
  LayerPropertiesWayPoints ();
  ~LayerPropertiesWayPoints ();

  enum ColorMapType { SolidColor = 0, HeatScale };

  enum ScalarType { ScalarLayer = 0, ScalarSet };
  
  enum PointSetType { WayPoints = 0, ControlPoints };

  double GetOpacity() const;
  void SetOpacity( double opacity );

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

  void SetRadius( double r );

  double GetSplineRadius()
  {
    return m_dSplineRadius;
  }

  void SetSplineRadius( double r );

  void SetColorMap( int nColorMap );
  int GetColorMap()
  {
    return m_nColorMap;
  }

  void SetHeatScaleMin( double iValue );
  double GetHeatScaleMin()
  {
    return m_dHeatScaleMin;
  }

  void SetHeatScaleMid ( double iValue );
  double GetHeatScaleMid ()
  {
    return m_dHeatScaleMid;
  }

  void SetHeatScaleMax ( double iValue );
  double GetHeatScaleMax ()
  {
    return m_dHeatScaleMax;
  }

  void SetHeatScaleOffset ( double iValue );
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

  bool LoadScalarsFromFile( const char* filename );

  int GetScalarType()
  {
    return m_nScalarType;
  }

  void SetScalarSet( int n );
  int GetScalarSet()
  {
    return m_nScalarSet;
  }
  ScalarValues GetActiveScalarSet()
  {
    return m_scalarSets[m_nScalarSet];
  }
  
  void SetSnapToVoxelCenter( bool bSnap );
  bool GetSnapToVoxelCenter()
  {
    return m_bSnapToVoxelCenter;
  }
  
  void SetShowSpline( bool bSpline );
  bool GetShowSpline()
  {
    return m_bShowSpline;
  }
  
  void SetType( int nType );
  int  GetType()
  {
    return m_nType;
  }

protected:
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );

private:
  void ColorMapChanged ();
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
