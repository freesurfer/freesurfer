/**
 * @file  LayerPropertyMRI.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
 * Reimplemented by: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/02/05 20:51:41 $
 *    $Revision: 1.11 $
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
 *
 */

#ifndef LayerPropertyMRI_h
#define LayerPropertyMRI_h

#include "LayerProperty.h"
#include "vtkSmartPointer.h"
#include <QVariantMap>
#include <QColor>

extern "C"
{
#include "colortab.h"
}

class vtkFreesurferLookupTable;
class vtkRGBAColorTransferFunction;
class vtkLookupTable;
class vtkScalarsToColors;
class FSVolume;
class QString;

class LayerPropertyMRI : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyMRI ( QObject* parent = NULL );
  ~LayerPropertyMRI ();

  // Description:
  // The color map types in which a volume can be drawn.
  enum ColorMapType
  {
    NoColorMap=-1, Grayscale, LUT, Heat, Jet, GEColor, NIH, DirectionCoded
  };

  enum VectorInversion
  {
    VI_None = 0, VI_X, VI_Y, VI_Z
  };

  enum VectorRepresentation
  {
    VR_Line = 0, VR_Bar
  };

  enum TensorRepresentation
  {
    TR_Boxoid = 0, TR_Ellipsoid
  };

  enum UpSampleMethod
  {
    UM_None = 0, UM_NearestNeighbor, UM_Linear, UM_Cubic
  };

  QVariantMap GetSettings();
  QVariantMap GetActiveSettings();
  QVariantMap GetFullSettings();
  void CopySettings  ( const LayerPropertyMRI* p );
  void RestoreSettings( const QVariantMap& map);
  void RestoreSettings( const QString& filename );
  void SaveSettings   ( const QString& filename );
  void RestoreFullSettings( const QVariantMap& map );

  void SetVolumeSource( FSVolume* source );

  double  GetOpacity() const;

  vtkRGBAColorTransferFunction*  GetLUTTable()  const;
  vtkRGBAColorTransferFunction*  GetGrayScaleTable () const;
  vtkRGBAColorTransferFunction*  GetHeatScaleTable () const;
  vtkRGBAColorTransferFunction*  GetColorMapTable () const;
  vtkLookupTable*           GetDirectionCodedTable () const;

  COLOR_TABLE*  GetLUTCTAB () const;
  void          SetLUTCTAB ( COLOR_TABLE* ct );

  virtual vtkScalarsToColors* GetActiveLookupTable();

  // Description:
  // Set the color map appropriately.
  //BTX
  void SetColorMap ( int iType );
  int GetColorMap () const;
  //ETX
  void SetColorMapToGrayScale ();
  void SetColorMapToHeatScale ();
  void SetColorMapToLUT ();
  virtual void OnColorMapChanged ();

  // Description.
  // Set the reslice object's interpolation.
  int   GetResliceInterpolation () const;

  int   GetTextureSmoothing () const;

  // Description:
  // This determines the levels above and below which voxels are drawn
  // as black, or no color.
  void    SetMinVisibleValue ( double iValue );
  double  GetMinVisibleValue ();
  void    SetMaxVisibleValue ( double iValue );
  double  GetMaxVisibleValue ();
  void    SetMinMaxVisibleValue ( double iMinValue, double iMaxValue );

  // Description:
  // These determine the visible value range. Values below the range
  // are drawn at the 'lowest' color (e.g. black) and values above are
  // drawn in the 'highest' color (e.g. white).
  double  GetWindow ();
  double  GetLevel ();
  void    SetWindowLevel( double iWindow, double iLevel );
  void    SetMinMaxGrayscaleWindow ( double iMin, double iMax );
  void    SetMinGrayscaleWindow ( double iMin );
  void    SetMaxGrayscaleWindow ( double iMax );

  void    SetMinMaxGenericThreshold ( double iMin, double iMax );
  double  GetMinGenericThreshold();
  void    SetMinGenericThreshold ( double iMin );
  double  GetMaxGenericThreshold();
  void    SetMaxGenericThreshold ( double iMax );

  // Description:
  // These determine the heatscale color map. The threshold is mirrored:
  // -> cyan -> blue -> trans_blue -> clear -> trans_orange -> orange -> red ->
  // -> -max -> -mid ->    -min    ->   0   ->     min      ->   mid  -> max ->
  void    SetHeatScaleMinThreshold ( double iValue );
  double  GetHeatScaleMinThreshold ();
  void    SetHeatScaleMidThreshold ( double iValue );
  double  GetHeatScaleMidThreshold ();
  void    SetHeatScaleMaxThreshold ( double iValue );
  double  GetHeatScaleMaxThreshold ();
  void    SetHeatScaleOffset ( double iValue );
  double  GetHeatScaleOffset ();
  void    SetHeatScale( double dMin, double dMid, double dMax );

  bool GetHeatScaleClearHigh()
  {
    return m_bHeatScaleClearHigh;
  }

  bool GetHeatScaleTruncate()
  {
    return m_bHeatScaleTruncate;
  }

  bool GetHeatScaleInvert()
  {
    return m_bHeatScaleInvert;
  }

  // Description:
  // Various heat scale characterists.
  void SetReverseHeatScale ( bool ib );
  bool GetReverseHeatScale ();
  void SetShowPositiveHeatScaleValues ( bool ib );
  bool GetShowPositiveHeatScaleValues ();
  void SetShowNegativeHeatScaleValues ( bool ib );
  bool GetShowNegativeHeatScaleValues ();

  bool GetClearZero();

  // void EditorChangedHeatScale ();

  // Description:
  // Load the LUT into the LUT color map.
  void LoadFreesurferLUT ( const QString& ifn );

  double GetMinValue();
  double GetMaxValue();

  bool UpdateValueRange( double dValue );

  double* GetWindowRange()
  {
    return mWindowRange;
  }

  double* GetLevelRange()
  {
    return mLevelRange;
  }

  bool GetDisplayVector()
  {
    return m_bDisplayVector;
  }

  bool GetDisplayTensor()
  {
    return m_bDisplayTensor;
  }

  int GetVectorInversion()
  {
    return m_nVectorInversion;
  }

  int GetTensorInversion()
  {
    return m_nTensorInversion;
  }

  bool GetVectorRepresentation()
  {
    return m_nVectorRepresentation;
  }

  bool GetTensorRepresentation()
  {
    return m_nTensorRepresentation;
  }

  bool GetShowAsContour()
  {
    return mbShowAsContour;
  }

  double GetContourMinThreshold()
  {
    return mMinContourThreshold;
  }

  double GetContourMaxThreshold()
  {
    return mMaxContourThreshold;
  }

  void SetContourThreshold( double dMin, double dMax );

  void SetContourColor( double r, double g, double b );
  void GetContourColor( double* c )
  {
    c[0] = m_rgbContour[0];
    c[1] = m_rgbContour[1];
    c[2] = m_rgbContour[2];
  }

  double* GetContourColor()
  {
    return m_rgbContour;
  }

  bool GetContourUseImageColorMap()
  {
    return m_bContourUseImageColorMap;
  }

  bool GetContourExtractAllRegions()
  {
    return m_bContourExtractAll;
  }

  bool GetShowLabelOutline()
  {
    return m_bShowLabelOutline;
  }

  int GetUpSampleMethod()
  {
    return m_nUpSampleMethod;
  }

  int GetContourSmoothIterations()
  {
    return m_nContourSmoothIterations;
  }

  bool GetShowProjectionMap()
  {
    return m_bShowProjectionMap;
  }

  bool GetRememberFrameSettings()
  {
    return m_bRememberFrameSettings;
  }

  bool GetShowAsLabelContour()
  {
    return this->m_bShowAsLabelContour;
  }

  bool GetContourUpsample()
  {
    return this->m_bContourUpsample;
  }

  void SetActiveFrame(int nFrame);

  void GetLabelContourRange(double* th1, double* th2)
  {
    *th1 = m_dLabelContourRange[0];
    *th2 = m_dLabelContourRange[1];
  }

public slots:
  void SetOpacity( double opacity );
  void SetUpSampleMethod( int nUpSampleMethod );
  void SetShowLabelOutline( bool bOutline );
  void SetContourSmoothIterations( int nIterations );
  void SetTextureSmoothing ( int iSmooth );
  void SetShowAsContour( bool bContour );
  void SetShowAsLabelContour(bool bLabelContour);
  void SetClearZero( bool bClear );
  void SetResliceInterpolation ( int iMode );
  void SetWindow( double iWindow );
  void SetLevel ( double iLevel );
  void SetDisplayVector( bool b );
  void SetDisplayTensor( bool b );
  void SetVectorInversion( int n );
  void SetVectorRepresentation( int n );
  void SetContourMinThreshold( double dValue );
  void SetContourMaxThreshold( double dValue );
  void SetHeatScaleClearHigh( bool bClear );
  void SetHeatScaleTruncate( bool bTruncate );
  void SetHeatScaleInvert( bool bInvert );
  void SetContourUseImageColorMap( bool bFlag );
  void SetContourUpsample( bool bFlag );
  void SetContourExtractAllRegions( bool bExtractAll );
  void SetContourColor(const QColor& c)
  {
    SetContourColor(c.redF(), c.greenF(), c.blueF());
  }
  void SetShowProjectionMap(bool bShow);
  void SetRememberFrameSettings(bool bFlag);

  void SetLabelContourRange(double dmin, double dmax);

signals:
  void ColorMapChanged();
  void ResliceInterpolationChanged();
  void TextureSmoothingChanged();
  void OpacityChanged( double );
  void WindowLevelChanged();
  void GenericWindowLevelChanged();
  void ContourShown( bool bShown );
  void ContourChanged();
  void ContourColorChanged();
  void ContourSmoothIterationChanged( int );
  void LabelOutlineChanged( bool bOutline );
  void UpSampleMethodChanged( int nMethod );
  void ProjectionMapShown(bool bShown);

private:

  void UpdateLUTTable();
  void BuildGenericLUT( const int colors[256][3] );

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  vtkSmartPointer<vtkRGBAColorTransferFunction> mGrayScaleTable;
  vtkSmartPointer<vtkRGBAColorTransferFunction> mHeatScaleTable;
  vtkSmartPointer<vtkRGBAColorTransferFunction> mColorMapTable;
  vtkSmartPointer<vtkLookupTable> mDirectionCodedTable;
  // ---------------------------------------------------------------------


  // Color map variables -------------------------------------------------
  int mColorMapType;

  // Applicable to all types.
  int     mResliceInterpolation;
  int     mTextureSmoothing;

  bool    mbClearZero;
  double  mMinVoxelValue, mMaxVoxelValue;

  double  mMinVisibleValue, mMaxVisibleValue;

  // For grayscale drawing.
  double  mMinGrayscaleWindow, mMaxGrayscaleWindow;

  // For heatScale drawing.
  double  mHeatScaleMinThreshold, mHeatScaleMidThreshold, mHeatScaleMaxThreshold;
  double  mHeatScaleOffset;
  bool    mbReverseHeatScale;
  bool    mbShowPositiveHeatScaleValues;
  bool    mbShowNegativeHeatScaleValues;
  bool    m_bHeatScaleClearHigh;
  bool    m_bHeatScaleTruncate;
  bool    m_bHeatScaleInvert;

  double  mMinGenericThreshold;
  double  mMaxGenericThreshold;

  double  mOpacity;

  double  mWindowRange[2];
  double  mLevelRange[2];

  bool    m_bRememberFrameSettings;
  QVariantMap m_frameSettings;

  // LUT drawing.
  COLOR_TABLE* mFreeSurferCTAB;

  bool    mbShowAsContour;
  double  mMinContourThreshold;
  double  mMaxContourThreshold;

  bool    m_bDisplayVector;
  int     m_nVectorInversion;
  int     m_nVectorRepresentation;

  bool    m_bDisplayTensor;
  int     m_nTensorInversion;
  int     m_nTensorRepresentation;

  double  m_rgbContour[3];
  bool    m_bContourUseImageColorMap;
  bool    m_bContourExtractAll;
  int     m_nContourSmoothIterations;
  bool    m_bContourUpsample;

  bool    m_bShowAsLabelContour;
  double  m_dLabelContourRange[2];

  bool    m_bShowLabelOutline;
  int     m_nUpSampleMethod;

  bool    m_bShowProjectionMap;

  // ---------------------------------------------------------------------

  FSVolume*   mSource;
  int     m_nActiveFrame;
  QString mfnVolume;
  //ETX

};

#endif
