/**
 * @file  LayerPropertiesMRI.h
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

#ifndef LayerPropertiesMRI_h
#define LayerPropertiesMRI_h

#include "vtkSmartPointer.h"
#include "LayerProperties.h"

extern "C"
{
#include "colortab.h"
}

class vtkFreesurferLookupTable;
class vtkRGBAColorTransferFunction;
class vtkLookupTable;
class vtkScalarsToColors;
class FSVolume;

class LayerPropertiesMRI : public LayerProperties
{

public:
  LayerPropertiesMRI ();
  ~LayerPropertiesMRI ();

  // Description:
  // The color map types in which a volume can be drawn.
  enum ColorMapType {
    NoColorMap=-1, Grayscale, LUT, Heat, Jet, GEColor, NIH, DirectionCoded 
  };
  
  enum VectorInversion {
    VI_None = 0, VI_X, VI_Y, VI_Z
  };

  enum VectorRepresentation {
    VR_Line = 0, VR_Bar 
  };
  
  enum TensorRepresentation {
    TR_Boxoid = 0, TR_Ellipsoid 
  };
  
  enum UpSampleMethod {
    UM_None = 0, UM_NearestNeighbor, UM_BiLinear
  };
  
  void CopySettings  ( const LayerPropertiesMRI* p );
  void RestoreSettings( const char* filename );
  void SaveSettings   ( const char* filename );
  
  void SetVolumeSource( FSVolume* source );

  double  GetOpacity() const;
  void    SetOpacity( double opacity );

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
  virtual void ColorMapChanged ();

  // Description.
  // Set the reslice object's interpolation.
  void  SetResliceInterpolation ( int iMode );
  int   GetResliceInterpolation () const;

  // Description.
  // Set the texture's interpolation (we call it smoothing to
  // differentiate from the reslice interpolation).
  void  SetTextureSmoothing ( int iSmooth );
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
  void    SetWindow ( double iWindow );
  double  GetWindow ();
  void    SetLevel ( double iLevel );
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
  void SetHeatScaleClearHigh( bool bClear );

  bool GetHeatScaleTruncate()
  {
    return m_bHeatScaleTruncate;
  }
  void SetHeatScaleTruncate( bool bTruncate );
   
  bool GetHeatScaleInvert()
  {
    return m_bHeatScaleInvert;
  }
  void SetHeatScaleInvert( bool bInvert );
  
  // Description:
  // Various heat scale characterists.
  void SetReverseHeatScale ( bool ib );
  bool GetReverseHeatScale ();
  void SetShowPositiveHeatScaleValues ( bool ib );
  bool GetShowPositiveHeatScaleValues ();
  void SetShowNegativeHeatScaleValues ( bool ib );
  bool GetShowNegativeHeatScaleValues ();

  bool GetClearZero();
  void SetClearZero( bool bClear );

  // void EditorChangedHeatScale ();

  // Description:
  // Pop up a dialog box asking for a file name, and pass the file
  // name to LoadFreesurferLUT.
  // void LoadLUTFromDlog ();

  // Description:
  // Load the LUT into the LUT color map.
  void LoadFreesurferLUT ( const char* ifn );

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

  void SetDisplayVector( bool b );
  bool GetDisplayVector()
  {
    return m_bDisplayVector;
  }
  
  void SetDisplayTensor( bool b );
  bool GetDisplayTensor()
  {
    return m_bDisplayTensor;
  }
  
  void SetVectorInversion( int n );
  int GetVectorInversion()
  {
    return m_nVectorInversion;
  }
  
  void SetTensorInversion( int n );
  int GetTensorInversion()
  {
    return m_nTensorInversion;
  }
  
  void SetVectorRepresentation( int n );
  bool GetVectorRepresentation()
  {
    return m_nVectorRepresentation;
  }
  
  void SetTensorRepresentation( int n );
  bool GetTensorRepresentation()
  {
    return m_nTensorRepresentation;
  }
  
  void SetShowAsContour( bool bContour );
  bool GetShowAsContour()
  {
    return mbShowAsContour;
  }

  void SetContourMinThreshold( double dValue );
  double GetContourMinThreshold()
  {
    return mMinContourThreshold;
  }
  
  void SetContourMaxThreshold( double dValue );
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
  void SetContourUseImageColorMap( bool bFlag );
  
  bool GetContourExtractAllRegions()
  {
    return m_bContourExtractAll;
  }
  void SetContourExtractAllRegions( bool bExtractAll );
  
  bool GetShowLabelOutline()
  {
    return m_bShowLabelOutline;
  }
  void SetShowLabelOutline( bool bOutline );
  
  int GetUpSampleMethod()
  {
    return m_nUpSampleMethod;
  }
  void SetUpSampleMethod( int nUpSampleMethod );
  
  int GetContourSmoothIterations()
  {
    return m_nContourSmoothIterations;
  }
  
  void SetContourSmoothIterations( int nIterations );
  
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

  bool    m_bShowLabelOutline;
  int     m_nUpSampleMethod;
  
  // ---------------------------------------------------------------------

  FSVolume*   mSource;
  std::string mfnVolume;
  //ETX

};

#endif
