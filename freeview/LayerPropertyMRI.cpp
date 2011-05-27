/**
 * @file  LayerPropertyMRI.cxx
 * @brief Implementation for MRI layer properties.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
 * Reimplemented by: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/05/27 14:16:52 $
 *    $Revision: 1.5 $
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


#include <assert.h>
#include "LayerPropertyMRI.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkImageReslice.h"
#include "vtkPointData.h"
#include "FSVolume.h"
#include "StockColorMap.h"
#include <QtGlobal>
#include <QSettings>
#include <QVariantMap>

//using namespace std;

LayerPropertyMRI::LayerPropertyMRI (QObject* parent) : LayerProperty( parent ),
  mColorMapType( LayerPropertyMRI::Grayscale ),
  mResliceInterpolation( 0 ),
  mTextureSmoothing( 0 ),
  mbClearZero( false ),
  mMinVoxelValue( 0 ),
  mMaxVoxelValue( 0 ),
  mMinVisibleValue( 0 ),
  mMaxVisibleValue( 0 ),
  mMinGrayscaleWindow( 0 ),
  mMaxGrayscaleWindow( 0 ),
  mHeatScaleMinThreshold( 0 ),
  mHeatScaleMidThreshold( 0 ),
  mHeatScaleMaxThreshold( 0 ),
  mHeatScaleOffset( 0 ),
  mbReverseHeatScale( false ),
  mbShowPositiveHeatScaleValues( true ),
  mbShowNegativeHeatScaleValues( true ),
  m_bHeatScaleClearHigh( false ),
  m_bHeatScaleTruncate( false ),
  m_bHeatScaleInvert( false ),
  mOpacity( 1 ),
  mFreeSurferCTAB( NULL ),
  mbShowAsContour( false ),
  mMinContourThreshold( 0 ),
  mMaxContourThreshold( 0 ),
  m_bDisplayVector( false ),
  m_nVectorInversion( VI_None ),
  m_nVectorRepresentation( VR_Line ),
  m_bDisplayTensor( false ),
  m_nTensorInversion( VI_None ),
  m_nTensorRepresentation( TR_Boxoid ),
  m_bContourUseImageColorMap( false ),
  m_bContourExtractAll( false ),
  m_bShowLabelOutline( false ),
  m_nUpSampleMethod( UM_None ),
  m_nContourSmoothIterations( 0 ),
  mSource( NULL )
{
  mGrayScaleTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mHeatScaleTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mColorMapTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mGrayScaleTable->ClampingOff();
  mHeatScaleTable->ClampingOn();
  m_rgbContour[0] = 0.92;
  m_rgbContour[1] = 0.78;
  m_rgbContour[2] = 0.54;

  mLUTTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourColorChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourShown(bool)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourSmoothIterationChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(DisplayModeChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(LabelOutlineChanged(bool)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ResliceInterpolationChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(TextureSmoothingChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(UpSampleMethodChanged(int)), this, SIGNAL(PropertyChanged()) );
}

LayerPropertyMRI::~LayerPropertyMRI ()
{}

void LayerPropertyMRI::CopySettings( const LayerPropertyMRI* p )
{
  blockSignals( true );
  mColorMapType           =   p->mColorMapType;
  mOpacity                =   p->mOpacity;
  mResliceInterpolation   =   p->mResliceInterpolation;
  mTextureSmoothing       =   p->mTextureSmoothing;
  mMinVisibleValue        =   p->mMinVisibleValue;
  mMaxVisibleValue        =   p->mMaxVisibleValue;
  mMinGrayscaleWindow     =   p->mMinGrayscaleWindow;
  mMaxGrayscaleWindow     =   p->mMaxGrayscaleWindow;
  mHeatScaleMinThreshold  =   p->mHeatScaleMinThreshold;
  mHeatScaleMidThreshold  =   p->mHeatScaleMidThreshold;
  mHeatScaleMaxThreshold  =   p->mHeatScaleMaxThreshold;
  mHeatScaleOffset        =   p->mHeatScaleOffset;
  mbReverseHeatScale      =   p->mbReverseHeatScale;
  mbShowPositiveHeatScaleValues =  p->mbShowPositiveHeatScaleValues;
  mbShowNegativeHeatScaleValues =  p->mbShowNegativeHeatScaleValues;
  mbClearZero             =   p->mbClearZero;
  mMinGenericThreshold    =   p->mMinGenericThreshold;
  mMaxGenericThreshold    =   p->mMaxGenericThreshold;
  m_bDisplayVector        =   p->m_bDisplayVector;
  m_nVectorInversion      =   p->m_nVectorInversion;
  m_nVectorRepresentation =   p->m_nVectorRepresentation;
  m_bHeatScaleClearHigh   =   p->m_bHeatScaleClearHigh;
  m_bHeatScaleTruncate    =   p->m_bHeatScaleTruncate;
  m_bHeatScaleInvert      =   p->m_bHeatScaleInvert;

  SetLUTCTAB  ( p->mFreeSurferCTAB );

  blockSignals( false );

  this->OnColorMapChanged();
}

void LayerPropertyMRI::RestoreSettings( const QString& filename )
{
  QSettings settings;
  QVariantMap rootmap = settings.value( "VolumeProperties" ).toMap();
  if ( !rootmap.contains( filename ) )
  {
    return;
  }

  QVariantMap map = rootmap[filename].toMap();
  RestoreSettings(map);
}

void LayerPropertyMRI::RestoreSettings(const QVariantMap& map)
{
  if ( map.contains("MinGrayscaleWindow") )
  {
    mMinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
  }

  if ( map.contains("MaxGrayscaleWindow") )
  {
    mMaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
  }

  if ( map.contains("HeatScaleMinThreshold") )
  {
    mHeatScaleMinThreshold = map["HeatScaleMinThreshold"].toDouble();
  }

  if ( map.contains("HeatScaleMidThreshold") )
  {
    mHeatScaleMidThreshold = map["HeatScaleMidThreshold"].toDouble();
  }

  if ( map.contains("HeatScaleMaxThreshold") )
  {
    mHeatScaleMaxThreshold = map["HeatScaleMaxThreshold"].toDouble();
  }

  if ( map.contains("HeatScaleOffset") )
  {
    mHeatScaleOffset = map["HeatScaleOffset"].toDouble();
  }

  if ( map.contains("MinGenericThreshold") )
  {
    mMinGenericThreshold = map["MinGenericThreshold"].toDouble();
  }

  if ( map.contains("MaxGenericThreshold") )
  {
    mMaxGenericThreshold = map["MaxGenericThreshold"].toDouble();
  }

  if ( map.contains("MinContourThreshold") )
  {
    mMinContourThreshold = map["MinContourThreshold"].toDouble();
  }

  if ( map.contains("MaxContourThreshold") )
  {
    mMaxContourThreshold = map["MaxContourThreshold"].toDouble();
  }

  this->OnColorMapChanged();
}

void LayerPropertyMRI::SaveSettings( const QString& filename )
{
  QSettings settings;
  QVariantMap rootmap = settings.value( "VolumeProperties" ).toMap();
  rootmap[filename] = GetSettings();
  settings.setValue( "VolumeProperties", rootmap );
}

QVariantMap LayerPropertyMRI::GetSettings()
{
  QVariantMap map;
  map["MinGrayscaleWindow"] = mMinGrayscaleWindow;
  map["MaxGrayscaleWindow"] = mMaxGrayscaleWindow;
  map["HeatScaleMinThreshold"] = mHeatScaleMinThreshold;
  map["HeatScaleMidThreshold"] = mHeatScaleMidThreshold;
  map["HeatScaleMaxThreshold"] = mHeatScaleMaxThreshold;
  map["HeatScaleOffset"] = mHeatScaleOffset;
  map["MinGenericThreshold"] = mMinGenericThreshold;
  map["MaxGenericThreshold"] = mMaxGenericThreshold;
  map["MinContourThreshold"] = mMinContourThreshold;
  map["MaxContourThreshold"] = mMaxContourThreshold;
  return map;
}

QVariantMap LayerPropertyMRI::GetActiveSettings()
{
  QVariantMap map;
  switch (mColorMapType)
  {
  case Grayscale:
    map["MinGrayscaleWindow"] = mMinGrayscaleWindow;
    map["MaxGrayscaleWindow"] = mMaxGrayscaleWindow;
    break;
  case Heat:
    map["HeatScaleMinThreshold"] = mHeatScaleMinThreshold;
    map["HeatScaleMidThreshold"] = mHeatScaleMidThreshold;
    map["HeatScaleMaxThreshold"] = mHeatScaleMaxThreshold;
    map["HeatScaleOffset"] = mHeatScaleOffset;
    break;
  case Jet:
  case GEColor:
  case NIH:
    map["MinGenericThreshold"] = mMinGenericThreshold;
    map["MaxGenericThreshold"] = mMaxGenericThreshold;
    break;
  }
  if (this->GetShowAsContour())
  {
      map["MinContourThreshold"] = mMinContourThreshold;
      map["MaxContourThreshold"] = mMaxContourThreshold;
  }
  return map;
}

vtkRGBAColorTransferFunction* LayerPropertyMRI::GetLUTTable () const
{
  return mLUTTable;
}

vtkRGBAColorTransferFunction* LayerPropertyMRI::GetGrayScaleTable () const
{
  return mGrayScaleTable;
}

vtkRGBAColorTransferFunction* LayerPropertyMRI::GetHeatScaleTable () const
{
  return mHeatScaleTable;
}

vtkRGBAColorTransferFunction* LayerPropertyMRI::GetColorMapTable () const
{
  return mColorMapTable;
}

vtkLookupTable* LayerPropertyMRI::GetDirectionCodedTable () const
{
  return mDirectionCodedTable;
}

vtkScalarsToColors* LayerPropertyMRI::GetActiveLookupTable()
{
  switch ( mColorMapType )
  {
  case NoColorMap:
    return NULL;
  case Grayscale:
    return mGrayScaleTable;
  case Heat:
    return mHeatScaleTable;
  case LUT:
    return mLUTTable;
  default:
    return mColorMapTable;
    break;
  }
  return NULL;
}

COLOR_TABLE* LayerPropertyMRI::GetLUTCTAB () const
{
  return mFreeSurferCTAB;
}

void LayerPropertyMRI::SetLUTCTAB( COLOR_TABLE* ct )
{
  if ( ct != mFreeSurferCTAB )
  {
    mFreeSurferCTAB = ct;
    UpdateLUTTable();
    if ( ct )
    {
      emit ColorMapChanged();
    }

  }
}

void LayerPropertyMRI::UpdateLUTTable()
{
  if ( mFreeSurferCTAB )
  {
    //  mLUTTable->BuildFromCTAB( mFreeSurferCTAB );
    mLUTTable->RemoveAllPoints();
    mLUTTable->AddRGBAPoint( 0, 0, 0, 0, 0 );
    int cEntries;
    CTABgetNumberOfTotalEntries( mFreeSurferCTAB, &cEntries );
    for ( int nEntry = 1; nEntry < cEntries; nEntry++ )
    {
      int bValid;
      CTABisEntryValid( mFreeSurferCTAB, nEntry, &bValid );
      if ( bValid )
      {
        float red, green, blue, alpha;
        CTABrgbaAtIndexf( mFreeSurferCTAB, nEntry, &red, &green, &blue, &alpha );
        mLUTTable->AddRGBAPoint( nEntry, red, green, blue, 1 );
      }
    }
    mLUTTable->Build();
  }
}

void LayerPropertyMRI::OnColorMapChanged ()
{
  if ( !mSource )
  {
    return;
  }

  double tiny_fraction = ( mMaxVoxelValue - mMinVoxelValue ) * 1e-12;
  switch ( mColorMapType )
  {
  case NoColorMap:
    break;

  case Grayscale:

    // Check the color map variables and update range sliders if
    // necessary.
    if ( mMinGrayscaleWindow < mMinVisibleValue )
    {
      mMinVisibleValue = mMinGrayscaleWindow;
    }
    if ( mMaxGrayscaleWindow > mMaxVisibleValue )
    {
      mMaxVisibleValue = mMaxGrayscaleWindow;
    }

    // Build our lookup table.
    assert( mGrayScaleTable.GetPointer() );
    mGrayScaleTable->RemoveAllPoints();
    mGrayScaleTable->AddRGBAPoint( mMinVisibleValue - tiny_fraction, 0, 0, 0, 0 );
    if ( mbClearZero )
    {
      mGrayScaleTable->AddRGBAPoint( mMinGrayscaleWindow - tiny_fraction, 0, 0, 0, 0 );
      mGrayScaleTable->AddRGBAPoint( mMinGrayscaleWindow,    0, 0, 0, 1 );
    }
    else
    {
      mGrayScaleTable->AddRGBAPoint( mMinVisibleValue,       0, 0, 0, 1 );
      mGrayScaleTable->AddRGBAPoint( mMinGrayscaleWindow,
                                     0, 0, 0, 1 );
    }
    mGrayScaleTable->AddRGBAPoint( mMinGrayscaleWindow,    0, 0, 0, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxGrayscaleWindow,    1, 1, 1, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue,       1, 1, 1, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue + tiny_fraction, 1, 1, 1, 0 );
    mGrayScaleTable->Build();

    break;

  case Heat:
    mHeatScaleTable->RemoveAllPoints();
    if ( m_bHeatScaleTruncate && m_bHeatScaleInvert )
    {
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold + mHeatScaleOffset - tiny_fraction, 1, 1, 0, 0 );
      }
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold + mHeatScaleOffset, 1, 1, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMidThreshold + mHeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMinThreshold + mHeatScaleOffset, 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  0 + mHeatScaleOffset, 0, 0, 0, 0 );
    }
    else if ( m_bHeatScaleTruncate )
    {
      mHeatScaleTable->AddRGBAPoint(  0 + mHeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold + mHeatScaleOffset, 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMidThreshold + mHeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset, 1, 1, 0, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
      }
    }
    else if ( m_bHeatScaleInvert )
    {
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold + mHeatScaleOffset, 1, 1, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMidThreshold + mHeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMinThreshold + mHeatScaleOffset, 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  0 + mHeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold + mHeatScaleOffset, 0, 0, 1, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMidThreshold + mHeatScaleOffset, 0, 0, 1, 1 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset, 0, 1, 1, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 0, 1, 1, 0 );
      }
    }
    else
    {
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold + mHeatScaleOffset, 0, 1, 1, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMidThreshold + mHeatScaleOffset, 0, 0, 1, 1 );
      mHeatScaleTable->AddRGBAPoint( -mHeatScaleMinThreshold + mHeatScaleOffset, 0, 0, 1, 0 );
      mHeatScaleTable->AddRGBAPoint(  0 + mHeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold + mHeatScaleOffset, 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMidThreshold + mHeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset, 1, 1, 0, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
      }
    }
    mHeatScaleTable->Build();
    break;

  case Jet:
    mColorMapTable->RemoveAllPoints();
    mColorMapTable->AddRGBAPoint( qMin( 0.0, mMinGenericThreshold), 0, 0, 0, 0 );
    mColorMapTable->AddRGBAPoint( mMinGenericThreshold, 0, 0, 1, 1 );
    mColorMapTable->AddRGBAPoint( mMinGenericThreshold + (mMaxGenericThreshold - mMinGenericThreshold) / 4, 0, 1, 1, 1 );
    mColorMapTable->AddRGBAPoint( (mMinGenericThreshold + mMaxGenericThreshold) / 2, 0, 1, 0, 1 );
    mColorMapTable->AddRGBAPoint( mMaxGenericThreshold - (mMaxGenericThreshold - mMinGenericThreshold) / 4, 1, 1, 0, 1 );
    mColorMapTable->AddRGBAPoint( mMaxGenericThreshold, 1, 0, 0, 1 );
    //  mColorMapTable->AddRGBAPoint( mMaxGenericThreshold + (mMaxGenericThreshold - mMinJGenericThreshold), 1, 0, 0, 1 );
    mColorMapTable->Build();
    break;

  case GEColor:
    BuildGenericLUT( stock_ge_color );
    break;

  case NIH:
    BuildGenericLUT( stock_nih );
    break;

  case LUT:
    break;

  default:
    break;
  }

  // Notify the layers that use the color map stuff.
  emit ColorMapChanged();

}

void LayerPropertyMRI::BuildGenericLUT( const int colors[256][3] )
{
  mColorMapTable->RemoveAllPoints();
  mColorMapTable->AddRGBAPoint( mMinGenericThreshold, 0, 0, 0, 0 );
  {
    double stepsize = ( mMaxGenericThreshold - mMinGenericThreshold ) / 256;
    for ( int i = 0; i < 256; i++ )
    {
      mColorMapTable->AddRGBAPoint( mMinGenericThreshold + stepsize*i,
                                    colors[i][0]/255.0,
                                    colors[i][1]/255.0,
                                    colors[i][2]/255.0,
                                    1 );
    }
  }
  mColorMapTable->Build();
}

void LayerPropertyMRI::SetDisplayVector( bool b )
{
  m_bDisplayVector = b;
  emit DisplayModeChanged();

}

void LayerPropertyMRI::SetVectorInversion( int nInversion )
{
  if ( m_bDisplayVector )
  {
    m_nVectorInversion = nInversion;
    emit DisplayModeChanged();

  }
  else if ( m_bDisplayTensor )
  {
    m_nTensorInversion = nInversion;
    emit DisplayModeChanged();

  }
}

void LayerPropertyMRI::SetVectorRepresentation( int n )
{
  if ( m_bDisplayVector )
  {
    m_nVectorRepresentation = n;
    emit DisplayModeChanged();

  }
  else if ( m_bDisplayTensor )
  {
    m_nTensorRepresentation = n;
    emit DisplayModeChanged();
  }
}

void LayerPropertyMRI::SetDisplayTensor( bool b )
{
  m_bDisplayTensor = b;
  emit DisplayModeChanged();

}

void LayerPropertyMRI::SetColorMap ( int iType )
{
  if ( mColorMapType != iType )
  {
    mColorMapType = iType;
    if ( iType == LUT )
    {
      SetUpSampleMethod( UM_None );
    }
    this->OnColorMapChanged();
  }
}

int LayerPropertyMRI::GetColorMap () const
{
  return mColorMapType;
}

void LayerPropertyMRI::SetColorMapToGrayScale ()
{
  SetColorMap( Grayscale );
}

void LayerPropertyMRI::SetColorMapToHeatScale ()
{
  SetColorMap( Heat );
}

void LayerPropertyMRI::SetColorMapToLUT ()
{
  SetColorMap( LUT );
}

void LayerPropertyMRI::SetResliceInterpolation ( int iMode )
{
  if ( mResliceInterpolation != iMode )
  {
    mResliceInterpolation = iMode;
    emit ResliceInterpolationChanged();

  }
}

int LayerPropertyMRI::GetResliceInterpolation () const
{
  return mResliceInterpolation;
}

void LayerPropertyMRI::SetTextureSmoothing ( int iSmooth )
{
  if ( mTextureSmoothing != iSmooth )
  {
    mTextureSmoothing = iSmooth;
    emit TextureSmoothingChanged();

  }
}

int LayerPropertyMRI::GetTextureSmoothing () const
{
  return mTextureSmoothing;
}

void LayerPropertyMRI::SetMinVisibleValue ( double iValue )
{
  if ( mMinVisibleValue != iValue )
  {
    mMinVisibleValue = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetMinVisibleValue ()
{
  return mMinVisibleValue;
}

void LayerPropertyMRI::SetMaxVisibleValue ( double iValue )
{
  if ( mMaxVisibleValue != iValue )
  {
    mMaxVisibleValue = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetMaxVisibleValue ()
{
  return mMaxVisibleValue;
}

void LayerPropertyMRI::SetMinMaxVisibleValue ( double iMinValue,
    double iMaxValue )
{
  if ( mMinVisibleValue != iMinValue ||
       mMaxVisibleValue != iMaxValue )
  {
    mMinVisibleValue = iMinValue;
    mMaxVisibleValue = iMaxValue;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetWindow ( double iWindow )
{
  double level = this->GetLevel();
  this->SetMinMaxGrayscaleWindow( level - iWindow/2.0,
                                  level + iWindow/2.0 );
}

double LayerPropertyMRI::GetWindow ()
{
  return ( mMaxGrayscaleWindow - mMinGrayscaleWindow );
}

void LayerPropertyMRI::SetLevel ( double iLevel )
{
  double window = this->GetWindow();
  this->SetMinMaxGrayscaleWindow( iLevel - window/2.0,
                                  iLevel + window/2.0 );
}

double LayerPropertyMRI::GetLevel ()
{
  return ((mMaxGrayscaleWindow - mMinGrayscaleWindow) / 2.0) +
         mMinGrayscaleWindow;
}

void LayerPropertyMRI::SetWindowLevel( double iWindow, double iLevel )
{
  this->SetMinMaxGrayscaleWindow( iLevel - iWindow/2.0, iLevel + iWindow/2.0 );
}

void LayerPropertyMRI::SetMinMaxGrayscaleWindow ( double iMin, double iMax )
{
  if ( mMinGrayscaleWindow != iMin || mMaxGrayscaleWindow != iMax )
  {
    mMinGrayscaleWindow = iMin;
    mMaxGrayscaleWindow = iMax;
    this->OnColorMapChanged();
    emit WindowLevelChanged();
  }
}

void LayerPropertyMRI::SetMinGrayscaleWindow ( double iMin )
{
  if ( mMinGrayscaleWindow != iMin )
  {
    mMinGrayscaleWindow = iMin;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetMaxGrayscaleWindow ( double iMax )
{
  if ( mMaxGrayscaleWindow != iMax )
  {
    mMaxGrayscaleWindow = iMax;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetHeatScaleMinThreshold ( double iValue )
{
  if ( mHeatScaleMinThreshold != iValue )
  {
    mHeatScaleMinThreshold = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleMinThreshold ()
{
  return mHeatScaleMinThreshold;
}

void LayerPropertyMRI::SetHeatScaleMidThreshold ( double iValue )
{
  if ( mHeatScaleMidThreshold != iValue )
  {
    mHeatScaleMidThreshold = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleMidThreshold ()
{
  return mHeatScaleMidThreshold;
}

void LayerPropertyMRI::SetHeatScaleMaxThreshold ( double iValue )
{
  if ( mHeatScaleMaxThreshold != iValue )
  {
    mHeatScaleMaxThreshold = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleMaxThreshold ()
{
  return mHeatScaleMaxThreshold;
}

void LayerPropertyMRI::SetHeatScale( double dMin, double dMid, double dMax )
{
  mHeatScaleMinThreshold = dMin;
  mHeatScaleMidThreshold = dMid;
  mHeatScaleMaxThreshold = dMax;
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetHeatScaleOffset ( double iValue )
{
  if ( mHeatScaleOffset != iValue )
  {
    mHeatScaleOffset = iValue;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleOffset ()
{
  return mHeatScaleOffset;
}

void LayerPropertyMRI::SetReverseHeatScale ( bool ib )
{
  if ( mbReverseHeatScale != ib )
  {
    mbReverseHeatScale = ib;
    this->OnColorMapChanged();
  }
}

bool LayerPropertyMRI::GetReverseHeatScale ()
{
  return mbReverseHeatScale;
}

void LayerPropertyMRI::SetHeatScaleClearHigh( bool bClear )
{
  m_bHeatScaleClearHigh = bClear;
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetHeatScaleTruncate( bool bTruncate )
{
  m_bHeatScaleTruncate = bTruncate;
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetHeatScaleInvert( bool bInvert )
{
  m_bHeatScaleInvert = bInvert;
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetShowPositiveHeatScaleValues ( bool ib )
{
  if ( mbShowPositiveHeatScaleValues != ib )
  {
    mbShowPositiveHeatScaleValues = ib;
    this->OnColorMapChanged();
  }
}

bool LayerPropertyMRI::GetShowPositiveHeatScaleValues ()
{
  return mbShowPositiveHeatScaleValues;
}

void LayerPropertyMRI::SetShowNegativeHeatScaleValues ( bool ib )
{
  if ( mbShowNegativeHeatScaleValues != ib )
  {
    mbShowNegativeHeatScaleValues = ib;
    this->OnColorMapChanged();
  }
}

bool LayerPropertyMRI::GetShowNegativeHeatScaleValues ()
{
  return mbShowNegativeHeatScaleValues;
}

double LayerPropertyMRI::GetOpacity() const
{
  return mOpacity;
}

void LayerPropertyMRI::SetOpacity( double opacity )
{
  if ( mOpacity != opacity && opacity >= 0 && opacity <= 1 )
  {
    mOpacity = opacity;
    emit OpacityChanged( opacity );
  }
}

bool LayerPropertyMRI::GetClearZero()
{
  return mbClearZero;
}

void LayerPropertyMRI::SetClearZero( bool bClear )
{
  if ( mbClearZero != bClear )
  {
    mbClearZero = bClear;
    // this->SendBroadcast( "ClearZeroChanged", this );
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetMinValue()
{
  return mMinVoxelValue;
}

double LayerPropertyMRI::GetMaxValue()
{
  return mMaxVoxelValue;
}


bool LayerPropertyMRI::UpdateValueRange( double dValue )
{
  if ( mMinVoxelValue <= dValue && mMaxVoxelValue >= dValue )
  {
    return false;
  }

  if ( mMinVoxelValue > dValue )
  {
    mMinVoxelValue = dValue;
  }
  else if ( mMaxVoxelValue < dValue )
  {
    mMaxVoxelValue = dValue;
  }

  mMinVisibleValue = mMinVoxelValue - ( mMaxVoxelValue - mMinVoxelValue );
  mMaxVisibleValue = mMaxVoxelValue + ( mMaxVoxelValue - mMinVoxelValue );
  mWindowRange[0] = 0;
  mWindowRange[1] = (mMaxVoxelValue - mMinVoxelValue) * 2;
  mLevelRange[0] = mMinVoxelValue;
  mLevelRange[1] = mMaxVoxelValue;

//  this->OnColorMapChanged();
  return true;
}

void LayerPropertyMRI::SetVolumeSource ( FSVolume* source )
{
  // Source object reads the volume and outputs a volume.
  mSource = source;

  // Init our color scale values.
  mMinVoxelValue = mSource->GetMinValue();
  mMaxVoxelValue = mSource->GetMaxValue();
  mMinGrayscaleWindow = mMinVoxelValue;
  mMaxGrayscaleWindow = mMaxVoxelValue;
  mMinVisibleValue = mMinVoxelValue - ( mMaxVoxelValue - mMinVoxelValue );
  mMaxVisibleValue = mMaxVoxelValue + ( mMaxVoxelValue - mMinVoxelValue );
  mWindowRange[0] = 0;
  mWindowRange[1] = (mMaxVoxelValue - mMinVoxelValue) * 2;
  mLevelRange[0] = mMinVoxelValue;
  mLevelRange[1] = mMaxVoxelValue;

  double oneTenth;
  double highestAbsValue;
  highestAbsValue = qMax( fabs(mMinVoxelValue), fabs(mMaxVoxelValue) );
  oneTenth = highestAbsValue / 10.0;
  mHeatScaleMinThreshold = mMinGrayscaleWindow;
  mHeatScaleMidThreshold = highestAbsValue / 2.0;
  mHeatScaleMaxThreshold = highestAbsValue - oneTenth;

  UpdateLUTTable();

  // Init the colors in the heat scale. The editor will handle the
  // editing from here.

  mMinGenericThreshold = mMinVoxelValue;
  mMaxGenericThreshold = mMaxVoxelValue;
  mColorMapTable->ClampingOn();

  mDirectionCodedTable = vtkSmartPointer<vtkLookupTable>::New();
  int ns = 64;
  mDirectionCodedTable->SetNumberOfTableValues(ns*ns*ns);
  mDirectionCodedTable->ForceBuild();
  mDirectionCodedTable->SetTableRange(0, ns*ns*ns-1);
  for (int i = 0; i < ns; i++)
  {
    for (int j = 0; j < ns; j++)
    {
      for (int k = 0; k < ns; k++)
      {
        mDirectionCodedTable->SetTableValue(i*ns*ns+j*ns+k, (float)i/ns, (float)j/ns, (float)k/ns);
      }
    }
  }

  mMinContourThreshold = mHeatScaleMidThreshold;
  mMaxContourThreshold = mSource->GetMaxValue();

  if ( source->GetEmbeddedColorTable() )
  {
    SetColorMap( LUT );
    SetLUTCTAB( source->GetEmbeddedColorTable() );
  }

  // Set up our initial tables.
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetMinMaxGenericThreshold ( double iMin, double iMax )
{
  if ( mMinGenericThreshold != iMin || mMaxGenericThreshold != iMax )
  {
    mMinGenericThreshold = iMin;
    mMaxGenericThreshold = iMax;
    this->OnColorMapChanged();
    emit GenericWindowLevelChanged();
  }
}

double LayerPropertyMRI::GetMinGenericThreshold()
{
  return mMinGenericThreshold;
}

void LayerPropertyMRI::SetMinGenericThreshold ( double iMin )
{
  if ( mMinGenericThreshold != iMin )
  {
    mMinGenericThreshold = iMin;
    this->OnColorMapChanged();
  }
}


double LayerPropertyMRI::GetMaxGenericThreshold()
{
  return mMaxGenericThreshold;
}

void LayerPropertyMRI::SetMaxGenericThreshold ( double iMax )
{
  if ( mMaxGenericThreshold != iMax )
  {
    mMaxGenericThreshold = iMax;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetShowAsContour( bool bContour )
{
  if ( mbShowAsContour != bContour )
  {
    mbShowAsContour = bContour;
    emit ContourShown( mbShowAsContour );

  }
}

void LayerPropertyMRI::SetContourMinThreshold( double dValue )
{
  if ( mMinContourThreshold != dValue )
  {
    mMinContourThreshold = dValue;
    emit ContourChanged();

  }
}

void LayerPropertyMRI::SetContourMaxThreshold( double dValue )
{
  if ( mMaxContourThreshold != dValue )
  {
    mMaxContourThreshold = dValue;
    emit ContourChanged();

  }
}

void LayerPropertyMRI::SetContourThreshold( double dMin, double dMax )
{
  if ( mMinContourThreshold != dMin || mMaxContourThreshold != dMax )
  {
    mMinContourThreshold = dMin;
    mMaxContourThreshold = dMax;
    emit ContourChanged();

  }
}

void LayerPropertyMRI::SetContourExtractAllRegions( bool bExtractAll )
{
  if ( m_bContourExtractAll != bExtractAll )
  {
    m_bContourExtractAll = bExtractAll;
    emit ContourChanged();

  }
}


void LayerPropertyMRI::SetContourColor( double r, double g, double b )
{
  m_rgbContour[0] = r;
  m_rgbContour[1] = g;
  m_rgbContour[2] = b;

  emit ContourColorChanged();
}

void LayerPropertyMRI::SetContourUseImageColorMap( bool bFlag )
{
  if ( m_bContourUseImageColorMap != bFlag )
  {
    m_bContourUseImageColorMap = bFlag;

    emit ContourColorChanged();
  }
}

void LayerPropertyMRI::SetShowLabelOutline( bool bOutline )
{
  if ( bOutline )
  {
    SetUpSampleMethod( 0 );
  }

  if ( m_bShowLabelOutline != bOutline )
  {
    m_bShowLabelOutline = bOutline;

    emit LabelOutlineChanged( bOutline );

  }
}

void LayerPropertyMRI::SetUpSampleMethod( int nSampleMethod )
{
  if ( nSampleMethod > 0 )
  {
    this->SetShowLabelOutline( false );
  }

  if ( nSampleMethod != m_nUpSampleMethod )
  {
    m_nUpSampleMethod = nSampleMethod;

    emit UpSampleMethodChanged( nSampleMethod );

  }
}

void LayerPropertyMRI::SetContourSmoothIterations( int nIterations )
{
  if ( m_nContourSmoothIterations != nIterations )
  {
    m_nContourSmoothIterations = nIterations;
    emit ContourSmoothIterationChanged( nIterations );

  }
}
