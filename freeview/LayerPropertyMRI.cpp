/**
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


#include <assert.h>
#include "LayerPropertyMRI.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkImageReslice.h"
#include "vtkPointData.h"
#include "vtkImageData.h"
#include "FSVolume.h"
#include "LayerMRI.h"
#include "StockColorMap.h"
#include "TurboColorMap.h"
#include <QtGlobal>
#include <QSettings>
#include <QVariantMap>
#include <QDebug>

//using namespace std;

LayerPropertyMRI::LayerPropertyMRI (QObject* parent) : LayerProperty( parent ),
  mColorMapType( LayerPropertyMRI::Grayscale ),
  mResliceInterpolation( 0 ),
  mTextureSmoothing( 0 ),
  mbClearBackground( false ),
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
  m_bShow2DContour( false ),
  m_bDisplayVector( false ),
  m_nVectorInversion( VI_None ),
  m_nVectorRepresentation( VR_Line ),
  m_bDisplayTensor( false ),
  m_nTensorInversion( VI_None ),
  m_nTensorRepresentation( TR_Boxoid ),
  m_bContourUseImageColorMap( false ),
  m_bContourExtractAll( true ),
  m_bContourDilateFirst(false),
  m_bShowLabelOutline( false ),
  m_nUpSampleMethod( UM_None ),
  m_nContourSmoothIterations( 5 ),
  mSource( NULL ),
  m_bRememberFrameSettings( false ),
  m_nActiveFrame( 0 ),
  m_bShowAsLabelContour( false ),
  m_bShowVoxelizedContour(false),
  m_bContourUpsample(false),
  m_dVectorScale(1.0),
  m_bUsePercentile(false),
  m_bAutoAdjustFrameLevel(false),
  m_bNormalizeVector(true),
  m_dVectorDisplayScale(1.0),
  m_bHeatScaleAutoMid(true),
  m_bHeatScaleSetMidToMin(false),
  m_nProjectionMapType(0),
  m_bDisplayRGB(false),
  m_dVectorLineWidth(1),
  m_nVectorSkip(0),
  m_dVectorNormThreshold(0),
  m_colorBinary(Qt::white)
{
  mGrayScaleTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mHeatScaleTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mColorMapTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();
  mGrayScaleTable->ClampingOn();
  mHeatScaleTable->ClampingOn();
  mColorMapTable->ClampingOn();
  m_rgbContour[0] = 0.92;
  m_rgbContour[1] = 0.78;
  m_rgbContour[2] = 0.54;

  for (int i = 0; i < 3; i++)
  {
    m_nProjectionMapRange[i*2] = 0;
    m_nProjectionMapRange[i*2+1] = -1;
  }

  mLUTTable = vtkSmartPointer<vtkRGBAColorTransferFunction>::New();

  connect( this, SIGNAL(ColorMapChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourColorChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourShown(bool)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ContourSmoothIterationChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(DisplayModeChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(LabelOutlineChanged(bool)), this, SIGNAL(PropertyChanged()) );
//  connect( this, SIGNAL(OpacityChanged(double)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(ResliceInterpolationChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(TextureSmoothingChanged()), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(UpSampleMethodChanged(int)), this, SIGNAL(PropertyChanged()) );
  connect( this, SIGNAL(VectorLineWidthChanged(double)), this, SIGNAL(PropertyChanged()) );
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
  mbClearBackground             =   p->mbClearBackground;
  mClearBackgroundValue   =   p->mClearBackgroundValue;
  mMinGenericThreshold    =   p->mMinGenericThreshold;
  mMaxGenericThreshold    =   p->mMaxGenericThreshold;
  m_bDisplayVector        =   p->m_bDisplayVector;
  m_nVectorInversion      =   p->m_nVectorInversion;
  m_nVectorRepresentation =   p->m_nVectorRepresentation;
  m_bHeatScaleClearHigh   =   p->m_bHeatScaleClearHigh;
  m_bHeatScaleTruncate    =   p->m_bHeatScaleTruncate;
  m_bHeatScaleInvert      =   p->m_bHeatScaleInvert;
  m_bRememberFrameSettings = p->m_bRememberFrameSettings;
  m_dVectorLineWidth      =   p->m_dVectorLineWidth;
  m_dVectorNormThreshold  =   p->m_dVectorNormThreshold;
  m_colorBinary           =   p->m_colorBinary;

  SetLUTCTAB  ( p->mFreeSurferCTAB );

  blockSignals( false );

  this->OnColorMapChanged();
}

void LayerPropertyMRI::CopyWindowLevelSettings( const LayerPropertyMRI* p )
{
  blockSignals( true );
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
  mbClearBackground             =   p->mbClearBackground;
  mClearBackgroundValue   =   p->mClearBackgroundValue;
  mMinGenericThreshold    =   p->mMinGenericThreshold;
  mMaxGenericThreshold    =   p->mMaxGenericThreshold;
  m_bHeatScaleClearHigh   =   p->m_bHeatScaleClearHigh;
  m_bHeatScaleTruncate    =   p->m_bHeatScaleTruncate;
  m_bHeatScaleInvert      =   p->m_bHeatScaleInvert;

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
  map.remove("RememberFrameSettings");
  RestoreSettings(map);
}

void LayerPropertyMRI::RestoreSettings(const QVariantMap& map)
{
  LayerMRI* mri = qobject_cast<LayerMRI*>(parent());
//  m_bUsePercentile = (map["UsePercentile"].toDouble() > 0);
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

  if ( map.contains("RememberFrameSettings"))
  {
    //  m_bRememberFrameSettings = map["RememberFrameSettings"].toBool();
  }

  if (map.contains("AutoAdjustFrameLevel") && map["AutoAdjustFrameLevel"].toBool())
  {
    SetAutoAdjustFrameLevel(true);
  }

  if ( map.contains("FrameSettings") && mri->GetNumberOfFrames() > 1)
  {
    m_frameSettings = map["FrameSettings"].toMap();
  }

  if (map.contains("ClearBackground"))
  {
    mbClearBackground = map["ClearBackground"].toBool();
  }

  if (map.contains("ClearBackgroundValue"))
  {
    mClearBackgroundValue = map["ClearBackgroundValue"].toDouble();
  }

  if (m_bRememberFrameSettings && mri->GetNumberOfFrames() > 1)
  {
    QVariantMap map2 = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    QStringList keys = map.keys();
    foreach (QString key, keys)
      map2[key] = map[key];
    m_frameSettings[QString::number(m_nActiveFrame)] = map2;
  }

  if (m_bUsePercentile)
  {
    if ( map.contains("MinGrayscaleWindowInPercentile") )
    {
      mMinGrayscaleWindow = mri->GetHistoValueFromPercentile(map["MinGrayscaleWindowInPercentile"].toDouble());
    }

    if ( map.contains("MaxGrayscaleWindowInPercentile") )
    {
      mMaxGrayscaleWindow = mri->GetHistoValueFromPercentile(map["MaxGrayscaleWindowInPercentile"].toDouble());
    }

    if ( map.contains("HeatScaleMinThresholdInPercentile") )
    {
      mHeatScaleMinThreshold = mri->GetHistoValueFromPercentile(map["HeatScaleMinThresholdInPercentile"].toDouble());
    }

    if ( map.contains("HeatScaleMidThresholdInPercentile") )
    {
      mHeatScaleMidThreshold = mri->GetHistoValueFromPercentile(map["HeatScaleMidThresholdInPercentile"].toDouble());
    }

    if ( map.contains("HeatScaleMaxThresholdInPercentile") )
    {
      mHeatScaleMaxThreshold = mri->GetHistoValueFromPercentile(map["HeatScaleMaxThresholdInPercentile"].toDouble());
    }

    if ( map.contains("MinGenericThresholdInPercentile") )
    {
      mMinGenericThreshold = mri->GetHistoValueFromPercentile(map["MinGenericThresholdInPercentile"].toDouble());
    }

    if ( map.contains("MaxGenericThresholdInPercentile") )
    {
      mMaxGenericThreshold = mri->GetHistoValueFromPercentile(map["MaxGenericThresholdInPercentile"].toDouble());
    }
  }

  this->OnColorMapChanged();
  emit OpacityChanged( mOpacity );
  emit DisplayModeChanged();
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
  map["RememberFrameSettings"] = m_bRememberFrameSettings;
  map["FrameSettings"] = m_frameSettings;
  map["ClearBackground"] = mbClearBackground;
  map["ClearBackgroundValue"] = mClearBackgroundValue;
  map["UsePercentile"] = m_bUsePercentile;
  map["AutoAdjustFrameLevel"] = m_bAutoAdjustFrameLevel;
  return map;
}

QVariantMap LayerPropertyMRI::GetFullSettings()
{
  QVariantMap map = GetSettings();
  map["ColorMapType"] = mColorMapType;
  map["MinContourThreshold"] = mMinContourThreshold;
  map["MaxContourThreshold"] = mMaxContourThreshold;

  map["DisplayVector"] = m_bDisplayVector;
  map["VectorInversion"] = m_nVectorInversion;
  map["VectorRepresentation"] = m_nVectorRepresentation;
  map["VectorLengthScale"] = m_dVectorDisplayScale;
  map["VectorWidthScale"] = m_dVectorLineWidth;
  map["VectorNormalize"] = m_bNormalizeVector;

  map["DisplayTensor"] = m_bDisplayTensor;
  map["TensorInversion"] = m_nTensorInversion;
  map["TensorRepresentation"] = m_nTensorRepresentation;

  map["Opacity"] = this->mOpacity;

  return map;
}

void LayerPropertyMRI::RestoreFullSettings(const QVariantMap &map)
{
  if (map.contains("ColorMapType"))
    mColorMapType = map["ColorMapType"].toInt();

  if (map.contains("MinContourThreshold"))
    mMinContourThreshold = map["MinContourThreshold"].toDouble();

  if (map.contains("MaxContourThreshold"))
    mMaxContourThreshold = map["MaxContourThreshold"].toDouble();

  if (map.contains("DisplayVector"))
    m_bDisplayVector = map["DisplayVector"].toBool();

  if (map.contains("VectorInversion"))
    m_nVectorInversion = map["VectorInversion"].toInt();

  if (map.contains("VectorRepresentation"))
    m_nVectorRepresentation = map["VectorRepresentation"].toInt();

  if (map.contains("DisplayTensor"))
    m_bDisplayTensor = map["DisplayTensor"].toBool();

  if (map.contains("TensorInversion"))
    m_nTensorInversion = map["TensorInversion"].toInt();

  if (map.contains("TensorRepresentation"))
    m_nTensorRepresentation = map["TensorRepresentation"].toInt();

  if (map.contains("Opacity"))
    mOpacity = map["Opacity"].toDouble();

  if (map.contains("VectorLengthScale"))
    m_dVectorDisplayScale = map["VectorLengthScale"].toDouble();

  if (map.contains("VectorWidthScale"))
    m_dVectorLineWidth = map["VectorWidthScale"].toDouble();

  if (map.contains("VectorNormalize"))
    m_bNormalizeVector = map["VectorNormalize"].toBool();

  RestoreSettings(map);
}

QVariantMap LayerPropertyMRI::GetActiveSettings()
{
  QVariantMap map;
  LayerMRI* mri = qobject_cast<LayerMRI*>(parent());
  switch (mColorMapType)
  {
  case Grayscale:
    map["MinGrayscaleWindow"] = mMinGrayscaleWindow;
    map["MaxGrayscaleWindow"] = mMaxGrayscaleWindow;
    if (m_bUsePercentile)
    {
      map["MinGrayscaleWindowInPercentile"] = mri->GetHistoPercentileFromValue(mMinGrayscaleWindow);
      map["MaxGrayscaleWindowInPercentile"] = mri->GetHistoPercentileFromValue(mMaxGrayscaleWindow);
    }
    break;
  case Heat:
    map["HeatScaleMinThreshold"] = mHeatScaleMinThreshold;
    map["HeatScaleMidThreshold"] = mHeatScaleMidThreshold;
    map["HeatScaleMaxThreshold"] = mHeatScaleMaxThreshold;
    map["HeatScaleOffset"] = mHeatScaleOffset;
    if (m_bUsePercentile)
    {
      map["HeatScaleMinThresholdInPercentile"] = mri->GetHistoPercentileFromValue(mHeatScaleMinThreshold);
      map["HeatScaleMidThresholdInPercentile"] = mri->GetHistoPercentileFromValue(mHeatScaleMidThreshold);
      map["HeatScaleMaxThresholdInPercentile"] = mri->GetHistoPercentileFromValue(mHeatScaleMaxThreshold);
    }
    break;
  case Jet:
  case Turbo:
  case GEColor:
  case NIH:
  case PET:
    map["MinGenericThreshold"] = mMinGenericThreshold;
    map["MaxGenericThreshold"] = mMaxGenericThreshold;
    if (m_bUsePercentile)
    {
      map["MinGenericThresholdInPercentile"] = mri->GetHistoPercentileFromValue(mMinGenericThreshold);
      map["MaxGenericThresholdInPercentile"] = mri->GetHistoPercentileFromValue(mMaxGenericThreshold);
    }
    break;
  }
  if (this->GetShowAsContour())
  {
    map["MinContourThreshold"] = mMinContourThreshold;
    map["MaxContourThreshold"] = mMaxContourThreshold;
  }
  map["UsePercentile"] = (m_bUsePercentile?1.0:0.0);
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
  case Binary:
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
  //  if ( ct != mFreeSurferCTAB )
  {
    mFreeSurferCTAB = ct;
    SetSelectAllLabels();
    UpdateLUTTable();
    if ( ct )
    {
      emit ColorMapChanged();
    }
  }
}

bool LayerPropertyMRI::IsValueInColorTable(double nVal)
{
  int bValid = 0;
  if (mFreeSurferCTAB)
    CTABisEntryValid( mFreeSurferCTAB, (int)nVal, &bValid );

  return bValid;
}

void LayerPropertyMRI::UpdateLUTTable()
{
  if ( mFreeSurferCTAB )
  {
    mLUTTable->RemoveAllPoints();
    mLUTTable->AddRGBAPoint( 0, 0, 0, 0, 0 );
    int cEntries;
    CTABgetNumberOfTotalEntries( mFreeSurferCTAB, &cEntries );
    bool last_is_valid = true;
    for ( int nEntry = 1; nEntry < cEntries; nEntry++ )
    {
      int bValid;
      CTABisEntryValid( mFreeSurferCTAB, nEntry, &bValid );
      if ( bValid )
      {
        if (!last_is_valid)
          last_is_valid = true;

        float red, green, blue, alpha;
        CTABrgbaAtIndexf( mFreeSurferCTAB, nEntry, &red, &green, &blue, &alpha );
        mLUTTable->AddRGBAPoint( nEntry, red, green, blue, m_listVisibleLabels.contains(nEntry)?1:0 );
      }
      else if (last_is_valid)
      {
        mLUTTable->AddRGBAPoint( nEntry, 0, 0, 0, 0 );
        last_is_valid = false;
      }
    }
    mLUTTable->AddRGBAPoint( cEntries, 0, 0, 0, 0 );
    mLUTTable->Build();
  }
}

void LayerPropertyMRI::SetCustomColors(const QMap<int, QColor> &colors)
{
  mLUTTable->RemoveAllPoints();
  QList<int> idx = colors.keys();
  foreach (int id, idx)
  {
    mLUTTable->AddRGBAPoint(id, colors[id].redF(), colors[id].greenF(), colors[id].blueF(), colors[id].alphaF());
  }
  mLUTTable->Build();
}

void LayerPropertyMRI::OnColorMapChanged ()
{
  if ( !mSource )
  {
    return;
  }

  double MinGrayscaleWindow = mMinGrayscaleWindow;
  double MaxGrayscaleWindow = mMaxGrayscaleWindow;
  double HeatScaleMinThreshold = mHeatScaleMinThreshold;
  double HeatScaleMidThreshold = mHeatScaleMidThreshold;
  double HeatScaleMaxThreshold = mHeatScaleMaxThreshold;
  double HeatScaleOffset = mHeatScaleOffset;
  double MinGenericThreshold = mMinGenericThreshold;
  double MaxGenericThreshold = mMaxGenericThreshold;
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGrayscaleWindow"))
      MinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
    if (map.contains("MaxGrayscaleWindow"))
      MaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
    if (map.contains("HeatScaleMinThreshold"))
      HeatScaleMinThreshold = map["HeatScaleMinThreshold"].toDouble();
    if (map.contains("HeatScaleMidThreshold"))
      HeatScaleMidThreshold = map["HeatScaleMidThreshold"].toDouble();
    if (map.contains("HeatScaleMaxThreshold"))
      HeatScaleMaxThreshold = map["HeatScaleMaxThreshold"].toDouble();
    if (map.contains("HeatScaleOffset"))
      HeatScaleOffset = map["HeatScaleOffset"].toDouble();
    if (map.contains("MinGenericThreshold"))
      MinGenericThreshold = map["MinGenericThreshold"].toDouble();
    if (map.contains("MaxGenericThreshold"))
      MaxGenericThreshold = map["MaxGenericThreshold"].toDouble();
  }

  double tiny_fraction = qMin(1e-12, ( mMaxVoxelValue - mMinVoxelValue ) * 1e-12);
  bool has_mask = (((LayerMRI*)parent())->GetMaskLayer());
  switch ( mColorMapType )
  {
  case NoColorMap:
    break;

  case Grayscale:
    // Check the color map variables and update range sliders if
    // necessary.
    if ( MinGrayscaleWindow < mMinVisibleValue )
    {
      mMinVisibleValue = MinGrayscaleWindow;
    }
    if ( MaxGrayscaleWindow >= mMaxVisibleValue )
    {
      mMaxVisibleValue = MaxGrayscaleWindow+1;
    }

    // Build our lookup table.
    assert( mGrayScaleTable.GetPointer() );
    mGrayScaleTable->RemoveAllPoints();
    if ( mbClearBackground || has_mask )
    {
      if (mClearBackgroundValue < MinGrayscaleWindow)
      {
        mGrayScaleTable->AddRGBAPoint( mClearBackgroundValue, 0, 0, 0, 0 );
        mGrayScaleTable->AddRGBAPoint( mClearBackgroundValue + tiny_fraction, 0, 0, 0, 1 );
        mGrayScaleTable->AddRGBAPoint( MinGrayscaleWindow,    0, 0, 0, 1 );
      }
      else
      {
        mGrayScaleTable->AddRGBAPoint( mClearBackgroundValue, 0, 0, 0, 0 );
        double val = (mClearBackgroundValue-mMinGrayscaleWindow) /
            (mMaxGrayscaleWindow == mMinGrayscaleWindow ? 1:(mMaxGrayscaleWindow - mMinGrayscaleWindow));
        mGrayScaleTable->AddRGBAPoint( mClearBackgroundValue + tiny_fraction, val, val, val, 1 );
      }
    }
    else
    {
      mGrayScaleTable->AddRGBAPoint( mMinVisibleValue - tiny_fraction, 0, 0, 0, 0 );
      mGrayScaleTable->AddRGBAPoint( mMinVisibleValue,       0, 0, 0, 1 );
      mGrayScaleTable->AddRGBAPoint( MinGrayscaleWindow,
                                     0, 0, 0, 1 );
    }
    mGrayScaleTable->AddRGBAPoint( MaxGrayscaleWindow,    1, 1, 1, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue,       1, 1, 1, 1 );
//    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue + tiny_fraction, 1, 1, 1, 0 );
    mGrayScaleTable->Build();
    break;

  case Heat:
    mHeatScaleTable->RemoveAllPoints();
    if ( m_bHeatScaleTruncate && m_bHeatScaleInvert )
    {
      if ( m_bHeatScaleClearHigh || has_mask )
      {
        mHeatScaleTable->AddRGBAPoint( -HeatScaleMaxThreshold + HeatScaleOffset - tiny_fraction, 1, 1, 0, 0 );
      }
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMaxThreshold + HeatScaleOffset, 1, 1, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMidThreshold + HeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMinThreshold + HeatScaleOffset + ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 1, 0, 0, 0 );
//      mHeatScaleTable->AddRGBAPoint(  0 + HeatScaleOffset, 0, 0, 0, 0 );
    }
    else if ( m_bHeatScaleTruncate )
    {
//      mHeatScaleTable->AddRGBAPoint(  0 + HeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMinThreshold + HeatScaleOffset - ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMidThreshold + HeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset, 1, 1, 0, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
      }
    }
    else if ( m_bHeatScaleInvert )
    {
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMaxThreshold + HeatScaleOffset, 1, 1, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMidThreshold + HeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMinThreshold + HeatScaleOffset + ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  0 + mHeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMinThreshold + HeatScaleOffset - ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 0, 0, 1, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMidThreshold + HeatScaleOffset, 0, 0, 1, 1 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset, 0, 1, 1, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset + tiny_fraction, 0, 1, 1, 0 );
      }
    }
    else
    {
      if ( has_mask )
      {
        mHeatScaleTable->AddRGBAPoint( -HeatScaleMaxThreshold - qAbs(HeatScaleOffset) - tiny_fraction, 0, 0, 0, 0 );
      }
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMaxThreshold + HeatScaleOffset, 0, 1, 1, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMidThreshold + HeatScaleOffset, 0, 0, 1, 1 );
      mHeatScaleTable->AddRGBAPoint( -HeatScaleMinThreshold + HeatScaleOffset + ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 0, 0, 1, 0 );
      mHeatScaleTable->AddRGBAPoint(  0 + HeatScaleOffset, 0, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMinThreshold + HeatScaleOffset - ((HeatScaleMinThreshold == HeatScaleMidThreshold)?tiny_fraction:0), 1, 0, 0, 0 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMidThreshold + HeatScaleOffset, 1, 0, 0, 1 );
      mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset, 1, 1, 0, 1 );
      if ( m_bHeatScaleClearHigh )
      {
        mHeatScaleTable->AddRGBAPoint(  HeatScaleMaxThreshold + HeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
      }
    }
    mHeatScaleTable->Build();
    break;

  case Jet:
    mColorMapTable->RemoveAllPoints();
//    qDebug() << qMin( 0.0, MinGenericThreshold) - tiny_fraction;
    mColorMapTable->AddRGBAPoint( qMin( 0.0, MinGenericThreshold) - tiny_fraction, 0, 0, 0, 0 );
    mColorMapTable->AddRGBAPoint( MinGenericThreshold, 0, 0, 1, 1 );
    mColorMapTable->AddRGBAPoint( MinGenericThreshold + (MaxGenericThreshold - MinGenericThreshold) / 4, 0, 1, 1, 1 );
    mColorMapTable->AddRGBAPoint( (MinGenericThreshold + MaxGenericThreshold) / 2, 0, 1, 0, 1 );
    mColorMapTable->AddRGBAPoint( MaxGenericThreshold - (MaxGenericThreshold - MinGenericThreshold) / 4, 1, 1, 0, 1 );
    mColorMapTable->AddRGBAPoint( MaxGenericThreshold, 1, 0, 0, 1 );
    //  mColorMapTable->AddRGBAPoint( mMaxGenericThreshold + (mMaxGenericThreshold - mMinJGenericThreshold), 1, 0, 0, 1 );
    mColorMapTable->Build();
    break;

  case Turbo:
    BuildGenericLUT( stock_turbo_color );
    break;

  case GEColor:
    BuildGenericLUT( stock_ge_color );
    break;

  case NIH:
    BuildGenericLUT( stock_nih );
    break;

  case PET:
    BuildGenericLUT( stock_pet );
    break;

  case LUT:
  case Binary:
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
    // ignore the first black color (0, 0, 0)
    for ( int i = 1; i < 256; i++ )
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
  m_bDisplayRGB = false;
  emit DisplayModeChanged();
}

void LayerPropertyMRI::SetDisplayRGB(bool b)
{
  m_bDisplayRGB = b;
  m_bDisplayVector = false;
  emit DisplayModeChanged();
}

void LayerPropertyMRI::SetNormalizeVector(bool b)
{
  m_bNormalizeVector = b;
  emit DisplayModeChanged();
}

void LayerPropertyMRI::SetVectorDisplayScale(double val)
{
  m_dVectorDisplayScale = val;
  emit DisplayModeChanged();
}

void LayerPropertyMRI::SetVectorNormThreshold(double dVal)
{
  m_dVectorNormThreshold = dVal;
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
    else if (iType == Binary)
      SetBinaryColor(m_colorBinary);
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
  double MinGrayscaleWindow = mMinGrayscaleWindow;
  double MaxGrayscaleWindow = mMaxGrayscaleWindow;
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGrayscaleWindow"))
      MinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
    if (map.contains("MaxGrayscaleWindow"))
      MaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
  }
  return ( MaxGrayscaleWindow - MinGrayscaleWindow );
}

void LayerPropertyMRI::SetLevel ( double iLevel )
{
  double window = this->GetWindow();
  this->SetMinMaxGrayscaleWindow( iLevel - window/2.0,
                                  iLevel + window/2.0 );
}

double LayerPropertyMRI::GetLevel ()
{
  double MinGrayscaleWindow = mMinGrayscaleWindow;
  double MaxGrayscaleWindow = mMaxGrayscaleWindow;
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGrayscaleWindow"))
      MinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
    if (map.contains("MaxGrayscaleWindow"))
      MaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
  }
  return ((MaxGrayscaleWindow - MinGrayscaleWindow) / 2.0) +
      MinGrayscaleWindow;
}

void LayerPropertyMRI::SetWindowLevel( double iWindow, double iLevel )
{
  this->SetMinMaxGrayscaleWindow( iLevel - iWindow/2.0, iLevel + iWindow/2.0 );
}

void LayerPropertyMRI::SetMinMaxGrayscaleWindow ( double iMin, double iMax )
{
  double MinGrayscaleWindow = mMinGrayscaleWindow;
  double MaxGrayscaleWindow = mMaxGrayscaleWindow;
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGrayscaleWindow"))
      MinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
    if (map.contains("MaxGrayscaleWindow"))
      MaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
  }

  if ( MinGrayscaleWindow != iMin || MaxGrayscaleWindow != iMax )
  {
    if (m_bRememberFrameSettings)
    {
      QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
      map["MinGrayscaleWindow"] = iMin;
      map["MaxGrayscaleWindow"] = iMax;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMinGrayscaleWindow = iMin;
      mMaxGrayscaleWindow = iMax;
    }
    if (!mbClearBackground)
      mClearBackgroundValue = mMinGrayscaleWindow;
    this->OnColorMapChanged();
    emit WindowLevelChanged();
  }
}

void LayerPropertyMRI::SetMinGrayscaleWindow ( double iMin )
{
  double MinGrayscaleWindow = mMinGrayscaleWindow;
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGrayscaleWindow"))
      MinGrayscaleWindow = map["MinGrayscaleWindow"].toDouble();
  }

  if ( MinGrayscaleWindow != iMin )
  {
    if (m_bRememberFrameSettings)
    {
      QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
      map["MinGrayscaleWindow"] = iMin;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMinGrayscaleWindow = iMin;
    }
    if (!mbClearBackground)
      mClearBackgroundValue = mMinGrayscaleWindow;
    this->OnColorMapChanged();
    emit WindowLevelChanged();
  }
}

void LayerPropertyMRI::SetMaxGrayscaleWindow ( double iMax )
{
  double MaxGrayscaleWindow = mMaxGrayscaleWindow;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MaxGrayscaleWindow"))
      MaxGrayscaleWindow = map["MaxGrayscaleWindow"].toDouble();
  }

  if ( MaxGrayscaleWindow != iMax )
  {
    if (m_bRememberFrameSettings)
    {
      map["MaxGrayscaleWindow"] = iMax;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMaxGrayscaleWindow = iMax;
    }
    this->OnColorMapChanged();
    emit WindowLevelChanged();
  }
}

void LayerPropertyMRI::SetHeatScaleAutoMid(bool bAutoMid)
{
  if (bAutoMid != m_bHeatScaleAutoMid)
  {
    m_bHeatScaleAutoMid = bAutoMid;
//    if (!bAutoMid)
//      m_bHeatScaleSetMidToMin = false;
    if (bAutoMid)
    {
      if (m_bHeatScaleSetMidToMin)
        mHeatScaleMidThreshold = mHeatScaleMinThreshold;
      else
        mHeatScaleMidThreshold = (mHeatScaleMinThreshold + mHeatScaleMaxThreshold)/2;
    }
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetHeatScaleSetMidToMin(bool bMidToMin)
{
  m_bHeatScaleSetMidToMin = bMidToMin;
  if (bMidToMin)
    mHeatScaleMidThreshold = mHeatScaleMinThreshold;
  else if (m_bHeatScaleAutoMid)
    mHeatScaleMidThreshold = (mHeatScaleMinThreshold + mHeatScaleMaxThreshold)/2;
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetHeatScaleMinThreshold( double iValue)
{
  double HeatScaleMinThreshold = mHeatScaleMinThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMinThreshold"))
      HeatScaleMinThreshold = map["HeatScaleMinThreshold"].toDouble();
  }
  if ( HeatScaleMinThreshold != iValue )
  {
    if (m_bRememberFrameSettings)
    {
      map["HeatScaleMinThreshold"] = iValue;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mHeatScaleMinThreshold = iValue;
    }
  }
  if (m_bHeatScaleSetMidToMin)
      mHeatScaleMidThreshold = mHeatScaleMinThreshold;
  else if (m_bHeatScaleAutoMid)
      mHeatScaleMidThreshold = (mHeatScaleMinThreshold + mHeatScaleMaxThreshold)/2;
  this->OnColorMapChanged();
}

double LayerPropertyMRI::GetHeatScaleMinThreshold ()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMinThreshold"))
      return map["HeatScaleMinThreshold"].toDouble();
  }
  return mHeatScaleMinThreshold;
}

void LayerPropertyMRI::SetHeatScaleMidThreshold ( double iValue )
{
  double HeatScaleMidThreshold = mHeatScaleMidThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMidThreshold"))
      HeatScaleMidThreshold = map["HeatScaleMidThreshold"].toDouble();
  }
  if ( HeatScaleMidThreshold != iValue )
  {
    if (m_bRememberFrameSettings)
    {
      map["HeatScaleMidThreshold"] = iValue;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mHeatScaleMidThreshold = iValue;
    }
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleMidThreshold ()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMidThreshold"))
      return map["HeatScaleMidThreshold"].toDouble();
  }
  return mHeatScaleMidThreshold;
}

void LayerPropertyMRI::SetHeatScaleMaxThreshold ( double iValue )
{
  double HeatScaleMaxThreshold = mHeatScaleMaxThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMaxThreshold"))
      HeatScaleMaxThreshold = map["HeatScaleMaxThreshold"].toDouble();
  }
  if ( HeatScaleMaxThreshold != iValue )
  {
    if (m_bRememberFrameSettings)
    {
      map["HeatScaleMaxThreshold"] = iValue;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mHeatScaleMaxThreshold = iValue;
    }
    if (m_bHeatScaleAutoMid && !m_bHeatScaleSetMidToMin)
      mHeatScaleMidThreshold = (mHeatScaleMinThreshold + mHeatScaleMaxThreshold)/2;
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleMaxThreshold ()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleMaxThreshold"))
      return map["HeatScaleMaxThreshold"].toDouble();
  }
  return mHeatScaleMaxThreshold;
}

void LayerPropertyMRI::SetHeatScale( double dMin, double dMid, double dMax )
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    map["HeatScaleMinThreshold"] = dMin;
    map["HeatScaleMidThreshold"] = dMid;
    map["HeatScaleMaxThreshold"] = dMax;
    m_frameSettings[QString::number(m_nActiveFrame)] = map;
  }
  else
  {
    mHeatScaleMinThreshold = dMin;
    mHeatScaleMidThreshold = dMid;
    mHeatScaleMaxThreshold = dMax;
  }
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetHeatScaleOffset ( double iValue )
{
  double HeatScaleOffset = mHeatScaleOffset;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleOffset"))
      HeatScaleOffset = map["HeatScaleOffset"].toDouble();
  }
  if ( HeatScaleOffset != iValue )
  {
    if (m_bRememberFrameSettings)
    {
      map["HeatScaleOffset"] = iValue;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mHeatScaleOffset = iValue;
    }
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetHeatScaleOffset ()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("HeatScaleOffset"))
      return map["HeatScaleOffset"].toDouble();
  }
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

bool LayerPropertyMRI::GetClearBackground()
{
  return mbClearBackground;
}

void LayerPropertyMRI::SetClearBackground( bool bClear )
{
  if ( mbClearBackground != bClear )
  {
    mbClearBackground = bClear;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetClearBackgroundValue(double val)
{
  if (mClearBackgroundValue != val)
  {
    mClearBackgroundValue = val;
    OnColorMapChanged();
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

  double voxel_size[3];
  source->GetImageOutput()->GetSpacing(voxel_size);
  m_dVectorScale = qMin( qMin( voxel_size[0], voxel_size[1] ), voxel_size[2] ) / 2.0;

  // Init our color scale values.
  UpdateMinMaxValues();

//  double dscale = qMax(fabs(mMaxVoxelValue), fabs(mMinVoxelValue));
//  if (dscale > 0)
//    m_dVectorDisplayScale = 1/dscale;

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

  UpdateLUTTable();
  if ( source->GetEmbeddedColorTable() )
  {
    SetColorMap( LUT );
    SetLUTCTAB( source->GetEmbeddedColorTable() );
  }

  MRI* mri = source->GetMRI();
  m_mapMeanValues[0] = MRImeanFrame(mri, 0);

  // Set up our initial tables.
  this->OnColorMapChanged();
}

void LayerPropertyMRI::UpdateMinMaxValues()
{
  if (!mSource)
    return;

//  if (m_bAutoAdjustFrameLevel)
//  {
//    if (m_mapMinMaxValues.contains(m_nActiveFrame))
//    {
//      QPair<double, double> pair = m_mapMinMaxValues[m_nActiveFrame];
//      mMinVoxelValue = pair.first;
//      mMaxVoxelValue = pair.second;
//    }
//    else
//    {
//      double range[2];
//      mSource->GetFrameValueRange(m_nActiveFrame, range);
//      mMinVoxelValue = range[0];
//      mMaxVoxelValue = range[1];
//    }
//  }
//  else
  {
    mMinVoxelValue = mSource->GetMinValue();
    mMaxVoxelValue = mSource->GetMaxValue();
  }

  mMinGrayscaleWindow = mMinVoxelValue;
  mMaxGrayscaleWindow = mMaxVoxelValue;
  if (mSource->HasValidHistogram())
  {
    if (m_mapMinMaxValues.contains(m_nActiveFrame))
    {
      QPair<double, double> pair = m_mapMinMaxValues[m_nActiveFrame];
      mMinGrayscaleWindow = pair.first;
      mMaxGrayscaleWindow = pair.second;
    }
    else
    {
      mMinGrayscaleWindow = mSource->GetHistoValueFromPercentile(0.05, m_nActiveFrame);
      mMaxGrayscaleWindow = mSource->GetHistoValueFromPercentile(0.99, m_nActiveFrame);
      mMaxGrayscaleWindow = mMaxGrayscaleWindow+(mMaxGrayscaleWindow-mMinGrayscaleWindow)/5;
      if (mMaxGrayscaleWindow == mMinGrayscaleWindow)
      {
        mMinGrayscaleWindow = mMinVoxelValue;
        mMaxGrayscaleWindow = mMaxVoxelValue;
      }
      QPair<double, double> pair;
      pair.first = mMinGrayscaleWindow;
      pair.second = mMaxGrayscaleWindow;
      m_mapMinMaxValues[m_nActiveFrame] = pair;
    }
  }

  mMinVisibleValue = mMinVoxelValue - ( mMaxVoxelValue - mMinVoxelValue );
  mMaxVisibleValue = mMaxVoxelValue + ( mMaxVoxelValue - mMinVoxelValue );
  mWindowRange[0] = 0;
  mWindowRange[1] = (mMaxVoxelValue - mMinVoxelValue) * 2;
  mLevelRange[0] = mMinVoxelValue;
  mLevelRange[1] = mMaxVoxelValue;

  double oneTenth;
  double highestAbsValue;
  highestAbsValue = qMax( fabs(mMinVoxelValue), fabs(mMaxVoxelValue) );
  if (mSource->HasValidHistogram())
  {
    highestAbsValue = mMaxGrayscaleWindow;
  }
  oneTenth = highestAbsValue / 10.0;
  mHeatScaleMinThreshold = mMinGrayscaleWindow;
  if (mHeatScaleMinThreshold < 0)
    mHeatScaleMinThreshold = 0;
  mHeatScaleMidThreshold = highestAbsValue / 2.0;
  mHeatScaleMaxThreshold = highestAbsValue - oneTenth;
  mHeatScaleOffset = 0;

  // Init the colors in the heat scale. The editor will handle the
  // editing from here.

  mMinGenericThreshold = mMinGrayscaleWindow;
  mMaxGenericThreshold = mMaxGrayscaleWindow;
  mMinContourThreshold = mHeatScaleMidThreshold;
  mMaxContourThreshold = mMaxGrayscaleWindow;
  mClearBackgroundValue = mMinGrayscaleWindow;
}

void LayerPropertyMRI::ResetWindowLevel()
{
  UpdateMinMaxValues();
  this->OnColorMapChanged();
}

void LayerPropertyMRI::SetMinMaxGenericThreshold ( double iMin, double iMax )
{
  double MinGenericThreshold = mMinGenericThreshold;
  double MaxGenericThreshold = mMaxGenericThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGenericThreshold"))
      MinGenericThreshold = map["MinGenericThreshold"].toDouble();
    if (map.contains("MaxGenericThreshold"))
      MaxGenericThreshold = map["MaxGenericThreshold"].toDouble();
  }

  if ( MinGenericThreshold != iMin || MaxGenericThreshold != iMax )
  {
    if (m_bRememberFrameSettings)
    {
      map["MinGenericThreshold"] = iMin;
      map["MaxGenericThreshold"] = iMax;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMinGenericThreshold = iMin;
      mMaxGenericThreshold = iMax;
    }
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetMinGenericThreshold()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGenericThreshold"))
      return map["MinGenericThreshold"].toDouble();
  }
  return mMinGenericThreshold;
}

void LayerPropertyMRI::SetMinGenericThreshold ( double iMin )
{
  double MinGenericThreshold = mMinGenericThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MinGenericThreshold"))
      MinGenericThreshold = map["MinGenericThreshold"].toDouble();
  }
  if ( MinGenericThreshold != iMin )
  {
    if (m_bRememberFrameSettings)
    {
      map["MinGenericThreshold"] = iMin;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMinGenericThreshold = iMin;
    }
    this->OnColorMapChanged();
  }
}

double LayerPropertyMRI::GetMaxGenericThreshold()
{
  if (m_bRememberFrameSettings)
  {
    QVariantMap map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MaxGenericThreshold"))
      return map["MaxGenericThreshold"].toDouble();
  }
  return mMaxGenericThreshold;
}

void LayerPropertyMRI::SetMaxGenericThreshold ( double iMax )
{
  double MaxGenericThreshold = mMaxGenericThreshold;
  QVariantMap map;
  if (m_bRememberFrameSettings)
  {
    map = m_frameSettings[QString::number(m_nActiveFrame)].toMap();
    if (map.contains("MaxGenericThreshold"))
      MaxGenericThreshold = map["MaxGenericThreshold"].toDouble();
  }
  if ( MaxGenericThreshold != iMax )
  {
    if (m_bRememberFrameSettings)
    {
      map["MaxGenericThreshold"] = iMax;
      m_frameSettings[QString::number(m_nActiveFrame)] = map;
    }
    else
    {
      mMaxGenericThreshold = iMax;
    }
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetShowAsContour( bool bContour )
{
  if ( mbShowAsContour != bContour )
  {
    mbShowAsContour = bContour;
    if (bContour && GetColorMap() == LUT)
    {
      m_bShowAsLabelContour = true;
    }
    emit ContourShown( mbShowAsContour );
  }
}

void LayerPropertyMRI::SetShowAsLabelContour(bool bLabelContour)
{
  if (m_bShowAsLabelContour != bLabelContour)
  {
    m_bShowAsLabelContour = bLabelContour;
    emit ContourChanged();
  }
}

void LayerPropertyMRI::SetShowVoxelizedContour(bool bVoxelize)
{
  if (m_bShowVoxelizedContour != bVoxelize)
  {
    m_bShowVoxelizedContour = bVoxelize;
    emit ContourNeedsRebuild();
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

void LayerPropertyMRI::SetContourShow2D(bool bShow)
{
  m_bShow2DContour = bShow;
  emit Contour2DShown(bShow);
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

void LayerPropertyMRI::SetContourDilateFirst(bool bDilate)
{
  if (m_bContourDilateFirst != bDilate)
  {
    m_bContourDilateFirst = bDilate;
    emit ContourNeedsRebuild();
  }
}

void LayerPropertyMRI::SetContourUpsample(bool bFlag)
{
  if (m_bContourUpsample != bFlag)
  {
    m_bContourUpsample = bFlag;
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

void LayerPropertyMRI::SetProjectionMapType(int nType)
{
  if (nType != m_nProjectionMapType)
  {
    m_nProjectionMapType = nType;
    emit ProjectionMapChanged();
  }
}

void LayerPropertyMRI::SetProjectionMapRange(int n, int start, int end)
{
  m_nProjectionMapRange[n*2] = start;
  m_nProjectionMapRange[n*2+1] = end;
  emit ProjectionMapChanged();
}

void LayerPropertyMRI::SetRememberFrameSettings(bool bFlag)
{
  if (bFlag != m_bRememberFrameSettings)
  {
    m_bRememberFrameSettings = bFlag;
    this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::UpdateActiveFrame(int nFrame)
{
  if (nFrame != m_nActiveFrame)
  {
    m_nActiveFrame = nFrame;
    if (m_bAutoAdjustFrameLevel)
      UpdateMinMaxValues();
    if (m_bRememberFrameSettings || m_bAutoAdjustFrameLevel)
      this->OnColorMapChanged();
  }
}

void LayerPropertyMRI::SetVectorScale(double dval)
{
  if (m_dVectorScale != dval)
  {
    m_dVectorScale = dval;
  }
}

void LayerPropertyMRI::SetVectorLineWidth(double val)
{
  if (m_dVectorLineWidth != val)
  {
    m_dVectorLineWidth = val;
    emit VectorLineWidthChanged(val);
  }
}

void LayerPropertyMRI::SetAutoAdjustFrameLevel(bool b)
{
  if (m_bAutoAdjustFrameLevel != b)
  {
    m_bAutoAdjustFrameLevel = b;
    if (b)
      UpdateMinMaxValues();
    this->OnColorMapChanged();
    emit AutoAdjustFrameContrastChanged(b);
  }
}

void LayerPropertyMRI::SetSelectAllLabels()
{
  m_listVisibleLabels.clear();
  if (mFreeSurferCTAB)
  {
    int cEntries;
    CTABgetNumberOfTotalEntries( mFreeSurferCTAB, &cEntries );
    for ( int nEntry = 0; nEntry < cEntries; nEntry++ )
    {
      int bValid;
      CTABisEntryValid( mFreeSurferCTAB, nEntry, &bValid );
      if ( bValid )
      {
        m_listVisibleLabels << nEntry;
      }
    }
  }
  UpdateLUTTable();
  this->OnColorMapChanged();
  emit LabelContourChanged();
}

void LayerPropertyMRI::SetUnselectAllLabels()
{
  m_listVisibleLabels.clear();
  UpdateLUTTable();
  this->OnColorMapChanged();
  emit LabelContourChanged();
}

void LayerPropertyMRI::SetSelectLabel(int nVal, bool bSelected)
{
  if (bSelected && !m_listVisibleLabels.contains(nVal))
    m_listVisibleLabels << nVal;
  else if (!bSelected)
    m_listVisibleLabels.removeOne(nVal);
  UpdateLUTTable();
//  this->OnColorMapChanged();
  emit LabelContourChanged(nVal);
}

void LayerPropertyMRI::SetVectorSkip(int nSkip)
{
  if (nSkip != m_nVectorSkip)
  {
    m_nVectorSkip = nSkip;
    emit VectorSkipChanged(nSkip);
  }
}

void LayerPropertyMRI::SetBinaryColor(const QColor &color)
{
  m_colorBinary = color;
  if (GetColorMap() == Binary)
  {
    QMap<int, QColor> colors;
    colors[0] = QColor(0,0,0,0);
    colors[1] = color;
    SetCustomColors(colors);
  }
}

double LayerPropertyMRI::GetFrameMeanValue(int frame)
{
  LayerMRI* layer = qobject_cast<LayerMRI*>(parent());
  if (frame < 0)
    frame = layer->GetActiveFrame();

  if (m_mapMeanValues.contains(frame))
    return m_mapMeanValues[frame];
  else
  {
    MRI* mri = layer->GetSourceVolume()->GetMRI();
    double val = MRImeanFrame(mri, frame);
    m_mapMeanValues[frame] = val;
    return val;
  }
}
