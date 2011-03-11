/**
 * @file  LayerPropertiesMRI.cxx
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


#include <assert.h>
#include "LayerPropertiesMRI.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkFSVolumeSource.h"
#include "vtkImageReslice.h"
#include "vtkPointData.h"
#include "FSVolume.h"
#include "StockColorMap.h"
#include <wx/filename.h>
#include <wx/config.h>

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

using namespace std;

LayerPropertiesMRI::LayerPropertiesMRI () : LayerProperties( ),
    mColorMapType( LayerPropertiesMRI::Grayscale ),
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
    m_nContourSmoothIterations( 0 ),
    m_bShowLabelOutline( false ),
    m_nUpSampleMethod( UM_None ),
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
}

LayerPropertiesMRI::~LayerPropertiesMRI ()
{}

void LayerPropertiesMRI::CopySettings( const LayerPropertiesMRI* p )
{
  BlockBroadcast( true );
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
  m_nContourSmoothIterations = p->m_nContourSmoothIterations;
  m_bContourExtractAll    = p->m_bContourExtractAll;

  SetLUTCTAB  ( p->mFreeSurferCTAB );

  BlockBroadcast( false );

  this->ColorMapChanged();
  this->SendBroadcast( "TextureSmoothingChanged", NULL );
  this->SendBroadcast( "ResliceInterpolationChanged", NULL );
  this->SendBroadcast( "WindowLevelChanged", NULL );
  this->SendBroadcast( "OpacityChanged", this );
}


void LayerPropertiesMRI::RestoreSettings( const char* filename )
{
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int n = 0;
    wxString groupname = wxString::Format( _T("/VolumeProperties/Volume%d"), n );
    wxString fn;
    wxFileName wxfn( filename );
    wxfn.Normalize();
    while ( ( config->Read( groupname + _T( "/Filename" ), &fn ) ) )
    {
      if ( wxfn.SameAs( fn ) )
      {
        double value;
        if ( config->Read( groupname + _T( "/MinGrayscaleWindow" ), &value ) )
          mMinGrayscaleWindow = value;
        
        if ( config->Read( groupname + _T( "/MaxGrayscaleWindow" ), &value ) )
          mMaxGrayscaleWindow = value;
        
        if ( config->Read( groupname + _T( "/HeatScaleMinThreshold" ), &value ) )
          mHeatScaleMinThreshold = value;
        
        if ( config->Read( groupname + _T( "/HeatScaleMidThreshold" ), &value ) )
          mHeatScaleMidThreshold = value;
        
        if ( config->Read( groupname + _T( "/HeatScaleMaxThreshold" ), &value ) )
          mHeatScaleMaxThreshold = value;
        
        if ( config->Read( groupname + _T( "/HeatScaleOffset" ), &value ) )
          mHeatScaleOffset = value;
        
        if ( config->Read( groupname + _T( "/MinGenericThreshold" ), &value ) )
          mMinGenericThreshold = value;
        
        if ( config->Read( groupname + _T( "/MaxGenericThreshold" ), &value ) )
          mMaxGenericThreshold = value;
        
        if ( config->Read( groupname + _T( "/MinContourThreshold" ), &value ) )
          mMinContourThreshold = value;
        
        if ( config->Read( groupname + _T( "/MaxContourThreshold" ), &value ) )
          mMaxContourThreshold = value;
        
        this->ColorMapChanged();
        
        break;
      }
      else
      {
        n++;
        groupname = wxString::Format( _T("/VolumeProperties/Volume%d"), n );
      }
    } 
  }
}
    
void LayerPropertiesMRI::SaveSettings( const char* filename )
{
  wxConfigBase* config = wxConfigBase::Get();
  if ( config )
  {
    int n = 0;
    wxString groupname = wxString::Format( _T("/VolumeProperties/Volume%d"), n );
    wxString fn;
    wxFileName wxfn( filename );
    wxfn.Normalize();
    while ( ( config->Read( groupname + _T( "/Filename" ), &fn ) ) )
    {
      if ( wxfn.SameAs( fn ) )
      {
        break;
      }
      else
      {
        n++;
        groupname = wxString::Format( _T("/VolumeProperties/Volume%d"), n );
      }
    }
    
    fn = filename;
    config->Write( groupname + _T( "/Filename" ),               wxfn.GetFullPath() );
    config->Write( groupname + _T( "/MinGrayscaleWindow" ),     mMinGrayscaleWindow );
    config->Write( groupname + _T( "/MaxGrayscaleWindow" ),     mMaxGrayscaleWindow );
    config->Write( groupname + _T( "/HeatScaleMinThreshold" ),  mHeatScaleMinThreshold );
    config->Write( groupname + _T( "/HeatScaleMidThreshold" ),  mHeatScaleMidThreshold );
    config->Write( groupname + _T( "/HeatScaleMaxThreshold" ),  mHeatScaleMaxThreshold );
    config->Write( groupname + _T( "/HeatScaleOffset" ),        mHeatScaleOffset );
    config->Write( groupname + _T( "/MinGenericThreshold" ),    mMinGenericThreshold );
    config->Write( groupname + _T( "/MaxGenericThreshold" ),    mMaxGenericThreshold );
    config->Write( groupname + _T( "/MinContourThreshold" ),    mMinContourThreshold );
    config->Write( groupname + _T( "/MaxContourThreshold" ),    mMaxContourThreshold );
    config->Flush();
  }
}

vtkRGBAColorTransferFunction* LayerPropertiesMRI::GetLUTTable () const
{
  return mLUTTable;
}

vtkRGBAColorTransferFunction* LayerPropertiesMRI::GetGrayScaleTable () const
{
  return mGrayScaleTable;
}

vtkRGBAColorTransferFunction* LayerPropertiesMRI::GetHeatScaleTable () const
{
  return mHeatScaleTable;
}

vtkRGBAColorTransferFunction* LayerPropertiesMRI::GetColorMapTable () const
{
  return mColorMapTable;
}

vtkLookupTable* LayerPropertiesMRI::GetDirectionCodedTable () const
{
  return mDirectionCodedTable;
}

vtkScalarsToColors* LayerPropertiesMRI::GetActiveLookupTable()
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

COLOR_TABLE* LayerPropertiesMRI::GetLUTCTAB () const
{
  return mFreeSurferCTAB;
}

void LayerPropertiesMRI::SetLUTCTAB( COLOR_TABLE* ct )
{
  if ( ct != mFreeSurferCTAB )
  {
    mFreeSurferCTAB = ct;
    UpdateLUTTable();
  }
}

void LayerPropertiesMRI::UpdateLUTTable()
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

void LayerPropertiesMRI::ColorMapChanged ()
{
  if ( !mSource )
    return;

  double tiny_fraction = ( mMaxVoxelValue - mMinVoxelValue ) * 1e-12; 
  switch ( mColorMapType )
  {
  case NoColorMap:
    break;

  case Grayscale:

    // Check the color map variables and update range sliders if
    // necessary.
    if ( mMinGrayscaleWindow < mMinVisibleValue )
      mMinVisibleValue = mMinGrayscaleWindow;
    if ( mMaxGrayscaleWindow > mMaxVisibleValue )
      mMaxVisibleValue = mMaxGrayscaleWindow;

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
        mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold + mHeatScaleOffset - tiny_fraction, 1, 1, 0, 0 );
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
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
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
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 0, 1, 1, 0 );
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
        mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold + mHeatScaleOffset + tiny_fraction, 1, 1, 0, 0 );
    }
    mHeatScaleTable->Build();
    break;   
       
  case Jet:
    mColorMapTable->RemoveAllPoints();
    mColorMapTable->AddRGBAPoint( min( 0, mMinGenericThreshold), 0, 0, 0, 0 );
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
  this->SendBroadcast( "ColorMapChanged", NULL );
}

void LayerPropertiesMRI::BuildGenericLUT( const int colors[256][3] )
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

void LayerPropertiesMRI::SetDisplayVector( bool b )
{
  m_bDisplayVector = b;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}

void LayerPropertiesMRI::SetVectorInversion( int nInversion )
{
  m_nVectorInversion = nInversion;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}

void LayerPropertiesMRI::SetVectorRepresentation( int n )
{
  m_nVectorRepresentation = n;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}


void LayerPropertiesMRI::SetDisplayTensor( bool b )
{
  m_bDisplayTensor = b;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}

void LayerPropertiesMRI::SetTensorInversion( int nInversion )
{
  m_nTensorInversion = nInversion;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}

void LayerPropertiesMRI::SetTensorRepresentation( int n )
{
  m_nTensorRepresentation = n;
  this->SendBroadcast( "DisplayModeChanged", NULL );
}

void LayerPropertiesMRI::SetColorMap ( int iType )
{
  if ( mColorMapType != iType )
  {
    mColorMapType = iType;
    if ( iType == LUT )
      SetUpSampleMethod( UM_None );
    this->ColorMapChanged();
  }
}

int LayerPropertiesMRI::GetColorMap () const
{
  return mColorMapType;
}

void LayerPropertiesMRI::SetColorMapToGrayScale ()
{
  SetColorMap( Grayscale );
}

void LayerPropertiesMRI::SetColorMapToHeatScale ()
{
  SetColorMap( Heat );
}

void LayerPropertiesMRI::SetColorMapToLUT ()
{
  SetColorMap( LUT );
}

void LayerPropertiesMRI::SetResliceInterpolation ( int iMode )
{
  if ( mResliceInterpolation != iMode )
  {
    mResliceInterpolation = iMode;
    this->SendBroadcast( "ResliceInterpolationChanged", NULL );
  }
}

int LayerPropertiesMRI::GetResliceInterpolation () const
{
  return mResliceInterpolation;
}

void LayerPropertiesMRI::SetTextureSmoothing ( int iSmooth )
{
  if ( mTextureSmoothing != iSmooth )
  {
    mTextureSmoothing = iSmooth;
    this->SendBroadcast( "TextureSmoothingChanged", NULL );
  }
}

int LayerPropertiesMRI::GetTextureSmoothing () const
{
  return mTextureSmoothing;
}

void LayerPropertiesMRI::SetMinVisibleValue ( double iValue )
{
  if ( mMinVisibleValue != iValue )
  {
    mMinVisibleValue = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetMinVisibleValue ()
{
  return mMinVisibleValue;
}

void LayerPropertiesMRI::SetMaxVisibleValue ( double iValue )
{
  if ( mMaxVisibleValue != iValue )
  {
    mMaxVisibleValue = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetMaxVisibleValue ()
{
  return mMaxVisibleValue;
}

void LayerPropertiesMRI::SetMinMaxVisibleValue ( double iMinValue,
    double iMaxValue )
{
  if ( mMinVisibleValue != iMinValue ||
       mMaxVisibleValue != iMaxValue )
  {
    mMinVisibleValue = iMinValue;
    mMaxVisibleValue = iMaxValue;
    this->ColorMapChanged();
  }
}

void LayerPropertiesMRI::SetWindow ( double iWindow )
{
  double level = this->GetLevel();
  this->SetMinMaxGrayscaleWindow( level - iWindow/2.0,
                                  level + iWindow/2.0 );
}

double LayerPropertiesMRI::GetWindow ()
{
  return ( mMaxGrayscaleWindow - mMinGrayscaleWindow );
}

void LayerPropertiesMRI::SetLevel ( double iLevel )
{
  double window = this->GetWindow();
  this->SetMinMaxGrayscaleWindow( iLevel - window/2.0,
                                  iLevel + window/2.0 );
}

double LayerPropertiesMRI::GetLevel ()
{
  return ((mMaxGrayscaleWindow - mMinGrayscaleWindow) / 2.0) +
         mMinGrayscaleWindow;
}

void LayerPropertiesMRI::SetWindowLevel( double iWindow, double iLevel )
{
  this->SetMinMaxGrayscaleWindow( iLevel - iWindow/2.0, iLevel + iWindow/2.0 );
}

void LayerPropertiesMRI::SetMinMaxGrayscaleWindow ( double iMin, double iMax )
{
  if ( mMinGrayscaleWindow != iMin || mMaxGrayscaleWindow != iMax )
  {
    mMinGrayscaleWindow = iMin;
    mMaxGrayscaleWindow = iMax;
    this->ColorMapChanged();
    this->SendBroadcast( "WindowLevelChanged", NULL );
  }
}

void LayerPropertiesMRI::SetMinGrayscaleWindow ( double iMin )
{
  if ( mMinGrayscaleWindow != iMin )
  {
    mMinGrayscaleWindow = iMin;
    this->ColorMapChanged();
  }
}

void LayerPropertiesMRI::SetMaxGrayscaleWindow ( double iMax )
{
  if ( mMaxGrayscaleWindow != iMax )
  {
    mMaxGrayscaleWindow = iMax;
    this->ColorMapChanged();
  }
}

void LayerPropertiesMRI::SetHeatScaleMinThreshold ( double iValue )
{
  if ( mHeatScaleMinThreshold != iValue )
  {
    mHeatScaleMinThreshold = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetHeatScaleMinThreshold ()
{
  return mHeatScaleMinThreshold;
}

void LayerPropertiesMRI::SetHeatScaleMidThreshold ( double iValue )
{
  if ( mHeatScaleMidThreshold != iValue )
  {
    mHeatScaleMidThreshold = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetHeatScaleMidThreshold ()
{
  return mHeatScaleMidThreshold;
}

void LayerPropertiesMRI::SetHeatScaleMaxThreshold ( double iValue )
{
  if ( mHeatScaleMaxThreshold != iValue )
  {
    mHeatScaleMaxThreshold = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetHeatScaleMaxThreshold ()
{
  return mHeatScaleMaxThreshold;
}

void LayerPropertiesMRI::SetHeatScale( double dMin, double dMid, double dMax )
{
  mHeatScaleMinThreshold = dMin;
  mHeatScaleMidThreshold = dMid;
  mHeatScaleMaxThreshold = dMax;
  this->ColorMapChanged();
}

void LayerPropertiesMRI::SetHeatScaleOffset ( double iValue )
{
  if ( mHeatScaleOffset != iValue )
  {
    mHeatScaleOffset = iValue;
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetHeatScaleOffset ()
{
  return mHeatScaleOffset;
}

void LayerPropertiesMRI::SetReverseHeatScale ( bool ib )
{
  if ( mbReverseHeatScale != ib )
  {
    mbReverseHeatScale = ib;
    this->ColorMapChanged();
  }
}

bool LayerPropertiesMRI::GetReverseHeatScale ()
{
  return mbReverseHeatScale;
}

void LayerPropertiesMRI::SetHeatScaleClearHigh( bool bClear )
{
  m_bHeatScaleClearHigh = bClear;
  this->ColorMapChanged();
}

void LayerPropertiesMRI::SetHeatScaleTruncate( bool bTruncate )
{
  m_bHeatScaleTruncate = bTruncate;
  this->ColorMapChanged();
}

void LayerPropertiesMRI::SetHeatScaleInvert( bool bInvert )
{
  m_bHeatScaleInvert = bInvert;
  this->ColorMapChanged();
}

void LayerPropertiesMRI::SetShowPositiveHeatScaleValues ( bool ib )
{
  if ( mbShowPositiveHeatScaleValues != ib )
  {
    mbShowPositiveHeatScaleValues = ib;
    this->ColorMapChanged();
  }
}

bool LayerPropertiesMRI::GetShowPositiveHeatScaleValues ()
{
  return mbShowPositiveHeatScaleValues;
}

void LayerPropertiesMRI::SetShowNegativeHeatScaleValues ( bool ib )
{
  if ( mbShowNegativeHeatScaleValues != ib )
  {
    mbShowNegativeHeatScaleValues = ib;
    this->ColorMapChanged();
  }
}

bool LayerPropertiesMRI::GetShowNegativeHeatScaleValues ()
{
  return mbShowNegativeHeatScaleValues;
}

double LayerPropertiesMRI::GetOpacity() const
{
  return mOpacity;
}

void LayerPropertiesMRI::SetOpacity( double opacity )
{
  if ( mOpacity != opacity && opacity >= 0 && opacity <= 1 )
  {
    mOpacity = opacity;
    this->SendBroadcast( "OpacityChanged", this );
  }
}

bool LayerPropertiesMRI::GetClearZero()
{
  return mbClearZero;
}

void LayerPropertiesMRI::SetClearZero( bool bClear )
{
  if ( mbClearZero != bClear )
  {
    mbClearZero = bClear;
    // this->SendBroadcast( "ClearZeroChanged", this );
    this->ColorMapChanged();
  }
}

double LayerPropertiesMRI::GetMinValue()
{
  return mMinVoxelValue;
}

double LayerPropertiesMRI::GetMaxValue()
{
  return mMaxVoxelValue;
}


bool LayerPropertiesMRI::UpdateValueRange( double dValue )
{
  if ( mMinVoxelValue <= dValue && mMaxVoxelValue >= dValue )
    return false;
  
  if ( mMinVoxelValue > dValue )
    mMinVoxelValue = dValue;
  else if ( mMaxVoxelValue < dValue )
    mMaxVoxelValue = dValue;
  
  mMinVisibleValue = mMinVoxelValue - ( mMaxVoxelValue - mMinVoxelValue );
  mMaxVisibleValue = mMaxVoxelValue + ( mMaxVoxelValue - mMinVoxelValue );
  mWindowRange[0] = 0;
  mWindowRange[1] = (mMaxVoxelValue - mMinVoxelValue) * 2;
  mLevelRange[0] = mMinVoxelValue;
  mLevelRange[1] = mMaxVoxelValue;
  
//  this->ColorMapChanged();
  return true;
}

void LayerPropertiesMRI::SetVolumeSource ( FSVolume* source )
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
  highestAbsValue = max( fabs(mMinVoxelValue), fabs(mMaxVoxelValue) );
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
  this->ColorMapChanged();
}

void LayerPropertiesMRI::SetMinMaxGenericThreshold ( double iMin, double iMax )
{
  if ( mMinGenericThreshold != iMin || mMaxGenericThreshold != iMax )
  {
    mMinGenericThreshold = iMin;
    mMaxGenericThreshold = iMax;
    this->ColorMapChanged();
    this->SendBroadcast( "GenericWindowLevelChanged", NULL );
  }
}

double LayerPropertiesMRI::GetMinGenericThreshold()
{
  return mMinGenericThreshold;
}

void LayerPropertiesMRI::SetMinGenericThreshold ( double iMin )
{
  if ( mMinGenericThreshold != iMin )
  {
    mMinGenericThreshold = iMin;
    this->ColorMapChanged();
  }
}


double LayerPropertiesMRI::GetMaxGenericThreshold()
{
  return mMaxGenericThreshold;
}

void LayerPropertiesMRI::SetMaxGenericThreshold ( double iMax )
{
  if ( mMaxGenericThreshold != iMax )
  {
    mMaxGenericThreshold = iMax;
    this->ColorMapChanged();
  }
}

void LayerPropertiesMRI::SetShowAsContour( bool bContour )
{
  if ( mbShowAsContour != bContour )
  {
    mbShowAsContour = bContour;
    this->SendBroadcast( "LayerContourShown", NULL );
  }
}

void LayerPropertiesMRI::SetContourMinThreshold( double dValue )
{
  if ( mMinContourThreshold != dValue )
  {
    mMinContourThreshold = dValue;
    this->SendBroadcast( "LayerContourChanged", NULL );
  }
}

void LayerPropertiesMRI::SetContourMaxThreshold( double dValue )
{
  if ( mMaxContourThreshold != dValue )
  {
    mMaxContourThreshold = dValue;
    this->SendBroadcast( "LayerContourChanged", NULL );
  }
}

void LayerPropertiesMRI::SetContourThreshold( double dMin, double dMax )
{
  if ( mMinContourThreshold != dMin || mMaxContourThreshold != dMax )
  {
    mMinContourThreshold = dMin;
    mMaxContourThreshold = dMax;
    this->SendBroadcast( "LayerContourChanged", NULL );
  }
}

void LayerPropertiesMRI::SetContourExtractAllRegions( bool bExtractAll )
{
  if ( m_bContourExtractAll != bExtractAll )
  {
    m_bContourExtractAll = bExtractAll;
    this->SendBroadcast( "LayerContourChanged", NULL );
  }
}

void LayerPropertiesMRI::SetContourSmoothIterations( int nIterations )
{
  if ( m_nContourSmoothIterations != nIterations )
  {
    m_nContourSmoothIterations = nIterations;
    this->SendBroadcast( "LayerContourChanged", NULL );
  }
}

void LayerPropertiesMRI::SetContourColor( double r, double g, double b )
{
  m_rgbContour[0] = r;
  m_rgbContour[1] = g;
  m_rgbContour[2] = b;
  
  this->SendBroadcast( "LayerContourColorChanged", NULL );
}

void LayerPropertiesMRI::SetContourUseImageColorMap( bool bFlag )
{
  if ( m_bContourUseImageColorMap != bFlag )
  {
    m_bContourUseImageColorMap = bFlag;
    
    this->SendBroadcast( "LayerContourColorChanged", NULL );
  }
}

void LayerPropertiesMRI::SetShowLabelOutline( bool bOutline )
{
  if ( bOutline )
    SetUpSampleMethod( 0 );
  
  if ( m_bShowLabelOutline != bOutline )
  {
    m_bShowLabelOutline = bOutline;
    
    this->SendBroadcast( "LayerLabelOutlineChanged", NULL );
  }
}

void LayerPropertiesMRI::SetUpSampleMethod( int nSampleMethod )
{
  if ( nSampleMethod > 0 )
    this->SetShowLabelOutline( false );
  
  if ( nSampleMethod != m_nUpSampleMethod )
  {
    m_nUpSampleMethod = nSampleMethod;
    
    this->SendBroadcast( "LayerUpSampleMethodChanged", NULL );
  }
}
