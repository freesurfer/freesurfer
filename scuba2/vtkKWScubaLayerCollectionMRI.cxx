/**
 * @file  vtkKWScubaLayerCollectionMRI.cxx
 * @brief Implementation for MRI viewers.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.6 $
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
#include "vtkKWScubaLayerCollectionMRI.h"
#include "vtkKWScubaLayer2DMRI.h"
#include "vtkKWScubaLayer3DMRI.h"
#include "vtkObjectFactory.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkRGBATransferFunction.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRange.h"
#include "vtkKWFrame.h"
#include "vtkKWRGBATransferFunctionEditor.h"
#include "vtkKWLabel.h"
#include "vtkFSVolumeSource.h"
#include "vtkKWRadioButtonSet.h"
#include "vtkImageReslice.h"
#include "vtkPointData.h"
#include "vtkKWEntry.h"
#include "vtkKWPushButton.h"
#include "vtkKWLoadSaveDialog.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayerCollectionMRI );
vtkCxxRevisionMacro( vtkKWScubaLayerCollectionMRI, "$Revision: 1.6 $" );

vtkKWScubaLayerCollectionMRI::vtkKWScubaLayerCollectionMRI () :
  mColorMapType( ScubaCollectionPropertiesMRI::GrayScale ),
  mResliceInterpolation( 0 ),
  mTextureSmoothing( 0 ),
  mbClearZero( false ),
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
  mFreeSurferCTAB( NULL ),
  mfnVolume("") {
}

vtkKWScubaLayerCollectionMRI::~vtkKWScubaLayerCollectionMRI () {

}

void
vtkKWScubaLayerCollectionMRI::SetVolumeFileName ( const char* ifnVolume ) {

  // Set our file name and load the volume.
  mfnVolume = ifnVolume;
  this->LoadVolumeFromFileName();

  // Find a good label based on the filename and set it in the
  // layer.
  string fnVolume = ifnVolume;
  string::size_type lastSlash = fnVolume.rfind( "/" );
  this->SetLabel( fnVolume.substr( lastSlash+1, string::npos ).c_str() );

}

void
vtkKWScubaLayerCollectionMRI::AddControls ( vtkKWWidget* iPanel ) {

  // Interpolate radio buttons ------------------------------------------
  vtkSmartPointer<vtkKWRadioButtonSet> radBtnSetInterpolation = 
    vtkSmartPointer<vtkKWRadioButtonSet>::New();
  radBtnSetInterpolation->SetParent( iPanel );
  radBtnSetInterpolation->Create();
  radBtnSetInterpolation->PackHorizontallyOff();

  char sTclCmd[1024];
  vtkSmartPointer<vtkKWRadioButton> radBtnNearestNeighbor;
  radBtnNearestNeighbor.TakeReference(radBtnSetInterpolation->AddWidget( 0 ) );
  radBtnNearestNeighbor->SetText( "Nearest" );
  sprintf( sTclCmd, "SetResliceInterpolation %d", VTK_RESLICE_NEAREST );
  radBtnNearestNeighbor->SetCommand( this, sTclCmd );
  if ( mResliceInterpolation == VTK_RESLICE_NEAREST )
    radBtnNearestNeighbor->SelectedStateOn();

  vtkSmartPointer<vtkKWRadioButton> radBtnLinear;
  radBtnLinear.TakeReference( radBtnSetInterpolation->AddWidget( 1 ) );
  radBtnLinear->SetText( "Linear" );
  sprintf( sTclCmd, "SetResliceInterpolation %d", VTK_RESLICE_LINEAR );
  radBtnLinear->SetCommand( this, sTclCmd );
  if ( mResliceInterpolation == VTK_RESLICE_LINEAR )
    radBtnLinear->SelectedStateOn();

  vtkSmartPointer<vtkKWRadioButton> radBtnCubic;
  radBtnCubic.TakeReference( radBtnSetInterpolation->AddWidget( 2 ) );
  radBtnCubic->SetText( "Cubic" );
  sprintf( sTclCmd, "SetResliceInterpolation %d", VTK_RESLICE_CUBIC );
  radBtnCubic->SetCommand( this, sTclCmd );
  if ( mResliceInterpolation == VTK_RESLICE_CUBIC )
    radBtnCubic->SelectedStateOn();
  // --------------------------------------------------------------------

  // Smooth check button ------------------------------------------------
  vtkSmartPointer<vtkKWCheckButton> chkBtnSmooth = 
    vtkSmartPointer<vtkKWCheckButton>::New();
  chkBtnSmooth->SetParent( iPanel );
  chkBtnSmooth->Create();
  chkBtnSmooth->SetAnchorToWest();
  chkBtnSmooth->SetText( "Smooth" );
  chkBtnSmooth->SetCommand( this, "SetTextureSmoothing" );
  if ( mTextureSmoothing )
    chkBtnSmooth->SelectedStateOn();
  // --------------------------------------------------------------------

  // Color map radio buttons --------------------------------------------
  vtkSmartPointer<vtkKWRadioButtonSet> radBtnSetColorMap = 
    vtkSmartPointer<vtkKWRadioButtonSet>::New();
  radBtnSetColorMap->SetParent( iPanel );
  radBtnSetColorMap->Create();
  radBtnSetColorMap->PackHorizontallyOn();

  for ( int nType = (int)GrayScale; nType < (int)cColorMapTypes;  nType++ ) {

    maRadBtnColorMap[nType].
      TakeReference( radBtnSetColorMap->AddWidget( nType ) );

    switch ( nType ) {
    case GrayScale:
      maRadBtnColorMap[nType]->SetText( "Gray Scale" );
      maRadBtnColorMap[nType]->SetCommand( this, "SetColorMapToGrayScale" );
      break;
    case HeatScale:
      maRadBtnColorMap[nType]->SetText( "Heat Scale" );
      maRadBtnColorMap[nType]->SetCommand( this, "SetColorMapToHeatScale" );
      break;
    case LUT:
      maRadBtnColorMap[nType]->SetText( "LUT" );
      maRadBtnColorMap[nType]->SetCommand( this, "SetColorMapToLUT" );
      break;
    }

    if ( (int)mColorMapType == nType )
      maRadBtnColorMap[nType]->SelectedStateOn();
  }
  // --------------------------------------------------------------------

  // Color map specfic settings page -----------------------------------
  for ( int nType = (int)GrayScale; nType < (int)cColorMapTypes;  nType++ ) {

    maFrameColorMapSettings[nType] = 
      vtkSmartPointer<vtkKWFrame>::New();
    maFrameColorMapSettings[nType]->SetParent( iPanel );
    maFrameColorMapSettings[nType]->Create();
    vtkKWFrame* frame = maFrameColorMapSettings[nType];

    switch ( nType ) {
    case GrayScale:

      mRangeVisibleValue = vtkSmartPointer<vtkKWRange>::New();
      mRangeVisibleValue->SetParent( frame );
      mRangeVisibleValue->Create();
      mRangeVisibleValue->SetLabelText( "Visible" );
      mRangeVisibleValue->SetOrientationToHorizontal();
      mRangeVisibleValue->GetEntry1()->SetWidth( 6 );
      mRangeVisibleValue->GetEntry2()->SetWidth( 6 );
      mRangeVisibleValue->SetWholeRange( mSource->GetMinValue(),
                                         mSource->GetMaxValue() );
      mRangeVisibleValue->SetRange( mMinVisibleValue, mMaxVisibleValue );
      mRangeVisibleValue->
      SetResolution( mSource->GetPreferredValueIncrement() );
      mRangeVisibleValue->SetCommand( this, "SetMinMaxVisibleValue" );


      mRangeWindowLevel = vtkSmartPointer<vtkKWRange>::New();
      mRangeWindowLevel->SetParent( frame );
      mRangeWindowLevel->Create();
      mRangeWindowLevel->SetLabelText( "Grayscale Window" );
      mRangeWindowLevel->SetOrientationToHorizontal();
      mRangeWindowLevel->GetEntry1()->SetWidth( 6 );
      mRangeWindowLevel->GetEntry2()->SetWidth( 6 );
      mRangeWindowLevel->SymmetricalInteractionOn();
      mRangeWindowLevel->SetWholeRange( mSource->GetMinValue(),
                                        mSource->GetMaxValue() );
      mRangeWindowLevel->SetRange( mMinGrayscaleWindow, mMaxGrayscaleWindow );
      mRangeWindowLevel->
      SetResolution( mSource->GetPreferredValueIncrement() );
      mRangeWindowLevel->SetCommand( this, "SetMinMaxGrayscaleWindow" );

      this->Script( "pack %s %s -side top -fill x",
                    mRangeVisibleValue->GetWidgetName(),
                    mRangeWindowLevel->GetWidgetName() );

      break;
    case HeatScale: {

      vtkSmartPointer<vtkKWHistogram> histogram = 
	vtkSmartPointer<vtkKWHistogram>::New();
      histogram->
      BuildHistogram( mSource->GetOutput()->GetPointData()->GetScalars(), 0);

      mRGBAEditorHeatScale = 
	vtkSmartPointer<vtkKWRGBATransferFunctionEditor>::New();
      mRGBAEditorHeatScale->SetParent( frame );
      mRGBAEditorHeatScale->Create();

      mRGBAEditorHeatScale->SetBorderWidth( 2 );
      mRGBAEditorHeatScale->SetReliefToGroove();
      mRGBAEditorHeatScale->SetPadX( 2 );
      mRGBAEditorHeatScale->SetPadY( 2 );
      mRGBAEditorHeatScale->ExpandCanvasWidthOn();
      //      mRGBAEditorHeatScale->SetCanvasWidth( 450 );
      mRGBAEditorHeatScale->SetCanvasHeight( 150 );
      mRGBAEditorHeatScale->SetLabelText("Heat Scale Editor");
      mRGBAEditorHeatScale->SetRangeLabelPositionToTop();
      mRGBAEditorHeatScale->ColorSpaceOptionMenuVisibilityOff();
      mRGBAEditorHeatScale->ValueEntriesVisibilityOff();
      mRGBAEditorHeatScale->SetPointPositionInValueRangeToTop();
      mRGBAEditorHeatScale->SetPointStyleToCursorDown();
      mRGBAEditorHeatScale->FunctionLineVisibilityOff();
      mRGBAEditorHeatScale->PointGuidelineVisibilityOn();
      mRGBAEditorHeatScale->PointIndexVisibilityOff();
      mRGBAEditorHeatScale->SelectedPointIndexVisibilityOn();
      mRGBAEditorHeatScale->MidPointEntryVisibilityOff();
      mRGBAEditorHeatScale->SharpnessEntryVisibilityOff();
      mRGBAEditorHeatScale->SetLabelPositionToTop();
      mRGBAEditorHeatScale->ColorRampVisibilityOff();
      mRGBAEditorHeatScale->
	SetFunctionChangingCommand( this, "EditorChangedHeatScale" );
      mRGBAEditorHeatScale->
	SetFunctionChangedCommand( this, "EditorChangedHeatScale" );

      mRGBAEditorHeatScale->SetRGBATransferFunction( mHeatScaleTable );
      mRGBAEditorHeatScale->
      SetWholeParameterRange( -max(fabs(mSource->GetMinValue()),
                                   fabs(mSource->GetMaxValue())),
                              max(fabs(mSource->GetMinValue()),
                                  fabs(mSource->GetMaxValue())) );
      mRGBAEditorHeatScale->SetVisibleParameterRangeToWholeParameterRange();

      mRGBAEditorHeatScale->SetHistogram( histogram );

      mRGBAEditorHeatScale->ComputeHistogramColorFromValueOff();

      mRGBAEditorHeatScale->ParameterTicksVisibilityOn();
      mRGBAEditorHeatScale->ComputeValueTicksFromHistogramOn();
      mRGBAEditorHeatScale->SetParameterTicksFormat("%-#6.0f");

      this->Script( "pack %s -side top -fill x -expand y",
                    mRGBAEditorHeatScale->GetWidgetName() );

      // Set up the alpha values and symmetry.
      mRGBAEditorHeatScale->SetPointCountMinimum( 6 );
      mRGBAEditorHeatScale->SetPointCountMaximum( 6 );
      mRGBAEditorHeatScale->SetPointSymmetry( 0, 5 );
      mRGBAEditorHeatScale->SetPointSymmetry( 1, 4 );
      mRGBAEditorHeatScale->SetPointSymmetry( 2, 3 );

      } break;
    case LUT: {

      mLabelLUTFileName = vtkSmartPointer<vtkKWLabel>::New();
      mLabelLUTFileName->SetParent( frame );
      mLabelLUTFileName->Create();

      mLabelLUTFileName->SetText( "No LUT loaded" );

      vtkSmartPointer<vtkKWPushButton> loadButton = 
	vtkSmartPointer<vtkKWPushButton>::New();
      loadButton->SetParent( frame );
      loadButton->Create();

      loadButton->SetText( "Load LUT" );
      loadButton->SetCommand( this, "LoadLUTFromDlog" );

      this->Script( "pack %s -side left -fill x -expand y",
		    mLabelLUTFileName->GetWidgetName() );
      this->Script( "pack %s -side left",
		    loadButton->GetWidgetName() );
      } break;
    }
  }

  // --------------------------------------------------------------------


  this->Script( "pack %s %s %s %s -side top -fill x -anchor nw",
                radBtnSetInterpolation->GetWidgetName(),
                chkBtnSmooth->GetWidgetName(),
                radBtnSetColorMap->GetWidgetName(),
                maFrameColorMapSettings[GrayScale]->GetWidgetName() );

  this->ColorMapChanged();
}

void
vtkKWScubaLayerCollectionMRI::RemoveControls () {

  for ( int nType = (int)GrayScale; nType < (int)cColorMapTypes;  nType++ ) {
    maFrameColorMapSettings[nType] = NULL;
    maRadBtnColorMap[nType] = NULL;
  }

  mRangeWindowLevel = NULL;
  mRangeVisibleValue = NULL;
  mRGBAEditorHeatScale = NULL;
  mLabelLUTFileName = NULL;
}

vtkFSVolumeSource* 
vtkKWScubaLayerCollectionMRI::GetSource () const {
  return mSource;
}

vtkFreesurferLookupTable*
vtkKWScubaLayerCollectionMRI::GetLUTTable () const {
  return mLUTTable;
}

vtkRGBATransferFunction*
vtkKWScubaLayerCollectionMRI::GetGrayScaleTable () const {
  return mGrayScaleTable;
}

vtkRGBATransferFunction*
vtkKWScubaLayerCollectionMRI::GetHeatScaleTable () const {
  return mHeatScaleTable;
}

COLOR_TABLE*
vtkKWScubaLayerCollectionMRI::GetLUTCTAB () const {
  return mFreeSurferCTAB;
}

void
vtkKWScubaLayerCollectionMRI::ColorMapChanged () {

  switch ( mColorMapType ) {
  case NoColorMap:
    break;

  case GrayScale:

    // Check the color map variables and update range sliders if
    // necessary.
    if ( mMinGrayscaleWindow < mMinVisibleValue )
      mMinVisibleValue = mMinGrayscaleWindow;
    if ( mMaxGrayscaleWindow > mMaxVisibleValue )
      mMaxVisibleValue = mMaxGrayscaleWindow;

    // Build our lookup table.
    assert( mGrayScaleTable.GetPointer() );
    mGrayScaleTable->RemoveAllPoints();
    mGrayScaleTable->AddRGBAPoint( mMinVisibleValue-0.001, 0, 0, 0, 0 );
    mGrayScaleTable->AddRGBAPoint( mMinVisibleValue,       0, 0, 0, 1 );
    mGrayScaleTable->AddRGBAPoint( mMinGrayscaleWindow,    0, 0, 0, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxGrayscaleWindow,    1, 1, 1, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue,       1, 1, 1, 1 );
    mGrayScaleTable->AddRGBAPoint( mMaxVisibleValue+0.001, 1, 1, 1, 0 );
    mGrayScaleTable->Build();
  break;

  case HeatScale:
    break;

  case LUT:
    break;

  default:
    break;
  }

  // Update the packed settings panel.
  if ( maFrameColorMapSettings[mColorMapType].GetPointer() ) {
    for ( int nType = (int)GrayScale; nType < (int)cColorMapTypes; nType++ ) {
      if ( nType == (int)mColorMapType ) {
        this->Script( "pack %s -side top -fill x",
                      maFrameColorMapSettings[nType]->GetWidgetName() );
      } else {
        this->Script( "pack forget %s",
                      maFrameColorMapSettings[nType]->GetWidgetName() );
      }
    }
  }

  // Update the controls.
  if ( mRangeVisibleValue.GetPointer() ) {
    mRangeVisibleValue->DisableCommandsOn();
    mRangeVisibleValue->SetRange( mMinVisibleValue, mMaxVisibleValue );
    mRangeVisibleValue->DisableCommandsOff();
  }
  if ( mRangeWindowLevel.GetPointer() ) {
    mRangeWindowLevel->DisableCommandsOn();
    mRangeWindowLevel->SetRange( mMinGrayscaleWindow, mMaxGrayscaleWindow );
    mRangeWindowLevel->DisableCommandsOff();
  }
  if ( mRGBAEditorHeatScale.GetPointer() ) {
    mRGBAEditorHeatScale->Update();
  }

  if ( mLabelLUTFileName .GetPointer()) {
    if ( mFreeSurferCTAB ) {
      char fnLUT[1024];
      CTABcopyFileName( mFreeSurferCTAB, fnLUT, sizeof(fnLUT) );
      mLabelLUTFileName->SetText( fnLUT );
    } else {
      mLabelLUTFileName->SetText( "Choose an LUT file" );
    }
  }

  // Notify the layers that use the color map stuff.
  this->SendBroadcast( "ColorMapChanged", NULL );
}

void
vtkKWScubaLayerCollectionMRI::SetColorMap ( ColorMapType iType ) {

  if ( mColorMapType != iType ) {
    mColorMapType = iType;
    this->ColorMapChanged();
  }
}

ScubaCollectionPropertiesMRI::ColorMapType
vtkKWScubaLayerCollectionMRI::GetColorMap () const {
  return mColorMapType;
}

void
vtkKWScubaLayerCollectionMRI::SetColorMapToGrayScale () {
  SetColorMap( GrayScale );
}

void
vtkKWScubaLayerCollectionMRI::SetColorMapToHeatScale () {
  SetColorMap( HeatScale );
}

void
vtkKWScubaLayerCollectionMRI::SetColorMapToLUT () {
  SetColorMap( LUT );
}

void
vtkKWScubaLayerCollectionMRI::SetResliceInterpolation ( int iMode ) {

  if( mResliceInterpolation != iMode ) {
    mResliceInterpolation = iMode;
    this->SendBroadcast( "ResliceInterpolationChanged", NULL );
  }
}

int
vtkKWScubaLayerCollectionMRI::GetResliceInterpolation () const {
  return mResliceInterpolation;
}

void
vtkKWScubaLayerCollectionMRI::SetTextureSmoothing ( int iSmooth ) {

  if( mTextureSmoothing != iSmooth ) {
    mTextureSmoothing = iSmooth;
    this->SendBroadcast( "TextureSmoothingChanged", NULL );
  }
}

int 
vtkKWScubaLayerCollectionMRI::GetTextureSmoothing () const {
  return mTextureSmoothing;
}

void
vtkKWScubaLayerCollectionMRI::SetMinVisibleValue ( float iValue ) {
  if ( mMinVisibleValue != iValue ) {
    mMinVisibleValue = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetMinVisibleValue () {
  return mMinVisibleValue;
}

void
vtkKWScubaLayerCollectionMRI::SetMaxVisibleValue ( float iValue ) {
  if ( mMaxVisibleValue != iValue ) {
    mMaxVisibleValue = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetMaxVisibleValue () {
  return mMaxVisibleValue;
}

void
vtkKWScubaLayerCollectionMRI::SetMinMaxVisibleValue ( float iMinValue,
    float iMaxValue ) {
  if ( mMinVisibleValue != iMinValue ||
       mMaxVisibleValue != iMaxValue ) {
    mMinVisibleValue = iMinValue;
    mMaxVisibleValue = iMaxValue;
  }
  this->ColorMapChanged();
}

void
vtkKWScubaLayerCollectionMRI::SetWindow ( float iWindow ) {
  float level = this->GetLevel();
  this->SetMinMaxGrayscaleWindow( level - iWindow/2.0,
                                  level + iWindow/2.0 );
}

float
vtkKWScubaLayerCollectionMRI::GetWindow () {
  return (mMaxGrayscaleWindow - mMinGrayscaleWindow);
}

void
vtkKWScubaLayerCollectionMRI::SetLevel ( float iLevel ) {
  float window = this->GetWindow();
  this->SetMinMaxGrayscaleWindow( iLevel - window/2.0,
                                  iLevel + window/2.0 );
}

float
vtkKWScubaLayerCollectionMRI::GetLevel () {
  return ((mMaxGrayscaleWindow - mMinGrayscaleWindow) / 2.0) +
         mMinGrayscaleWindow;
}

void
vtkKWScubaLayerCollectionMRI::SetMinMaxGrayscaleWindow ( float iMin, float iMax ) {
  this->SetMinGrayscaleWindow( iMin );
  this->SetMaxGrayscaleWindow( iMax );
}

void
vtkKWScubaLayerCollectionMRI::SetMinGrayscaleWindow ( float iMin ) {
  if ( mMinGrayscaleWindow != iMin ) {
    mMinGrayscaleWindow = iMin;
    this->ColorMapChanged();
  }
}

void
vtkKWScubaLayerCollectionMRI::SetMaxGrayscaleWindow ( float iMax ) {
  if ( mMaxGrayscaleWindow != iMax ) {
    mMaxGrayscaleWindow = iMax;
    this->ColorMapChanged();
  }
}

void
vtkKWScubaLayerCollectionMRI::SetHeatScaleMinThreshold ( float iValue ) {
  if ( mHeatScaleMinThreshold != iValue ) {
    mHeatScaleMinThreshold = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetHeatScaleMinThreshold () {
  return mHeatScaleMinThreshold;
}

void
vtkKWScubaLayerCollectionMRI::SetHeatScaleMidThreshold ( float iValue ) {
  if ( mHeatScaleMidThreshold != iValue ) {
    mHeatScaleMidThreshold = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetHeatScaleMidThreshold () {
  return mHeatScaleMidThreshold;
}

void
vtkKWScubaLayerCollectionMRI::SetHeatScaleMaxThreshold ( float iValue ) {
  if ( mHeatScaleMaxThreshold != iValue ) {
    mHeatScaleMaxThreshold = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetHeatScaleMaxThreshold () {
  return mHeatScaleMaxThreshold;
}

void
vtkKWScubaLayerCollectionMRI::SetHeatScaleOffset ( float iValue ) {
  if ( mHeatScaleOffset != iValue ) {
    mHeatScaleOffset = iValue;
    this->ColorMapChanged();
  }
}

float
vtkKWScubaLayerCollectionMRI::GetHeatScaleOffset () {
  return mHeatScaleOffset;
}

void
vtkKWScubaLayerCollectionMRI::SetReverseHeatScale ( bool ib ) {
  if ( mbReverseHeatScale != ib ) {
    mbReverseHeatScale = ib;
    this->ColorMapChanged();
  }
}

bool 
vtkKWScubaLayerCollectionMRI::GetReverseHeatScale () {
  return mbReverseHeatScale;
}

void
vtkKWScubaLayerCollectionMRI::SetShowPositiveHeatScaleValues ( bool ib ) {
  if ( mbShowPositiveHeatScaleValues != ib ) {
    mbShowPositiveHeatScaleValues = ib;
    this->ColorMapChanged();
  }
}

bool 
vtkKWScubaLayerCollectionMRI::GetShowPositiveHeatScaleValues () {
  return mbShowPositiveHeatScaleValues;
}

void
vtkKWScubaLayerCollectionMRI::SetShowNegativeHeatScaleValues ( bool ib ) {
  if ( mbShowNegativeHeatScaleValues != ib ) {
    mbShowNegativeHeatScaleValues = ib;
    this->ColorMapChanged();
  }
}

bool
vtkKWScubaLayerCollectionMRI::GetShowNegativeHeatScaleValues () {
  return mbShowNegativeHeatScaleValues;
}

void
vtkKWScubaLayerCollectionMRI::EditorChangedHeatScale () {
  this->ColorMapChanged();
}

void
vtkKWScubaLayerCollectionMRI::LoadLUTFromDlog () { 

  // Create a Load dialog and set it up.
  vtkKWLoadSaveDialog* dialog = vtkKWLoadSaveDialog::New();
  dialog->SetApplication( this->GetApplication() );
  dialog->Create();
  dialog->SetTitle( "Load an LUT" );
  dialog->SetFileTypes( "{LUT {.txt}} "
                        "{All {*}}" );
  dialog->RetrieveLastPathFromRegistry( "LoadLUT" );
  dialog->SetDefaultExtension( ".txt" );

  // Show the dialog, and when it returns, Invoke() will be true if
  // they clicked OK and gave us a filename.
  if ( dialog->Invoke() ) {
    dialog->SaveLastPathToRegistry( "LoadLUT" );
    this->LoadFreesurferLUT( dialog->GetFileName() );
  }
}

void
vtkKWScubaLayerCollectionMRI::LoadFreesurferLUT ( const char* ifn ) {

  // Try to load the table.
  char fnLUT[1024];
  strncpy( fnLUT, ifn, sizeof(fnLUT) );
  COLOR_TABLE* ctab = CTABreadASCII( fnLUT );
  if ( NULL == ctab ) {
    throw runtime_error( string("Couldn't open color table file ") + ifn );
  }

  // If we got it, delete the existing table.
  if( NULL != mFreeSurferCTAB ) 
    CTABfree( &mFreeSurferCTAB );

  // Save the table.
  mFreeSurferCTAB = ctab;

  // Build our VTK object from the table.
  mLUTTable->BuildFromCTAB( mFreeSurferCTAB );

  this->ColorMapChanged();
}

vtkKWScubaLayer*
vtkKWScubaLayerCollectionMRI::MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode ) {

  vtkKWScubaLayer* layer = NULL;
  if( vtkKWScubaView::TwoDee == iMode ) {

    vtkKWScubaLayer2DMRI* layer2DMRI = vtkKWScubaLayer2DMRI::New();
    layer2DMRI->SetMRIProperties( this );
    
    layer = (vtkKWScubaLayer*)layer2DMRI;

  } else if( vtkKWScubaView::ThreeDee == iMode ) {

    vtkKWScubaLayer3DMRI* layer3DMRI = vtkKWScubaLayer3DMRI::New();
    layer3DMRI->SetMRIProperties( this );
    
    layer = (vtkKWScubaLayer*)layer3DMRI;

  }

  return layer;
}

void
vtkKWScubaLayerCollectionMRI::LoadVolumeFromFileName () {

  // Source object reads the volume and outputs a volume.
  mSource = vtkSmartPointer<vtkFSVolumeSource>::New();
  mSource->MRIRead( mfnVolume.c_str() );
  mSource->Update();

  // Init our color scale values.
  mMinVisibleValue = mSource->GetMinValue();
  mMaxVisibleValue = mSource->GetMaxValue();
  mMinGrayscaleWindow = mMinVisibleValue;
  mMaxGrayscaleWindow = mMaxVisibleValue;

  float oneTenth;
  float highestAbsValue;
  highestAbsValue = max( fabs(mMinVisibleValue), fabs(mMaxVisibleValue) );
  oneTenth = highestAbsValue / 10.0;
  mHeatScaleMinThreshold = oneTenth;
  mHeatScaleMidThreshold = highestAbsValue / 2.0;
  mHeatScaleMaxThreshold = highestAbsValue - oneTenth;

  //
  // Color tables.
  //
  mGrayScaleTable = vtkSmartPointer<vtkRGBATransferFunction>::New();
  mHeatScaleTable = vtkSmartPointer<vtkRGBATransferFunction>::New();
  mLUTTable = vtkSmartPointer<vtkFreesurferLookupTable>::New();
  mGrayScaleTable->ClampingOff();
  mHeatScaleTable->ClampingOn();

  // Init the colors in the heat scale. The editor will handle the
  // editing from here.
  mHeatScaleTable->AddRGBAPoint( -mHeatScaleMaxThreshold, 0, 1, 1, 1 );
  mHeatScaleTable->AddRGBAPoint( -mHeatScaleMidThreshold, 0, 0, 1, 1 );
  mHeatScaleTable->AddRGBAPoint( -mHeatScaleMinThreshold, 0, 0, 1, 0 );
  mHeatScaleTable->AddRGBAPoint(  mHeatScaleMinThreshold, 1, 0, 0, 0 );
  mHeatScaleTable->AddRGBAPoint(  mHeatScaleMidThreshold, 1, 0, 0, 1 );
  mHeatScaleTable->AddRGBAPoint(  mHeatScaleMaxThreshold, 1, 1, 0, 1 );
  mHeatScaleTable->Build();

  // Set up our initial tables.
  this->ColorMapChanged();
}
