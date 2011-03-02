/**
 * @file  ScubaLayer2DMRI.cpp
 * @brief Layer for displaying MRI volumes
 *
 * Displays MRI volumes in grayscale, heatscale, and with an
 * LUT. Handles multiframe volumes and can interpret them as having
 * time points and conditions. Contains the logic for the Edge Path
 * and other path tools.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:37 $
 *    $Revision: 1.166 $
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

#include <limits>
#include <fstream>
#include <iomanip>
#include "ScubaLayer2DMRI.h"
#include "ViewState.h"
#include "ProgressDisplayManager.h"
#include "Utilities.h"
#include "Array2.h"
#include "PathManager.h"
#include "VectorOps.h"
#include "ScubaView.h"
#include "DataManager.h"
extern "C" {
  // Have to do this define trick here because these header files use
  // 'this' as a variable name in the declaration, which messes up the
  // c++ compiler.
#define this xthis
#include "talairachex.h"
#include "mri2.h"
#include "error.h"
#undef this
}

using namespace std;

int const ScubaLayer2DMRI::cGrayscaleLUTEntries = 256;
int const ScubaLayer2DMRI::kMaxPixelComponentValue = 255;
float const ScubaLayer2DMRI::kMaxPixelComponentValueFloat = 255.0;
int const ScubaLayer2DMRI::kcTimersBetweenAutosaves = 60000;
char* const ScubaLayer2DMRI::kaReportableInfo[ScubaLayer2DMRI::kcReportableInfo] = 
{ (char*)"Value", (char*)"Index", (char*)"Talairach" };


/* When setting the initial window/level, use a histogram and cut this
   much off either end to find the initial range. */
float const kPercentValuesToCut = 0.05;
int const zGrayscaleHistogramBins = 200;



ScubaLayer2DMRI::ScubaLayer2DMRI () :
    mTimersSinceLastAutosave(0),
    mVolume(NULL),
    mSampleMethod(nearest),
    mColorMapMethod(grayscale),
    mbClearZero(false),
    mMinVisibleValue(-99999),
    mMaxVisibleValue(99999),
    mOldMinValue(-99999),
    mOldMaxValue(99999),
    mCurrentFrame(0),
    mBrightness(0.25),
    mContrast(12.0),
    mWindow(1.0),
    mLevel(0.5),
    mHeatScaleMinThreshold(0),
    mHeatScaleMidThreshold(0),
    mHeatScaleMaxThreshold(0),
    mHeatScaleOffset(0),
    mbReverseHeatScale(false),
    mbShowPositiveHeatScaleValues(true),
    mbShowNegativeHeatScaleValues(true),
    mColorLUT(NULL),
    mMaskVolume(NULL),
    mROIOpacity(0.7),
    mbEditableROI(true),
    mbDrawMIP(false),
    mbDrawEditingLine(false),
    mCurrentPath(false),
    mRowStartRAS(NULL),
    mColIncrementRAS(NULL),
    mCurrentDrawingOperationActionID(-1),
    mEdgePathFinder( NULL ) {

  SetOutputStreamToCerr();

  mBufferIncSize[0] = mBufferIncSize[1] = -1;
  mLastMouseUpRAS.Set( 0, 0, 0 );

  // Try setting our initial color LUT to the default LUT with
  // id 0. If it's not there, create it.
  try {
    mColorLUT = &(ScubaColorLUT::FindByID( 0 ));
  } catch (...) {

    ScubaColorLUT* lut = new ScubaColorLUT();
    lut->SetLabel( "Default" );

    try {
      mColorLUT = &(ScubaColorLUT::FindByID( 0 ));
    } catch (...) {
      DebugOutput( << "Couldn't make default lut!" );
    }
  }

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "Set2DMRILayerVolumeCollection", 2,
                         "layerID collectionID",
                         "Sets the volume collection for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerVolumeCollection", 1,
                         "layerID",
                         "Returns the volume collection for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerCurrentFrame", 2,
                         "layerID frame",
                         "Sets the current frame for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerCurrentFrame", 1, "layerID",
                         "Returns the current frame for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerCurrentCondition", 2,
                         "layerID condition",
                         "Sets the current condition for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerCurrentCondition", 1, "layerID",
                         "Returns the current condition for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerCurrentTimePoint", 2,
                         "layerID timePoint",
                         "Sets the current time point for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerCurrentTimePoint", 1, "layerID",
                         "Returns the current time point for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerColorMapMethod", 2,
                         "layerID method",
                         "Sets the color map method for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerColorMapMethod", 1, "layerID",
                         "Returns the color map method for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerSampleMethod", 2,
                         "layerID method",
                         "Sets the sample method for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerSampleMethod", 1, "layerID",
                         "Returns the sample method for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerBrightness", 2,
                         "layerID brightness",
                         "Sets the brightness for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerBrightness", 1, "layerID",
                         "Returns the brightness for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerContrast", 2,
                         "layerID contrast",
                         "Sets the contrast for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerContrast", 1, "layerID",
                         "Returns the contrast for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerWindow", 2,
                         "layerID window",
                         "Sets the window for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerWindow", 1, "layerID",
                         "Returns the window for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerLevel", 2,
                         "layerID level",
                         "Sets the level for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerLevel", 1, "layerID",
                         "Returns the level for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerColorLUT", 2, "layerID lutID",
                         "Sets the LUT  for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerColorLUT", 1, "layerID",
                         "Returns the LUT id for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerDrawZeroClear", 2,
                         "layerID drawClear", "Sets property for drawing"
                         "values of zero clear." );
  commandMgr.AddCommand( *this, "Get2DMRILayerDrawZeroClear", 1, "layerID",
                         "Returns the value of the property for drawing"
                         "values of zero clear." );
  commandMgr.AddCommand( *this, "Set2DMRILayerDrawMIP", 2,
                         "layerID drawMIP", "Sets property for drawing"
                         "the maximum intensity projection." );
  commandMgr.AddCommand( *this, "Get2DMRILayerDrawMIP", 1, "layerID",
                         "Returns the value of the property for drawing"
                         "the maximum intensity projection." );
  commandMgr.AddCommand( *this, "Set2DMRILayerMinVisibleValue", 2,
                         "layerID value", "Sets the minimum value to be drawn."
                         "values of zero clear." );
  commandMgr.AddCommand( *this, "Get2DMRILayerMinVisibleValue", 1, "layerID",
                         "Returns the minimum value to be drawn." );
  commandMgr.AddCommand( *this, "Set2DMRILayerMaxVisibleValue", 2,
                         "layerID value", "Sets the maximum value to be drawn."
                         "values of zero clear." );
  commandMgr.AddCommand( *this, "Get2DMRILayerMaxVisibleValue", 1, "layerID",
                         "Returns the maximum value to be drawn." );
  commandMgr.AddCommand( *this, "Get2DMRILayerMinValue", 1, "layerID",
                         "Returns the minimum value of the volume." );
  commandMgr.AddCommand( *this, "Get2DMRILayerMaxValue", 1, "layerID",
                         "Returns the maximum value of the volume." );
  commandMgr.AddCommand( *this, "Get2DMRILayerHeatScaleMin", 1, "layerID",
                         "Returns the heat scale min value for the layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerHeatScaleMin", 2,
                         "layerID value", "Sets the heat scale min value "
                         "for the layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerHeatScaleMid", 1, "layerID",
                         "Returns the heat scale mid value for the layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerHeatScaleMid", 2,
                         "layerID value", "Sets the heat scale mid value "
                         "for the layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerHeatScaleMax", 1, "layerID",
                         "Returns the heat scale max value for the layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerHeatScaleMax", 2,
                         "layerID value", "Sets the heat scale max value "
                         "for the layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerHeatScaleOffset", 1, "layerID",
                         "Returns the heat scale offset for the layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerHeatScaleOffset", 2,
                         "layerID value", "Sets the heat scale offset "
                         "for the layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerReverseHeatScale", 1,
                         "layerID", "Returns whether or not the heat "
                         "scale colors are reversed." );
  commandMgr.AddCommand( *this, "Set2DMRILayerReverseHeatScale", 2,
                         "layerID reverse", "Set whether or not the heat "
                         "scale colors are reversed." );
  commandMgr.AddCommand( *this, "Get2DMRILayerShowPositiveHeatScaleValues", 1,
                         "layerID", "Returns whether or not positive heat "
                         "scale values are being drawn." );
  commandMgr.AddCommand( *this, "Set2DMRILayerShowPositiveHeatScaleValues", 2,
                         "layerID show", "Set whether or not positive heat "
                         "scale values are being drawn." );
  commandMgr.AddCommand( *this, "Get2DMRILayerShowNegativeHeatScaleValues", 1,
                         "layerID", "Returns whether or not negative heat "
                         "scale values are being drawn." );
  commandMgr.AddCommand( *this, "Set2DMRILayerShowNegativeHeatScaleValues", 2,
                         "layerID show", "Set whether or not negative heat "
                         "scale values are being drawn." );
  commandMgr.AddCommand( *this, "Set2DMRILayerHeatScaleUsingFDR", 3,
                         "layerID rate maskVolume", "Uses a false "
                         "discovery rate algorithm to find a heat scale "
                         "threshold for the bolume. rate should be 0-1. "
                         "If maskVolume is not empty, will load "
                         "a volume to use as a mask." );
  commandMgr.AddCommand( *this, "Set2DMRILayerMaskVolume", 2,
                         "layerID volumeID", "Sets the ID of the mask volume "
                         "for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerMaskVolume", 1, "layerID",
                         "Gets the ID of the mask volume for this layer." );
  commandMgr.AddCommand( *this, "Remove2DMRILayerMaskVolume", 1, "layerID",
                         "Stops using a mask volume for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerROIOpacity", 2,"layerID opacity",
                         "Sets the opacity of the ROI for a layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerROIOpacity", 1, "layerID",
                         "Returns the opacity of the ROI for a layer." );
  commandMgr.AddCommand( *this, "Set2DMRILayerEditableROI", 2,
                         "layerID editable", "Specify whether or not this "
                         "layer's ROI is editable." );
  commandMgr.AddCommand( *this, "Get2DMRILayerEditableROI", 1, "layerID",
                         "Returns whether or not this layer's ROI is "
                         "editable." );
  commandMgr.AddCommand( *this, "Get2DMRIRASCoordsFromIndex", 4,
                         "layerID x y z", "Returns a list of RAS coords "
                         "converted from the input index coords." );
  commandMgr.AddCommand( *this, "Flood2DMRIVolume", 7,
                         "layerID x y z toolID viewID type",
                         "Floods the volume in a 2DMRIS layer at the input "
                         "RAS point x,y,z, using the input tool and view "
                         "for settings. Floods are of the following types: "
                         "voxelEditingNew, voxelEditingErase, "
                         "roiEditingSelect, roiEditingUnselect." );

  // Add the items we can report. 
  for( int nInfo = 0; nInfo < kcReportableInfo; nInfo++ )
    AddReportableInfo( kaReportableInfo[nInfo] );

  // Init our color opacity cache for our initial value.
  InitColorOpacityCache();
}

ScubaLayer2DMRI::~ScubaLayer2DMRI () {
  if ( NULL != mVolume )
    mVolume->RemoveListener( *this );

  delete mEdgePathFinder;
}

void
ScubaLayer2DMRI::SetVolumeCollection ( VolumeCollection& iVolume ) {

  mVolume = &iVolume;

  // This makes sure the volume is loaded here instead of later.
  mVolume->GetMRI();

  // Set our color scale settings to default values. Min/max visible
  // go to the min/max values, level is set to the middle of the value
  // range, and window is set to entire range.
  SetMinMaxVisibleValue( mVolume->GetMRIMinValue(),
                         mVolume->GetMRIMaxValue() );

  // MRIhistogram doesn't work on vols with more than one frame. :(
  if ( mVolume->GetNumberOfFrames() > 1 ) {

    SetLevel( ((mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue()) / 2.0) +
              mVolume->GetMRIMinValue() );
    SetWindow( mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue() );

  } else {

    // We're going to get a historgram of the volume and use it to
    // remove any outliers when determining our initial
    // window/level. Create the histogram.
    HISTOGRAM* histo = MRIhistogram( const_cast<MRI*>(mVolume->GetMRI()),
                                     zGrayscaleHistogramBins );
    if ( NULL != histo ) {

      // Calculate the number of values in the volume.
      int range[3];
      mVolume->GetMRIIndexRange( range );
      int zValues = range[0] * range[1] * range[2];

      // Get the number of values to cut based on the percentage.
      int zCutValues = (int)((float)zValues * (kPercentValuesToCut/100.0));

      // Get the index of the bin with the value zero in it. We'll skip
      // this bin because we want to ignore zeroes in our histogram.
      int nZeroBin = HISTOfindBin( histo, 0.0 );

      // Start at the bottom and go up. Sum up the number of the values
      // in the bins. When we get to the number of values we want to
      // cut, stop.
      int cCurValues = 0;
      int nBin = 0;
      while ( cCurValues < zCutValues &&
              nBin < zGrayscaleHistogramBins ) {
        if ( nBin != nZeroBin )
          cCurValues += (int)histo->counts[nBin];
        nBin++;
      }
      // Our bin values are the number at the top range of the bin, so
      // use the bin underneath us.
      float minValue = histo->bins[nBin>0?nBin-1:0];

      // Do the same for the top end, going down.
      cCurValues = 0;
      nBin = histo->nbins-1;
      while ( cCurValues < zCutValues &&
              nBin >= 0 ) {
        if ( nBin != nZeroBin )
          cCurValues += (int)histo->counts[nBin];
        nBin--;
      }
      float maxValue = histo->bins[nBin];

      // In the special case of the fixed constant volume, minValue and
      // MaxValue will be the same.
      if ( fabs( minValue - maxValue ) < 0.0001 ) {

        // Just set the level window appropriately.
        SetLevel( minValue );
        SetWindow( 1 );

        // If max < min, then it's a binary volume and we crossed each
        // other, so just set it to the normal full range.
      } else if ( maxValue < minValue ) {

        SetLevel( ((mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue())
                   / 2.0) +
                  mVolume->GetMRIMinValue() );
        SetWindow( mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue() );

      } else {

        // Use this max and min to define a range and set our level and
        // window.
        SetLevel( (maxValue - minValue)/2.0 + minValue );
        SetWindow( (maxValue - minValue) );
      }

      // Free the histogram.
      HISTOfree( &histo );

    } else {

      // Couldn't histogram for some reason, so use defaults.
      SetLevel( ((mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue())
                 / 2.0) +
                mVolume->GetMRIMinValue() );
      SetWindow( mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue() );

    }
  }

  // Calc one tenth of the abs range.
  float oneTenth;
  float highestAbsValue;
  highestAbsValue =
    max( fabs(mVolume->GetMRIMinValue()), fabs(mVolume->GetMRIMaxValue()) );
  oneTenth = highestAbsValue / 10.0;

  // Min heat goes to 1/10 above 0, max goes to 1/10 below the
  // highestAbsValue, and mid is in the middle.
  SetHeatScaleMaxThreshold( highestAbsValue - oneTenth );
  SetHeatScaleMidThreshold( highestAbsValue / 2.0 );
  SetHeatScaleMinThreshold( oneTenth );

  SetHeatScaleOffset( 0 );

  // Save these min/max values.
  mOldMinValue = mVolume->GetMRIMinValue();
  mOldMaxValue = mVolume->GetMRIMaxValue();

  // Listen to the volume for dataChanged messages.
  mVolume->AddListener( *this );

  // Init our grayscale.
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                                  ViewState& iViewState,
                                  ScubaWindowToRASTranslator& iTranslator ) {

  if ( NULL == mVolume ) {
    DebugOutput( << "No volume to draw" );
    return;
  }

  if ( !mbVisible ) {
    return;
  }

  if ( mbDrawMIP ) {
    DrawMIPIntoBuffer( iBuffer, iWidth, iHeight, iViewState, iTranslator );
    return;
  }

  GLubyte* dest = iBuffer;
  int window[2], window2[2];
  float RAS[3], RAS2[3];
  float value = 0;
  int color[3];

  // Init our buffers if necessary.
  if ( iWidth != mBufferIncSize[0] || iHeight != mBufferIncSize[1] ) {
    InitBufferCoordCache( iWidth, iHeight );
  }

  // Find the increments. To do this, take two window coords on this
  // row, one at col 0 and one at col 1. Then we convert both to RAS
  // and find the difference. We use that difference to increment RAS
  // coords across this row for all columns.
  for ( int nRow = 0; nRow < iHeight; nRow++ ) {

    window[0] = 0;
    window[1] = nRow;
    iTranslator.TranslateWindowToRAS( window, mRowStartRAS[nRow] );

    window2[0] = 1;
    window2[1] = nRow;
    iTranslator.TranslateWindowToRAS( window2, RAS2 );
    mColIncrementRAS[nRow][0] = RAS2[0] - mRowStartRAS[nRow][0];
    mColIncrementRAS[nRow][1] = RAS2[1] - mRowStartRAS[nRow][1];
    mColIncrementRAS[nRow][2] = RAS2[2] - mRowStartRAS[nRow][2];
  }

  // Get the update bounds. We'll only iterate over those.
  int windowUpdateBounds[4];
  iViewState.CopyUpdateRect( windowUpdateBounds );

  // Create a dummy location, we'll change it soon. We'll create a
  // location for the mask volume too, if it exists.
  RAS[0] = RAS[1] = RAS[2] = 0;
  VolumeLocation 
    loc( mVolume->MakeVolumeLocationFromRAS( RAS, mCurrentFrame ) );
  VolumeLocation* maskLoc = NULL;
  if( mMaskVolume )
    maskLoc = 
      new VolumeLocation( mMaskVolume->MakeVolumeLocationFromRAS( RAS, 0 ) );

  for ( window[1] = windowUpdateBounds[1];
        window[1] <= windowUpdateBounds[3]; window[1]++ ) {

    // Grab the RAS beginning for this row and column.
    RAS[0] = mRowStartRAS[window[1]][0] +
             windowUpdateBounds[0]*mColIncrementRAS[window[1]][0];
    RAS[1] = mRowStartRAS[window[1]][1] +
             windowUpdateBounds[0]*mColIncrementRAS[window[1]][1];
    RAS[2] = mRowStartRAS[window[1]][2] +
             windowUpdateBounds[0]*mColIncrementRAS[window[1]][2];

    // Find this row start.
    dest = iBuffer +
           (((iWidth * window[1]) + windowUpdateBounds[0]) * mBytesPerPixel);
    for ( window[0] = windowUpdateBounds[0];
          window[0] <= windowUpdateBounds[2]; window[0]++ ) {

      // Set the location from this RAS.
      loc.SetFromRAS( RAS );
      if ( maskLoc ) maskLoc->SetFromRAS( RAS );

      // If the location is in bounds, and if there's a mask volume, if
      // this value in the mask volume is not 0...
      int selectColor[3];
      if ( mVolume->IsInBounds( loc ) &&
           (NULL == mMaskVolume ||
            (mMaskVolume->IsInBounds( *maskLoc ) &&
             mMaskVolume->GetMRINearestValue( *maskLoc ) != 0)) ) {

        switch ( mSampleMethod ) {
        case nearest:
          value = mVolume->GetMRINearestValue( loc );
          break;
        case trilinear:
          value = mVolume->GetMRITrilinearValue( loc );
          break;
        case sinc:
          value = mVolume->GetMRISincValue( loc );
          break;
        case magnitude:
          value = mVolume->GetMRIMagnitudeValue( loc );
          break;
        }

        switch ( mColorMapMethod ) {
        case grayscale:
          GetGrayscaleColorForValue( value, dest, color );
          break;
        case heatScale:
          GetHeatscaleColorForValue( value, dest, color );
          break;
        case LUT:
          GetColorLUTColorForValue( value, dest, color );
          break;
        }

        dest[0] = mColorTimesOneMinusOpacity[dest[0]] +
                  mColorTimesOpacity[color[0]];
        dest[1] = mColorTimesOneMinusOpacity[dest[1]] +
                  mColorTimesOpacity[color[1]];
        dest[2] = mColorTimesOneMinusOpacity[dest[2]] +
                  mColorTimesOpacity[color[2]];

        if ( mVolume->IsSelected( loc, selectColor ) ) {
          // Write the RGB value to the buffer. Write a 255 in the
          // alpha byte.
          dest[0] = (GLubyte) (((float)dest[0] * (1.0 - mROIOpacity)) +
                               ((float)selectColor[0] * mROIOpacity));
          dest[1] = (GLubyte) (((float)dest[1] * (1.0 - mROIOpacity)) +
                               ((float)selectColor[1] * mROIOpacity));
          dest[2] = (GLubyte) (((float)dest[2] * (1.0 - mROIOpacity)) +
                               ((float)selectColor[2] * mROIOpacity));
        }
      }

      // Increment the dest buffer pointer.
      dest += mBytesPerPixel;

      // Increment the RAS point.
      RAS[0] += mColIncrementRAS[window[1]][0];
      RAS[1] += mColIncrementRAS[window[1]][1];
      RAS[2] += mColIncrementRAS[window[1]][2];
    }
  }

  delete maskLoc;

  if ( mbDrawEditingLine ) {

    Point2<int> windowA, windowB;
    int color[3] = { 255, 0, 255 };
    iTranslator.TranslateRASToWindow( mLastMouseUpRAS.xyz(), windowA.xy() );
    iTranslator.TranslateRASToWindow( mCurrentMouseRAS.xyz(), windowB.xy() );
    DrawLineIntoBuffer( iBuffer, iWidth, iHeight,
                        windowA.xy(), windowB.xy(), color, 2, 0.5 );
  }

}

void
ScubaLayer2DMRI::DrawIntoGL ( ViewState& iViewState,
                              ScubaWindowToRASTranslator& iTranslator ) {}

void
ScubaLayer2DMRI::DrawMIPIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                                     ViewState& iViewState,
                                     ScubaWindowToRASTranslator& iTranslator ) {

  // This functino is similar to DrawIntoBuffer except we iterate over
  // a series of coordinates, drawing all 'planes' of the volume, and
  // only draw pixel values that are lighter than what is already
  // there.

  if ( NULL == mVolume ) {
    DebugOutput( << "No volume to draw" );
    return;
  }

  if ( !mbVisible ) {
    return;
  }

  GLubyte* dest = iBuffer;
  int window[2], window2[2];
  float RAS[3], RAS2[3];
  float value = 0;
  int color[3];

  // Init our buffers if necessary.
  if ( iWidth != mBufferIncSize[0] || iHeight != mBufferIncSize[1] ) {
    InitBufferCoordCache( iWidth, iHeight );
  }

  // We need to iterate over all the planes in our current view
  // state. To do this, we switch on the current in plane and iterate
  // over the bounds of the volume in that plane, moving by our in
  // plane increment. First select the in plane and get the proper
  // bounds.
  float volumeRASRange[6];
  mVolume->GetDataRASBounds( volumeRASRange );
  float increments[3];
  GetPreferredThroughPlaneIncrements( increments );
  float minZ = 0, maxZ = 0, incZ = 0;
  switch ( iViewState.GetInPlane() ) {
  case ViewState::X:
    minZ = min( volumeRASRange[0], volumeRASRange[1] );
    maxZ = max( volumeRASRange[0], volumeRASRange[1] );
    incZ = fabs(increments[0]);
    break;
  case ViewState::Y:
    minZ = min( volumeRASRange[2], volumeRASRange[3] );
    maxZ = max( volumeRASRange[2], volumeRASRange[3] );
    incZ = fabs(increments[1]);
    break;
  case ViewState::Z:
    minZ = min( volumeRASRange[4], volumeRASRange[5] );
    maxZ = max( volumeRASRange[4], volumeRASRange[5] );
    incZ = fabs(increments[2]);
    break;
  }

  // Create a dummy location, we'll change it soon. Note to self:
  // learn how to use C++ references properly.
  RAS[0] = RAS[1] = RAS[2] = 0;
  VolumeLocation 
    loc( mVolume->MakeVolumeLocationFromRAS( RAS, mCurrentFrame ) );

  // Start a progress bar for the MIP since it takes a while.
  ProgressDisplayManager& progMgr = ProgressDisplayManager::GetManager();
  list<string> lButtons;
  lButtons.push_back( "Stop" );
  progMgr.NewTask( "Building MIP",
                   "Building the maximum intensity projection",
                   true, lButtons );

  // For each of the z values in our range...
  bool bCancel = false;
  for ( float z = minZ; z <= maxZ && !bCancel; z += incZ ) {

    // Update the progress bar.
    progMgr.UpdateTask( "Building the maximum intensity projection",
                        ((z - minZ) / (maxZ - minZ)) * 100.0 );

    // Calculate the vector adjustment for the point for this
    // iteration by scaling our plane normal.
    Point3<float> adjust;
    iViewState.GetPlaneNormal( adjust.xyz() );
    adjust = adjust * z;

    // Find the increments. To do this, take two window coords on this
    // row, one at col 0 and one at col 1. Then we convert both to RAS
    // and find the difference. We use that difference to increment
    // RAS coords across this row for all columns.
    for ( int nRow = 0; nRow < iHeight; nRow++ ) {

      window[0] = 0;
      window[1] = nRow;
      iTranslator.TranslateWindowToRAS( window, mRowStartRAS[nRow] );

      // Adjust the coord by the vector we calced earlier.
      Point3<float> coord;
      coord.Set( mRowStartRAS[nRow] );
      coord = coord + Point3<float>(adjust);
      mRowStartRAS[nRow][0] = coord[0];
      mRowStartRAS[nRow][1] = coord[1];
      mRowStartRAS[nRow][2] = coord[2];

      window2[0] = 1;
      window2[1] = nRow;
      iTranslator.TranslateWindowToRAS( window2, RAS2 );

      // Adjust the coord by the vector we calced earlier.
      coord.Set( RAS2 );
      coord = coord + Point3<float>(adjust);
      RAS2[0] = coord[0];
      RAS2[1] = coord[1];
      RAS2[2] = coord[2];

      mColIncrementRAS[nRow][0] = RAS2[0] - mRowStartRAS[nRow][0];
      mColIncrementRAS[nRow][1] = RAS2[1] - mRowStartRAS[nRow][1];
      mColIncrementRAS[nRow][2] = RAS2[2] - mRowStartRAS[nRow][2];
    }

    // Get the update bounds. We'll only iterate over those.
    int windowUpdateBounds[4];
    iViewState.CopyUpdateRect( windowUpdateBounds );

    for ( window[1] = windowUpdateBounds[1];
          window[1] <= windowUpdateBounds[3]; window[1]++ ) {

      // Check for cancel.
      int nButton = progMgr.CheckTaskForButton();
      if ( nButton == 0 ) {
        bCancel = true;
        break;
      }

      // Grab the RAS beginning for this row and column.
      RAS[0] = mRowStartRAS[window[1]][0] +
               windowUpdateBounds[0]*mColIncrementRAS[window[1]][0];
      RAS[1] = mRowStartRAS[window[1]][1] +
               windowUpdateBounds[0]*mColIncrementRAS[window[1]][1];
      RAS[2] = mRowStartRAS[window[1]][2] +
               windowUpdateBounds[0]*mColIncrementRAS[window[1]][2];

      // Find this row start.
      dest = iBuffer +
             (((iWidth * window[1]) + windowUpdateBounds[0]) * mBytesPerPixel);

      for ( window[0] = windowUpdateBounds[0];
            window[0] <= windowUpdateBounds[2]; window[0]++ ) {

        // Set the location from this RAS.
        loc.SetFromRAS( RAS );

        int selectColor[3];
        if ( mVolume->IsInBounds( loc ) ) {

          switch ( mSampleMethod ) {
          case nearest:
            value = mVolume->GetMRINearestValue( loc );
            break;
          case trilinear:
            value = mVolume->GetMRITrilinearValue( loc );
            break;
          case sinc:
            value = mVolume->GetMRISincValue( loc );
            break;
          case magnitude:
            value = mVolume->GetMRIMagnitudeValue( loc );
            break;
          }

          switch ( mColorMapMethod ) {
          case grayscale:
            GetGrayscaleColorForValue( value, dest, color );
            break;
          case heatScale:
            GetHeatscaleColorForValue( value, dest, color );
            break;
          case LUT:
            GetColorLUTColorForValue( value, dest, color );
            break;
          }

          int newColor = mColorTimesOneMinusOpacity[dest[0]] +
                         mColorTimesOpacity[color[0]];
          if ( newColor > dest[0] ) {
            dest[0] = newColor;
          }
          newColor = mColorTimesOneMinusOpacity[dest[1]] +
                     mColorTimesOpacity[color[1]];
          if ( newColor > dest[1] ) {
            dest[1] = newColor;
          }
          newColor = mColorTimesOneMinusOpacity[dest[2]] +
                     mColorTimesOpacity[color[2]];
          if ( newColor > dest[2] ) {
            dest[2] = newColor;
          }

          if ( mVolume->IsSelected( loc, selectColor ) ) {
            // Write the RGB value to the buffer.
            dest[0] = (GLubyte) (((float)dest[0] * (1.0 - mROIOpacity)) +
                                 ((float)selectColor[0] * mROIOpacity));
            dest[1] = (GLubyte) (((float)dest[1] * (1.0 - mROIOpacity)) +
                                 ((float)selectColor[1] * mROIOpacity));
            dest[2] = (GLubyte) (((float)dest[2] * (1.0 - mROIOpacity)) +
                                 ((float)selectColor[2] * mROIOpacity));
          }

        }

        // Increment the dest buffer pointer.
        dest += mBytesPerPixel;

        // Increment the RAS point.
        RAS[0] += mColIncrementRAS[window[1]][0];
        RAS[1] += mColIncrementRAS[window[1]][1];
        RAS[2] += mColIncrementRAS[window[1]][2];
      }
    }
  }

  progMgr.EndTask();

  if ( mbDrawEditingLine ) {

    Point2<int> windowA, windowB;
    int color[3] = { 255, 0, 255 };
    iTranslator.TranslateRASToWindow( mLastMouseUpRAS.xyz(), windowA.xy() );
    iTranslator.TranslateRASToWindow( mCurrentMouseRAS.xyz(), windowB.xy() );
    DrawLineIntoBuffer( iBuffer, iWidth, iHeight,
                        windowA.xy(), windowB.xy(), color, 2, 0.5 );
  }

}

void
ScubaLayer2DMRI::InitBufferCoordCache ( int iWidth, int iHeight ) {

  // Delete existing buffer if necessary.
  if ( mRowStartRAS != NULL ) {
    for ( int nRow = 0; nRow < mBufferIncSize[1]; nRow++ )
      free( mRowStartRAS[nRow] );
    free( mRowStartRAS );
  }

  if ( mColIncrementRAS != NULL ) {
    for ( int nRow = 0; nRow < mBufferIncSize[1]; nRow++ )
      free( mColIncrementRAS[nRow] );
    free( mColIncrementRAS );
  }

  // Init new buffer arrays.
  mRowStartRAS =     (float**) calloc(iHeight, sizeof(float*));
  mColIncrementRAS = (float**) calloc(iHeight, sizeof(float*));
  for ( int nRow = 0; nRow < iHeight; nRow++ ) {
    mRowStartRAS[nRow] =     (float*) calloc( 3, sizeof(float) );
    mColIncrementRAS[nRow] = (float*) calloc( 3, sizeof(float) );
  }

  // Save the size of our new buffer.
  mBufferIncSize[0] = iWidth;
  mBufferIncSize[1] = iHeight;

}


void
ScubaLayer2DMRI::SetOpacity ( float iOpacity ) {

  // Call superclass first.
  Layer::SetOpacity( iOpacity );

  InitColorOpacityCache();
}

void
ScubaLayer2DMRI::InitColorOpacityCache () {

  // We calculate values from 0 to 255 for drawing into the
  // buffer. Here is our cache of the values affected by
  // opacity. We'll use them in the DrawIntoBuffer functions.
  for ( int i = 0; i < 256; i++ ) {
    mColorTimesOpacity[i] = (GLubyte)( (float)i * mOpacity );
    mColorTimesOneMinusOpacity[i] = (GLubyte)( (float)i * (1.0 - mOpacity) );
  }
}

void
ScubaLayer2DMRI::GetGrayscaleColorForValue ( float iValue,GLubyte* const iBase,
    int* oColor ) {

  if ( (!mbClearZero || (mbClearZero && iValue != 0)) &&
       (iValue >= mMinVisibleValue && iValue <= mMaxVisibleValue) &&
       (fabs(mMinVisibleValue - mMaxVisibleValue) > 0.0001) ) {

    int nLUT = (int) floor( (cGrayscaleLUTEntries-1) *
                            ((iValue - mMinVisibleValue) /
                             (mMaxVisibleValue - mMinVisibleValue)) );

    oColor[0] = mGrayscaleLUT[nLUT];
    oColor[1] = mGrayscaleLUT[nLUT];
    oColor[2] = mGrayscaleLUT[nLUT];

  } else {

    oColor[0] = (int)iBase[0];
    oColor[1] = (int)iBase[1];
    oColor[2] = (int)iBase[2];
  }
}

void
ScubaLayer2DMRI::GetHeatscaleColorForValue ( float iValue,GLubyte* const iBase,
    int* oColor ) {

  // Apply the offset.
  iValue -= mHeatScaleOffset;

  if ( (!mbClearZero || (mbClearZero && iValue != 0)) &&
       (iValue >= (mMinVisibleValue-mHeatScaleOffset) && 
	iValue <= (mMaxVisibleValue-mHeatScaleOffset)) &&
       (fabs(mMinVisibleValue - mMaxVisibleValue) > 0.0001) &&
       (iValue >= mHeatScaleMinThreshold || iValue <= -mHeatScaleMinThreshold)&&
       ((iValue > 0 && mbShowPositiveHeatScaleValues) ||
        (iValue < 0 && mbShowNegativeHeatScaleValues)) ) {

    float minValue = mHeatScaleMinThreshold;
    float midValue = mHeatScaleMidThreshold;
    float maxValue = mHeatScaleMaxThreshold;

    float fValue = iValue;

    float tmp;
    if ( fabs(fValue) > minValue &&
         fabs(fValue) < midValue ) {
      tmp = fabs(fValue);
      tmp = (1.0/(midValue-minValue)) * (tmp-minValue)*(tmp-minValue) +
            minValue;
      fValue = (fValue<0) ? -tmp : tmp;
    }

    // Reverse the value if we have that flag set.
    if ( mbReverseHeatScale )
      fValue = -fValue;

    float background_component;
    float br, bg, bb;
    float red = 0, green = 0, blue = 0;
    if ( fValue >= 0 ) {

      // First we calculate our background color. This will let us
      // blend the functional overlay color appropriately, so that a
      // value close to the minimum threshold will be a faint color
      // instead of a flat black. background_component is 1 if the
      // value is < min (so all the background is used), 1->0 from
      // min->mid, and 0 if > mid.
      background_component = (fValue<minValue) ? 1.0 :
                             (fValue<midValue) ? 1.0 - (fValue-minValue)/(midValue-minValue) : 0.0;
      br = ((float)iBase[0] / 255.0) * background_component;
      bg = ((float)iBase[1] / 255.0) * background_component;
      bb = ((float)iBase[2] / 255.0) * background_component;

      // Now calculate the overlay value color, which is red/yellow
      // for positive values and green/blue for negative.
      red   = br + ((fValue<minValue) ? 0.0 :
                    (fValue<midValue) ? (fValue-minValue)/(midValue-minValue) :
                    1.0);
      green = bg + ((fValue<midValue) ? 0.0 :
                    (fValue<maxValue) ? (fValue-midValue)/(maxValue-midValue) :
                    1.0);
      blue  = bb;

    } else {
      fValue = -fValue;

      // Do this after reversing fValue so we can compare it to the
      // positive threshold values.
      background_component = (fValue<minValue) ? 1.0 :
                             (fValue<midValue) ? 1.0 - (fValue-minValue)/(midValue-minValue) : 0.0;
      br = ((float)iBase[0] / 255.0) * background_component;
      bg = ((float)iBase[1] / 255.0) * background_component;
      bb = ((float)iBase[2] / 255.0) * background_component;

      red   = br;
      green = bg + ((fValue<midValue) ? 0.0 :
                    (fValue<maxValue) ? (fValue-midValue)/(maxValue-midValue) :
                    1.0);
      blue  = bb + ((fValue<minValue) ? 0.0 :
                    (fValue<midValue) ? (fValue-minValue)/(midValue-minValue) :
                    1.0);
    }

    if ( red > 1.0 )   red = 1.0;
    if ( green > 1.0 ) green = 1.0;
    if ( blue > 1.0 )  blue = 1.0;

    oColor[0] = (int) (red * (float)kMaxPixelComponentValue);
    oColor[1] = (int) (green * (float)kMaxPixelComponentValue);
    oColor[2] = (int) (blue * (float)kMaxPixelComponentValue);

  } else {
    oColor[0] = (int)iBase[0];
    oColor[1] = (int)iBase[1];
    oColor[2] = (int)iBase[2];
  }
}

void
ScubaLayer2DMRI::GetColorLUTColorForValue ( float iValue, GLubyte* const iBase,
    int* oColor ) {

  if ( (NULL != mColorLUT) &&
       (!mbClearZero || (mbClearZero && iValue != 0)) &&
       (iValue >= mMinVisibleValue && iValue <= mMaxVisibleValue) &&
       (fabs(mMinVisibleValue - mMaxVisibleValue) > 0.0001) ) {

    mColorLUT->GetColorAtIndex( (int)iValue, oColor );

  } else {
    oColor[0] = (int)iBase[0];
    oColor[1] = (int)iBase[1];
    oColor[2] = (int)iBase[2];
  }
}



void
ScubaLayer2DMRI::GetInfoAtRAS ( float iRAS[3],
                                std::list<InfoAtRAS>& ioInfo ) {

  if ( !mbReportInfoAtRAS )
    return;

  if ( NULL == mVolume ) {
    return;
  }

  // Look up the value of the volume at this point.
  InfoAtRAS info;
  VolumeLocation 
    loc( mVolume->MakeVolumeLocationFromRAS( iRAS, mCurrentFrame ) );
  if ( mVolume->IsInBounds( loc ) ) {

    if( GetReportInfo( kaReportableInfo[Value] ) ) {
      
      float value;
      value = mVolume->GetMRINearestValue( loc );
      
      // If this is a LUT volume, use the label from the lookup file and
      // set the shorten hint, otherwise just display the value.
      stringstream ssValue;
      if ( mColorMapMethod == LUT && NULL != mColorLUT ) {
	ssValue << mColorLUT->GetLabelAtIndex((int)value);
	info.SetShortenHint( true );
      } else {
	ssValue << value;
      }

      info.SetLabel( mVolume->GetLabel() + ",value" );
      info.SetValue( ssValue.str() );
      ioInfo.push_back( info );
      info.Clear();
    }

    if( GetReportInfo( kaReportableInfo[Index] ) ) {
      
      int index[3];
      mVolume->RASToMRIIndex( iRAS, index );
      
      stringstream ssIndex;
      ssIndex << index[0] << " " << index[1] << " " << index[2];
      
      stringstream ssCallback;
      ssCallback << "SetCursorFromVolumeIndexCoords " << mVolume->GetID();
      
      info.SetLabel( mVolume->GetLabel() + ",index" );
      info.SetInputFilter( "3ui" );
      info.SetTclCallback( ssCallback.str() );
      info.SetValue( ssIndex.str() );
      ioInfo.push_back( info );
      info.Clear();
    }
	
    if( GetReportInfo( kaReportableInfo[Talairach] ) &&
	mVolume->IsTalTransformPresent() ) {
      
      float tal[3];
      mVolume->RASToTal( iRAS, tal );
      
      stringstream ssTal;
      ssTal << setprecision(2) << tal[0] << " " << tal[1] << " " << tal[2];
      
      info.SetLabel( mVolume->GetLabel() + ",Talairach" );
      info.SetInputFilter( "3f" );
      info.SetValue( ssTal.str() );
      ioInfo.push_back( info );
      info.Clear();
    }

  } else {

    // Even if we're out of bounds, report it.
    if( GetReportInfo( "Value" ) ) {
      info.SetLabel( mVolume->GetLabel() + ",value" );
      info.SetValue( "OOB" );
      ioInfo.push_back( info );
      info.Clear();
    }
    if( GetReportInfo( "Index" ) ) {
      stringstream ssCallback;
      ssCallback << "SetCursorFromVolumeIndexCoords " << mVolume->GetID();
      
      info.SetLabel( mVolume->GetLabel() + ",index" );
      info.SetValue( "OOB" );
      info.SetInputFilter( "3ui" );
      info.SetTclCallback( ssCallback.str() );
      ioInfo.push_back( info );
      info.Clear();
    }
  }
}

TclCommandListener::TclCommandResult
ScubaLayer2DMRI::DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv ) {

  // Set2DMRILayerVolumeCollection <layerID> <collectionID>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerVolumeCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int collectionID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad collection ID";
        return error;
      }

      try {
        DataCollection& data = DataCollection::FindByID( collectionID );
        VolumeCollection& volume = (VolumeCollection&)data;
        // VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);

        SetVolumeCollection( volume );
      } catch (...) {
        sResult = "bad collection ID, collection not found";
        return error;
      }
    }
  }

  // Get2DMRILayerVolumeCollection <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerVolumeCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      if ( NULL != mVolume ) {
        ssReturnValues << (int) (mVolume->GetID());
      } else {
        ssReturnValues << -1;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Get2DMRILayerCurrentFrame <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerCurrentFrame" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetCurrentFrame();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerCurrentFrame <layerID> <frame>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerCurrentFrame" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int frame = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      SetCurrentFrame( frame );
      return ok;
    }
  }

  // Get2DMRILayerCurrentCondition <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerCurrentCondition" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetCurrentCondition();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerCurrentCondition <layerID> <condition>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerCurrentCondition" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int condition = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      SetCurrentCondition( condition );
      return ok;
    }
  }

  // Get2DMRILayerCurrentTimePoint <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerCurrentTimePoint" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetCurrentTimePoint();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerCurrentTimePoint <layerID> <timePoint>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerCurrentTimePoint" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int timePoint = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      SetCurrentTimePoint( timePoint );
      return ok;
    }
  }

  // Set2DMRILayerColorMapMethod <layerID> <method>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerColorMapMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      ColorMapMethod method;
      if ( 0 == strcmp( iasArgv[2], "grayscale" ) ) {
        method = grayscale;
      } else if ( 0 == strcmp( iasArgv[2], "heatScale" ) ) {
        method = heatScale;
      } else if ( 0 == strcmp( iasArgv[2], "lut" ) ) {
        method = LUT;
      } else {
        sResult = "bad method \"" + string(iasArgv[2]) +
                  "\", should be grayscale, heatScale or LUT";
        return error;
      }

      SetColorMapMethod( method );
    }
  }

  // Get2DMRILayerColorMapMethod <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerColorMapMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      switch ( mColorMapMethod ) {
      case grayscale:
        sReturnValues = "grayscale";
        break;
      case heatScale:
        sReturnValues = "heatScale";
        break;
      case LUT:
        sReturnValues = "lut";
        break;
      }
      sReturnFormat = "s";
    }
  }

  // Set2DMRILayerSampleMethod <layerID> <sampleMethod>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerSampleMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      SampleMethod sampleMethod;
      if ( 0 == strcmp( iasArgv[2], "nearest" ) ) {
        sampleMethod = nearest;
      } else if ( 0 == strcmp( iasArgv[2], "trilinear" ) ) {
        sampleMethod = trilinear;
      } else if ( 0 == strcmp( iasArgv[2], "sinc" ) ) {
        sampleMethod = sinc;
      } else if ( 0 == strcmp( iasArgv[2], "magnitude" ) ) {
        sampleMethod = magnitude;
      } else {
        sResult = "bad sampleMethod \"" + string(iasArgv[2]) +
                  "\", should be nearest, trilinear, or sinc";
        return error;
      }

      SetSampleMethod( sampleMethod );
    }
  }

  // Get2DMRILayerSampleMethod <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerSampleMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      switch ( mSampleMethod ) {
      case nearest:
        sReturnValues = "nearest";
        break;
      case trilinear:
        sReturnValues = "trilinear";
        break;
      case sinc:
        sReturnValues = "sinc";
        break;
      case magnitude:
        sReturnValues = "magnitude";
        break;
      }
      sReturnFormat = "s";
    }
  }

  // Set2DMRILayerBrightness <layerID> <brightness>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerBrightness" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float brightness = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad brightness";
        return error;
      }

      if ( brightness > 0 && brightness < 1 ) {
        SetBrightness( brightness );
        BuildGrayscaleLUT();
      }
    }
  }

  // Get2DMRILayerBrightness <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerBrightness" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << mBrightness;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerContrast <layerID> <contrast>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerContrast" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float contrast = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad contrast";
        return error;
      }

      if ( contrast > 0 && contrast < 30 ) {
        SetContrast( contrast );
        BuildGrayscaleLUT();
      }
    }
  }

  // Get2DMRILayerContrast <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerContrast" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << mContrast;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerWindow <layerID> <window>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerWindow" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float window = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad window";
        return error;
      }

      SetWindow( window );
    }
  }

  // Get2DMRILayerWindow <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerWindow" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetWindow();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerLevel <layerID> <level>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerLevel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float level = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad level";
        return error;
      }

      SetLevel( level );
    }
  }

  // Get2DMRILayerLevel <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerLevel" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetLevel();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerColorLUT <layerID> <lutID>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerColorLUT" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int lutID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad lut ID";
        return error;
      }

      SetColorLUT( lutID );
    }
  }

  // Get2DMRILayerColorLUT <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerColorLUT" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      if ( NULL != mColorLUT ) {
        ssReturnValues << mColorLUT->GetID();
      } else {
        ssReturnValues << -1;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";

    }
  }

  // Set2DMRIDrawZeroClear <layerID> <drawClear>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerDrawZeroClear" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        mbClearZero = true;
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
                  0 == strcmp( iasArgv[2], "0" ) ) {
        mbClearZero = false;
      } else {
        sResult = "bad drawClear \"" + string(iasArgv[2]) +
                  "\", should be true, 1, false, or 0";
        return error;
      }
    }
  }

  // Get2DMRIDrawZeroClear <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerDrawZeroClear" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << (int)mbClearZero;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRIDrawMIP <layerID> <drawMIP>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerDrawMIP" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetDrawMIP( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
                  0 == strcmp( iasArgv[2], "0" ) ) {
        SetDrawMIP( false );
      } else {
        sResult = "bad drawMIP \"" + string(iasArgv[2]) +
                  "\", should be true, 1, false, or 0";
        return error;
      }
    }
  }

  // Get2DMRIDrawMIP <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerDrawMIP" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << (int)GetDrawMIP();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Get2DMRILayerMinVisibleValue <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerMinVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetMinVisibleValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerMinVisibleValue <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerMinVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetMinVisibleValue( value );
      BuildGrayscaleLUT();
    }
  }

  // Get2DMRILayerMaxVisibleValue <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerMaxVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetMaxVisibleValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerMaxVisibleValue <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerMaxVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetMaxVisibleValue( value );
      BuildGrayscaleLUT();
    }
  }

  // Get2DMRILayerMinValue <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerMinValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << mVolume->GetMRIMinValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Get2DMRILayerMaxValue <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerMaxValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << mVolume->GetMRIMaxValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Get2DMRILayerHeatScaleMin <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerHeatScaleMin" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetHeatScaleMinThreshold();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerHeatScaleMin <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerHeatScaleMin" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetHeatScaleMinThreshold( value );
    }
  }

  // Get2DMRILayerHeatScaleMid <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerHeatScaleMid" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetHeatScaleMidThreshold();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerHeatScaleMid <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerHeatScaleMid" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetHeatScaleMidThreshold( value );
    }
  }

  // Get2DMRILayerHeatScaleMax <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerHeatScaleMax" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetHeatScaleMaxThreshold();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerHeatScaleMax <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerHeatScaleMax" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetHeatScaleMaxThreshold( value );
    }
  }

  // Get2DMRILayerHeatScaleOffset <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerHeatScaleOffset" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetHeatScaleOffset();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerHeatScaleOffset <layerID> <value>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerHeatScaleOffset" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad value";
        return error;
      }

      SetHeatScaleOffset( value );
    }
  }

  // Get2DMRILayerReverseHeatScale <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerReverseHeatScale" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      bool bReverse = GetReverseHeatScale();
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( bReverse );
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerReverseHeatScale <layerID> <reverse>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerReverseHeatScale" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bReverse =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetReverseHeatScale( bReverse );
      } catch ( runtime_error& e ) {
        sResult = "bad reverse \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // Get2DMRILayerShowPositiveHeatScaleValues <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerShowPositiveHeatScaleValues" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      bool bShow = GetShowPositiveHeatScaleValues();
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( bShow );
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerShowPositiveHeatScaleValues <layerID> <show>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerShowPositiveHeatScaleValues" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bShow =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetShowPositiveHeatScaleValues( bShow );
      } catch ( runtime_error& e ) {
        sResult = "bad show \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // Get2DMRILayerShowNegativeHeatScaleValues <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerShowNegativeHeatScaleValues" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      bool bShow = GetShowNegativeHeatScaleValues();
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( bShow );
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerShowNegativeHeatScaleValues <layerID> <show>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerShowNegativeHeatScaleValues" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bShow =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetShowNegativeHeatScaleValues( bShow );
      } catch ( runtime_error& e ) {
        sResult = "bad show \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // Set2DMRILayerHeatScaleUsingFDR <layerID> <rate> <maskVolumeName>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerHeatScaleUsingFDR" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      float rate;
      try {
        rate = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = "bad rate \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }

      // This will be an empty string if no mask.
      string fnMaskVolume =  iasArgv[3];

      SetHeatScaleThresholdUsingFDR( rate, fnMaskVolume );
    }
  }

  // Set2DMRILayerMaskVolume <layerID> <volumdID>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerMaskVolume" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int collectionID;
      try {
        collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad volumeID: ") + e.what();
        return error;
      }

      try {
        DataCollection& data = DataCollection::FindByID( collectionID );
        VolumeCollection& volume = (VolumeCollection&)data;
        // VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);

        SetMaskVolume( volume );
      } catch (...) {
        sResult = "bad collection ID, collection not found";
        return error;
      }
    }
  }

  // Get2DMRILayerMaskVolume <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerMaskVolume" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << GetMaskVolume();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Remove2DMRILayerMaskVolume <layerID>
  if ( 0 == strcmp( isCommand, "Remove2DMRILayerMaskVolume" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      RemoveMaskVolume();
    }
  }

  // Get2DMRILayerROIOpacity <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerROIOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetROIOpacity();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerROIOpacity <layerID> <opacity>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerROIOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      float opacity = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad opacity";
        return error;
      }

      SetROIOpacity( opacity );
    }
  }

  // Set2DMRILayerEditableROI <layerID> <editable>
  if ( 0 == strcmp( isCommand, "Set2DMRILayerEditableROI" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        mbEditableROI =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = "bad editable \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // Get2DMRILayerEditableROI <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRILayerEditableROI" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mbEditableROI );
      sReturnFormat = "i";
    }
  }

  // Get2DMRIRASCoordsFromIndex <layerID x y z>
  if ( 0 == strcmp( isCommand, "Get2DMRIRASCoordsFromIndex" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int index[3];
      index[0] = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      index[1] = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      index[2] = TclCommandManager::ConvertArgumentToInt( iasArgv[4] );

      float ras[3];
      mVolume->MRIIndexToRAS( index, ras );

      stringstream ssReturnValues;
      ssReturnValues << ras[0] << " " << ras[1] << " " << ras[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }

  // Flood2DMRISVolume <layerID x y z toolID type>
  if ( 0 == strcmp( isCommand, "Flood2DMRIVolume" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      float ras[3];
      try {
        ras[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
        ras[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
        ras[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );
      } catch ( runtime_error& e ) {
        sResult = string("bad RAS coord: ") + e.what();
        return error;
      }

      int toolID;
      try {
        toolID = TclCommandManager::ConvertArgumentToInt( iasArgv[5] );
        ScubaToolState::FindByID( toolID );
      } catch ( runtime_error& e ) {
        sResult = string("bad toolID: ") + e.what();
        return error;
      }

      int viewID;
      try {
        viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[6] );
        View::FindByID( viewID );
      } catch ( runtime_error& e ) {
        sResult = string("bad viewID: ") + e.what();
        return error;
      }

      // Get the view and tool and use them to set the params.
      View& genericView = View::FindByID( viewID );
      ScubaView& view = (ScubaView&) genericView;
      ScubaToolState& tool = ScubaToolState::FindByID( toolID );
      VolumeCollectionFlooder::Params params;
      SetFloodParams( tool, view.GetViewState(), params );

      // Make sure we have a source volume.
      if ( tool.GetFloodSourceCollection() < 0 ) {
        sResult = "Specify a fill source data collection.";
        return error;
      }

      // Make the right kind of flooder.
      VolumeCollectionFlooder* flooder = NULL;
      if ( 0 == strcmp( iasArgv[7], "voxelEditingNew" ) ) {
        flooder = new ScubaLayer2DMRIFloodVoxelEdit( tool.GetNewValue() );
      } else if ( 0 == strcmp( iasArgv[7], "voxelEditingErase" ) ) {
        flooder = new ScubaLayer2DMRIFloodVoxelEdit( tool.GetEraseValue() );
      } else if ( 0 == strcmp( iasArgv[7], "roiEditingSelect" ) ) {
        flooder = new ScubaLayer2DMRIFloodSelect( true );
      } else if ( 0 == strcmp( iasArgv[7], "roiEditingUnselect" ) ) {
        flooder = new ScubaLayer2DMRIFloodSelect( false );
      } else {
        sResult = "bad type \"" + string(iasArgv[7]) +
                  "\", should be voxelEditingNew, voxelEditingErase, "
                  "roiEditingSelect, roiEditingUnselect.";
        return error;
      }

      // Do the flood.
      flooder->Flood( *mVolume, ras, mCurrentFrame, params );
      delete flooder;

      return ok;
    }
  }

  return Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
ScubaLayer2DMRI::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "DataDeleted" ) {
    mVolume = NULL;
  }

  return Layer::DoListenToMessage( isMessage, iData );
}

void
ScubaLayer2DMRI::DataChanged () {

  // Our data has changed. We want to check the new min and max and
  // see if they've changed, so we can update our min/max visible
  // levels and window/level ranges.
  float newMinValue, newMaxValue;
  newMinValue = mVolume->GetMRIMinValue();
  newMaxValue = mVolume->GetMRIMaxValue();
  if ( newMinValue < mOldMinValue ||
       newMaxValue > mOldMaxValue ) {

    // Set the min and max visible value to include this new range, if
    // necessary.
    if ( newMinValue < mOldMinValue ) {
      SetMinVisibleValue( newMinValue );
    }
    if ( newMaxValue > mOldMaxValue ) {
      SetMaxVisibleValue( newMaxValue );
    }

    try {
      stringstream ssCommand;
      ssCommand << "LayerSettingsChanged " << GetID();
      TclCommandManager& mgr = TclCommandManager::GetManager();
      mgr.SendCommand( ssCommand.str() );
    } catch (...) {}

    // Rebuild our grayscale table.
    BuildGrayscaleLUT();
  }

  // Save these values.
  mOldMinValue = newMinValue;
  mOldMaxValue = newMaxValue;

  RequestRedisplay();
}

void
ScubaLayer2DMRI::HandleTool ( float iRAS[3], ViewState& iViewState,
                              ScubaWindowToRASTranslator& iTranslator,
                              ScubaToolState& iTool, InputState& iInput ) {

  // Only do this if we're the target layer.
  if ( iTool.GetTargetLayer() != GetID() )
    return;

#if 0
  if ( iInput.IsButtonUpEvent() ) {
    VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( iRAS ) );
    cerr << "Clicked " << Point3<float>(iRAS) << " -> "
    << Point3<int>(loc.Index()) << ", "
    << Point3<float>(loc.IndexF()) << endl;
    Point3<int> MRIIdx( loc.Index() );
    mVolume->PrintVoxelCornerCoords( cerr, MRIIdx );
    Point3<int> MRIIdx2( (int)loc.IndexF()[0], (int)loc.IndexF()[1],
                         (int)loc.IndexF()[2] );
    cerr << "And " << MRIIdx2 << endl;
    mVolume->PrintVoxelCornerCoords( cerr, MRIIdx2 );

  }
#endif

  // If this is our contrast/brightness key combo..
  if ( iInput.IsShiftKeyDown() &&
       iInput.Button () == 1 ) {

    if ( iInput.IsButtonDownEvent() ) {

      mOriginalLevel = GetLevel();
      mOriginalWindow = GetWindow();

    } else if ( iInput.IsButtonDragEvent() ) {

      float deltaLevel =
        ( (float)iInput.GetTotalButtonDeltaX() *
          (mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue()) ) /
        (float)iViewState.GetBufferWidth();

      float deltaWindow =
        ( (float)iInput.GetTotalButtonDeltaY() *
          (mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue()) ) /
        (float)iViewState.GetBufferHeight();

      try {
        SetLevel( mOriginalLevel + deltaLevel );
      } catch (...) {}
      try {
        SetWindow( mOriginalWindow + deltaWindow );
      } catch (...) {}
      RequestRedisplay();

      try {
        stringstream ssCommand;
        ssCommand << "LayerSettingsChanged " << GetID();
        TclCommandManager& mgr = TclCommandManager::GetManager();
        mgr.SendCommand( ssCommand.str() );
      } catch (...) {}

    }

    return;
  }


  switch ( iTool.GetMode() ) {
  case ScubaToolState::voxelEditing:
  case ScubaToolState::roiEditing: {

      // If roiEditing, check editable ROI flag.
      if ( ScubaToolState::roiEditing == iTool.GetMode() &&
           !mbEditableROI )
        return;

      // These tools have basically three sets of bahavior. First is
      // straight up painting. In this case, on a mouse down, we do undo
      // begin stuff. On a mouse drag, we get points in the brush shape
      // and paint/select them. On mouse up, we do end undo stuff.
      bool bBrush = false, bLine = false, bEyedropper = false;
      if ( !iInput.IsShiftKeyDown() &&
           !iInput.IsControlKeyDown() &&
           (2 == iInput.Button() || 3 == iInput.Button()) )
        bBrush = true;

      // Second is the shift-line. In this case, on a mouse move with
      // the shift key down, we draw we set flags to draw a line from
      // the last mouse down to the current lock. On a mouse down, we do
      // undo begin stuff, and on mouse up, paint/select the line and
      // end undo stuff.
      if ( iInput.IsShiftKeyDown() &&
           !iInput.IsControlKeyDown() )
        bLine = true;

      // Finally there is just a straight eyedropper tool for shift+ctrl
      // mouse up.
      if ( iInput.IsShiftKeyDown() &&
           iInput.IsControlKeyDown() &&
           2 == iInput.Button() )
        bEyedropper = true;


      // Eyedropper the color and set the tool.
      if ( bEyedropper && iInput.IsButtonUpEvent() ) {

        VolumeLocation 
	  loc( mVolume->MakeVolumeLocationFromRAS( iRAS, mCurrentFrame ) );
        if ( mVolume->IsInBounds( loc ) ) {

          float value = mVolume->GetMRINearestValue( loc );
          iTool.SetNewValue( value );

          stringstream ssCommand;
          ssCommand << "ToolSettingsChanged " << iTool.GetID();
          TclCommandManager& mgr = TclCommandManager::GetManager();
          mgr.SendCommand( ssCommand.str() );
        }
      }

      // Line tool. This will tell the DrawIntoBuffer to draw line from
      // last point to here.
      if ( bLine ) {

        mCurrentMouseRAS.Set( iRAS );
        mbDrawEditingLine = true;
        RequestRedisplay();

      } else {

        if ( mbDrawEditingLine ) {
          mbDrawEditingLine = false;
          RequestRedisplay();
        }
      }


      // Handle the brushing/selecting stuff here. The only difference
      // is that they find their points differently, and work on
      // different events.
      if ( (bBrush || bLine) &&
           (2 == iInput.Button() || 3 == iInput.Button()) ) {

        // If this is a mouse down event, open up an undo action and
        // enter batch mode.
        UndoManager& undoList = UndoManager::GetManager();

        if ( iInput.IsButtonDownEvent() ) {
          if ( ScubaToolState::voxelEditing == iTool.GetMode() ) {
            if ( iInput.Button() == 2 ) {
              mCurrentDrawingOperationActionID = 
		undoList.BeginAction( "Edit Voxel" );
            } else if ( iInput.Button() == 3 ) {
              mCurrentDrawingOperationActionID = 
		undoList.BeginAction( "Erase Voxel" );
            }
          } else if ( ScubaToolState::roiEditing == iTool.GetMode() ) {
            if ( iInput.Button() == 2 ) {
              mCurrentDrawingOperationActionID =
		undoList.BeginAction( "Selection Brush" );
            } else if ( iInput.Button() == 3 ) {
              mCurrentDrawingOperationActionID = 
		undoList.BeginAction( "Unselection Brush" );
            }
          }

          // Begin batch change mode.
          switch ( iTool.GetMode() ) {
          case ScubaToolState::voxelEditing:
            mVolume->BeginBatchChanges();
            break;
          case ScubaToolState::roiEditing:
            mVolume->BeginBatchROIChanges();
            break;
          default:
            break;
          }
        }

        // We'll find voxels. In brush mode, it's only for button
        // drag and up. For line, it's only on mouse up.
        if ( (bBrush &&
              (iInput.IsButtonDragEvent() || iInput.IsButtonUpEvent()))  ||
             (bLine && iInput.IsButtonUpEvent()) ) {

          list<Point3<float> > points;

          if ( bBrush ) {

            // Brush. Find the radius according to the brush shape. Calc
            // a rect for with that radius around the clicked point.
            float radiusRAS;
            ScubaToolState::Shape shape = iTool.GetBrushShape();
            if ( ScubaToolState::voxel == shape ) {
              radiusRAS = MAX( MAX( mVolume->GetVoxelXSize(),
                                    mVolume->GetVoxelYSize() ),
                               mVolume->GetVoxelZSize() );
            } else {
              radiusRAS = iTool.GetBrushRadius();
            }

            // Calc a square in the view plane.
            float squareRAS[4][3];
            CalcRASSquareInViewPlane( iRAS, radiusRAS, iTranslator, iViewState,
                                      squareRAS );

            // Now get the RAS points in this square or circle.
            switch ( shape ) {
            case ScubaToolState::voxel:
              points.push_back( Point3<float>(iRAS) );
              break;
            case ScubaToolState::square: {
                mVolume->FindRASPointsInSquare( iRAS,
                                                squareRAS[0], squareRAS[1],
                                                squareRAS[2], squareRAS[3],
                                                0,
                                                points );
              }
              break;
            case ScubaToolState::circle:
              mVolume->FindRASPointsInCircle( squareRAS[0], squareRAS[1],
                                              squareRAS[2], squareRAS[3],
                                              0, iRAS, radiusRAS,
                                              points );
              break;
            }

            // Go ahead and add our update square now while we have the info.
            CalcAndAddUpdateSquare( iRAS, radiusRAS, iTranslator, iViewState );


          } else {

            // Line. Just find the points on the segment between last
            // mouse up and where we clicked here.
            mVolume->FindRASPointsOnSegment( mLastMouseUpRAS.xyz(), iRAS,
                                             points );
          }


          // For each point we got...
          list<Point3<float> >::iterator tPoints;
          for ( tPoints = points.begin(); tPoints != points.end(); ++tPoints ) {

            // If the point is in bounds...
            Point3<float> point = *tPoints;
            VolumeLocation 
	      loc( mVolume->MakeVolumeLocationFromRAS( point.xyz(),
						       mCurrentFrame ));

            if ( mVolume->IsInBounds( loc ) ) {

              // Depending on whether we're editing or selecting, and
              // whether we're actioning or unactioning, perform the
              // proper action. Also create the proper undo item.
              UndoAction* action = NULL;
              if ( ScubaToolState::voxelEditing == iTool.GetMode() ) {

                // Editing. Get the original value for the undo item.
                float origValue = mVolume->GetMRINearestValue( loc );

                // If only brushing zero and are brushing, skip if not
                // zero.
                if ( iTool.GetOnlyBrushZero() &&
                     iInput.Button() == 2 &&
                     origValue != 0 )
                  continue;

                // If we're using a threhold and we're not in LUT mode,
                // and this value doesn't fall in it, skip.
                if ( iTool.GetUseEditThreshold() &&
                     mColorMapMethod != LUT ) {
                  if ( iInput.Button() == 2 &&
                       (origValue < iTool.GetNewValueMinThreshold() ||
                        origValue > iTool.GetNewValueMaxThreshold() ))
                    continue;
                  if ( iInput.Button() == 3 &&
                       (origValue < iTool.GetEraseValueMinThreshold() ||
                        origValue > iTool.GetEraseValueMaxThreshold() ))
                    continue;
                }

                // New value depends on voxel button.
                float newValue = iInput.Button() == 2 ?
                                 iTool.GetNewValue() : iTool.GetEraseValue();

                // Set value and make undo item.
                mVolume->SetMRIValue( loc, newValue );
                action =
                  new UndoVoxelEditAction( mVolume, newValue, origValue,
                                           point.xyz(), mCurrentFrame );


                // Selecting.
              } else if ( ScubaToolState::roiEditing == iTool.GetMode() ) {

                // New value depends on voxel button.
                bool bSelect = iInput.Button() == 2 ? true : false;

                // Select or unselect.
                if ( iInput.Button() == 2 ) {
                  mVolume->Select( loc );
                } else if ( iInput.Button() == 3 ) {
                  mVolume->Unselect( loc );
                }

                action = new UndoSelectionAction( mVolume, bSelect, point.xyz());
              }

              // Add the undo item.
              undoList.AddAction( mCurrentDrawingOperationActionID, action );
            }
          }

          RequestRedisplay();

        }

        // If this is a mouse up event, close the undo stuff and end
        // batch mode.
        if ( iInput.IsButtonUpEvent() ) {
          undoList.EndAction( mCurrentDrawingOperationActionID );
	  mCurrentDrawingOperationActionID = -1;

          // End batch change mode.
          switch ( iTool.GetMode() ) {
          case ScubaToolState::voxelEditing:
            mVolume->EndBatchChanges();
            break;
          case ScubaToolState::roiEditing:
            mVolume->EndBatchROIChanges();
            break;
          default:
            break;
          }
        }
      }

      // Also request a redisplay if this is a mouse up.
      if ( iInput.IsButtonUpEvent() ) {
        RequestRedisplay();
      }

    }
    break;

  case ScubaToolState::voxelFilling:
  case ScubaToolState::roiFilling:

    // If roiFilling, check editable ROI flag.
    if ( ScubaToolState::roiFilling == iTool.GetMode() &&
         !mbEditableROI )
      return;

    // Eyedropper the color and set the tool.
    if ( iInput.IsShiftKeyDown() && iInput.IsControlKeyDown() &&
         iInput.IsButtonDownEvent() && 2 == iInput.Button() ) {

      VolumeLocation
	loc( mVolume->MakeVolumeLocationFromRAS( iRAS, mCurrentFrame ) );
      if ( mVolume->IsInBounds( loc ) ) {

        float value = mVolume->GetMRINearestValue( loc );
        iTool.SetNewValue( value );

        stringstream ssCommand;
        ssCommand << "ToolSettingsChanged " << iTool.GetID();
        TclCommandManager& mgr = TclCommandManager::GetManager();
        mgr.SendCommand( ssCommand.str() );
      }
    }

    // Make a flood params object and fill it out, then make a flood
    // select object, specifying select or unselect in the ctor. Then
    // run the flood object with the params.  Hack: Added a 'f' key
    // shortcut for Jean. Definitely a better way of doing this but
    // don't want to do it now.
    if ( !iInput.IsShiftKeyDown() && !iInput.IsControlKeyDown() &&
         (( iInput.IsButtonDownEvent() &&
            (2 == iInput.Button() || 3 == iInput.Button()) )      ||
          ( iInput.Key().GetKeyCode() == ScubaKeyCombo::Key_F )) ) {

      VolumeLocation
	loc( mVolume->MakeVolumeLocationFromRAS( iRAS, mCurrentFrame ) );
      if ( mVolume->IsInBounds( loc ) ) {

        VolumeCollectionFlooder::Params params;
        SetFloodParams( iTool, iViewState, params );

        // Create and run the flood object.
        if ( ScubaToolState::voxelFilling == iTool.GetMode() ) {
          if ( iInput.Button() == 2 ||
               (iInput.Key().GetKeyCode() == ScubaKeyCombo::Key_F &&
                !iInput.IsControlKeyDown()) ) {
            ScubaLayer2DMRIFloodVoxelEdit flooder( iTool.GetNewValue() );
            flooder.Flood( *mVolume, iRAS, mCurrentFrame, params );
          } else if ( iInput.Button() == 3||
                      (iInput.Key().GetKeyCode() == ScubaKeyCombo::Key_F &&
                       iInput.IsControlKeyDown()) ) {
            ScubaLayer2DMRIFloodVoxelEdit flooder(iTool.GetEraseValue());
            flooder.Flood( *mVolume, iRAS, mCurrentFrame, params );
          }
        } else if ( ScubaToolState::roiFilling == iTool.GetMode() ) {
          if ( iInput.Button() == 2||
               (iInput.Key().GetKeyCode() == ScubaKeyCombo::Key_F &&
                !iInput.IsControlKeyDown()) ) {
            ScubaLayer2DMRIFloodSelect flooder( true );
            flooder.Flood( *mVolume, iRAS, mCurrentFrame, params );
          } else if ( iInput.Button() == 3||
                      (iInput.Key().GetKeyCode() == ScubaKeyCombo::Key_F &&
                       iInput.IsControlKeyDown()) ) {
            ScubaLayer2DMRIFloodSelect flooder( false );
            flooder.Flood( *mVolume, iRAS, mCurrentFrame, params );
          }
        }

        RequestRedisplay();
      }
    }



    break;

  case ScubaToolState::straightPath:
  case ScubaToolState::edgePath: {

      PathManager& pathMgr = PathManager::GetManager();
      UndoManager& undoList = UndoManager::GetManager();

      if ( iInput.IsButtonDownEvent() ) {

        Point3<float> ras( iRAS );

        // Button down, no current path. If button 1, create a new path
        // and set the first vertex at that point. If 2 or 3, find a
        // path near the click and make it the current path. Select the
        // path we made or found.
        if ( NULL == mCurrentPath ) {

          switch ( iInput.Button() ) {
          case 1: {
            mFirstPathRAS.Set( iRAS );
            mCurrentPath = new Path<float>;
            pathMgr.ManagePath( *mCurrentPath );
            mCurrentPath->AddVertex( ras );
            mCurrentPath->MarkEndOfSegment();
            mCurrentPath->SetSelected( true );

            // Make a undo list entry.
            int actionListID = undoList.BeginAction( "New Path" );
            undoList.AddAction( actionListID, 
				new UndoNewPathAction( mCurrentPath ) );
            undoList.EndAction( actionListID );

	  } break;
          case 2:
          case 3:
            if ( !iInput.IsShiftKeyDown() ) {

              mCurrentPath = FindClosestPathInPlane( iRAS, iViewState );
              if ( mCurrentPath ) {
                mCurrentPath->SetSelected( true );
                mLastPathMoveRAS.Set( iRAS );
              }

              // If shift is down, we're going to just select this path
              // in the ROI and then clear the current path. If button
              // 2, select, if button 3, unselect.
            } else {

              Path<float>* path = FindClosestPathInPlane( iRAS, iViewState );

              if ( path ) {

                bool bSelect;
                if ( iInput.Button() == 2 )
                  bSelect = true;
                else
                  bSelect = false;

                SelectVoxelsOnPath( *path, bSelect );
                RequestRedisplay();
              }
            }

            break;
          }

          // Button down, current path. If 1, add a new vertex at this
          // click. If 2, end the path at the click. If button 3,
          // stretch the path back to the first point and end it
          // there. In these cases, unselect the path.
        } else {

          switch ( iInput.Button() ) {
          case 1:
            // Button 1, add a new vertex and segment to the path.
            mCurrentPath->AddVertex( ras );
            mCurrentPath->MarkEndOfSegment();
            break;
          case 2:
            // Button 2, add the last vertex, tell the path manager to
            // manage this path, and 'let it go' by setting the current
            // path to NULL.
            mCurrentPath->AddVertex( ras );
            mCurrentPath->SetSelected( false );
            mCurrentPath = NULL;
            break;
          case 3:
            // Button 3, same as button 2, but stretch the path back to
            // the first point first.
            mCurrentPath->AddVertex( ras );

            if ( iTool.GetMode() == ScubaToolState::straightPath ) {
              StretchPathStraight( *mCurrentPath, ras.xyz(),mFirstPathRAS.xyz());
            } else if ( iTool.GetMode() == ScubaToolState::edgePath ) {
              StretchPathAsEdge( *mCurrentPath, ras.xyz(),
                                 mFirstPathRAS.xyz(), iViewState, iTranslator,
                                 iTool.GetEdgePathEdgeBias() );
            }

            mCurrentPath->SetSelected( false );
            mCurrentPath = NULL;

            break;
          }
        }

        RequestRedisplay();

      } else if ( iInput.IsButtonDragEvent() ) {

        /* Drag event, current line. This only happens with button
        two. We should already have a current path so move that. */
        if ( 2 == iInput.Button() ) {

          if ( NULL != mCurrentPath ) {
            Point3<float> deltaRAS( iRAS[0] - mLastPathMoveRAS.x(),
                                    iRAS[1] - mLastPathMoveRAS.y(),
                                    iRAS[2] - mLastPathMoveRAS.z() );
            mCurrentPath->Move( deltaRAS );
            mLastPathMoveRAS.Set( iRAS );

            RequestRedisplay();
          }
        }

      } else if ( iInput.IsButtonUpEvent() ) {

        /* Button up, current path. If we were dragging with button 2,
        don't do anything, but if button 3, find the closest path to
        the mouse up point, and if it's still the same path we
        clicked before, delete it. Then unselect the path. */
        if ( NULL != mCurrentPath ) {

          if ( 3 == iInput.Button() ) {
            Path<float>* delPath = FindClosestPathInPlane( iRAS, iViewState );
            if ( delPath == mCurrentPath ) {

              pathMgr.UnmanagePath( *delPath );

              // The UndoAction will keep a pointer to the path until it
              // goes off the undo list. It will be deleted then.
              int actionListID = undoList.BeginAction( "Delete Path" );
              undoList.AddAction( actionListID,
				  new UndoDeletePathAction( delPath ) );
              undoList.EndAction( actionListID );
            }
          }

          if ( 3 == iInput.Button() ||
               2 == iInput.Button() ) {

            mCurrentPath->SetSelected( false );
            mCurrentPath = NULL;
          }

          RequestRedisplay();
        }

      } else {

        /* No mouse event, current path. Stretch the path out to our
        current mouse position. */
        if ( NULL != mCurrentPath ) {

          Point3<float> end = mCurrentPath->GetPointAtEndOfLastSegment();
          if ( iTool.GetMode() == ScubaToolState::straightPath ) {
            StretchPathStraight( *mCurrentPath, end.xyz(), iRAS );
          } else if ( iTool.GetMode() == ScubaToolState::edgePath ) {
            StretchPathAsEdge( *mCurrentPath, end.xyz(),
                               iRAS, iViewState, iTranslator,
                               iTool.GetEdgePathEdgeBias() );
          }

          RequestRedisplay();
        }
      }


    }
    break;
  default:
    break;
  }

  // Save the mouse up location.
  if ( iInput.IsButtonUpEvent() ) {
    mLastMouseUpRAS.Set( iRAS );
  }

}

void
ScubaLayer2DMRI::CalcRASSquareInViewPlane ( float iRAS[3], float iRadiusRAS,
    ScubaWindowToRASTranslator& iTranslator,
    ViewState& iViewState,
    float oSquare[4][3] ) {

  // Find a square centered on the point we clicked with the
  // radius of the brush radius. We get the window point, offset
  // by the radius (* zoom level), and get an RAS point from
  // that. These are the corners of our square.
  Point2<int> window;

  // Calc radius in window terms.
  int windowBrushRad = (int)( iViewState.GetZoomLevel() * iRadiusRAS );

  // Get our four plane points.
  iTranslator.TranslateRASToWindow( iRAS, window.xy() );
  window[0] -= windowBrushRad;
  window[1] -= windowBrushRad;
  iTranslator.TranslateWindowToRAS( window.xy(), oSquare[0] );
  window[0] += 2 * windowBrushRad;
  iTranslator.TranslateWindowToRAS( window.xy(), oSquare[1] );
  window[1] += 2 * windowBrushRad;
  iTranslator.TranslateWindowToRAS( window.xy(), oSquare[2] );
  window[0] -= 2* windowBrushRad;
  iTranslator.TranslateWindowToRAS( window.xy(), oSquare[3] );
}

void
ScubaLayer2DMRI::CalcAndAddUpdateSquare ( float iRAS[3], float iRadiusRAS,
    ScubaWindowToRASTranslator& iTranslator,
    ViewState& iViewState ) {

  // Calc radius in window terms.
  int windowBrushRad = (int)( iViewState.GetZoomLevel() * iRadiusRAS );

  // Now we need to calculate an update square in window
  // coordinates. When we brush, this will be the area that is
  // updated on the screen. But we want to make this at least as
  // big enough to cover whole voxels. So we'll increase the
  // square by the size of the voxel in window coords.
  Point2<int> updateSquareWindow[2];

  // Get a 0,0,0 voxel and a voxel-size voxel, convert both to
  // window, and get the difference. This is the voxel size in
  // the window.
  Point3<float> voxelRAS;
  Point3<float> curRAS;
  voxelRAS[0] = iRAS[0] + mVolume->GetVoxelXSize();
  voxelRAS[1] = iRAS[1] + mVolume->GetVoxelYSize();
  voxelRAS[2] = iRAS[2] + mVolume->GetVoxelZSize();
  curRAS.Set( iRAS );
  Point2<int> voxelWindow;
  Point2<int> curWindow;
  iTranslator.TranslateRASToWindow( voxelRAS.xyz(), voxelWindow.xy() );
  iTranslator.TranslateRASToWindow( curRAS.xyz(), curWindow.xy() );
  Point2<int> voxelSizeWindow;
  voxelSizeWindow.Set( abs( voxelWindow[0] - curWindow[0] ),
                       abs( voxelWindow[1] - curWindow[1] ) );

  Point2<int> window;
  iTranslator.TranslateRASToWindow( iRAS, window.xy() );
  updateSquareWindow[0][0] =
    window[0] - windowBrushRad - voxelSizeWindow[0];
  updateSquareWindow[0][1] =
    window[1] - windowBrushRad - voxelSizeWindow[1];
  updateSquareWindow[1][0] =
    window[0] + windowBrushRad + voxelSizeWindow[0];
  updateSquareWindow[1][1] =
    window[1] + windowBrushRad + voxelSizeWindow[1];

  iViewState.AddUpdateRect( updateSquareWindow[0][0], updateSquareWindow[0][1], updateSquareWindow[1][0], updateSquareWindow[1][1] );
}

void
ScubaLayer2DMRI::ProcessOption ( string isOption, string isValue ) {

  char sValue[1024];
  strcpy( sValue, isValue.c_str() );

  if ( 0 == isOption.compare( "colormap" ) ) {
    if ( 0 == isValue.compare( "grayscale" ) ) {
      SetColorMapMethod( grayscale );
    } else if ( 0 == isValue.compare( "heatscale" ) ) {
      SetColorMapMethod( heatScale );
    } else if ( 0 == isValue.compare( "lut" ) ) {
      SetColorMapMethod( LUT );
    } else {
      throw runtime_error( "Bad colormap value" );
    }

  } else if ( 0 == isOption.compare( "samplemethod" ) ) {
    if ( 0 == isValue.compare( "nearest" ) ) {
      SetSampleMethod( nearest );
    } else if ( 0 == isValue.compare( "trilinear" ) ) {
      SetSampleMethod( trilinear );
    } else if ( 0 == isValue.compare( "sinc" ) ) {
      SetSampleMethod( sinc );
    } else {
      throw runtime_error( "Bad colormap value" );
    }

  } else if ( 0 == isOption.compare( "lut" ) ) {
    list<int> lutIDList;
    ScubaColorLUT::GetIDList( lutIDList );
    list<int>::iterator tID;
    for ( tID = lutIDList.begin(); tID != lutIDList.end(); ++tID ) {
      int id = *tID;
      ScubaColorLUT& lut = ScubaColorLUT::FindByID( id );
      if ( 0 == isValue.compare( lut.GetLabel() ) ) {
        SetColorMapMethod( LUT );
        SetColorLUT( id );
        return;
      }
    }
    if ( ERANGE == errno ) {
      throw runtime_error( "LUT not found" );
    }

  } else if ( 0 == isOption.compare( "drawzeroclear" ) ||
              /* Alternate Bruce spellings */
              0 == isOption.compare( "drawzerosclear" ) ||
              0 == isOption.compare( "drawzeroesclear" ) ) {
    int bDrawZeroClear = strtol( sValue, (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad drawzeroclear value" );
    }
    SetDrawZeroClear( bDrawZeroClear );

  } else if ( 0 == isOption.compare( "brightness" ) ) {
    float brightness = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad brightness value" );
    }
    SetBrightness( brightness );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "contrast" ) ) {
    float contrast = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad contrast value" );
    }
    SetContrast( contrast );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "window" ) ) {
    float window = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad window value" );
    }
    SetWindow( window );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "level" ) ) {
    float level = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad level value" );
    }
    SetLevel( level );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "visiblemin" ) ) {
    float min = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad visible min value" );
    }
    SetMinVisibleValue( min );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "visiblemax" ) ) {
    float max = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad visible max value" );
    }
    SetMaxVisibleValue( max );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "visible" ) ) {
    // Value will be "min,max" so we need to separate out the three
    // values.
    vector<string> lResults;
    Utilities::SplitString( isValue, ",", lResults );
    if ( 2 != lResults.size() ) {
      throw runtime_error( "Couldn't parse two values from value string" );
    }

    // Covert each to a float.
    float min = strtol( lResults[0].c_str(), (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert first value" );
    }

    float max = strtol( lResults[1].c_str(), (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert second value" );
    }

    // Set the visible range.
    SetMinMaxVisibleValue( min, max );
    BuildGrayscaleLUT();

  } else if ( 0 == isOption.compare( "heatscalemin" ) ) {
    float min = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad heat scale min value" );
    }
    SetHeatScaleMinThreshold( min );

  } else if ( 0 == isOption.compare( "heatscalemid" ) ) {
    float mid = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad heat scale mid value" );
    }
    SetHeatScaleMidThreshold( mid );

  } else if ( 0 == isOption.compare( "heatscalemax" ) ) {
    float max = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad heat scale max value" );
    }
    SetHeatScaleMaxThreshold( max );

  } else if ( 0 == isOption.compare( "heatscaleoffset" ) ) {
    float offset = (float) strtod( sValue, (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad heat scale offset" );
    }
    SetHeatScaleOffset( offset );

  } else if ( 0 == isOption.compare( "heatscale" ) ) {
    // Value will be "min,mid,max" so we need to separate out the
    // three values.
    vector<string> lResults;
    Utilities::SplitString( isValue, ",", lResults );
    if ( 3 != lResults.size() ) {
      throw runtime_error( "Couldn't parse three values from value string" );
    }

    // Covert each to a float.
    float min = strtod( lResults[0].c_str(), (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert first value" );
    }

    float mid = strtod( lResults[1].c_str(), (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert second value" );
    }

    float max = strtod( lResults[2].c_str(), (char**)NULL );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert third value" );
    }

    // Set the heatscale thresholds.
    SetHeatScaleMinThreshold( min );
    SetHeatScaleMidThreshold( mid );
    SetHeatScaleMaxThreshold( max );

    // Current frame
  } else if ( 0 == isOption.compare( "frame" ) ) {
    int frame = (int) strtol( sValue, (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad frame number" );
    }
    SetCurrentFrame( frame );

  } else {

    return Layer::ProcessOption( isOption, isValue );
  }
}

int
ScubaLayer2DMRI::GetCurrentFrame () {
  return mCurrentFrame;
}

void
ScubaLayer2DMRI::SetCurrentFrame ( int iFrame ) {

  if ( NULL == mVolume ) {
    throw runtime_error( "No volume" );
  }

  if ( iFrame >= 0 && iFrame < mVolume->GetNumberOfFrames() ) {
    mCurrentFrame = iFrame;

    // Send the layerChanged message so the view will know to request
    // info from us again, since the value under the cursor will have
    // changed.
    int id = GetID();
    SendBroadcast( "layerChanged", (void*)&id );

    // Draw with the new frame.
    RequestRedisplay();
  } else {
    throw runtime_error( "Invalid frame" );
  }
}

int
ScubaLayer2DMRI::GetCurrentTimePoint () {
  if ( NULL == mVolume ) {
    throw runtime_error( "No volume" );
  }

  return mVolume->ExtractTimePointFromFrame( mCurrentFrame );
}

int
ScubaLayer2DMRI::GetCurrentCondition () {
  if ( NULL == mVolume ) {
    throw runtime_error( "No volume" );
  }

  return mVolume->ExtractConditionFromFrame( mCurrentFrame );
}

void
ScubaLayer2DMRI::SetCurrentTimePoint ( int iTimePoint ) {
  if ( NULL == mVolume ) {
    throw runtime_error( "No volume" );
  }

  int condition = mVolume->ExtractConditionFromFrame( mCurrentFrame );
  int frame =
    mVolume->ConvertConditionAndTimePointToFrame( condition, iTimePoint );

  SetCurrentFrame( frame );
}

void
ScubaLayer2DMRI::SetCurrentCondition ( int iCondition ) {
  if ( NULL == mVolume ) {
    throw runtime_error( "No volume" );
  }

  int timePoint = mVolume->ExtractTimePointFromFrame( mCurrentFrame );
  int frame =
    mVolume->ConvertConditionAndTimePointToFrame( iCondition, timePoint );

  SetCurrentFrame( frame );
}

void
ScubaLayer2DMRI::SetColorMapMethod ( ColorMapMethod iMethod ) {
  mColorMapMethod = iMethod;

  /* Automatically set draw zero clear to true if we're switching to
     LUT. */
  if ( mColorMapMethod == LUT ) {
    SetDrawZeroClear( true );
  }
}

string
ScubaLayer2DMRI::GetColorMapMethodAsString () {

  switch ( mColorMapMethod ) {
  case grayscale:
    return "grayscale";
    break;
  case heatScale:
    return "heatScale";
    break;
  case LUT:
    return "lut";
    break;
  }
  return "Unknown";
}

string
ScubaLayer2DMRI::GetSampleMethodAsString () {

  switch ( mSampleMethod ) {
  case nearest:
    return "nearest";
    break;
  case trilinear:
    return "trilinear";
    break;
  case sinc:
    return "sinc";
    break;
  case magnitude:
    return "magnitude";
    break;
  }
  return "Unknown";
}

void
ScubaLayer2DMRI::SetColorLUT ( int iLUTID ) {

  try {
    mColorLUT = &(ScubaColorLUT::FindByID( iLUTID ));
  } catch (...) {
    DebugOutput( << "Couldn't find color LUT " << iLUTID );
  }

}

int
ScubaLayer2DMRI::GetColorLUT () {

  return mColorLUT->GetID();
}

void
ScubaLayer2DMRI::SetBrightness ( float iBrightness ) {

  if ( iBrightness < 0 || iBrightness > 1 ) {
    throw runtime_error( "Invalid brightness" );
  }
  mBrightness = iBrightness;
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::SetContrast ( float iContrast ) {

  if ( iContrast < 0 || iContrast > 30 ) {
    throw runtime_error( "Invalid contrast" );
  }
  mContrast = iContrast;
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::SetWindow ( float iWindow ) {

  mWindow = iWindow;
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::SetLevel ( float iLevel ) {

  mLevel = iLevel;
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::BuildGrayscaleLUT () {

  // Calculate the window so that it is centered on the level and
  // extends window/2 in either direction.
  float window[2];
  window[0] = mLevel - mWindow*0.5;
  window[1] = mLevel + mWindow*0.5;

#if 0
  cerr << "BuildGrayscaleLUT: " << endl
  << "\tVisible: " << mMinVisibleValue << " -> "
  << mMaxVisibleValue << endl
  << "\tWindow: " << window[0] << " -> " << window[1] << endl
  << "\tBrightness: " << mBrightness << " Contrast: " << mContrast
  << endl;
#endif

  if ( mMinVisibleValue == mMaxVisibleValue ||
       (mLevel == 0 && mWindow == 0) ) {

    // If same min and max visible values, all values should be 0.
    memset( mGrayscaleLUT, 0, sizeof(mGrayscaleLUT) );

  } else {

    for ( float nEntry = 0; nEntry < cGrayscaleLUTEntries; nEntry+=1 ) {

      // Get the value using the visible min/max to get highest
      // granularity within the 0 - cGrayscaleLUTEntries range.
      float value = ((nEntry * (mMaxVisibleValue-mMinVisibleValue)) /
                     (cGrayscaleLUTEntries-1)) + mMinVisibleValue;

      // Get an intensity from 0-1 based on our window.
      float intensity = (value - window[0]) / (window[1] - window[0]);

      // Cap the intensity.
      if ( intensity < 0.0 ) intensity = 0.0;
      if ( intensity > 1.0 ) intensity = 1.0;

      // Now apply a sigmoid function with the brightness and contrast
      // values, with the brightness adjusting the x shift and the
      // contrast adjusting the sharpness of the curve.
      float bcdValue =
        1.0 / (1.0 + exp( (intensity-mBrightness) * -mContrast) );

      // Set the value.
      mGrayscaleLUT[(int)nEntry] =
        (int) floor( bcdValue * kMaxPixelComponentValueFloat );
    }
  }
}

void
ScubaLayer2DMRI::SetMinVisibleValue ( float iValue ) {

  mMinVisibleValue = iValue;
  if ( mMinVisibleValue > mMaxVisibleValue ) {
    mMinVisibleValue = mMaxVisibleValue - 0.00001;
  }
}

void
ScubaLayer2DMRI::SetMaxVisibleValue ( float iValue ) {

  mMaxVisibleValue = iValue;
  if ( mMaxVisibleValue < mMinVisibleValue ) {
    mMaxVisibleValue = mMinVisibleValue + 0.00001;
  }
}

void
ScubaLayer2DMRI::SetMinMaxVisibleValue ( float iMinValue, float iMaxValue ) {

  mMinVisibleValue = iMinValue;
  mMaxVisibleValue = iMaxValue;
  if ( mMinVisibleValue > mMaxVisibleValue ) {
    mMinVisibleValue = mMaxVisibleValue - 0.00001;
  }
  if ( mMaxVisibleValue < mMinVisibleValue ) {
    mMaxVisibleValue = mMinVisibleValue + 0.00001;
  }
}



void
ScubaLayer2DMRI::SetHeatScaleMinThreshold ( float iValue ) {

  // Since the threshold is mirrored, it has to be positive, and the
  // range is up to fabs(min) or max, whichever is greater.
  if ( iValue < 0 || iValue > MAX(fabs(mVolume->GetMRIMinValue()),
                                  fabs(mVolume->GetMRIMaxValue())) )
    return;

  mHeatScaleMinThreshold = iValue;
  if ( mHeatScaleMinThreshold >= mHeatScaleMidThreshold ) {
    mHeatScaleMinThreshold = mHeatScaleMidThreshold - 0.0001;
  }
}

void
ScubaLayer2DMRI::SetHeatScaleMidThreshold ( float iValue ) {

  if ( iValue < 0 || iValue > MAX(fabs(mVolume->GetMRIMinValue()),
                                  fabs(mVolume->GetMRIMaxValue())) )
    return;

  mHeatScaleMidThreshold = iValue;
  if ( mHeatScaleMidThreshold <= mHeatScaleMinThreshold ) {
    mHeatScaleMidThreshold = mHeatScaleMinThreshold + 0.0001;
  }
  if ( mHeatScaleMidThreshold >= mHeatScaleMaxThreshold ) {
    mHeatScaleMidThreshold = mHeatScaleMaxThreshold - 0.0001;
  }
}

void
ScubaLayer2DMRI::SetHeatScaleMaxThreshold ( float iValue ) {

  if ( iValue < 0 || iValue > MAX(fabs(mVolume->GetMRIMinValue()),
                                  fabs(mVolume->GetMRIMaxValue())) )
    return;

  mHeatScaleMaxThreshold = iValue;
  if ( mHeatScaleMaxThreshold <= mHeatScaleMidThreshold ) {
    mHeatScaleMaxThreshold = mHeatScaleMidThreshold + 0.0001;
  }
}

void
ScubaLayer2DMRI::SetHeatScaleOffset ( float iValue ) {

  mHeatScaleOffset = iValue;
}

void
ScubaLayer2DMRI::SetHeatScaleThresholdUsingFDR ( float  iRate,
    string ifnMask ) {

  if ( NULL == mVolume ) {
    return;
  }

  MRI* maskVol = NULL;
  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();

  // If we got a mask file name, try to load it.
  if ( ifnMask != "" ) {
    try {
      maskVol = mriLoader.GetData( ifnMask );
    } catch (...) {
      // Couldn't load it, return an error.
      throw runtime_error( "Couldn't load the mask volume." );
    }
  }

  // Find the sign for the calculation based on our existing settings.
  int sign = 0;
  if ( !GetShowNegativeHeatScaleValues() )
    sign = 1;
  if ( !GetShowPositiveHeatScaleValues() )
    sign = -1;
  if ( GetReverseHeatScale() )
    sign = -sign;

  // Run the function.
  double newHeatScaleMin;
  int eMRI = MRIfdr2vwth( const_cast<MRI*>(mVolume->GetMRI()),
                          mCurrentFrame,
                          iRate,
                          sign,
                          TRUE,
                          maskVol,
                          &newHeatScaleMin,
                          NULL );

  // If we got a mask volume, delete it.
  if ( NULL != maskVol ) {
    mriLoader.ReleaseData( &maskVol );
  }

  // If we got an error, return.
  if ( ERROR_NONE != eMRI ) {
    throw runtime_error( "Error running the FDR function: check shell output for details." );
  }

  // Set our threshold. The mid is 1.5 above the min, and the max is
  // based on a slope of 0.66 --- (0.5 / slope) + mid
  SetHeatScaleMinThreshold( newHeatScaleMin );
  SetHeatScaleMidThreshold( newHeatScaleMin + 1.5 );
  SetHeatScaleMaxThreshold( (0.5 / 0.66) + newHeatScaleMin + 1.5 );
}

void
ScubaLayer2DMRI::SetMaskVolume ( VolumeCollection& iVolume ) {

  mMaskVolume = &iVolume;
}

void
ScubaLayer2DMRI::RemoveMaskVolume () {

  mMaskVolume = NULL;
}

int
ScubaLayer2DMRI::GetMaskVolume () {

  if ( NULL != mMaskVolume ) {
    return mMaskVolume->GetID();
  } else {
    return -1;
  }
}

// PATHS =================================================================

void
ScubaLayer2DMRI::StretchPathStraight ( Path<float>& iPath,
                                       float iRASBegin[3],
                                       float iRASEnd[3] ) {

  // Just clear the path and add the begin and end point.
  iPath.ClearLastSegment ();
  Point3<float> begin( iRASBegin );
  Point3<float> end( iRASEnd );
  iPath.AddVertex( begin );
  iPath.AddVertex( end );
}


void
ScubaLayer2DMRI::StretchPathAsEdge ( Path<float>& iPath,
                                     float iRASBegin[3],
                                     float iRASEnd[3],
                                     ViewState& iViewState,
                                     ScubaWindowToRASTranslator& iTranslator,
                                     float iEdgeBias ) {

  // Make an edge path finder if we don't have one already.
  if( NULL == mEdgePathFinder ) {
    mEdgePathFinder = 
      new EdgePathFinder( iViewState.GetBufferWidth(),
			  iViewState.GetBufferHeight(),
			  iTranslator, *mVolume );
  } else {
    // If we already have one, make sure the dimensions are what we
    // need, and set them if not.
    mEdgePathFinder->SetDimensions( iViewState.GetBufferWidth(),
				    iViewState.GetBufferHeight() );
  }

  // Set the bias.
  mEdgePathFinder->SetEdgeBias( iEdgeBias );

  // Get the first point from the path and the last point as passed
  // in. Convert to window points. Then find the path between them.
  Point3<float> beginRAS( iRASBegin );
  Point3<float> endRAS( iRASEnd );

  Point2<int> beginWindow;
  Point2<int> endWindow;
  iTranslator.TranslateRASToWindow( beginRAS.xyz(), beginWindow.xy() );
  iTranslator.TranslateRASToWindow( endRAS.xyz(), endWindow.xy() );
  list<Point2<int> > windowPoints;

  // Set the start point and find the path.
  mEdgePathFinder->SetStartPoint( beginWindow ),
  mEdgePathFinder->FindPathToEndPoint( endWindow, windowPoints );

  // Clear the path points and add the points we got from the path,
  // converting them to RAS one the way.
  iPath.ClearLastSegment();
  list<Point2<int> >::iterator tWindowPoint;
  for ( tWindowPoint = windowPoints.begin();
        tWindowPoint != windowPoints.end();
        ++tWindowPoint ) {
    Point3<float> currentRAS;
    iTranslator.TranslateWindowToRAS( (*tWindowPoint).xy(), currentRAS.xyz() );
    iPath.AddVertex( currentRAS );
  }
}

void
ScubaLayer2DMRI::SelectVoxelsOnPath( Path<float>& iPath, bool ibSelect ) {

  UndoManager& undoList = UndoManager::GetManager();
  int actionListID;
  if ( ibSelect )
    actionListID = undoList.BeginAction( "Select Path" );
  else
    actionListID = undoList.BeginAction( "Unselect Path" );

  int cVertices = iPath.GetNumVertices();
  for ( int nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {

    int nBackVertex = nCurVertex - 1;

    Point3<float> curVertex = iPath.GetVertexAtIndex( nCurVertex );
    Point3<float> backVertex = iPath.GetVertexAtIndex( nBackVertex );

    list<Point3<float> > rasPoints;
    mVolume->FindRASPointsOnSegment( backVertex.xyz(),
                                     curVertex.xyz(),
                                     rasPoints );

    list<Point3<float> >::iterator tRASPoint;
    for ( tRASPoint = rasPoints.begin();
          tRASPoint != rasPoints.end(); ++tRASPoint ) {
      Point3<float> rasPoint = *tRASPoint;
      VolumeLocation
	loc( mVolume->MakeVolumeLocationFromRAS( rasPoint.xyz(),
						 mCurrentFrame ) );

      if ( ibSelect )
        mVolume->Select( loc );
      else
        mVolume->Unselect( loc );

      UndoAction* action =
        new UndoSelectionAction( mVolume, ibSelect, rasPoint.xyz() );
      undoList.AddAction( actionListID, action );
    }
  }

  undoList.EndAction( actionListID );
}



Path<float>*
ScubaLayer2DMRI::FindClosestPathInPlane ( float iRAS[3],
    ViewState& iViewState ) {

  float minDistance = mWidth * mHeight;
  Path<float>* closestPath = NULL;
  Point3<float> whereRAS( iRAS );

  float range = 0;
  switch ( iViewState.GetInPlane() ) {
  case 0:
    range = mVolume->GetVoxelXSize() / 2.0;
    break;
  case 1:
    range = mVolume->GetVoxelYSize() / 2.0;
    break;
  case 2:
    range = mVolume->GetVoxelZSize() / 2.0;
    break;
  }

  PathManager& pathMgr = PathManager::GetManager();
  vector<Path<float>*>::iterator tPath;
  vector<Path<float>*>& paths = pathMgr.GetPathList();
  for ( tPath = paths.begin(); tPath != paths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    Point3<float> beginRAS = path->GetVertexAtIndex( 0 );
    if ( iViewState.IsRASVisibleInPlane( beginRAS.xyz(), range ) ) {

      float minDistanceInPath = 999999;

      int cVertices = path->GetNumVertices();
      int nCurVertex = 1;
      for ( nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {

        int nBackVertex = nCurVertex - 1;

        Point3<float> curVertex  = path->GetVertexAtIndex( nCurVertex );
        Point3<float> backVertex = path->GetVertexAtIndex( nBackVertex );

        float distance =
          Utilities::DistanceFromSegmentToPoint3f( curVertex, backVertex,
              whereRAS );

        if ( distance < minDistanceInPath ) {
          minDistanceInPath = distance;
        }
      }

      if ( minDistanceInPath < minDistance ) {
        minDistance = minDistanceInPath;
        closestPath = path;
      }
    }
  }

  return closestPath;
}

void
ScubaLayer2DMRI::DrawRASPathIntoBuffer ( GLubyte* iBuffer,
    int iWidth, int iHeight,
    int iColor[3],
    ViewState&,
    ScubaWindowToRASTranslator& iTranslator,
    Path<float>& iPath ) {


  // For every two RAS vertices on our path, translate them to window
  // points, and draw the path.
  int cVertices = iPath.GetNumVertices();
  for ( int nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {

    int nBackVertex = nCurVertex - 1;

    Point3<float> curVertex  = iPath.GetVertexAtIndex( nCurVertex );
    Point3<float> backVertex = iPath.GetVertexAtIndex( nBackVertex );

    int curWindow[2], backWindow[2];
    iTranslator.TranslateRASToWindow( curVertex.xyz(), curWindow );
    iTranslator.TranslateRASToWindow( backVertex.xyz(), backWindow );

    DrawLineIntoBuffer( iBuffer, iWidth, iHeight, backWindow, curWindow,
                        iColor, 1, 1 );
  }
}


void
ScubaLayer2DMRI::GetPreferredThroughPlaneIncrements ( float oIncrements[3] ) {

  oIncrements[0] = mVolume->GetVoxelXSize();
  oIncrements[1] = mVolume->GetVoxelYSize();
  oIncrements[2] = mVolume->GetVoxelZSize();
}

float
ScubaLayer2DMRI::GetPreferredBrushRadiusIncrement () {

  float smallestVoxelSize =
    MIN( mVolume->GetVoxelXSize(),
         MIN ( mVolume->GetVoxelYSize(), mVolume->GetVoxelZSize() ) );
  return (smallestVoxelSize / 2.0);
}

float
ScubaLayer2DMRI::GetPreferredValueIncrement () {

  return mVolume->GetPreferredValueIncrement();
}

void
ScubaLayer2DMRI::SetFloodParams ( ScubaToolState& iTool, ViewState& iViewState,
                                  VolumeCollectionFlooder::Params& ioParams ) {

  ioParams.mSourceCollection  = iTool.GetFloodSourceCollection();
  ioParams.mbStopAtPaths      = iTool.GetFloodStopAtPaths();
  ioParams.mbStopAtROIs       = iTool.GetFloodStopAtROIs();
  ioParams.mb3D               = iTool.GetFlood3D();
  ioParams.mFuzziness         = iTool.GetFloodFuzziness();
  ioParams.mMaxDistance       = iTool.GetFloodMaxDistance();
  iViewState.GetPlaneNormal( ioParams.mViewNormal );
  ioParams.mbOnlyZero         = iTool.GetOnlyFloodZero();
  ioParams.mFuzzinessType     =
    (VolumeCollectionFlooder::Params::FuzzinessType) iTool.GetFuzzinessType();
  if ( !iTool.GetFlood3D() ) {
    ioParams.mbWorkPlaneX     = (iViewState.GetInPlane() == ViewState::X);
    ioParams.mbWorkPlaneY     = (iViewState.GetInPlane() == ViewState::Y);
    ioParams.mbWorkPlaneZ     = (iViewState.GetInPlane() == ViewState::Z);
  }
}


void
ScubaLayer2DMRI::DoTimer () {

  mTimersSinceLastAutosave++;
  if ( mTimersSinceLastAutosave > kcTimersBetweenAutosaves ) {

    if ( mVolume->IsAutosaveDirty() ) {

      TclCommandManager& mgr = TclCommandManager::GetManager();
      mgr.SendCommand( "SetStatusBarText \"Autosaving...\"" );

      mVolume->AutosaveIfDirty();

      mgr.SendCommand( "SetStatusBarText \"Autosaving... done.\"" );
    }
    mTimersSinceLastAutosave = 0;
  }
}

// ======================================================================

ScubaLayer2DMRIFloodVoxelEdit::ScubaLayer2DMRIFloodVoxelEdit ( float iValue ) {
  mValue = iValue;
}


void
ScubaLayer2DMRIFloodVoxelEdit::DoBegin () {

  // Create a task in the progress display manager.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();

  list<string> lButtons;
  lButtons.push_back( "Stop" );

  manager.NewTask( "Filling", "Filling voxels", false, lButtons );

  // Start our undo action.
  UndoManager& undoList = UndoManager::GetManager();

  mActionListID = undoList.BeginAction( "Voxel Fill" );

  mVolume->BeginBatchChanges();
}

void
ScubaLayer2DMRIFloodVoxelEdit::DoEnd () {

  // End the task.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();
  manager.EndTask();

  // End our undo action.
  UndoManager& undoList = UndoManager::GetManager();
  undoList.EndAction( mActionListID );

  mVolume->EndBatchChanges();
}

bool
ScubaLayer2DMRIFloodVoxelEdit::DoStopRequested () {

  // Check for the stop button.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();
  int nButton = manager.CheckTaskForButton();
  if ( nButton == 0 ) {
    return true;
  }

  return false;
}

bool
ScubaLayer2DMRIFloodVoxelEdit::CompareVoxel ( float[3], int iFrame ) {

  // Always return true.
  return true;
}

void
ScubaLayer2DMRIFloodVoxelEdit::DoVoxel ( float iRAS[3], int iFrame ) {
  UndoManager& undoList = UndoManager::GetManager();

  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( iRAS, iFrame ) );

  // Save the original value. Set the new value.
  float origValue = mVolume->GetMRINearestValue( loc );
  mVolume->SetMRIValue( loc, mValue );

  // Make an undo item with the old value and add it to the list.
  UndoVoxelEditAction* action =
    new UndoVoxelEditAction( mVolume, mValue, origValue, iRAS, iFrame );

  undoList.AddAction( mActionListID, action );
}

UndoVoxelEditAction::UndoVoxelEditAction ( VolumeCollection* iVolume,
    float iNewValue, float iOrigValue,
    float iRAS[3], int iFrame ) {
  mVolume = iVolume;
  mNewValue = iNewValue;
  mOrigValue = iOrigValue;
  mRAS[0] = iRAS[0];
  mRAS[1] = iRAS[1];
  mRAS[2] = iRAS[2];
  mFrame = iFrame;
}

void
UndoVoxelEditAction::Undo () {

  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( mRAS, mFrame ) );
  mVolume->SetMRIValue( loc, mOrigValue );
}

void
UndoVoxelEditAction::Redo () {

  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( mRAS, mFrame ) );
  mVolume->SetMRIValue( loc, mNewValue );
}

// ============================================================

ScubaLayer2DMRIFloodSelect::ScubaLayer2DMRIFloodSelect ( bool ibSelect ) {
  mbSelect = ibSelect;
}


void
ScubaLayer2DMRIFloodSelect::DoBegin () {

  // Create a task in the progress display manager.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();

  list<string> lButtons;
  lButtons.push_back( "Stop" );

  if ( mbSelect ) {
    manager.NewTask( "Selecting", "Selecting voxels", false, lButtons );
  } else {
    manager.NewTask( "Unselecting", "Unselecting voxels", false, lButtons );
  }

  // Start our undo action.
  UndoManager& undoList = UndoManager::GetManager();

  if ( mbSelect ) {
    mActionListID = undoList.BeginAction( "Selection Fill" );
  } else {
    mActionListID = undoList.BeginAction( "Unselection Fill" );
  }

  mVolume->BeginBatchROIChanges();
}

void
ScubaLayer2DMRIFloodSelect::DoEnd () {

  // End the task.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();
  manager.EndTask();

  // End our undo action.
  UndoManager& undoList = UndoManager::GetManager();
  undoList.EndAction( mActionListID );

  mVolume->EndBatchROIChanges();
}

bool
ScubaLayer2DMRIFloodSelect::DoStopRequested () {

  // Check for the stop button.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();
  int nButton = manager.CheckTaskForButton();
  if ( nButton == 0 ) {
    return true;
  }

  return false;
}

bool
ScubaLayer2DMRIFloodSelect::CompareVoxel ( float[3], int iFrame ) {

  // Always return true.
  return true;
}

void
ScubaLayer2DMRIFloodSelect::DoVoxel ( float iRAS[3], int iFrame ) {
  UndoManager& undoList = UndoManager::GetManager();

  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( iRAS, iFrame ) );

  if ( mbSelect ) {
    mVolume->Select( loc );
    UndoSelectionAction* action =
      new UndoSelectionAction( mVolume, true, iRAS );
    undoList.AddAction( mActionListID, action );
  } else {
    mVolume->Unselect( loc );
    UndoSelectionAction* action =
      new UndoSelectionAction( mVolume, false, iRAS );
    undoList.AddAction( mActionListID, action );
  }
}

UndoSelectionAction::UndoSelectionAction ( VolumeCollection* iVolume,
    bool ibSelect, float iRAS[3] ) {
  mVolume = iVolume;
  mbSelect = ibSelect;
  mRAS[0] = iRAS[0];
  mRAS[1] = iRAS[1];
  mRAS[2] = iRAS[2];
}

void
UndoSelectionAction::Undo () {
  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( mRAS ) );

  if ( mbSelect ) {
    mVolume->Unselect( loc );
  } else {
    mVolume->Select( loc );
  }
}

void
UndoSelectionAction::Redo () {
  VolumeLocation loc( mVolume->MakeVolumeLocationFromRAS( mRAS ) );

  if ( mbSelect ) {
    mVolume->Select( loc );
  } else {
    mVolume->Unselect( loc );
  }
}


UndoPathAction::UndoPathAction ( Path<float>* iPath ) {
  mPath = iPath;
}

UndoPathAction::~UndoPathAction () {}

UndoNewPathAction::UndoNewPathAction ( Path<float>* iPath ) :
    UndoPathAction( iPath ) {}

void
UndoNewPathAction::Undo () {

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.UnmanagePath( *mPath );
}

void
UndoNewPathAction::Redo () {

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.ManagePath( *mPath );
}

UndoDeletePathAction::UndoDeletePathAction ( Path<float>* iPath ) :
    UndoPathAction( iPath ) {}

UndoDeletePathAction::~UndoDeletePathAction () {
  delete mPath;
}

void
UndoDeletePathAction::Undo () {

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.ManagePath( *mPath );
}

void
UndoDeletePathAction::Redo () {

  PathManager& pathMgr = PathManager::GetManager();
  pathMgr.UnmanagePath( *mPath );
}


// ============================================================

EdgePathFinder::EdgePathFinder ( int iViewWidth, int iViewHeight,
                                 ScubaWindowToRASTranslator& iTranslator,
                                 VolumeCollection& iVolume ) :
  mVolume( iVolume ),
  mTranslator( iTranslator ) {

  SetDimensions( iViewWidth, iViewHeight );
}

float
EdgePathFinder::GetEdgeCost ( Point2<int> const& iPoint ) const {


  // Get the magnitude value at this point.
  float cost = 0.0;
  float RAS[3];
  mTranslator.TranslateWindowToRAS( iPoint.xy(), RAS );
  VolumeLocation loc( mVolume.MakeVolumeLocationFromRAS( RAS ) );

  if ( mVolume.IsInBounds( loc ) ) {
    cost = 1.0 / (mVolume.GetMRIMagnitudeValue( loc ) + 0.0001);
  } else {
    cost = numeric_limits<float>::max();
  }
  return cost;
}
