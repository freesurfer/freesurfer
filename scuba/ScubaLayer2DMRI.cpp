#include <fstream>
#include "ScubaLayer2DMRI.h"
#include "ViewState.h"
#include "talairachex.h"
#include "ProgressDisplayManager.h"
#include "Utilities.h"
#include "Array2.h"
#include "PathManager.h"
#include "VectorOps.h"
#include "ScubaView.h"

using namespace std;

int const ScubaLayer2DMRI::cGrayscaleLUTEntries = 256;
int const ScubaLayer2DMRI::kMaxPixelComponentValue = 255;
float const ScubaLayer2DMRI::kMaxPixelComponentValueFloat = 255.0;
int const ScubaLayer2DMRI::kcTimersBetweenAutosaves = 60000;

ScubaLayer2DMRI::ScubaLayer2DMRI () {

  SetOutputStreamToCerr();

  mVolume = NULL;
  mSampleMethod = nearest;
  mColorMapMethod = grayscale;
  mBrightness = 0.25;
  mContrast = 12.0; 
  mCurrentPath = NULL;
  mROIOpacity = 0.7;
  mbEditableROI = true;
  mbClearZero = false;
  mTimersSinceLastAutosave = 0;
  mMinVisibleValue = 0;
  mMaxVisibleValue = 0;
  mBufferIncSize[0] = mBufferIncSize[1] = -1;
  mRowStartRAS = mColIncrementRAS = NULL;

  // Try setting our initial color LUT to the default LUT with
  // id 0. If it's not there, create it.
  try { 
    mColorLUT = &(ScubaColorLUT::FindByID( 0 ));
  }
  catch(...) {

    ScubaColorLUT* lut = new ScubaColorLUT();
    lut->SetLabel( "Default" );
    
    try {
      mColorLUT = &(ScubaColorLUT::FindByID( 0 ));
    }
    catch(...) {
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
}

ScubaLayer2DMRI::~ScubaLayer2DMRI () {

}

void
ScubaLayer2DMRI::SetVolumeCollection ( VolumeCollection& iVolume ) {

  mVolume = &iVolume;

  float oneTenth;
  oneTenth = (mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue()) / 10.0;

  mVolume->GetMRI();
  SetMinVisibleValue( mVolume->GetMRIMinValue() );
  SetHeatScaleMaxThreshold( mVolume->GetMRIMaxValue() - oneTenth );
  SetHeatScaleMidThreshold
    ( ((mVolume->GetMRIMaxValue() - mVolume->GetMRIMinValue()) / 2.0) + 
      mVolume->GetMRIMinValue() );
  SetHeatScaleMinThreshold( mVolume->GetMRIMinValue() + oneTenth );
  SetMaxVisibleValue( mVolume->GetMRIMaxValue() );

  mVolume->AddListener( this );

  BuildGrayscaleLUT();
}

void 
ScubaLayer2DMRI::DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				  ViewState& iViewState,
				  ScubaWindowToRASTranslator& iTranslator ) {

  if( NULL == mVolume ) {
    DebugOutput( << "No volume to draw" );
    return;
  }

  if( !mbVisible ) {
    return;
  }

  // Precalc our color * opacity.
  GLubyte aColorTimesOpacity[256];
  GLubyte aColorTimesOneMinusOpacity[256];
  for( int i = 0; i < 256; i++ ) {
    aColorTimesOpacity[i] = (GLubyte)( (float)i * mOpacity );
    aColorTimesOneMinusOpacity[i] = (GLubyte)( (float)i * (1.0 - mOpacity) );
  }

  GLubyte* dest = iBuffer;
  int window[2], window2[2];
  float RAS[3], RAS2[3];
  float value = 0;
  int color[3];

  // Init our buffers if necessary.
  if( iWidth != mBufferIncSize[0] || iHeight != mBufferIncSize[1] ) {

    // Delete existing buffer if necessary.
    if( mRowStartRAS != NULL ) {
      for( int nRow = 0; nRow < mBufferIncSize[1]; nRow++ )
	free( mRowStartRAS[nRow] );
      free( mRowStartRAS );
    }

    if( mColIncrementRAS != NULL ) {
      for( int nRow = 0; nRow < mBufferIncSize[1]; nRow++ )
	free( mColIncrementRAS[nRow] );
      free( mColIncrementRAS );
    }

    // Init new buffer arrays.
    mRowStartRAS =     (float**) calloc(iHeight, sizeof(float*));
    mColIncrementRAS = (float**) calloc(iHeight, sizeof(float*));
    for( int nRow = 0; nRow < iHeight; nRow++ ) {
      mRowStartRAS[nRow] =     (float*) calloc( 3, sizeof(float) );
      mColIncrementRAS[nRow] = (float*) calloc( 3, sizeof(float) );
    }

    // Save the size of our new buffer.
    mBufferIncSize[0] = iWidth;
    mBufferIncSize[1] = iHeight;
  }

  // Find the increments.
  for( int nRow = 0; nRow < iHeight; nRow++ ) {

    window[0] = 0; window[1] = nRow;
    iTranslator.TranslateWindowToRAS( window, mRowStartRAS[nRow] );
    
    window2[0] = 1; window2[1] = nRow;
    iTranslator.TranslateWindowToRAS( window2, RAS2 );
    mColIncrementRAS[nRow][0] = RAS2[0] - mRowStartRAS[nRow][0];
    mColIncrementRAS[nRow][1] = RAS2[1] - mRowStartRAS[nRow][1];
    mColIncrementRAS[nRow][2] = RAS2[2] - mRowStartRAS[nRow][2];
  }

  // Get the update bounds. We'll only iterate over those.
  int windowUpdateBounds[4];
  iViewState.CopyUpdateRect( windowUpdateBounds );

  // Create a dummy location, we'll change it soon. Note to self:
  // learn how to use C++ references properly.
  RAS[0] = RAS[1] = RAS[2] = 0;
  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( RAS );

  for( window[1] = windowUpdateBounds[1];
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

    for( window[0] = windowUpdateBounds[0]; 
	 window[0] <= windowUpdateBounds[2]; window[0]++ ) {

      // Set the location from this RAS.
      loc.SetFromRAS( RAS );

      int selectColor[3];
      if( mVolume->IsInBounds( loc ) ) {

	switch( mSampleMethod ) {
	case nearest:  value = mVolume->GetMRINearestValue( loc );   break;
	case trilinear:value = mVolume->GetMRITrilinearValue( loc ); break;
	case sinc:     value = mVolume->GetMRISincValue( loc );      break;
	case magnitude:value = mVolume->GetMRIMagnitudeValue( loc ); break;
	}
	
	switch( mColorMapMethod ) { 
	case grayscale: GetGrayscaleColorForValue( value, dest, color );break;
	case heatScale: GetHeatscaleColorForValue( value, dest, color );break;
	case LUT:       GetColorLUTColorForValue( value, dest, color ); break;
	}
	
	dest[0] = aColorTimesOneMinusOpacity[dest[0]] + 
	  aColorTimesOpacity[color[0]];
	dest[1] = aColorTimesOneMinusOpacity[dest[1]] + 
	  aColorTimesOpacity[color[1]];
	dest[2] = aColorTimesOneMinusOpacity[dest[2]] + 
	  aColorTimesOpacity[color[2]];
	
	if( mVolume->IsSelected( loc, selectColor ) ) {
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

  delete &loc;

}

void
ScubaLayer2DMRI::GetGrayscaleColorForValue ( float iValue,GLubyte* const iBase,
					     int* oColor ) {

  if( (!mbClearZero || (mbClearZero && iValue != 0)) &&
       (iValue >= mMinVisibleValue && iValue <= mMaxVisibleValue) ) {

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

  if( (!mbClearZero || (mbClearZero && iValue != 0)) &&
       (iValue >= mMinVisibleValue && iValue <= mMaxVisibleValue) ) {
    
    float minValue = mHeatScaleMinThreshold;
    float midValue = mHeatScaleMidThreshold;
    float maxValue = mHeatScaleMaxThreshold;
    
    float tmp;
    if ( fabs(iValue) > minValue &&
	 fabs(iValue) < midValue ) {
      tmp = fabs(iValue);
      tmp = (1.0/(midValue-minValue)) * (tmp-minValue)*(tmp-minValue) + 
	minValue;
      iValue = (iValue<0) ? -tmp : tmp;
    }
    
    /* calc the color */
    float red, green, blue;
    if( iValue >= 0 ) {
      red = ((iValue<minValue) ? 0.0 : 
	     (iValue<midValue) ? 
	     (iValue-minValue)/
	     (midValue-minValue) :
	     1.0);
      green = ((iValue<midValue) ? 0.0 :
	       (iValue<maxValue) ? 
	       (iValue-midValue)/
	       (maxValue-midValue) : 1.0);
      blue = 0.0; 
    } else {
      iValue = -iValue;
      red = 0.0;
      green = ((iValue<midValue) ? 0.0 :
	       (iValue<maxValue) ? 
	       (iValue-midValue)/
	       (maxValue-midValue) : 1.0);
      blue = ((iValue<minValue) ? 0.0 :
	      (iValue<midValue) ? 
	      (iValue-minValue)/
	      (midValue-minValue) : 
	      1.0);
    }
    
    if( red > 1.0 )   red = 1.0;
    if( green > 1.0 ) green = 1.0;
    if( blue > 1.0 )  blue = 1.0;
    
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
  
  if( (NULL != mColorLUT) && 
      (!mbClearZero || (mbClearZero && iValue != 0)) &&
      (iValue >= mMinVisibleValue && iValue <= mMaxVisibleValue) ) {

    mColorLUT->GetColorAtIndex( (int)iValue, oColor );

  } else {
    oColor[0] = (int)iBase[0]; 
    oColor[1] = (int)iBase[1]; 
    oColor[2] = (int)iBase[2];
  }
}
  

  
void 
ScubaLayer2DMRI::GetInfoAtRAS ( float iRAS[3], 
				map<string,string>& iLabelValues ) {

  if( NULL == mVolume ) {
    return;
  }

  // Look up the value of the volume at this point.
  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( iRAS );
  if ( mVolume->IsInBounds( loc ) ) {
    
    float value;
    value = mVolume->GetMRINearestValue( loc ); 

    // If this is a LUT volume, use the label from the lookup file,
    // otherwise just display the value.
    stringstream ssValue;
    if( mColorMapMethod == LUT && NULL != mColorLUT ) {
      ssValue << mColorLUT->GetLabelAtIndex((int)value);
    } else {
      ssValue << value;
    }

    iLabelValues[mVolume->GetLabel() + ",value"] = ssValue.str();

    int index[3];
    mVolume->RASToMRIIndex( iRAS, index );

    stringstream ssIndex;
    ssIndex << index[0] << " " << index[1] << " " << index[2];
    iLabelValues[mVolume->GetLabel() + ",index"] = ssIndex.str();
  }

  delete &loc;
}
  
TclCommandListener::TclCommandResult 
ScubaLayer2DMRI::DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv ) {

  // Set2DMRILayerVolumeCollection <layerID> <collectionID>
  if( 0 == strcmp( isCommand, "Set2DMRILayerVolumeCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      int collectionID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad collection ID";
	return error;
      }
	
      try { 
	DataCollection& data = DataCollection::FindByID( collectionID );
	VolumeCollection& volume = (VolumeCollection&)data;
	// VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);
	
	SetVolumeCollection( volume );
      }
      catch(...) {
	sResult = "bad collection ID, collection not found";
	return error;
      }
    }
  }

  // Get2DMRILayerVolumeCollection <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerVolumeCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {

      stringstream ssReturnValues;
      if( NULL != mVolume ) {
	ssReturnValues << (int) (mVolume->GetID());
      } else {
	ssReturnValues << -1;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRILayerColorMapMethod <layerID> <method>
  if( 0 == strcmp( isCommand, "Set2DMRILayerColorMapMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      ColorMapMethod method;
      if( 0 == strcmp( iasArgv[2], "grayscale" ) ) {
	method = grayscale;
      } else if( 0 == strcmp( iasArgv[2], "heatScale" ) ) {
	method = heatScale;
      } else if( 0 == strcmp( iasArgv[2], "lut" ) ) {
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
  if( 0 == strcmp( isCommand, "Get2DMRILayerColorMapMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
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
  if( 0 == strcmp( isCommand, "Set2DMRILayerSampleMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      SampleMethod sampleMethod;
      if( 0 == strcmp( iasArgv[2], "nearest" ) ) {
	sampleMethod = nearest;
      } else if( 0 == strcmp( iasArgv[2], "trilinear" ) ) {
	sampleMethod = trilinear;
      } else if( 0 == strcmp( iasArgv[2], "sinc" ) ) {
	sampleMethod = sinc;
      } else if( 0 == strcmp( iasArgv[2], "magnitude" ) ) {
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
  if( 0 == strcmp( isCommand, "Get2DMRILayerSampleMethod" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      switch( mSampleMethod ) {
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
  if( 0 == strcmp( isCommand, "Set2DMRILayerBrightness" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float brightness = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad brightness";
	return error;
      }

      if( brightness > 0 && brightness < 1 ) {
	SetBrightness( brightness );
	BuildGrayscaleLUT();
      }
    }
  }

  // Get2DMRILayerBrightness <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerBrightness" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << mBrightness;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerContrast <layerID> <contrast>
  if( 0 == strcmp( isCommand, "Set2DMRILayerContrast" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float contrast = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad contrast";
	return error;
      }

      if( contrast > 0 && contrast < 30 ) {
	SetContrast( contrast );
	BuildGrayscaleLUT();
      }
    }
  }

  // Get2DMRILayerContrast <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerContrast" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << mContrast;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // Set2DMRILayerColorLUT <layerID> <lutID>
  if( 0 == strcmp( isCommand, "Set2DMRILayerColorLUT" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      int lutID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad lut ID";
	return error;
      }
    
      SetColorLUT( lutID );
    }
  }

  // Get2DMRILayerColorLUT <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerColorLUT" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      stringstream ssReturnValues;
      if( NULL != mColorLUT ) {
	ssReturnValues << mColorLUT->GetID();
      } else {
	ssReturnValues << -1;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";

    }
  }

  // Set2DMRIDrawZeroClear <layerID> <drawClear>
  if( 0 == strcmp( isCommand, "Set2DMRILayerDrawZeroClear" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	mbClearZero = true;
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
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
  if( 0 == strcmp( isCommand, "Get2DMRILayerDrawZeroClear" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << (int)mbClearZero;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Get2DMRILayerMinVisibleValue <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerMinVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetMinVisibleValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerMinVisibleValue <layerID> <value>
  if( 0 == strcmp( isCommand, "Set2DMRILayerMinVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad value";
	return error;
      }

      SetMinVisibleValue( value );
      BuildGrayscaleLUT();
    }
  }

  // Get2DMRILayerMaxVisibleValue <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerMaxVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetMaxVisibleValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerMaxVisibleValue <layerID> <value>
  if( 0 == strcmp( isCommand, "Set2DMRILayerMaxVisibleValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float value = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad value";
	return error;
      }
      
      SetMaxVisibleValue( value );
      BuildGrayscaleLUT();
    }
  }
  
  // Get2DMRILayerMinValue <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerMinValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << mVolume->GetMRIMinValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Get2DMRILayerMaxValue <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerMaxValue" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << mVolume->GetMRIMaxValue();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Get2DMRILayerROIOpacity <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerROIOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetROIOpacity();
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRILayerROIOpacity <layerID> <opacity>
  if( 0 == strcmp( isCommand, "Set2DMRILayerROIOpacity" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      float opacity = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad opacity";
	return error;
      }
      
      SetROIOpacity( opacity );
    }
  }
  
  // Set2DMRILayerEditableROI <layerID> <editable>
  if( 0 == strcmp( isCommand, "Set2DMRILayerEditableROI" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }
    
    if( mID == layerID ) {
      
      try {
	mbEditableROI =
	  TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
      }
      catch( runtime_error e ) {
	sResult = "bad editable \"" + string(iasArgv[2]) + "\"," + e.what();
	return error;	
      }
    }
  }

  // Get2DMRILayerEditableROI <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerEditableROI" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }
    
    if( mID == layerID ) {

      sReturnValues =
	TclCommandManager::ConvertBooleanToReturnValue( mbEditableROI );
      sReturnFormat = "i";
    }
  }

  // Get2DMRIRASCoordsFromIndex <layerID x y z>
  if( 0 == strcmp( isCommand, "Get2DMRIRASCoordsFromIndex" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }
    
    if( mID == layerID ) {

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
  if( 0 == strcmp( isCommand, "Flood2DMRIVolume" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }
    
    if( mID == layerID ) {

      float ras[3];
      try { 
	ras[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
	ras[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
	ras[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );
      }
      catch( runtime_error e ) {
	sResult = string("bad RAS coord: ") + e.what();
	return error;
      }

      int toolID;
      try {
	toolID = TclCommandManager::ConvertArgumentToInt( iasArgv[5] );
	ScubaToolState::FindByID( toolID );
      }
      catch( runtime_error e ) {
	sResult = string("bad toolID: ") + e.what();
	return error;
      }
      
      int viewID;
      try {
	viewID = TclCommandManager::ConvertArgumentToInt( iasArgv[6] );
	View::FindByID( viewID );
      }
      catch( runtime_error e ) {
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
      if( tool.GetFloodSourceCollection() < 0 ) {
	sResult = "Specify a fill source data collection.";
	return error;
      }

      // Make the right kind of flooder.
      VolumeCollectionFlooder* flooder = NULL;
      if( 0 == strcmp( iasArgv[7], "voxelEditingNew" ) ) {
	flooder = new ScubaLayer2DMRIFloodVoxelEdit( tool.GetNewValue() );
      } else if( 0 == strcmp( iasArgv[7], "voxelEditingErase" ) ) {
	flooder = new ScubaLayer2DMRIFloodVoxelEdit( tool.GetEraseValue() );
      } else if( 0 == strcmp( iasArgv[7], "roiEditingSelect" ) ) {
	flooder = new ScubaLayer2DMRIFloodSelect( true );
      } else if( 0 == strcmp( iasArgv[7], "roiEditingUnselect" ) ) {
	flooder = new ScubaLayer2DMRIFloodSelect( false );
      } else {
	sResult = "bad type \"" + string(iasArgv[7]) +
	  "\", should be voxelEditingNew, voxelEditingErase, "
	  "roiEditingSelect, roiEditingUnselect.";
	return error;	
      }

      // Do the flood.
      flooder->Flood( *mVolume, ras, params );
      delete flooder;

      return ok;
    }
  }

  return Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
ScubaLayer2DMRI::DoListenToMessage ( string isMessage, void* iData ) {

  if( isMessage == "DataDeleted" ) {
    mVolume = NULL;
  }

  return Layer::DoListenToMessage( isMessage, iData );
}

void
ScubaLayer2DMRI::DataChanged () {

  float newMinValue, newMaxValue;
  newMinValue = mVolume->GetMRIMinValue();
  newMaxValue = mVolume->GetMRIMaxValue();
  if ( newMinValue < mMinVisibleValue ||
       newMaxValue > mMaxVisibleValue ) {

    SetMinVisibleValue( mVolume->GetMRIMinValue() );
    SetMaxVisibleValue( mVolume->GetMRIMaxValue() );

    stringstream ssCommand;
    ssCommand << "2DMRILayerMinMaxValueChanged " << GetID();
    TclCommandManager& mgr = TclCommandManager::GetManager();
    mgr.SendCommand( ssCommand.str() );
  }

  RequestRedisplay();
}

void
ScubaLayer2DMRI::HandleTool ( float iRAS[3], ViewState& iViewState,
			      ScubaWindowToRASTranslator& iTranslator,
			      ScubaToolState& iTool, InputState& iInput ) {
  
  // Only do this if we're the target layer.
  if( iTool.GetTargetLayer() != GetID() )
    return;

  switch( iTool.GetMode() ) {
  case ScubaToolState::voxelEditing:
  case ScubaToolState::roiEditing:

    // If roiEditing, check editable ROI flag.
    if( ScubaToolState::roiEditing == iTool.GetMode() &&
	!mbEditableROI ) 
      return;

    // If shift key is down, we're filling. Make a flood params object
    // and fill it out, then make a flood select object, specifying
    // select or unselect in the ctor. Then run the flood object with
    // the params.  Hack: Added a 'f' key shortcut for
    // Jean. Definitely a better way of doing this but don't want to
    // do it now.
    if( iInput.IsShiftKeyDown() && !iInput.IsControlKeyDown() &&
	(( iInput.IsButtonDownEvent() && 
	   (2 == iInput.Button() || 3 == iInput.Button()) )      ||
	 ( iInput.Key()[0] == 'f' )) ) {

      VolumeLocation& loc =
	(VolumeLocation&) mVolume->MakeLocationFromRAS( iRAS );
      if( mVolume->IsInBounds( loc ) ) {

	VolumeCollectionFlooder::Params params;
	SetFloodParams( iTool, iViewState, params );

	// Create and run the flood object.
	if( ScubaToolState::voxelEditing == iTool.GetMode() ) {
	  if( iInput.Button() == 2 || 
	      (iInput.Key()[0] == 'f' && !iInput.IsControlKeyDown()) ) {
	    ScubaLayer2DMRIFloodVoxelEdit flooder( iTool.GetNewValue() );
	    flooder.Flood( *mVolume, iRAS, params );
	  } else if( iInput.Button() == 3|| 
	      (iInput.Key()[0] == 'f' && iInput.IsControlKeyDown()) ) {
	    ScubaLayer2DMRIFloodVoxelEdit flooder(iTool.GetEraseValue());
	    flooder.Flood( *mVolume, iRAS, params );
	  }
	} else if( ScubaToolState::roiEditing == iTool.GetMode() ) {
	  if( iInput.Button() == 2|| 
	      (iInput.Key()[0] == 'f' && !iInput.IsControlKeyDown()) ) {
	    ScubaLayer2DMRIFloodSelect flooder( true );
	    flooder.Flood( *mVolume, iRAS, params );
	  } else if( iInput.Button() == 3|| 
	      (iInput.Key()[0] == 'f' && iInput.IsControlKeyDown()) ) {
	    ScubaLayer2DMRIFloodSelect flooder( false );
	    flooder.Flood( *mVolume, iRAS, params );
	  }
	}
	  
	RequestRedisplay();
      }      
      delete &loc;
    }

    // Eyedropper the color and set the tool.
    if( iInput.IsShiftKeyDown() && iInput.IsControlKeyDown() &&
	iInput.IsButtonDownEvent() && 2 == iInput.Button() ) {

      VolumeLocation& loc =
	(VolumeLocation&) mVolume->MakeLocationFromRAS( iRAS );
      if( mVolume->IsInBounds( loc ) ) {

	float value = mVolume->GetMRINearestValue( loc );
	iTool.SetNewValue( value );

	stringstream ssCommand;
	ssCommand << "ToolSettingsChanged " << iTool.GetID();
	TclCommandManager& mgr = TclCommandManager::GetManager();
	mgr.SendCommand( ssCommand.str() );
      }
    }

    if( (iInput.IsButtonDownEvent() || iInput.IsButtonUpEvent() ||
	 iInput.IsButtonDragEvent()) &&
	!iInput.IsShiftKeyDown() && !iInput.IsControlKeyDown() &&
	(2 == iInput.Button() || 3 == iInput.Button()) ) {

      // Otherwise we're just brushing. If this is a mouse down event,
      // open up an undo action, and if it's a mouse up event, close
      // it. Then if the mouse is down switch on the brush shape and
      // call a GetRASPointsIn{Shape} function to get the points we
      // need to brush.

      UndoManager& undoList = UndoManager::GetManager();
      if( iInput.IsButtonDownEvent() ) {
	if( ScubaToolState::voxelEditing == iTool.GetMode() ) {
	  if( iInput.Button() == 2 ) {
	    undoList.BeginAction( "Edit Voxel" );
	  } else if( iInput.Button() == 3 ) {
	    undoList.BeginAction( "Erase Voxel" );
	  }
	} else if( ScubaToolState::roiEditing == iTool.GetMode() ) {
	  if( iInput.Button() == 2 ) {
	    undoList.BeginAction( "Selection Brush" );
	  } else if( iInput.Button() == 3 ) {
	    undoList.BeginAction( "Unselection Brush" );
	  }
	}
      }

      if( iInput.IsButtonUpEvent() ) {
	undoList.EndAction();
      }

      if( iInput.IsButtonDown() ) {

	// Find a square centered on the point we clicked with the
	// radius of the brush radius. We get the window point, offset
	// by the radius (* zoom level), and get an RAS point from
	// that. These are the corners of our plane.
	Point2<int> window;
	Point3<float> sq[4];
	float rad;
	ScubaToolState::Shape shape = iTool.GetBrushShape();
	if( ScubaToolState::voxel == shape ) {
	  rad = MAX( MAX( mVolume->GetVoxelXSize(), mVolume->GetVoxelYSize() ),
		     mVolume->GetVoxelZSize() );
	} else {
	  rad = iTool.GetBrushRadius();
	}

	int windowBrushRad = (int)( iViewState.mZoomLevel * rad );

	// Get our four plane points.
	iTranslator.TranslateRASToWindow( iRAS, window.xy() );
	window[0] -= windowBrushRad;
	window[1] -= windowBrushRad;
	iTranslator.TranslateWindowToRAS( window.xy(), sq[0].xyz() );
	window[0] += 2 * windowBrushRad;
	iTranslator.TranslateWindowToRAS( window.xy(), sq[1].xyz() );
	window[1] += 2 * windowBrushRad;
	iTranslator.TranslateWindowToRAS( window.xy(), sq[2].xyz() );
	window[0] -= 2* windowBrushRad;
	iTranslator.TranslateWindowToRAS( window.xy(), sq[3].xyz() );

	// Now we need to calculate an update rect in window
	// coordinates. When we brush, this will be the area that is
	// updated on the screen. But we want to make this at least as
	// big enough to cover whole voxels. So we'll increase the
	// rect by the size of the voxel in window coords.
	Point2<int> updateRectWindow[2];

	// Get a 0,0,0 voxel and a voxel-size voxel, convert both to
	// window, and get the difference. This is the voxel size in
	// the window.
	Point3<float> voxelRAS;
	Point3<float> originRAS;
	voxelRAS[0] = mVolume->GetVoxelXSize();
	voxelRAS[1] = mVolume->GetVoxelYSize();
	voxelRAS[2] = mVolume->GetVoxelZSize();
	originRAS.Set( 0, 0, 0 );
	Point2<int> voxelWindow;
	Point2<int> originWindow;
	iTranslator.TranslateRASToWindow( voxelRAS.xyz(), voxelWindow.xy() );
	iTranslator.TranslateRASToWindow( originRAS.xyz(), originWindow.xy() );
	Point2<int> voxelSizeWindow;
	voxelSizeWindow.Set( abs( voxelWindow[0] - originWindow[0] ),
			     abs( voxelWindow[1] - originWindow[1] ) );
					  
	iTranslator.TranslateRASToWindow( iRAS, window.xy() );
	updateRectWindow[0][0] = 
	  window[0] - windowBrushRad - voxelSizeWindow[0];
	updateRectWindow[0][1] = 
	  window[1] - windowBrushRad - voxelSizeWindow[1];
	updateRectWindow[1][0] = 
	  window[0] + windowBrushRad + voxelSizeWindow[0];
	updateRectWindow[1][1] =
	  window[1] + windowBrushRad + voxelSizeWindow[1];

	iViewState.AddUpdateRect( updateRectWindow[0][0], updateRectWindow[0][1], updateRectWindow[1][0], updateRectWindow[1][1] );

	// Now get the RAS points in this square or circle.
	list<Point3<float> > points;
	switch( shape ) {
	case ScubaToolState::voxel:
	  points.push_back( Point3<float>(iRAS) );
	  break;
	case ScubaToolState::square: {
	  mVolume->FindRASPointsInSquare( sq[0].xyz(), sq[1].xyz(),
					  sq[2].xyz(), sq[3].xyz(),
					  0,
					  points );
	  } break;
	case ScubaToolState::circle:
	  mVolume->FindRASPointsInCircle( sq[0].xyz(), sq[1].xyz(),
					  sq[2].xyz(), sq[3].xyz(),
					  0, iRAS, rad,
					  points );
	  break;
	}
	
	// For each point we got...
	list<Point3<float> >::iterator tPoints;
	for( tPoints = points.begin(); tPoints != points.end(); ++tPoints ) {

	  // If the point is in bounds...
	  Point3<float> point = *tPoints;
	  VolumeLocation& loc =
	    (VolumeLocation&) mVolume->MakeLocationFromRAS( point.xyz() );
	  
	  if( mVolume->IsInBounds( loc ) ) {

	    // Depending on whether we're editing or selecting, and
	    // whether we're actioning or unactioning, perform the
	    // proper action. Also create the proper undo item.
	    UndoAction* action = NULL;
	    if( ScubaToolState::voxelEditing == iTool.GetMode() ) {

	      // Editing. Get the original value for the undo item.
	      float origValue = mVolume->GetMRINearestValue( loc );

	      // If only brushing zero and are brushing, skip if not
	      // zero.
	      if( iTool.GetOnlyBrushZero() && 
		  iInput.Button() == 2 &&
		  origValue != 0 )
		continue; 

	      // New value depends on voxel button.
	      float newValue = iInput.Button() == 2 ? 
		iTool.GetNewValue() : iTool.GetEraseValue();

	      // Set value and make undo item.
	      mVolume->SetMRIValue( loc, newValue );
	      action = new UndoVoxelEditAction( mVolume, newValue,
						origValue, point.xyz() );


	      // Selecting. 
	    } else if( ScubaToolState::roiEditing == iTool.GetMode() ) {

	      // New value depends on voxel button.
	      bool bSelect = iInput.Button() == 2 ? true : false;

	      // Select or unselect.
	      if( iInput.Button() == 2 ) {
		mVolume->Select( loc );
	      } else if( iInput.Button() == 3 ) {
		mVolume->Unselect( loc );
	      }

	      action = new UndoSelectionAction( mVolume, bSelect, point.xyz());
	    }
	    
	    // Add the undo item.
	    undoList.AddAction( action );
	  }

	  delete &loc;
	}
	RequestRedisplay();
      }
    }

    // Also request a redisplay if this is a mouse up.
    if( iInput.IsButtonUpEvent() ) {
      RequestRedisplay();
    }

    if( ScubaToolState::voxelEditing == iTool.GetMode() ) {
      
      // Adjust the min/max visible value.
      if( iTool.GetMode() == ScubaToolState::voxelEditing ) {
	float newValue = iInput.Button() == 2 ? 
	  iTool.GetNewValue() : iTool.GetEraseValue();
	if( newValue < mMinVisibleValue )
	  mMinVisibleValue = newValue;
	if( newValue > mMaxVisibleValue )
	  mMaxVisibleValue = newValue;
      }
    }

    break;

  case ScubaToolState::straightPath:
  case ScubaToolState::edgePath: {
    
    PathManager& pathMgr = PathManager::GetManager();
    UndoManager& undoList = UndoManager::GetManager();

    if( iInput.IsButtonDownEvent() ) {

      Point3<float> ras( iRAS );
	  
      /* Button down, no current path. If button 1, create a new path
	 and set the first vertex at that point. If 2 or 3, find a
	 path near the click and make it the current path. Select the
	 path we made or found. */
      if( NULL == mCurrentPath ) {

	switch( iInput.Button() ) {
	case 1:
	  mFirstPathRAS.Set( iRAS );
	  mCurrentPath = new Path<float>;
	  pathMgr.ManagePath( *mCurrentPath );
	  mCurrentPath->AddVertex( ras );
	  mCurrentPath->MarkEndOfSegment();
	  mCurrentPath->SetSelected( true );

	  // Make a undo list entry.
	  undoList.BeginAction( "New Path" );
	  undoList.AddAction( new UndoNewPathAction( mCurrentPath ) );
	  undoList.EndAction();

	  break;
	case 2:
	case 3:
	  mCurrentPath = FindClosestPathInPlane( iRAS, iViewState );
	  if( mCurrentPath ) {
	    mCurrentPath->SetSelected( true );
	    mLastPathMoveRAS.Set( iRAS );
	  }
	  break;
	}

	/* Button down, current path. If 1, add a new vertex at this
	   click. If 2, end the path at the click. If button 3,
	   stretch the path back to the first point and end it
	   there. In these cases, unselect the path. */
      } else {

	switch( iInput.Button() ) {
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

	  if( iTool.GetMode() == ScubaToolState::straightPath ) {
	    StretchPathStraight( *mCurrentPath, ras.xyz(),mFirstPathRAS.xyz());
	  } else if( iTool.GetMode() == ScubaToolState::edgePath ) {
	    StretchPathAsEdge( *mCurrentPath, ras.xyz(), 
			       mFirstPathRAS.xyz(), iViewState, iTranslator,
			       iTool.GetEdgePathStraightBias(),
			       iTool.GetEdgePathEdgeBias() );
	  }

	  // TEST CODE: If Shift is down, fill this path.
	  if( iInput.IsShiftKeyDown() ) {
	    Point3<float> seed;
	    float offsets[6][3];
	    bzero( offsets, sizeof(float) * 6 * 3 );
	    offsets[0][0] = GetPreferredBrushRadiusIncrement();
	    offsets[1][1] = GetPreferredBrushRadiusIncrement();
	    offsets[2][2] = GetPreferredBrushRadiusIncrement();
	    offsets[3][0] = -GetPreferredBrushRadiusIncrement();
	    offsets[4][1] = -GetPreferredBrushRadiusIncrement();
	    offsets[5][2] = -GetPreferredBrushRadiusIncrement();
	    for( int nOffset = 0; nOffset < 6; nOffset++ ) {
	      try {
		seed.Set( ras[0] + offsets[nOffset][0], 
			  ras[1] + offsets[nOffset][1], 
			  ras[2] + offsets[nOffset][2] );
		if( mCurrentPath->PointInPath( seed ) ) {
		  cerr << "filling at " << seed << endl;
		  VolumeCollectionFlooder::Params params;
		  SetFloodParams( iTool, iViewState, params );
		  params.mMaxDistance = 10;
		  params.mb3D = false;
		  params.mbStopAtPaths = true;
		  ScubaLayer2DMRIFloodVoxelEdit flooder( iTool.GetNewValue() );
		  flooder.Flood( *mVolume, seed.xyz(), params );
		  break;
		}
	      }
	      catch(...) {}
	    }
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
      if( 2 == iInput.Button() ) {
	
	if( NULL != mCurrentPath ) {
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
      if( NULL != mCurrentPath ) {
	
	if( 3 == iInput.Button() ) {
	  Path<float>* delPath = FindClosestPathInPlane( iRAS, iViewState );
	  if( delPath == mCurrentPath ) {

	    pathMgr.UnmanagePath( *delPath );

	    // The UndoAction will keep a pointer to the path until it
	    // goes off the undo list. It will be deleted then.
	    undoList.BeginAction( "Delete Path" );
	    undoList.AddAction( new UndoDeletePathAction( delPath ) );
	    undoList.EndAction();
	  }
	}
	
	if( 3 == iInput.Button() ||
	    2 == iInput.Button() ) {
	  
	  mCurrentPath->SetSelected( false );
	  mCurrentPath = NULL;
	}

	RequestRedisplay();
      }
      
    } else {

      /* No mouse event, current path. Stretch the path out to our
	 current mouse position. */
      if( NULL != mCurrentPath ) {

	Point3<float>& end = mCurrentPath->GetPointAtEndOfLastSegment();
	if( iTool.GetMode() == ScubaToolState::straightPath ) {
	  StretchPathStraight( *mCurrentPath, end.xyz(), iRAS );
	} else if( iTool.GetMode() == ScubaToolState::edgePath ) {
	  StretchPathAsEdge( *mCurrentPath, end.xyz(), 
			     iRAS, iViewState, iTranslator,
			     iTool.GetEdgePathStraightBias(),
			     iTool.GetEdgePathEdgeBias() );
	}

	RequestRedisplay();
      }
    }
    
    
  } break;
  default:
    break;
  }
}

void
ScubaLayer2DMRI::ProcessOption ( string isOption, string isValue ) {

  char sValue[1024];
  strcpy( sValue, isValue.c_str() );

  if( 0 == isOption.compare( "colormap" ) ) {
    if( 0 == isValue.compare( "grayscale" ) ) {
      SetColorMapMethod( grayscale );
    } else if( 0 == isValue.compare( "heatscale" ) ) {
      SetColorMapMethod( heatScale );
    } else if( 0 == isValue.compare( "lut" ) ) {
      SetColorMapMethod( LUT );
    } else {
      throw runtime_error( "Bad colormap value" );
    }

  } else if( 0 == isOption.compare( "samplemethod" ) ) {
    if( 0 == isValue.compare( "nearest" ) ) {
      SetSampleMethod( nearest );
    } else if( 0 == isValue.compare( "trilinear" ) ) {
      SetSampleMethod( trilinear );
    } else if( 0 == isValue.compare( "sinc" ) ) {
      SetSampleMethod( sinc );
    } else {
      throw runtime_error( "Bad colormap value" );
    }

  } else if( 0 == isOption.compare( "lut" ) ) {
    list<int> lutIDList;
    ScubaColorLUT::GetIDList( lutIDList );
    list<int>::iterator tID;
    for( tID = lutIDList.begin(); tID != lutIDList.end(); ++tID ) {
      int id = *tID;
      ScubaColorLUT& lut = ScubaColorLUT::FindByID( id );
      if( 0 == isValue.compare( lut.GetLabel() ) ) {
	SetColorMapMethod( LUT );
	SetColorLUT( id );
	return;
      }
    }
    if( ERANGE == errno ) {
      throw runtime_error( "LUT not found" );
    }

  } else if( 0 == isOption.compare( "drawzeroclear" ) ) {
    int bDrawZeroClear = strtol( sValue, (char**)NULL, 10 );
    if( ERANGE == errno ) {
      throw runtime_error( "Bad drawzeroclear value" );
    }
    SetDrawZeroClear( bDrawZeroClear );

  } else if( 0 == isOption.compare( "brightness" ) ) {
    float brightness = (float) strtod( sValue, (char**)NULL );
    if( ERANGE == errno ) {
      throw runtime_error( "Bad brightness value" );
    }
    SetBrightness( brightness );

  } else if( 0 == isOption.compare( "contrast" ) ) {
    float contrast = (float) strtod( sValue, (char**)NULL );
    if( ERANGE == errno ) {
      throw runtime_error( "Bad contrast value" );
    }
    SetContrast( contrast );

  } else {
    
    return Layer::ProcessOption( isOption, isValue );
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
  
  switch( mSampleMethod ) {
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
  }
  catch(...) {
    DebugOutput( << "Couldn't find color LUT " << iLUTID );
  }
  
}

int
ScubaLayer2DMRI::GetColorLUT () {

  return mColorLUT->GetID();
}

void
ScubaLayer2DMRI::SetBrightness ( float iBrightness ) { 

  if( iBrightness < 0 || iBrightness > 1 ) {
    throw runtime_error( "Invalid brightness" );
  }
  mBrightness = iBrightness; 
  BuildGrayscaleLUT();
}

void
ScubaLayer2DMRI::SetContrast ( float iContrast ) { 

  if( iContrast < 0 || iContrast > 30 ) {
    throw runtime_error( "Invalid contrast" );
  }
  mContrast = iContrast; 
  BuildGrayscaleLUT();
}



void
ScubaLayer2DMRI::BuildGrayscaleLUT () {

  for( float nEntry = 0; nEntry < cGrayscaleLUTEntries; nEntry+=1 ) {

    // Get the value using the actual min/max to get highest
    // granularity within the 0 - cGrayscaleLUTEntries range.
    float value = ((nEntry * (mMaxVisibleValue-mMinVisibleValue)) / 
		   cGrayscaleLUTEntries) + mMinVisibleValue;

    // Use sigmoid to apply brightness/contrast. Gets 0-1 value.
    float bcdValue;
    if( value != 0 ) {
      bcdValue = (1.0 / (1.0 + exp( (((value-mMinVisibleValue)/(mMaxVisibleValue-mMinVisibleValue))-mBrightness) * -mContrast)));
    } else {
      bcdValue = 0;
    }

    // Normalize back to pixel component value.
    float normValue = bcdValue * kMaxPixelComponentValueFloat;

    // Assign in table.
    mGrayscaleLUT[(int)nEntry] = (int) floor( normValue );
  }
}

void
ScubaLayer2DMRI::SetMinVisibleValue ( float iValue ) {

  mMinVisibleValue = iValue; 
  if( mMinVisibleValue >= mMaxVisibleValue ) {
    mMinVisibleValue = mMinVisibleValue - 0.00001;
  }
}

void
ScubaLayer2DMRI::SetMaxVisibleValue ( float iValue ) { 
  mMaxVisibleValue = iValue; 
  if( mMaxVisibleValue <= mMinVisibleValue ) {
    mMaxVisibleValue = mMinVisibleValue + 0.00001;
  }
}
 
void
ScubaLayer2DMRI::SetHeatScaleMinThreshold ( float iValue ) {

  mHeatScaleMinThreshold = iValue; 
  if( mHeatScaleMinThreshold >= mHeatScaleMidThreshold ) {
    mHeatScaleMinThreshold = mHeatScaleMidThreshold - 0.0001;
  }
}

void
ScubaLayer2DMRI::SetHeatScaleMidThreshold ( float iValue ) {

  mHeatScaleMidThreshold = iValue; 
  if( mHeatScaleMidThreshold <= mHeatScaleMinThreshold ) {
    mHeatScaleMidThreshold = mHeatScaleMinThreshold + 0.0001;
  }
  if( mHeatScaleMidThreshold >= mHeatScaleMaxThreshold ) {
    mHeatScaleMidThreshold = mHeatScaleMaxThreshold - 0.0001;
  }
}

void
ScubaLayer2DMRI::SetHeatScaleMaxThreshold ( float iValue ) {

  mHeatScaleMaxThreshold = iValue; 
  if( mHeatScaleMaxThreshold <= mHeatScaleMidThreshold ) {
    mHeatScaleMaxThreshold = mHeatScaleMidThreshold + 0.0001;
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
				     float iStraightBias, float iEdgeBias ){

  // Make an edge path finder.
  EdgePathFinder finder( iViewState.mBufferWidth, iViewState.mBufferHeight,
			 (int)mVolume->GetMRIMagnitudeMaxValue(),
			 &iTranslator, mVolume );
  finder.DisableOutput();
  finder.SetStraightBias( iStraightBias );
  finder.SetEdgeBias( iEdgeBias );

  // Get the first point from the path and the last point as passed
  // in. Convert to window points. Then find the path between them.
  Point3<float> beginRAS( iRASBegin );
  Point3<float> endRAS( iRASEnd );

  Point2<int> beginWindow;
  Point2<int> endWindow;
  iTranslator.TranslateRASToWindow( beginRAS.xyz(), beginWindow.xy() );
  iTranslator.TranslateRASToWindow( endRAS.xyz(), endWindow.xy() );
  list<Point2<int> > windowPoints;
  
  finder.FindPath( beginWindow, endWindow, windowPoints );

  // Clear the path points and add the points we got from the path,
  // converting them to RAS one the way.
  iPath.ClearLastSegment();
  list<Point2<int> >::iterator tWindowPoint;
  for( tWindowPoint = windowPoints.begin();
       tWindowPoint != windowPoints.end();
       ++tWindowPoint ) {
    Point3<float> currentRAS;
    iTranslator.TranslateWindowToRAS( (*tWindowPoint).xy(), currentRAS.xyz() );
    iPath.AddVertex( currentRAS );
  }
}


Path<float>*
ScubaLayer2DMRI::FindClosestPathInPlane ( float iRAS[3],
					  ViewState& iViewState ) {

  float minDistance = mWidth * mHeight;
  Path<float>* closestPath = NULL;
  Point3<float> whereRAS( iRAS );

  float range = 0;
  switch( iViewState.mInPlane ) {
  case 0: range = mVolume->GetVoxelXSize() / 2.0; break;
  case 1: range = mVolume->GetVoxelYSize() / 2.0; break;
  case 2: range = mVolume->GetVoxelZSize() / 2.0; break;
  }

  PathManager& pathMgr = PathManager::GetManager();
  list<Path<float>*>::iterator tPath;
  list<Path<float>*>& paths = pathMgr.GetPathList();  
  for( tPath = paths.begin(); tPath != paths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    Point3<float>& beginRAS = path->GetVertexAtIndex( 0 );
    if( iViewState.IsRASVisibleInPlane( beginRAS.xyz(), range ) ) {
      
      float minDistanceInPath = 999999;
      
      int cVertices = path->GetNumVertices();
      int nCurVertex = 1;
      for( nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {
	
	int nBackVertex = nCurVertex - 1;
	
	Point3<float>& curVertex  = path->GetVertexAtIndex( nCurVertex );
	Point3<float>& backVertex = path->GetVertexAtIndex( nBackVertex );
	
	float distance = 
	  Utilities::DistanceFromSegmentToPoint3f( curVertex, backVertex,
						   whereRAS );

	if( distance < minDistanceInPath ) {
	  minDistanceInPath = distance;
	}
      }
      
      if( minDistanceInPath < minDistance ) {
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
  for( int nCurVertex = 1; nCurVertex < cVertices; nCurVertex++ ) {
    
    int nBackVertex = nCurVertex - 1;
    
    Point3<float>& curVertex  = iPath.GetVertexAtIndex( nCurVertex );
    Point3<float>& backVertex = iPath.GetVertexAtIndex( nBackVertex );
    
    int curWindow[2], backWindow[2];
    iTranslator.TranslateRASToWindow( curVertex.xyz(), curWindow );
    iTranslator.TranslateRASToWindow( backVertex.xyz(), backWindow );
    
    DrawLineIntoBuffer( iBuffer, iWidth, iHeight, backWindow, curWindow,
			iColor, 1, 1 );
  }
}


void
ScubaLayer2DMRI::GetPreferredInPlaneIncrements ( float oIncrements[3] ) {
  
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

void
ScubaLayer2DMRI::SetFloodParams ( ScubaToolState& iTool, ViewState& iViewState,
				  VolumeCollectionFlooder::Params& ioParams ) {

  ioParams.mSourceCollection  = iTool.GetFloodSourceCollection();
  ioParams.mbStopAtPaths      = iTool.GetFloodStopAtPaths();
  ioParams.mbStopAtROIs       = iTool.GetFloodStopAtROIs();
  ioParams.mb3D               = iTool.GetFlood3D();
  ioParams.mFuzziness         = iTool.GetFloodFuzziness();
  ioParams.mMaxDistance       = iTool.GetFloodMaxDistance();
  ioParams.mViewNormal[0]     = iViewState.mPlaneNormal[0];
  ioParams.mViewNormal[1]     = iViewState.mPlaneNormal[1];
  ioParams.mViewNormal[2]     = iViewState.mPlaneNormal[2];
  ioParams.mbOnlyZero         = iTool.GetOnlyFloodZero();
  ioParams.mFuzzinessType     = 
    (VolumeCollectionFlooder::Params::FuzzinessType) iTool.GetFuzzinessType();
  if( !iTool.GetFlood3D() ) {
    ioParams.mbWorkPlaneX     = (iViewState.mInPlane == 0);
    ioParams.mbWorkPlaneY     = (iViewState.mInPlane == 1);
    ioParams.mbWorkPlaneZ     = (iViewState.mInPlane == 2);
  }
}


void
ScubaLayer2DMRI::DoTimer () {

  mTimersSinceLastAutosave++;
  if( mTimersSinceLastAutosave > kcTimersBetweenAutosaves ) {

    if( mVolume->IsAutosaveDirty() ) {
      
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

  undoList.BeginAction( "Voxel Fill" );

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
  undoList.EndAction();

  mVolume->EndBatchChanges();
}

bool
ScubaLayer2DMRIFloodVoxelEdit::DoStopRequested () {

  // Check for the stop button.
  ProgressDisplayManager& manager = 
    ProgressDisplayManager::GetManager();
  int nButton = manager.CheckTaskForButton();
  if( nButton == 0 ) {
    return true;
  } 

  return false;
}

bool
ScubaLayer2DMRIFloodVoxelEdit::CompareVoxel ( float[3] ) {

  // Always return true.
  return true;
}

void
ScubaLayer2DMRIFloodVoxelEdit::DoVoxel ( float iRAS[3] ) {
  UndoManager& undoList = UndoManager::GetManager();

  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( iRAS );

  // Save the original value. Set the new value. 
  float origValue = mVolume->GetMRINearestValue( loc );
  mVolume->SetMRIValue( loc, mValue );

  // Make an undo item with the old value and add it to the list.
  UndoVoxelEditAction* action = 
    new UndoVoxelEditAction( mVolume, mValue, origValue, iRAS );

  undoList.AddAction( action );

  delete &loc;
}

UndoVoxelEditAction::UndoVoxelEditAction ( VolumeCollection* iVolume,
					   float iNewValue, float iOrigValue, 
					   float iRAS[3] ) {
  mVolume = iVolume;
  mNewValue = iNewValue;
  mOrigValue = iOrigValue;
  mRAS[0] = iRAS[0];
  mRAS[1] = iRAS[1];
  mRAS[2] = iRAS[2];
}

void
UndoVoxelEditAction::Undo () {

  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( mRAS );
  mVolume->SetMRIValue( loc, mOrigValue );
  delete &loc;
}

void
UndoVoxelEditAction::Redo () {

  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( mRAS );
  mVolume->SetMRIValue( loc, mNewValue );
  delete &loc;
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
  
  if( mbSelect ) {
    manager.NewTask( "Selecting", "Selecting voxels", false, lButtons );
  } else {
    manager.NewTask( "Unselecting", "Unselecting voxels", false, lButtons );
  }

  // Start our undo action.
  UndoManager& undoList = UndoManager::GetManager();

  if( mbSelect ) {
    undoList.BeginAction( "Selection Fill" );
  } else {
    undoList.BeginAction( "Unselection Fill" );
  }

}

void 
ScubaLayer2DMRIFloodSelect::DoEnd () {

  // End the task.
  ProgressDisplayManager& manager =
    ProgressDisplayManager::GetManager();
  manager.EndTask();

  // End our undo action.
  UndoManager& undoList = UndoManager::GetManager();
  undoList.EndAction();
}

bool
ScubaLayer2DMRIFloodSelect::DoStopRequested () {

  // Check for the stop button.
  ProgressDisplayManager& manager = 
    ProgressDisplayManager::GetManager();
  int nButton = manager.CheckTaskForButton();
  if( nButton == 0 ) {
    return true;
  } 

  return false;
}

bool
ScubaLayer2DMRIFloodSelect::CompareVoxel ( float[3] ) {

  // Always return true.
  return true;
}

void
ScubaLayer2DMRIFloodSelect::DoVoxel ( float iRAS[3] ) {
  UndoManager& undoList = UndoManager::GetManager();

  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( iRAS );

  if( mbSelect ) {
    mVolume->Select( loc );
    UndoSelectionAction* action = 
      new UndoSelectionAction( mVolume, true, iRAS );
    undoList.AddAction( action );
  } else {
    mVolume->Unselect( loc );
    UndoSelectionAction* action =
      new UndoSelectionAction( mVolume, false, iRAS );
    undoList.AddAction( action );
  }

  delete &loc;
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
  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( mRAS );

  if( mbSelect ) {
    mVolume->Unselect( loc );
  } else {
    mVolume->Select( loc );
  }

  delete &loc;
}

void
UndoSelectionAction::Redo () {
  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( mRAS );

  if( mbSelect ) {
    mVolume->Select( loc );
  } else {
    mVolume->Unselect( loc );
  }

  delete &loc;
}


UndoPathAction::UndoPathAction ( Path<float>* iPath ) {
  mPath = iPath;
}

UndoPathAction::~UndoPathAction () {

}

UndoNewPathAction::UndoNewPathAction ( Path<float>* iPath ) :
  UndoPathAction( iPath ) {
}

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
  UndoPathAction( iPath ) {
}

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
				 int iLongestEdge,
				 ScubaWindowToRASTranslator* iTranslator,
				 VolumeCollection* iVolume ) {

  SetDimensions( iViewWidth, iViewHeight, iLongestEdge );
  mVolume = iVolume;
  mTranslator = iTranslator;
}
		  
float 
EdgePathFinder::GetEdgeCost ( Point2<int>& iPoint ) {

  // Get the magnitude value at this point.
  float RAS[3];
  mTranslator->TranslateWindowToRAS( iPoint.xy(), RAS );
  VolumeLocation& loc = (VolumeLocation&) mVolume->MakeLocationFromRAS( RAS );

  if( mVolume->IsInBounds( loc ) ) {
    return 1.0 / (mVolume->GetMRIMagnitudeValue( loc ) + 0.0001);
  } else {
    return mLongestEdge;
  }

  delete &loc;
}
