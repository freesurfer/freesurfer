#include <fstream>
#include "ScubaLayer2DMRI.h"
#include "ViewState.h"
#include "talairachex.h"

using namespace std;

int const ScubaLayer2DMRI::cGrayscaleLUTEntries = 256;
int const ScubaLayer2DMRI::kMaxPixelComponentValue = 255;
int const ScubaLayer2DMRI::cDefaultFileLUTEntries = 500;

ScubaLayer2DMRI::ScubaLayer2DMRI () {

  SetOutputStreamToCerr();

  mVolume = NULL;
  mSampleMethod = nearest;
  mColorMapMethod = grayscale;
  mBrightness = 0.25;
  mContrast = 12.0;
  mfnLUT = "";
  mcFileLUTEntries = cDefaultFileLUTEntries;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "Set2DMRILayerVolumeCollection", 2, 
			 "layerID collectionID",
			 "Sets the volume collection for this layer." );
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
  commandMgr.AddCommand( *this, "Set2DMRILayerFileLUTFileName", 2, 
			 "layerID fileName",
			 "Sets the LUT file name for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRILayerFileLUTFileName", 1, "layerID",
			 "Returns the LUT file name for this layer." );
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
  
}

ScubaLayer2DMRI::~ScubaLayer2DMRI () {

}

void
ScubaLayer2DMRI::SetVolumeCollection ( VolumeCollection& iVolume ) {

  mVolume = &iVolume;
  BuildGrayscaleLUT();
  BuildLUTFromFile();

  mMinVisibleValue = mVolume->GetMRIMinValue();
  mMaxVisibleValue = mVolume->GetMRIMaxValue();
}

void 
ScubaLayer2DMRI::DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				  ViewState& iViewState,
				  ScubaWindowToRASTranslator& iTranslator ) {

  if( NULL == mVolume ) {
    DebugOutput( << "No volume to draw" );
    return;
  }

  float midValue;
  midValue = (mMaxVisibleValue - mMinVisibleValue) / 2.0 + mMinVisibleValue;

  GLubyte* dest = iBuffer;

  float halfWidth  = ((float)iWidth / 2.0);
  float halfHeight = ((float)iHeight / 2.0);

  // Point to the beginning of the buffer. For each pixel in the
  // buffer...
  int window[2];
  for( window[1] = 0; window[1] < iHeight; window[1]++ ) {
    for( window[0] = 0; window[0] < iWidth; window[0]++ ) {

      // Use our translator to get an RAS point.
      float RAS[3];
      iTranslator.TranslateWindowToRAS( window, RAS );

      // Make sure this is within the bounds. If it is...
      GLubyte anaColor[3];
      if( mVolume->IsRASInMRIBounds( RAS ) ) {
	
	float value;
	switch( mSampleMethod ) {
	case nearest:
	  value = mVolume->GetMRINearestValueAtRAS( RAS ); 
 	  break;
	case trilinear:
	  value = mVolume->GetMRITrilinearValueAtRAS( RAS ); 
	  break;
	case sinc:
	  value = mVolume->GetMRISincValueAtRAS( RAS ); 
	  break;
	}
	
	if( (mColorMapMethod == heatScale || value >= mMinVisibleValue) &&
	    (mColorMapMethod == heatScale || value <= mMaxVisibleValue) &&
	    ((value != 0 && mbClearZero) || !mbClearZero) ) {
	  
	  switch( mColorMapMethod ) { 
	  case grayscale: {
	    int nLUT = (int) floor( (cGrayscaleLUTEntries-1) * 
				    ((value - mMinVisibleValue) /
				     (mMaxVisibleValue - mMinVisibleValue)) );
	    if( nLUT < 0 ) nLUT = 0;
	    if( nLUT > cGrayscaleLUTEntries ) nLUT = cGrayscaleLUTEntries-1;
	    anaColor[0] = (GLubyte) mGrayscaleLUT[nLUT];
	    anaColor[1] = anaColor[2] = anaColor[0];
	  }
	    break;
	  case heatScale:
	    
	    if( fabs(value) >= mMinVisibleValue &&
		fabs(value) <= mMaxVisibleValue ) {
	      
	      float tmp;
	      if ( fabs(value) > mMinVisibleValue &&
		   fabs(value) < midValue ) {
		tmp = fabs(value);
		tmp = (1.0/(midValue-mMinVisibleValue)) *
		  (tmp-mMinVisibleValue)*(tmp-mMinVisibleValue) + 
		  mMinVisibleValue;
		value = (value<0) ? -tmp : tmp;
	      }
	      
	      /* calc the color */
	      float red, green, blue;
	      if( value >= 0 ) {
		red = ((value<mMinVisibleValue) ? 0.0 : 
		       (value<midValue) ? 
		       (value-mMinVisibleValue)/(midValue-mMinVisibleValue) :
		       1.0);
		green = ((value<midValue) ? 0.0 :
			 (value<mMaxVisibleValue) ? 
			 (value-midValue)/(mMaxVisibleValue-midValue) : 1.0);
		blue = 0.0; 
	      } else {
		value = -value;
		red = 0.0;
		green = ((value<midValue) ? 0.0 :
			 (value<mMaxVisibleValue) ? 
			 (value-midValue)/(mMaxVisibleValue-midValue) : 1.0);
		blue = ((value<mMinVisibleValue) ? 0.0 :
			(value<midValue) ? 
			(value-mMinVisibleValue)/(midValue-mMinVisibleValue) : 
			1.0);
	      }
	      
	      if( red > 1.0 )   red = 1.0;
	      if( green > 1.0 ) green = 1.0;
	      if( blue > 1.0 )  blue = 1.0;
	      
	      anaColor[0] = (GLubyte)
		(red * (float)kMaxPixelComponentValue);
	      anaColor[1] = (GLubyte) 
		(green * (float)kMaxPixelComponentValue);
	      anaColor[2] = (GLubyte) 
		(blue * (float)kMaxPixelComponentValue);
	    } else {
	      anaColor[0] = dest[0];
	      anaColor[1] = dest[1];
	      anaColor[2] = dest[2];
	    }
	    break;
	  case LUT:
	    if( value > 0 && value < mcFileLUTEntries ) {
	      int nLUT = (int)value;
	      anaColor[0] = (GLubyte) mFileLUT[nLUT].r;
	      anaColor[1] = (GLubyte) mFileLUT[nLUT].g;
	      anaColor[2] = (GLubyte) mFileLUT[nLUT].b;
	    } else  {
	      anaColor[0] = anaColor[1] = anaColor[2] = 0;
	    }
	    break;
	  }
	  
	  // Write the RGB value to the buffer. Write a 255 in the
	  // alpha byte.
	  dest[0] = (GLubyte) (((float)dest[0] * (1.0 - mOpacity)) +
			       ((float)anaColor[0] * mOpacity));
	  dest[1] = (GLubyte) (((float)dest[1] * (1.0 - mOpacity)) +
			       ((float)anaColor[1] * mOpacity));
	  dest[2] = (GLubyte) (((float)dest[2] * (1.0 - mOpacity)) +
			       ((float)anaColor[2] * mOpacity));
	  dest[3] = (GLubyte)255;
	}
      }

      // Advance our pixel buffer pointer.
      dest += 4;
      
    }
  }

}
  
void 
ScubaLayer2DMRI::GetInfoAtRAS ( float iRAS[3], 
				map<string,string>& iLabelValues ) {

  if( NULL == mVolume ) {
    return;
  }

  // Look up the value of the volume at this point.
  if ( mVolume->IsRASInMRIBounds( iRAS ) ) {
    
    float value;
    value = mVolume->GetMRINearestValueAtRAS( iRAS ); 

    // If this is a LUT volume, use the label from the lookup file,
    // otherwise just display the value.
    stringstream ssValue;
    if( mColorMapMethod == LUT ) {
      ssValue << mFileLUT[(int)value].msLabel;
    } else {
      ssValue << value;
    }

    iLabelValues[mVolume->GetLabel()] = ssValue.str();
  }
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
	if( data.GetID() != collectionID ) {
	  cerr << "IDs didn't match" << endl;
	}
	VolumeCollection& volume = (VolumeCollection&)data;
	// VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);
	
	SetVolumeCollection( volume );
      }
      catch( std::bad_cast& e ) {
	DebugOutput( << "Bad cast from DataCollection" );
	sResult = "bad collection ID, collection not a volume collection";
	return error;
      }
      catch(...) {
	sResult = "bad collection ID, collection not found";
	return error;
      }
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
      if( LUT == mColorMapMethod ) {
	BuildLUTFromFile();
      }
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
      
      float brightness = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad brightness";
	return error;
      }

      SetBrightness( brightness );
      BuildGrayscaleLUT();
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
      
      float contrast = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad contrast";
	return error;
      }

      SetContrast( contrast );
      BuildGrayscaleLUT();
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

  // Set2DMRIFileLUTFileName <layerID> <fileName>
  if( 0 == strcmp( isCommand, "Set2DMRILayerFileLUTFileName" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      SetLUTFileName( iasArgv[2] );
      BuildLUTFromFile();
    }
  }

  // Get2DMRIFileLUTFileName <layerID>
  if( 0 == strcmp( isCommand, "Get2DMRILayerFileLUTFileName" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      stringstream ssReturnValues;
      ssReturnValues << "\"" << mfnLUT << "\"";
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "s";
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
      
      float value = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad opacity";
	return error;
      }
      
      SetMinVisibleValue( value );
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
      
      float value = strtof( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad opacity";
	return error;
      }
      
      SetMaxVisibleValue( value );
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

  return Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
ScubaLayer2DMRI::HandleTool ( float iRAS[3],
			      ScubaToolState& iTool, InputState& iInput ) {
  
  switch( iTool.GetMode() ) {
  case ScubaToolState::voxelEditing:

    switch( iInput.Button() ) {
    case 1: 
      mVolume->SetMRIValueAtRAS( iRAS, 255 );
      RequestRedisplay();
      break;
    }
  }

}


void
ScubaLayer2DMRI::BuildGrayscaleLUT () {

  MRI* mri = mVolume->GetMRI();
  if( NULL == mri ) return;

  for( float nEntry = 0; nEntry < cGrayscaleLUTEntries; nEntry+=1 ) {

    // Get the value using the actual min/max to get highest
    // granularity within the 0 - cGrayscaleLUTEntries range.
    float value = ((nEntry * (mMaxVisibleValue-mMinVisibleValue)) / 
		   cGrayscaleLUTEntries) + mMinVisibleValue;

    // Use sigmoid to apply brightness/contrast. Gets 0-1 value.
    float bcdValue = (1.0 / (1.0 + exp( (((value-mMinVisibleValue)/(mMaxVisibleValue-mMinVisibleValue))-mBrightness) * -mContrast)));

    // Normalize back to pixel component value.
    float normValue = bcdValue * (float)kMaxPixelComponentValue;

    // Assign in table.
    mGrayscaleLUT[(int)nEntry] = normValue;
  }
}

void
ScubaLayer2DMRI::BuildLUTFromFile () {

  ifstream fLUT( mfnLUT.c_str() );
  if( fLUT.is_open() ) {

    int nLine = 0;
    while( !fLUT.eof() ) {

      string sLine;
      getline( fLUT, sLine );
      
      int nEntry;
      char sLabel[1024];
      int red, green, blue;
      int eRead = sscanf( sLine.c_str(), "%d %s %d %d %d %*s",
			  &nEntry, sLabel, &red, &green, &blue );
      if( 5 != eRead &&
	  -1 != eRead ) {
	DebugOutput(  "Error parsing " << mfnLUT << ": Malformed line " 
		      << nLine );
	continue;
      }

      mFileLUT[nEntry].r  = red;
      mFileLUT[nEntry].g  = green;
      mFileLUT[nEntry].b  = blue;
      mFileLUT[nEntry].msLabel = sLabel;
    
      nLine++;
    }

  } else {

    for( int nEntry = 0; nEntry < mcFileLUTEntries; nEntry+=1 ) {
      // Assign random colors.
      mFileLUT[nEntry].r = rand() % kMaxPixelComponentValue;
      mFileLUT[nEntry].g = rand() % kMaxPixelComponentValue;
      mFileLUT[nEntry].b = rand() % kMaxPixelComponentValue;
      mFileLUT[nEntry].msLabel = "";
    }
  }
}
