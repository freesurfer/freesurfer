#include "ScubaLayer2DMRI.h"
#include "ViewState.h"
#include "talairachex.h"

using namespace std;

ScubaLayer2DMRI::ScubaLayer2DMRI () {

  mVolume = NULL;
  mWorldToIndexMatrix = NULL;
  mWorldCoord = VectorAlloc( 4, MATRIX_REAL );
  mIndexCoord = VectorAlloc( 4, MATRIX_REAL );

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollection" );
}

ScubaLayer2DMRI::~ScubaLayer2DMRI () {

  if( NULL != mWorldCoord ) {
    VectorFree( &mWorldCoord );
  }
  if( NULL != mIndexCoord ) {
    VectorFree( &mIndexCoord );
  }
  if( NULL != mWorldToIndexMatrix ) {
    MatrixFree( &mWorldToIndexMatrix );
  }
}

void
ScubaLayer2DMRI::SetVolumeCollection ( VolumeCollection& iVolume ) {

  mVolume = &iVolume;
}

void 
ScubaLayer2DMRI::DrawIntoBuffer ( GLubyte* iBuffer, ViewState& iViewState,
				  ScubaWindowToRASTranslator& iTranslator ) {

  if( NULL == mVolume ) {
    DebugOutput( << "No volume to draw" );
    return;
  }

  MRI* mri = mVolume->GetMRI();

  if( NULL == mWorldToIndexMatrix ) {
    mWorldToIndexMatrix = extract_r_to_i( mri );
  }

  float minValue, maxValue;
  MRIvalRange( mri, &minValue, &maxValue );
  GLubyte* dest = iBuffer;

  float halfWidth  = ((float)mWidth / 2.0);
  float halfHeight = ((float)mHeight / 2.0);

  // Point to the beginning of the buffer. For each pixel in the
  // buffer...
  for( int nY = 0; nY < mHeight; nY++ ) {
    for( int nX = 0; nX < mWidth; nX++ ) {

      // Use our translator to get an RAS point.
      float world[3];
      iTranslator.TranslateWindowToRAS( nX, nY, world );
      
      // Make sure this is within the bounds of the volume using
      // {x,y,z}{start,end}. If it is...
      GLubyte anaColor[3];
      if ( world[0] > mri->xstart && world[0] < mri->xend &&
	   world[1] > mri->ystart && world[1] < mri->yend &&
	   world[2] > mri->zstart && world[2] < mri->zend ) {
	
	// Convert to an index. Look up the value of the volume at this index.
	Real index[3];
	VECTOR_ELT( mWorldCoord, 1 ) = world[0];
	VECTOR_ELT( mWorldCoord, 2 ) = world[1];
	VECTOR_ELT( mWorldCoord, 3 ) = world[2];
	VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
	MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
	index[0] = VECTOR_ELT( mIndexCoord, 1 );
	index[1] = VECTOR_ELT( mIndexCoord, 2 );
	index[2] = VECTOR_ELT( mIndexCoord, 3 );
	
	if( index[0] >= 0 && index[0] < mri->width &&
	    index[1] >= 0 && index[1] < mri->height &&
	    index[2] >= 0 && index[2] < mri->depth ) {

	  float anaValue;
	  switch( mri->type ) {
	  case MRI_UCHAR:
	    anaValue = 
	      MRIvox( mri, (int)index[0], (int)index[1], (int)index[2] );
	    break;
	  case MRI_INT:
	    anaValue = 
	      MRIIvox( mri, (int)index[0], (int)index[1], (int)index[2] );
	    break;
	  case MRI_LONG:
	    anaValue = 
	      MRILvox( mri, (int)index[0], (int)index[1], (int)index[2] );
	    break;
	  case MRI_FLOAT:
	    anaValue = 
	      MRIFvox( mri, (int)index[0], (int)index[1], (int)index[2] );
	    break;
	  case MRI_SHORT:
	    anaValue = 
	      MRISvox( mri, (int)index[0], (int)index[1], (int)index[2] );
	    break;
	  default:
	    anaValue = 0;
	    break ;
	  }

	  // If we have a lookup table, use that to get the color. Otherwise
	  // make a grayscale color.
	  int color = (int) floor( 255.0 * 
    		            ((anaValue - minValue) / (maxValue - minValue)) );
	  if( color < 0 ) color = 0;
	  if( color > 255 ) color =  255;
	  anaColor[0] = (GLubyte) color;
	  anaColor[1] = anaColor[2] = anaColor[0];

	} else {

	  // In RAS bounds but out of index bounds.
	  anaColor[0] = (GLubyte) 255;
	  anaColor[1] = (GLubyte) 0;
	  anaColor[2] = (GLubyte) 0;
	}
	
      } else {
	
	// Out of bounds, so fill with background color.
	anaColor[0] = (GLubyte) 0;
	anaColor[1] = (GLubyte) 0;
	anaColor[2] = (GLubyte) 0;
      }

      // Write the RGB value to the buffer. Write a 255 in the alpha byte.
      dest[0] = anaColor[0];
      dest[1] = anaColor[1];
      dest[2] = anaColor[2];
      dest[3] = (GLubyte)255;
    
      // Advance our pixel buffer pointer.
      dest += 4;
      
    }
  }

}
  
void 
ScubaLayer2DMRI::GetInfoAtRAS ( float inX, float inY, float inZ,
				map<string,string>& iLabelValues ) {

  if( NULL == mVolume ) {
    return;
  }

  MRI* mri = mVolume->GetMRI();

 // Convert to an index. Look up the value of the volume at this index.
  Real index[3];
  VECTOR_ELT( mWorldCoord, 1 ) = inX;
  VECTOR_ELT( mWorldCoord, 2 ) = inY;
  VECTOR_ELT( mWorldCoord, 3 ) = inZ;
  VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
  MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
  index[0] = VECTOR_ELT( mIndexCoord, 1 );
  index[1] = VECTOR_ELT( mIndexCoord, 2 );
  index[2] = VECTOR_ELT( mIndexCoord, 3 );
  
  if( index[0] >= 0 && index[0] < mri->width &&
      index[1] >= 0 && index[1] < mri->height &&
      index[2] >= 0 && index[2] < mri->depth ) {
    
    float anaValue;
    switch( mri->type ) {
    case MRI_UCHAR:
      anaValue = 
	MRIvox( mri, (int)index[0], (int)index[1], (int)index[2] );
      break;
    case MRI_INT:
      anaValue = 
	MRIIvox( mri, (int)index[0], (int)index[1], (int)index[2] );
      break;
    case MRI_LONG:
      anaValue = 
	MRILvox( mri, (int)index[0], (int)index[1], (int)index[2] );
      break;
    case MRI_FLOAT:
      anaValue = 
	MRIFvox( mri, (int)index[0], (int)index[1], (int)index[2] );
      break;
    case MRI_SHORT:
      anaValue = 
	MRISvox( mri, (int)index[0], (int)index[1], (int)index[2] );
      break;
    default:
      anaValue = 0;
      break ;
    }

    stringstream ssValue;
    ssValue << anaValue << " ";
    
    iLabelValues[msLabel] = ssValue.str();
  }
}
  
void 
ScubaLayer2DMRI::DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv ) {

  // SetVolumeCollection <layerID> <collectionID>
  if( 0 == strcmp( isCommand, "SetVolumeCollection" ) ) {
    if( 3 == iArgc ) {
      int layerID = strtol(iasArgv[1], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad layer ID";
	return;
      }

      if( mID == layerID ) {
	
	int collectionID = strtol(iasArgv[2], (char**)NULL, 10);
	if( ERANGE == errno ) {
	  sResult = "bad collection ID";
	  return;
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
	  return;
	}
	catch(...) {
	  sResult = "bad collection ID, collection not found";
	  return;
	}
      }
    } else {
      sResult = "wrong # args: should be \"SetVolumeCollection "
	"layerID collectionID\"";
      DebugOutput( << sResult );
      return;
    }
  }

  Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

