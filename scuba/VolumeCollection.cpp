#include <stdexcept>
#include <vector>
#include "VolumeCollection.h"
#include "DataManager.h"
#include "Point3.h"

using namespace std;


VolumeCollection::VolumeCollection () :
  DataCollection() {
  mMRI = NULL;
  mWorldToIndexMatrix = NULL;
  mIndexToWorldMatrix = NULL;
  mWorldCoord = VectorAlloc( 4, MATRIX_REAL );
  mIndexCoord = VectorAlloc( 4, MATRIX_REAL );
  mEdgeVoxels = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollectionFileName", 2, 
			 "collectionID fileName", 
			 "Sets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "GetVolumeCollectionFileName", 1, 
			 "collectionID", 
			 "Gets the file name for a given volume collection.");
}

VolumeCollection::~VolumeCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();
  try { 
    mriLoader.ReleaseData( &mMRI );
  } 
  catch(...) {
    cerr << "Couldn't release data"  << endl;
  }

  if( NULL != mWorldCoord ) {
    VectorFree( &mWorldCoord );
  }
  if( NULL != mIndexCoord ) {
    VectorFree( &mIndexCoord );
  }
  if( NULL != mWorldToIndexMatrix ) {
    MatrixFree( &mWorldToIndexMatrix );
  }
  if( NULL != mIndexToWorldMatrix ) {
    MatrixFree( &mIndexToWorldMatrix );
  }
}

void
VolumeCollection::SetFileName ( string& ifnMRI ) {

  mfnMRI = ifnMRI;
}

MRI*
VolumeCollection::GetMRI() { 

  if( NULL == mMRI ) {
    
    DataManager dataMgr = DataManager::GetManager();
    MRILoader mriLoader = dataMgr.GetMRILoader();

    try { 
      mMRI = mriLoader.GetData( mfnMRI );
    }
    catch( exception e ) {
      throw logic_error( "Couldn't load MRI" );
    }

    if( msLabel == "" ) {
      SetLabel( mfnMRI );
    }

    mWorldToIndexMatrix = extract_r_to_i( mMRI );
    mIndexToWorldMatrix = extract_i_to_r( mMRI );
    UpdateMRIValueRange();

    // Size all the rois we may have.
    int bounds[3];
    bounds[0] = mMRI->width;
    bounds[1] = mMRI->height;
    bounds[2] = mMRI->depth;
    map<int,ScubaROI*>::iterator tIDROI;
    for( tIDROI = mROIMap.begin();
	 tIDROI != mROIMap.end(); ++tIDROI ) {
	ScubaROIVolume* roi = (ScubaROIVolume*)(*tIDROI).second;
	roi->SetROIBounds( bounds );
    }

    // Init the edge volume.
    InitEdgeVolume();
  }

  return mMRI; 
}

void
VolumeCollection::UpdateMRIValueRange () {

  if( NULL != mMRI ) {
    MRIvalRange( mMRI, &mMRIMinValue, &mMRIMaxValue );
  }
}

float 
VolumeCollection::GetVoxelXSize () {

  if( NULL != mMRI ) {
    return mMRI->xsize;
  }
}

float 
VolumeCollection::GetVoxelYSize () {

  if( NULL != mMRI ) {
    return mMRI->ysize;
  }
}

float 
VolumeCollection::GetVoxelZSize () {

  if( NULL != mMRI ) {
    return mMRI->zsize;
  }
}


void
VolumeCollection::RASToMRIIndex ( float iRAS[3], int oIndex[3] ) {
  
  VECTOR_ELT( mWorldCoord, 1 ) = iRAS[0];
  VECTOR_ELT( mWorldCoord, 2 ) = iRAS[1];
  VECTOR_ELT( mWorldCoord, 3 ) = iRAS[2];
  VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
  MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
  oIndex[0] = (int) rint( VECTOR_ELT( mIndexCoord, 1 ) );
  oIndex[1] = (int) rint( VECTOR_ELT( mIndexCoord, 2 ) );
  oIndex[2] = (int) rint( VECTOR_ELT( mIndexCoord, 3 ) );
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], float oIndex[3] ) {
  
  VECTOR_ELT( mWorldCoord, 1 ) = iRAS[0];
  VECTOR_ELT( mWorldCoord, 2 ) = iRAS[1];
  VECTOR_ELT( mWorldCoord, 3 ) = iRAS[2];
  VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
  MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
  oIndex[0] = VECTOR_ELT( mIndexCoord, 1 );
  oIndex[1] = VECTOR_ELT( mIndexCoord, 2 );
  oIndex[2] = VECTOR_ELT( mIndexCoord, 3 );
}

void
VolumeCollection::MRIIndexToRAS ( int iIndex[3], float oRAS[3] ) {
  
  VECTOR_ELT( mIndexCoord, 1 ) = iIndex[0];
  VECTOR_ELT( mIndexCoord, 2 ) = iIndex[1];
  VECTOR_ELT( mIndexCoord, 3 ) = iIndex[2];
  VECTOR_ELT( mIndexCoord, 4 ) = 1.0;
  MatrixMultiply( mIndexToWorldMatrix, mIndexCoord, mWorldCoord );
  oRAS[0] = VECTOR_ELT( mWorldCoord, 1 );
  oRAS[1] = VECTOR_ELT( mWorldCoord, 2 );
  oRAS[2] = VECTOR_ELT( mWorldCoord, 3 );
}

void
VolumeCollection::MRIIndexToRAS ( float iIndex[3], float oRAS[3] ) {
  
  VECTOR_ELT( mIndexCoord, 1 ) = iIndex[0];
  VECTOR_ELT( mIndexCoord, 2 ) = iIndex[1];
  VECTOR_ELT( mIndexCoord, 3 ) = iIndex[2];
  VECTOR_ELT( mIndexCoord, 4 ) = 1.0;
  MatrixMultiply( mIndexToWorldMatrix, mIndexCoord, mWorldCoord );
  oRAS[0] = VECTOR_ELT( mWorldCoord, 1 );
  oRAS[1] = VECTOR_ELT( mWorldCoord, 2 );
  oRAS[2] = VECTOR_ELT( mWorldCoord, 3 );
}

bool 
VolumeCollection::IsRASInMRIBounds ( float iRAS[3] ) {

  if( NULL != mMRI ) {
      return ( iRAS[0] > mMRI->xstart && iRAS[0] < mMRI->xend &&
	       iRAS[1] > mMRI->ystart && iRAS[1] < mMRI->yend &&
	       iRAS[2] > mMRI->zstart && iRAS[2] < mMRI->zend );
  } else {
    return false;
  }
}

bool 
VolumeCollection::IsMRIIndexInMRIBounds ( int iIndex[3] ) {

  if( NULL != mMRI ) {
      return ( iIndex[0] >= 0 && iIndex[0] < mMRI->width &&
	       iIndex[1] >= 0 && iIndex[1] < mMRI->height &&
	       iIndex[2] >= 0 && iIndex[2] < mMRI->depth );
  } else {
    return false;
  }
}

float 
VolumeCollection::GetMRINearestValueAtRAS ( float iRAS[3] ) {

  Real value = 0;
  if( NULL != mMRI ) {
    float index[3];
    RASToMRIIndex( iRAS, index );
    MRIsampleVolumeType( mMRI, index[0], index[1], index[2],
			 &value, SAMPLE_NEAREST );
  }
  return (float)value;
}

float 
VolumeCollection::GetMRITrilinearValueAtRAS ( float iRAS[3] ) {

  Real value = 0;
  if( NULL != mMRI ) {
    float index[3];
    RASToMRIIndex( iRAS, index );
    MRIsampleVolumeType( mMRI, index[0], index[1], index[2],
			 &value, SAMPLE_TRILINEAR );
  }
  return (float)value;
}

float 
VolumeCollection::GetMRISincValueAtRAS ( float iRAS[3] ) {
  
  Real value = 0;
  if( NULL != mMRI ) {
    float index[3];
    RASToMRIIndex( iRAS, index );
    MRIsampleVolumeType( mMRI, index[0], index[1], index[2],
			 &value, SAMPLE_SINC );
  }
  return (float)value;
}

void
VolumeCollection::SetMRIValueAtRAS ( float iRAS[3], float iValue ) {
  Real value = 0;
  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    switch( mMRI->type ) {
    default:
      break ;
    case MRI_UCHAR:
      MRIvox( mMRI, index[0], index[1], index[2] ) = (BUFTYPE) iValue;
      break ;
    case MRI_SHORT:
      MRISvox( mMRI, index[0], index[1], index[2] ) = (short) iValue;
      break ;
    case MRI_FLOAT:
      MRIFvox( mMRI, index[0], index[1], index[2] ) = (float) iValue;
      break ;
    case MRI_LONG:
      MRILvox( mMRI, index[0], index[1], index[2] ) = (long) iValue;
      break ;
    case MRI_INT:
      MRIIvox( mMRI, index[0], index[1], index[2] ) = (int) iValue;
      break ;
    }
  }
}

TclCommandListener::TclCommandResult 
VolumeCollection::DoListenToTclCommand ( char* isCommand, int iArgc, char** iasArgv ) {

  // SetVolumeCollectionFileName <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "SetVolumeCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      string fnVolume = iasArgv[2];
      SetFileName( fnVolume );
    }
  }
  
  // GetVolumeCollectionFileName <collectionID>
  if( 0 == strcmp( isCommand, "GetVolumeCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      sReturnFormat = "s";
      sReturnValues = mfnMRI;
    }
  }
  
  return DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

ScubaROI*
VolumeCollection::DoNewROI () {

  ScubaROIVolume* roi = new ScubaROIVolume();

  if( NULL != mMRI ) {
    int bounds[3];
    bounds[0] = mMRI->width;
    bounds[1] = mMRI->height;
    bounds[2] = mMRI->depth;
    
    roi->SetROIBounds( bounds );
  }
  
  return roi;
}

void 
VolumeCollection::SelectRAS ( float iRAS[3] ) {

  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
    volumeROI->SelectVoxel( index );
  }
}

void 
VolumeCollection::UnselectRAS ( float iRAS[3] ) {
  
  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
    volumeROI->UnselectVoxel( index );
  }
}

bool 
VolumeCollection::IsRASSelected ( float iRAS[3], int oColor[3] ) {

  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );

    bool bSelected = false;
    bool bFirstColor = true;

    map<int,ScubaROI*>::iterator tIDROI;
    for( tIDROI = mROIMap.begin();
	 tIDROI != mROIMap.end(); ++tIDROI ) {
      int roiID = (*tIDROI).first;
      
      ScubaROI* roi = &ScubaROI::FindByID( roiID );
      //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
      ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
      if( volumeROI->IsVoxelSelected( index ) ) {
	bSelected = true;
	int color[3];
	volumeROI->GetDrawColor( color );
	if( bFirstColor ) {
	  oColor[0] = color[0];
	  oColor[1] = color[1];
	  oColor[2] = color[2];
	  bFirstColor = false;
	} else {
	  oColor[0] = (int) (((float)color[0] * 0.5) +
			     ((float)oColor[0] * 0.5));
	  oColor[1] = (int) (((float)color[1] * 0.5) +
			     ((float)oColor[1] * 0.5));
	  oColor[2] = (int) (((float)color[2] * 0.5) +
			     ((float)oColor[2] * 0.5));
	}
      }
    }
    
    return bSelected;

  } else {
    return false;
  }
}

void 
VolumeCollection::InitEdgeVolume () {

  if( NULL != mMRI ) {

    if( NULL != mEdgeVoxels ) {
      delete mEdgeVoxels;
    }
    
    mEdgeVoxels = 
      new Volume3<bool>( mMRI->width, mMRI->height, mMRI->depth, false );
  }
}

void 
VolumeCollection::MarkRASEdge ( float iRAS[3] ) {

  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    mEdgeVoxels->Set( index[2], index[1], index[0], true );
  }
}

void 
VolumeCollection::UnmarkRASEdge ( float iRAS[3] ) {

  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    mEdgeVoxels->Set( index[2], index[1], index[0], false );
  }
}

bool 
VolumeCollection::IsRASEdge ( float iRAS[3] ) {

  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    return mEdgeVoxels->Get( index[2], index[1], index[0] );
  } else {
    return false;
  }
}

VolumeCollectionFlooder::VolumeCollectionFlooder () {
  mVolume = NULL;
  mParams = NULL;
}

VolumeCollectionFlooder::~VolumeCollectionFlooder () {
}

VolumeCollectionFlooder::Params::Params () {
  mbStopAtEdges = true;
  mbStopAtROIs  = true;
  mb3D          = true;
  mFuzziness    = 1;
  mbDiagonal    = false;
}



void 
VolumeCollectionFlooder::DoVoxel ( float iRAS[3] ) {

}

bool 
VolumeCollectionFlooder::CompareVoxel ( float iRAS[3] ) {
return true;
}

void
VolumeCollectionFlooder::Flood ( VolumeCollection& iVolume, 
				 float iRASSeed[3], Params& iParams ) {
  
  mVolume = &iVolume;
  mParams = &iParams;

  Volume3<bool> bVisited( iVolume.mMRI->width, 
			  iVolume.mMRI->height, 
			  iVolume.mMRI->depth, false );

  // Push the seed onto the list. 
  int index[3];
  iVolume.RASToMRIIndex( iRASSeed, index );
  Point3<int> seed( index );
  vector<Point3<int> > points;
  points.push_back( seed );
  while( points.size() > 0 ) {
    
    Point3<int> point = points.back();
    points.pop_back();

    int index[3];
    index[0] = point.x();
    index[1] = point.y();
    index[2] = point.z();

    if( !iVolume.IsMRIIndexInMRIBounds( index ) ) {
      continue;
    }

    if( bVisited.Get( point.z(), point.y(), point.x() ) ) {
      continue;
    }
    bVisited.Set( point.z(), point.y(), point.x(), true );

    // Get RAS.
    float ras[3];
    iVolume.MRIIndexToRAS( index, ras );

    // Check if this is an edge or an ROI. If so, and our params say
    // not go to here, continue.
    if( iParams.mbStopAtEdges ) {
      if( iVolume.IsRASEdge( ras ) ) {
	continue;
      }
    }
    if( iParams.mbStopAtROIs ) {
      int color[3];
      if( iVolume.IsRASSelected( ras, color ) ) {
	continue;
      }
    }
    
    // Check max distance.

    // Call the user compare function to give them a chance to bail.
    if( !this->CompareVoxel( ras ) ) {
      continue;
    }
    
    // Call the user function.
    this->DoVoxel( ras );

    // Add adjacent voxels.
    int nBeginX = MAX( index[0] - 1, 0 );
    int nEndX   = MIN( index[0] + 1, iVolume.mMRI->width );
    int nBeginY = MAX( index[1] - 1, 0 );
    int nEndY   = MIN( index[1] + 1, iVolume.mMRI->height );
    int nBeginZ = MAX( index[2] - 1, 0 );
    int nEndZ   = MIN( index[2] + 1, iVolume.mMRI->depth );
    Point3<int> newPoint;
    if( iParams.mbDiagonal ) {
      for( int nZ = nBeginZ; nZ <= nEndZ; nZ++ ) {
	for( int nY = nBeginY; nY <= nEndY; nY++ ) {
	  for( int nX = nBeginX; nX <= nEndX; nX++ ) {
	    newPoint.Set( nX, nY, nZ );
	    points.push_back( newPoint );
	  }
	}
      }
    } else {
      newPoint.Set( nBeginX, index[1], index[2] );
      points.push_back( newPoint );
      newPoint.Set( nEndX, index[1], index[2] );
      points.push_back( newPoint );
      newPoint.Set( index[0], nBeginY, index[2] );
      points.push_back( newPoint );
      newPoint.Set( index[0], nEndY, index[2] );
      points.push_back( newPoint );
      newPoint.Set( index[0], index[1], nBeginZ );
      points.push_back( newPoint );
      newPoint.Set( index[0], index[1], nEndZ );
      points.push_back( newPoint );
    }
  }

  mVolume = NULL;
  mParams = NULL;
}
