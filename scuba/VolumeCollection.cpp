#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include <vector>
#include "VolumeCollection.h"
#include "DataManager.h"
#include "Point3.h"
#include "error.h"

using namespace std;


VolumeCollection::VolumeCollection () :
  DataCollection() {
  mMRI = NULL;
  mMagnitudeMRI = NULL;
  mWorldToIndexMatrix = NULL;
  mIndexToWorldMatrix = NULL;
  mWorldCoord = VectorAlloc( 4, MATRIX_REAL );
  mIndexCoord = VectorAlloc( 4, MATRIX_REAL );
  mEdgeVoxels = NULL;
  mSelectedVoxels = NULL;
  mWorldToIndexCache = NULL;
  mVoxelSize[0] = mVoxelSize[1] = mVoxelSize[2] = 0;
  mOneOverVoxelSize[0] = mOneOverVoxelSize[1] = mOneOverVoxelSize[2] = 0;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollectionFileName", 2, 
			 "collectionID fileName", 
			 "Sets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "GetVolumeCollectionFileName", 1, 
			 "collectionID", 
			 "Gets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "WriteVolumeROIToLabel", 3, 
			 "collectionID roiID fileName", 
			 "Writes an ROI to a label file." );
  commandMgr.AddCommand( *this, "NewVolumeROIFromLabel", 2, 
			 "collectionID fileName", 
			 "Creates an ROI from a label file and returns the "
			 "ID of the new ROI." );
  commandMgr.AddCommand( *this, "WriteVolumeROIsToSegmentation", 2, 
			 "collectionID fileName", 
			 "Writes a series of structure ROIs to a "
			 "segmentation volume." );
  
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

    mWorldToIndexMatrix = voxelFromSurfaceRAS_( mMRI );
    mIndexToWorldMatrix = surfaceRASFromVoxel_( mMRI );

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

    // Init the edge and selection volume.
    InitEdgeVolume();
    InitSelectionVolume();

    UpdateRASBounds();

    mVoxelSize[0] = mMRI->xsize;
    mVoxelSize[1] = mMRI->ysize;
    mVoxelSize[2] = mMRI->zsize;
    mOneOverVoxelSize[0] = 1.0 / mVoxelSize[0];
    mOneOverVoxelSize[1] = 1.0 / mVoxelSize[1];
    mOneOverVoxelSize[2] = 1.0 / mVoxelSize[2];

    try { 
      // CalcWorldToIndexCache(); 
    }
    catch(...) { 
      DebugOutput( << "Failed while calcing world to index cache  " );
    }
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
VolumeCollection::GetMRIMagnitudeMinValue () { 
  if( NULL == mMagnitudeMRI ) {
    MakeMagnitudeVolume();
  }
  return mMRIMagMinValue; 
}

float
VolumeCollection::GetMRIMagnitudeMaxValue () { 
  if( NULL == mMagnitudeMRI ) {
    MakeMagnitudeVolume();
  }
  return mMRIMagMaxValue; 
}


void
VolumeCollection::UpdateRASBounds () {

  if( NULL != mMRI ) {
    
    int minIndex[3];
    minIndex[0] = 0;
    minIndex[1] = 0;
    minIndex[2] = 0;
    int maxIndex[3];
    maxIndex[0] = mMRI->width - 1;
    maxIndex[1] = mMRI->height - 1;
    maxIndex[2] = mMRI->depth - 1;
    float rasCorner1[3];
    float rasCorner2[3];
    MRIIndexToRAS( minIndex, rasCorner1 );
    MRIIndexToRAS( maxIndex, rasCorner2 );
    mMinRASBounds[0] = MIN( rasCorner1[0], rasCorner2[0] );
    mMinRASBounds[1] = MIN( rasCorner1[1], rasCorner2[1] );
    mMinRASBounds[2] = MIN( rasCorner1[2], rasCorner2[2] );
    mMaxRASBounds[0] = MAX( rasCorner1[0], rasCorner2[0] );
    mMaxRASBounds[1] = MAX( rasCorner1[1], rasCorner2[1] );
    mMaxRASBounds[2] = MAX( rasCorner1[2], rasCorner2[2] );
  }
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], int oIndex[3] ) {
  

  float rasX = iRAS[0];
  float rasY = iRAS[1];
  float rasZ = iRAS[2];

  float m11 = *MATRIX_RELT(mWorldToIndexMatrix,1,1);
  float m12 = *MATRIX_RELT(mWorldToIndexMatrix,1,2);
  float m13 = *MATRIX_RELT(mWorldToIndexMatrix,1,3);
  float m14 = *MATRIX_RELT(mWorldToIndexMatrix,1,4);

  float m21 = *MATRIX_RELT(mWorldToIndexMatrix,2,1);
  float m22 = *MATRIX_RELT(mWorldToIndexMatrix,2,2);
  float m23 = *MATRIX_RELT(mWorldToIndexMatrix,2,3);
  float m24 = *MATRIX_RELT(mWorldToIndexMatrix,2,4);
  
  float m31 = *MATRIX_RELT(mWorldToIndexMatrix,3,1);
  float m32 = *MATRIX_RELT(mWorldToIndexMatrix,3,2);
  float m33 = *MATRIX_RELT(mWorldToIndexMatrix,3,3);
  float m34 = *MATRIX_RELT(mWorldToIndexMatrix,3,4);
  
  float a = m11 * rasX;
  float b = m12 * rasY;
  float c = m13 * rasZ;
  float sum0 = a + b + c + m14;
  int sumI0 = (int) sum0;
  
  float d = m21 * rasX;
  float e = m22 * rasY;
  float f = m23 * rasZ;
  float sum1 = d + e + f + m24;
  int sumI1 = (int) sum1;

  float g = m31 * rasX;
  float h = m32 * rasY;
  float i = m33 * rasZ;
  float sum2 = g + h + i + m34;
  int sumI2 = (int) sum2;


  oIndex[0] = sumI0;
  oIndex[1] = sumI1;
  oIndex[2] = sumI2;


#if 0

  oIndex[0] = (int) (
		     *MATRIX_RELT(mWorldToIndexMatrix,1,1) * iRAS[0] +
		     *MATRIX_RELT(mWorldToIndexMatrix,1,2) * iRAS[1] +
		     *MATRIX_RELT(mWorldToIndexMatrix,1,3) * iRAS[2] +
		     *MATRIX_RELT(mWorldToIndexMatrix,1,4) );
  oIndex[1] = (int) (
		     *MATRIX_RELT(mWorldToIndexMatrix,2,1) * iRAS[0] +
		     *MATRIX_RELT(mWorldToIndexMatrix,2,2) * iRAS[1] +
		     *MATRIX_RELT(mWorldToIndexMatrix,2,3) * iRAS[2] +
		     *MATRIX_RELT(mWorldToIndexMatrix,2,4) );
  oIndex[2] = (int) (
		     *MATRIX_RELT(mWorldToIndexMatrix,3,1) * iRAS[0] +
		     *MATRIX_RELT(mWorldToIndexMatrix,3,2) * iRAS[1] +
		     *MATRIX_RELT(mWorldToIndexMatrix,3,3) * iRAS[2] +
		     *MATRIX_RELT(mWorldToIndexMatrix,3,4) );
#endif
  
#if 0
  int cacheIndex[3];
  WorldToIndexCacheIndex( iRAS, cacheIndex );
  Point3<int> index = 
    mWorldToIndexCache->Get( cacheIndex[0], cacheIndex[1], cacheIndex[2] );

  oIndex[0] = index.x();
  oIndex[1] = index.y();
  oIndex[2] = index.z();

#endif
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], float oIndex[3] ) {

#if 0  
  VECTOR_ELT( mWorldCoord, 1 ) = iRAS[0];
  VECTOR_ELT( mWorldCoord, 2 ) = iRAS[1];
  VECTOR_ELT( mWorldCoord, 3 ) = iRAS[2];
  VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
  MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
  oIndex[0] = VECTOR_ELT( mIndexCoord, 1 );
  oIndex[1] = VECTOR_ELT( mIndexCoord, 2 );
  oIndex[2] = VECTOR_ELT( mIndexCoord, 3 );
#else
  oIndex[0] = 
    *MATRIX_RELT(mWorldToIndexMatrix,1,1) * iRAS[0] +
    *MATRIX_RELT(mWorldToIndexMatrix,1,2) * iRAS[1] +
    *MATRIX_RELT(mWorldToIndexMatrix,1,3) * iRAS[2] +
    *MATRIX_RELT(mWorldToIndexMatrix,1,4);
  oIndex[1] =
    *MATRIX_RELT(mWorldToIndexMatrix,2,1) * iRAS[0] +
    *MATRIX_RELT(mWorldToIndexMatrix,2,2) * iRAS[1] +
    *MATRIX_RELT(mWorldToIndexMatrix,2,3) * iRAS[2] +
    *MATRIX_RELT(mWorldToIndexMatrix,2,4);
  oIndex[2] =
    *MATRIX_RELT(mWorldToIndexMatrix,3,1) * iRAS[0] +
    *MATRIX_RELT(mWorldToIndexMatrix,3,2) * iRAS[1] +
    *MATRIX_RELT(mWorldToIndexMatrix,3,3) * iRAS[2] +
    *MATRIX_RELT(mWorldToIndexMatrix,3,4);

#endif
}

void
VolumeCollection::MRIIndexToRAS ( int iIndex[3], float oRAS[3] ) {
  
#if 0
  VECTOR_ELT( mIndexCoord, 1 ) = iIndex[0];
  VECTOR_ELT( mIndexCoord, 2 ) = iIndex[1];
  VECTOR_ELT( mIndexCoord, 3 ) = iIndex[2];
  VECTOR_ELT( mIndexCoord, 4 ) = 1.0;
  MatrixMultiply( mIndexToWorldMatrix, mIndexCoord, mWorldCoord );
  oRAS[0] = VECTOR_ELT( mWorldCoord, 1 );
  oRAS[1] = VECTOR_ELT( mWorldCoord, 2 );
  oRAS[2] = VECTOR_ELT( mWorldCoord, 3 );
#else
  oRAS[0] = 
    *MATRIX_RELT(mIndexToWorldMatrix,1,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,4);
  oRAS[1] =
    *MATRIX_RELT(mIndexToWorldMatrix,2,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,4);
  oRAS[2] =
    *MATRIX_RELT(mIndexToWorldMatrix,3,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,4);

#endif
}

void
VolumeCollection::MRIIndexToRAS ( float iIndex[3], float oRAS[3] ) {
  
#if 0
  VECTOR_ELT( mIndexCoord, 1 ) = iIndex[0];
  VECTOR_ELT( mIndexCoord, 2 ) = iIndex[1];
  VECTOR_ELT( mIndexCoord, 3 ) = iIndex[2];
  VECTOR_ELT( mIndexCoord, 4 ) = 1.0;
  MatrixMultiply( mIndexToWorldMatrix, mIndexCoord, mWorldCoord );
  oRAS[0] = VECTOR_ELT( mWorldCoord, 1 );
  oRAS[1] = VECTOR_ELT( mWorldCoord, 2 );
  oRAS[2] = VECTOR_ELT( mWorldCoord, 3 );
#else
  oRAS[0] = 
    *MATRIX_RELT(mIndexToWorldMatrix,1,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,1,4);
  oRAS[1] =
    *MATRIX_RELT(mIndexToWorldMatrix,2,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,2,4);
  oRAS[2] =
    *MATRIX_RELT(mIndexToWorldMatrix,3,1) * iIndex[0] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,2) * iIndex[1] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,3) * iIndex[2] +
    *MATRIX_RELT(mIndexToWorldMatrix,3,4);

#endif
}

bool 
VolumeCollection::IsRASInMRIBounds ( float iRAS[3] ) {

  if( NULL != mMRI ) {
      return ( iRAS[0] >= mMinRASBounds[0] && iRAS[0] <= mMaxRASBounds[0] &&
	       iRAS[1] >= mMinRASBounds[1] && iRAS[1] <= mMaxRASBounds[1] &&
	       iRAS[2] >= mMinRASBounds[2] && iRAS[2] <= mMaxRASBounds[2] );
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

    int index[3];
    RASToMRIIndex( iRAS, index );

    switch( mMRI->type ) {
    case MRI_UCHAR:
      value = (float)MRIvox(mMRI, index[0], index[1], index[2] );
      break ;
    case MRI_SHORT:
      value = (float)MRISvox(mMRI, index[0], index[1], index[2] );
      break ;
    case MRI_INT:
      value = (float)MRIIvox(mMRI, index[0], index[1], index[2] );
      break ;
    case MRI_FLOAT:
      value = MRIFvox(mMRI, index[0], index[1], index[2] );
      break ;
    default:
      value = 0;
    }
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

void
VolumeCollection::MakeMagnitudeVolume () {

  if( NULL != mMRI ) {
    mMagnitudeMRI = 
      MRIallocSequence( mMRI->width, mMRI->height, mMRI->depth,
			MRI_FLOAT, mMRI->nframes );
    MRI* gradMRI = MRIsobel( mMRI, NULL, mMagnitudeMRI );
    MRIfree( &gradMRI );
    
    MRIvalRange( mMagnitudeMRI, &mMRIMagMinValue, &mMRIMagMaxValue );
  }
}

float 
VolumeCollection::GetMRIMagnitudeValueAtRAS ( float iRAS[3] ) {

  Real value = 0;

  // If we don't have the magnitude volume, calculate it.
  if( NULL == mMagnitudeMRI ) {
    MakeMagnitudeVolume();
  }

  // Get the value.
  if( NULL != mMagnitudeMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    value = MRIFvox( mMagnitudeMRI, index[0], index[1], index[2] );
  }
  return (float)value;
}

TclCommandListener::TclCommandResult 
VolumeCollection::DoListenToTclCommand ( char* isCommand, 
					 int iArgc, char** iasArgv ) {

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

  // WriteVolumeROIToLabel <collectionID> <roiID> <fileName>
  if( 0 == strcmp( isCommand, "WriteVolumeROIToLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
     int roiID = strtol(iasArgv[2], (char**)NULL, 10);
     if( ERANGE == errno ) {
       sResult = "bad roi ID";
       return error;
     }
     
     try {
       WriteROIToLabel( roiID, string(iasArgv[3]) );
     }
     catch(...) {
       sResult = "That ROI doesn't belong to this collection";
       return error;
     }
    }
  }
  
  // NewVolumeROIFromLabel <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "NewVolumeROIFromLabel" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
       int roiID = NewROIFromLabel( string(iasArgv[2]) );
       stringstream ssReturnValues;
       ssReturnValues << roiID;
       sReturnValues = ssReturnValues.str();
       sReturnFormat = "i";
    }
  }
  
  // WriteVolumeROIsToSegmentation <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "WriteVolumeROIsToSegmentation" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
      WriteROIsToSegmentation( string(iasArgv[2]) );
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
VolumeCollection::InitSelectionVolume () {

  if( NULL != mMRI ) {

    if( NULL != mSelectedVoxels ) {
      delete mSelectedVoxels;
    }
    
    mSelectedVoxels = 
      new Volume3<bool>( mMRI->width, mMRI->height, mMRI->depth, false );
  }
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

    // Also mark this in the selection voxel.
    mSelectedVoxels->Set_Unsafe( index[0], index[1], index[2], true );
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


    // If there are no more ROIs with this voxel selected, unselect it
    // in the selection volume.
    bool bSelected = false;
    map<int,ScubaROI*>::iterator tIDROI;
    for( tIDROI = mROIMap.begin();
	 tIDROI != mROIMap.end(); ++tIDROI ) {
      int roiID = (*tIDROI).first;
      
      ScubaROI* roi = &ScubaROI::FindByID( roiID );
      //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
      ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
      if( volumeROI->IsVoxelSelected( index ) ) {
	bSelected = true;
	break;
      }
    }
    if( !bSelected ) {
      mSelectedVoxels->Set_Unsafe( index[0], index[1], index[2], false );
    }
  }
}

bool 
VolumeCollection::IsRASSelected ( float iRAS[3], int oColor[3] ) {

  // Check the selection volume cache first.
  int index[3];
  RASToMRIIndex( iRAS, index );
  if( !(mSelectedVoxels->Get_Unsafe( index[0], index[1], index[2] )) )
    return false;

  try {
    
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
	  oColor[0] = (int) (((float)color[0] * 0.5) + ((float)oColor[0]*0.5));
	  oColor[1] = (int) (((float)color[1] * 0.5) + ((float)oColor[1]*0.5));
	  oColor[2] = (int) (((float)color[2] * 0.5) + ((float)oColor[2]*0.5));
	}
      }
    }
    
    return bSelected;
  }
  catch(...) {
    return false;
  }
  
}


bool 
VolumeCollection::IsOtherRASSelected ( float iRAS[3], int iThisROIID ) {

  // Check the selectin volume cache first.
  int index[3];
  RASToMRIIndex( iRAS, index );
  if( !(mSelectedVoxels->Get_Unsafe( index[0], index[1], index[2] )) )
    return false;

  bool bSelected = false;
  
  map<int,ScubaROI*>::iterator tIDROI;
  for( tIDROI = mROIMap.begin();
       tIDROI != mROIMap.end(); ++tIDROI ) {
    int roiID = (*tIDROI).first;
    
    ScubaROI* roi = &ScubaROI::FindByID( roiID );
    if( roiID == iThisROIID ) {
	continue;
    }
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
    if( volumeROI->IsVoxelSelected( index ) ) {
      bSelected = true;
      
    }
  }
  
  return bSelected;
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
    mEdgeVoxels->Set_Unsafe( index[0], index[1], index[2], true );
  }
}

void 
VolumeCollection::UnmarkRASEdge ( float iRAS[3] ) {

  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    mEdgeVoxels->Set_Unsafe( index[0], index[1], index[2], false );
  }
}

bool 
VolumeCollection::IsRASEdge ( float iRAS[3] ) {

  if( NULL != mMRI ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    return mEdgeVoxels->Get_Unsafe( index[0], index[1], index[2] );
  } else {
    return false;
  }
}

void
VolumeCollection::GetRASPointsInCube ( float iCenterRAS[3], float iRadius,
				       bool ibBrushX, bool ibBrushY,
				       bool ibBrushZ,
				       list<Point3<float> >& oPoints ) {
  
  // Find out the RAS.bounds.
  float beginX = MAX( mMRI->xstart, iCenterRAS[0] - iRadius );
  float endX   = MIN( mMRI->xend,   iCenterRAS[0] + iRadius );
  float beginY = MAX( mMRI->ystart, iCenterRAS[1] - iRadius );
  float endY   = MIN( mMRI->yend,   iCenterRAS[1] + iRadius );
  float beginZ = MAX( mMRI->zstart, iCenterRAS[2] - iRadius );
  float endZ   = MIN( mMRI->zend,   iCenterRAS[2] + iRadius );

  // Limit according to our dimensions.
  if( !ibBrushX ) {
    beginX = endX = iCenterRAS[0];
  }
  if( !ibBrushY ) {
    beginY = endY = iCenterRAS[1];
  }
  if( !ibBrushZ ) {
    beginZ = endZ = iCenterRAS[2];
  }

  // Go through the RAS coords and step in half the voxel size
  // amount. Then add each to the list.
  for( float nZ = beginZ; nZ <= endZ; nZ += GetVoxelXSize()/2.0 ) {
    for( float nY = beginY; nY <= endY; nY += GetVoxelYSize()/2.0 ) {
      for( float nX = beginX; nX <= endX; nX += GetVoxelZSize()/2.0 ) {
	Point3<float> ras( nX, nY, nZ );
	oPoints.push_back( ras );
      }
    }
  }
}

void
VolumeCollection::GetRASPointsInSphere ( float iCenterRAS[3], float iRadius,
					 bool ibBrushX, bool ibBrushY,
					 bool ibBrushZ,
					 list<Point3<float> >& oPoints ) {
  // Find out the RAS.bounds.
  float beginX = MAX( mMRI->xstart, iCenterRAS[0] - iRadius );
  float endX   = MIN( mMRI->xend,   iCenterRAS[0] + iRadius );
  float beginY = MAX( mMRI->ystart, iCenterRAS[1] - iRadius );
  float endY   = MIN( mMRI->yend,   iCenterRAS[1] + iRadius );
  float beginZ = MAX( mMRI->zstart, iCenterRAS[2] - iRadius );
  float endZ   = MIN( mMRI->zend,   iCenterRAS[2] + iRadius );

  // Limit according to our dimensions.
  if( !ibBrushX ) {
    beginX = endX = iCenterRAS[0];
  }
  if( !ibBrushY ) {
    beginY = endY = iCenterRAS[1];
  }
  if( !ibBrushZ ) {
    beginZ = endZ = iCenterRAS[2];
  }

  // Go through the RAS coords and step in half the voxel size
  // amount. Check the distance for the sphere shape. Then add each to
  // the list.
  for( float nZ = beginZ; nZ <= endZ; nZ += GetVoxelXSize()/2.0 ) {
    for( float nY = beginY; nY <= endY; nY += GetVoxelYSize()/2.0 ) {
      for( float nX = beginX; nX <= endX; nX += GetVoxelZSize()/2.0 ) {

	float distance = sqrt( ((nX-iCenterRAS[0]) * (nX-iCenterRAS[0])) + 
			       ((nY-iCenterRAS[1]) * (nY-iCenterRAS[1])) + 
			       ((nZ-iCenterRAS[2]) * (nZ-iCenterRAS[2])) );

	if( distance > iRadius ) {
	  continue;
	}

	Point3<float> ras( nX, nY, nZ );
	oPoints.push_back( ras );
      }
    }
  }

}

void
VolumeCollection::WriteROIToLabel ( int iROIID, string ifnLabel ) {
  
  map<int,ScubaROI*>::iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if( tIDROI != mROIMap.end() ) {
    ScubaROIVolume* roi = (ScubaROIVolume*)(*tIDROI).second;

    int bounds[3];
    roi->GetROIBounds( bounds );

    int cSelectedVoxels = roi->NumSelectedVoxels();
    if( 0 == cSelectedVoxels ) {
      throw runtime_error( "No selected voxels." );
    }

    char* fnLabel = strdup( ifnLabel.c_str() );
    LABEL* label = LabelAlloc( cSelectedVoxels, NULL, fnLabel );
    if( NULL == label ) {
      throw runtime_error( "Couldn't allocate label" );
    }
    label->n_points = cSelectedVoxels;

    int nPoint = 0;
    int voxel[3];
    for( voxel[2] = 0; voxel[2] < bounds[2]; voxel[2]++ ) {
      for( voxel[1] = 0; voxel[1] < bounds[1]; voxel[1]++ ) {
	for( voxel[0] = 0; voxel[0] < bounds[0]; voxel[0]++ ) {
	  
	  if( roi->IsVoxelSelected( voxel ) ) {

	    float ras[3];
	    MRIIndexToRAS( voxel, ras );

	    label->lv[nPoint].x = ras[0];
	    label->lv[nPoint].y = ras[1];
	    label->lv[nPoint].z = ras[2];
	    label->lv[nPoint].stat = GetMRINearestValueAtRAS( ras );
	    label->lv[nPoint].vno = -1;
	    label->lv[nPoint].deleted = false;

	    nPoint++;
	  }
	}
      }
    }

    int error = LabelWrite( label, fnLabel );
    if( NO_ERROR != error ) {
      throw runtime_error( "Couldn't write label" );
    }

    free( fnLabel );

  } else {
    throw runtime_error( "ROI doesn't belong to this collection" );
  }
}

int 
VolumeCollection::NewROIFromLabel ( string ifnLabel ) {

  char* fnLabel = strdup( ifnLabel.c_str() );
  LABEL* label = LabelRead( NULL, fnLabel );
  free( fnLabel );
  if( NULL == label ) {
    throw runtime_error( "Couldn't read label" );
  }

  ScubaROIVolume* volumeROI = NULL;
  try { 
    int roiID = NewROI();
    ScubaROI* roi = &ScubaROI::FindByID( roiID );
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    volumeROI = (ScubaROIVolume*)roi;
  }
  catch(...) {
    throw runtime_error( "Couldn't make ROI" );
  }

  for( int nPoint = 0; nPoint < label->n_points; nPoint++ ) {

    float ras[3];
    ras[0] = label->lv[nPoint].x;
    ras[1] = label->lv[nPoint].y;
    ras[2] = label->lv[nPoint].z;

    int index[3];
    RASToMRIIndex( ras, index );

    volumeROI->SelectVoxel( index );
  }
 
  LabelFree( &label );

  return volumeROI->GetID();
}

void
VolumeCollection::WriteROIsToSegmentation ( string ifnVolume ) {
  

  // Create a volume of the same size as our own.
  MRI* segVolume = MRIallocSequence( mMRI->width, mMRI->height, mMRI->depth, 
				     MRI_UCHAR, mMRI->nframes );
  if( NULL == segVolume ) {
    throw runtime_error( "Couldn't create seg volume" );
  }

  // Go through the volume...
  int index[3];
  for( index[2] = 0; index[2] < mMRI->depth; index[2]++ ) {
    for( index[1] = 0; index[1] < mMRI->height; index[1]++ ) {
      for( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {
	
	// For each of our ROIs, if one is selected here and if it's a
	// structure ROI, set the value of the seg volume to the
	// structure index.
	map<int,ScubaROI*>::iterator tIDROI;
	for( tIDROI = mROIMap.begin();
	     tIDROI != mROIMap.end(); ++tIDROI ) {
	  int roiID = (*tIDROI).first;
	  
	  ScubaROI* roi = &ScubaROI::FindByID( roiID );
	  //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
	  ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;
	  if( volumeROI->GetType() == ScubaROI::Structure &&
	      volumeROI->IsVoxelSelected( index ) ) {

	    MRIvox( segVolume, index[0], index[1], index[2] ) = 
	      (BUFTYPE) volumeROI->GetStructure();

	  }
	}
      }
    }
  }

  // Write the volume.
  char* fnVolume = strdup( ifnVolume.c_str() );
  int error = MRIwrite( segVolume, fnVolume );
  free( fnVolume );
  if( NO_ERROR != error ) {
    throw runtime_error( "Couldn't write segmentation." );
  }

  MRIfree( &segVolume );
}

void
VolumeCollection::CalcWorldToIndexCache () {

  float* minBounds = GetMinRASBounds();
  float* maxBounds = GetMaxRASBounds();

  int zX = (int) ceil((maxBounds[0]-minBounds[0]+1.0) * (1.0/GetVoxelXSize()));
  int zY = (int) ceil((maxBounds[1]-minBounds[1]+1.0) * (1.0/GetVoxelYSize()));
  int zZ = (int) ceil((maxBounds[2]-minBounds[2]+1.0) * (1.0/GetVoxelZSize()));

  if( NULL != mWorldToIndexCache ) {
    delete( mWorldToIndexCache );
  }

  mWorldToIndexCache = 
    new Volume3<Point3<int> >( zX, zY, zZ, Point3<int>(0,0,0) );

  DebugOutput( << "Voxel size: " << mVoxelSize[0] << ", " << mVoxelSize[1]
	       << ", " << mVoxelSize[2] );
  DebugOutput( << "Cachel volume size: " << zX << ", " << zY << ", " << zZ );

  Point3<int> index;
  float RAS[3];
  int cacheIndex[3];
  for( RAS[2] = minBounds[2]; RAS[2] < maxBounds[2]; 
       RAS[2] += GetVoxelZSize() ) {
    for( RAS[1] = minBounds[1]; RAS[1] < maxBounds[1]; 
	 RAS[1] += GetVoxelYSize() ) {
      for( RAS[0] = minBounds[0]; RAS[0] < maxBounds[0]; 
	   RAS[0] += GetVoxelXSize() ) {
	
	index.Set( (int) rint(*MATRIX_RELT(mWorldToIndexMatrix,1,1) * RAS[0] +
			      *MATRIX_RELT(mWorldToIndexMatrix,1,2) * RAS[1] +
			      *MATRIX_RELT(mWorldToIndexMatrix,1,3) * RAS[2] +
			      *MATRIX_RELT(mWorldToIndexMatrix,1,4) ),
		   (int) rint(*MATRIX_RELT(mWorldToIndexMatrix,2,1) * RAS[0] +
			      *MATRIX_RELT(mWorldToIndexMatrix,2,2) * RAS[1] +
			      *MATRIX_RELT(mWorldToIndexMatrix,2,3) * RAS[2] +
			      *MATRIX_RELT(mWorldToIndexMatrix,2,4) ),
		   (int) rint(*MATRIX_RELT(mWorldToIndexMatrix,3,1) * RAS[0] +
			      *MATRIX_RELT(mWorldToIndexMatrix,3,2) * RAS[1] +
			      *MATRIX_RELT(mWorldToIndexMatrix,3,3) * RAS[2] +
			      *MATRIX_RELT(mWorldToIndexMatrix,3,4) ) );

	WorldToIndexCacheIndex( RAS, cacheIndex );

	mWorldToIndexCache->Set_Unsafe( cacheIndex[0], cacheIndex[1], 
					cacheIndex[2], index );
	

      }
    }
  }
}

void
VolumeCollection::WorldToIndexCacheIndex ( float const iRAS[3],
					   int oCacheIndex[3] ) const {

  float iRASX = iRAS[0];
  float iRASY = iRAS[1];
  float iRASZ = iRAS[2];

  float oneOverVoxelX = mOneOverVoxelSize[0];
  float oneOverVoxelY = mOneOverVoxelSize[1];
  float oneOverVoxelZ = mOneOverVoxelSize[2];

  float mMinRASX = mMinRASBounds[0];
  float mMinRASY = mMinRASBounds[1];
  float mMinRASZ = mMinRASBounds[2];

  float cacheX = (iRASX - mMinRASX) * oneOverVoxelX;
  float cacheY = (iRASY - mMinRASY) * oneOverVoxelY;
  float cacheZ = (iRASZ - mMinRASZ) * oneOverVoxelZ;

  oCacheIndex[0] = (int) cacheX;
  oCacheIndex[1] = (int) cacheY;
  oCacheIndex[2] = (int) cacheZ;

#if 0
  oCacheIndex[0] = (int)((iRAS[0] - mMinRASBounds[0]) * mOneOverVoxelSize[0]);
  oCacheIndex[1] = (int)((iRAS[1] - mMinRASBounds[1]) * mOneOverVoxelSize[1]);
  oCacheIndex[2] = (int)((iRAS[2] - mMinRASBounds[2]) * mOneOverVoxelSize[2]);
#endif
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
  mbWorkPlaneX  = true;
  mbWorkPlaneY  = true;
  mbWorkPlaneZ  = true;
  mFuzziness    = 1;
  mbDiagonal    = false;
}


void
VolumeCollectionFlooder::DoBegin () {

}

void 
VolumeCollectionFlooder::DoEnd () {

}

bool
VolumeCollectionFlooder::DoStopRequested () {
  return false;
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

  this->DoBegin();

  Volume3<bool>* bVisited =
    new Volume3<bool>( iVolume.mMRI->width, 
		       iVolume.mMRI->height, 
		       iVolume.mMRI->depth, false );

  // Save the initial value.
  float seedValue = iVolume.GetMRINearestValueAtRAS( iRASSeed );

  // Push the seed onto the list. 
  Point3<int> seed;
  iVolume.RASToMRIIndex( iRASSeed, seed.xyz() );
  vector<Point3<int> > points;
  points.push_back( seed );
  while( points.size() > 0 &&
	 !this->DoStopRequested() ) {
    
    Point3<int> point = points.back();
    points.pop_back();


    if( !iVolume.IsMRIIndexInMRIBounds( point.xyz() ) ) {
      continue;
    }

    if( bVisited->Get_Unsafe( point.x(), point.y(), point.z() ) ) {
      continue;
    }
    bVisited->Set_Unsafe( point.x(), point.y(), point.z(), true );

    // Get RAS.
    float ras[3];
    iVolume.MRIIndexToRAS( point.xyz(), ras );

    // Check if this is an edge or an ROI. If so, and our params say
    // not go to here, continue.
    if( iParams.mbStopAtEdges ) {
      if( iVolume.IsRASEdge( ras ) ) {
	continue;
      }
    }
    if( iParams.mbStopAtROIs ) {
      if( iVolume.IsOtherRASSelected( ras, iVolume.GetSelectedROI() ) ) {
	continue;
      }
    }
    
    // Check max distance.
    if( iParams.mMaxDistance > 0 ) {
      float distance = sqrt( ((ras[0]-iRASSeed[0]) * (ras[0]-iRASSeed[0])) + 
			     ((ras[1]-iRASSeed[1]) * (ras[1]-iRASSeed[1])) + 
			     ((ras[2]-iRASSeed[2]) * (ras[2]-iRASSeed[2])) );
      if( distance > iParams.mMaxDistance ) {
	continue;
      }
    }

    // Check fuzziness.
    if( iParams.mFuzziness > 0 ) {
      float value = iVolume.GetMRINearestValueAtRAS( ras );
      if( fabs( value - seedValue ) > iParams.mFuzziness ) {
	continue;
      }
    }

    // Call the user compare function to give them a chance to bail.
    if( !this->CompareVoxel( ras ) ) {
      continue;
    }
    
    // Call the user function.
    this->DoVoxel( ras );

    // Add adjacent voxels.
    int beginX = MAX( point.x() - 1, 0 );
    int endX   = MIN( point.x() + 1, iVolume.mMRI->width );
    int beginY = MAX( point.y() - 1, 0 );
    int endY   = MIN( point.y() + 1, iVolume.mMRI->height );
    int beginZ = MAX( point.z() - 1, 0 );
    int endZ   = MIN( point.z() + 1, iVolume.mMRI->depth );
    if( !iParams.mb3D && iParams.mbWorkPlaneX ) {
      beginX = endX = point.x();
    }
    if( !iParams.mb3D && iParams.mbWorkPlaneY ) {
      beginY = endY = point.y();
    }
    if( !iParams.mb3D && iParams.mbWorkPlaneZ ) {
      beginZ = endZ = point.z();
    }
    Point3<int> newPoint;
    if( iParams.mbDiagonal ) {
      for( int nZ = beginZ; nZ <= endZ; nZ++ ) {
	for( int nY = beginY; nY <= endY; nY++ ) {
	  for( int nX = beginX; nX <= endX; nX++ ) {
	    newPoint.Set( nX, nY, nZ );
	    points.push_back( newPoint );
	  }
	}
      }
    } else {
      newPoint.Set( beginX, point.y(), point.z() );
      points.push_back( newPoint );
      newPoint.Set( endX, point.y(), point.z() );
      points.push_back( newPoint );
      newPoint.Set( point.x(), beginY, point.z() );
      points.push_back( newPoint );
      newPoint.Set( point.x(), endY, point.z() );
      points.push_back( newPoint );
      newPoint.Set( point.x(), point.y(), beginZ );
      points.push_back( newPoint );
      newPoint.Set( point.x(), point.y(), endZ );
      points.push_back( newPoint );
    }

  }

  mVolume = NULL;
  mParams = NULL;

  this->DoEnd();
}
