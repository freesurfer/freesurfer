#include <stdexcept>
#include "VolumeCollection.h"
#include "DataManager.h"

using namespace std;


VolumeCollection::VolumeCollection () :
  DataCollection() {
  mMRI = NULL;
  mWorldToIndexMatrix = NULL;
  mWorldCoord = VectorAlloc( 4, MATRIX_REAL );
  mIndexCoord = VectorAlloc( 4, MATRIX_REAL );

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
  }

  return mMRI; 
}

void
VolumeCollection::UpdateMRIValueRange () {

  if( NULL != mMRI ) {
    MRIvalRange( mMRI, &mMRIMinValue, &mMRIMaxValue );
  }
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], int oIndex[3] ) {
  
  Real index[3];
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
  
  Real index[3];
  VECTOR_ELT( mWorldCoord, 1 ) = iRAS[0];
  VECTOR_ELT( mWorldCoord, 2 ) = iRAS[1];
  VECTOR_ELT( mWorldCoord, 3 ) = iRAS[2];
  VECTOR_ELT( mWorldCoord, 4 ) = 1.0;
  MatrixMultiply( mWorldToIndexMatrix, mWorldCoord, mIndexCoord );
  oIndex[0] = VECTOR_ELT( mIndexCoord, 1 );
  oIndex[1] = VECTOR_ELT( mIndexCoord, 2 );
  oIndex[2] = VECTOR_ELT( mIndexCoord, 3 );
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
    ScubaROIVolume* roi = 
      dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( mSelectedROIID ));
    roi->SelectVoxel( index );
  }
}

void 
VolumeCollection::UnselectRAS ( float iRAS[3] ) {
  
  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    ScubaROIVolume* roi = 
      dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( mSelectedROIID ));
    roi->UnselectVoxel( index );
  }
}

bool 
VolumeCollection::IsRASSelected ( float iRAS[3] ) {

  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );
    ScubaROIVolume* roi = 
      dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( mSelectedROIID ));
    return roi->IsVoxelSelected( index );
  } else {
    return false;
  }
}

