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

#if 0
    mWorldToIndexMatrix = extract_r_to_i( mMRI );
    mIndexToWorldMatrix = extract_i_to_r( mMRI );
#else 
    mWorldToIndexMatrix = voxelFromSurfaceRAS_( mMRI );
    mIndexToWorldMatrix = surfaceRASFromVoxel_( mMRI );
#endif
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
  } else {
    return 0;
  }
}

float 
VolumeCollection::GetVoxelYSize () {

  if( NULL != mMRI ) {
    return mMRI->ysize;
  } else {
    return 0;
  }
}

float 
VolumeCollection::GetVoxelZSize () {

  if( NULL != mMRI ) {
    return mMRI->zsize;
  } else {
    return 0;
  }
}


#if 1
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

#else

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], int oIndex[3] ) {
  oIndex[0] = iRAS[0]; oIndex[1] = iRAS[1]; oIndex[2] = iRAS[2];
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], float oIndex[3] ) {
  oIndex[0] = iRAS[0]; oIndex[1] = iRAS[1]; oIndex[2] = iRAS[2];
}

void
VolumeCollection::MRIIndexToRAS ( int iIndex[3], float oRAS[3] ) {
  oRAS[0] = iIndex[0]; oRAS[1] = iIndex[1]; oRAS[2] = iIndex[2];
}

void
VolumeCollection::MRIIndexToRAS ( float iIndex[3], float oRAS[3] ) {
  oRAS[0] = iIndex[0]; oRAS[1] = iIndex[1]; oRAS[2] = iIndex[2];
}
#endif

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

    try {

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
    }
    catch(...) {
      return false;
    }

  } else {
    return false;
  }
}


bool 
VolumeCollection::IsOtherRASSelected ( float iRAS[3], int iThisROIID ) {

  if( mSelectedROIID >= 0 ) {
    int index[3];
    RASToMRIIndex( iRAS, index );

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

void
VolumeCollection::GetRASPointsInCube ( float iCenterRAS[3], int iRadius,
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
VolumeCollection::GetRASPointsInSphere ( float iCenterRAS[3], int iRadius,
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

	if( distance > (float)iRadius ) {
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

  Volume3<bool> bVisited( iVolume.mMRI->width, 
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

    if( bVisited.Get( point.z(), point.y(), point.x() ) ) {
      continue;
    }
    bVisited.Set( point.z(), point.y(), point.x(), true );

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
