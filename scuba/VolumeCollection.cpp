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
  mEdgeVoxels = NULL;
  mSelectedVoxels = NULL;
  mVoxelSize[0] = mVoxelSize[1] = mVoxelSize[2] = 0;

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


    // Get our surfaceRAS -> index transform.
    Matrix44 m;
    //    MATRIX* voxelFromSurfaceRAS = voxelFromSurfaceRAS_( mMRI );
    MATRIX* voxelFromSurfaceRAS = extract_r_to_i( mMRI );
    m.SetMatrix( voxelFromSurfaceRAS );
    MatrixFree( &voxelFromSurfaceRAS );

    mDataToIndexTransform.SetMainTransform( m );

    CalcWorldToIndexTransform();



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


    mVoxelSize[0] = mMRI->xsize;
    mVoxelSize[1] = mMRI->ysize;
    mVoxelSize[2] = mMRI->zsize;

    try { 
      // CalcDataToIndexCache(); 
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
VolumeCollection::RASToMRIIndex ( float iRAS[3], int oIndex[3] ) {
  

#if 0
  int cacheIndex[3];
  DataToIndexCacheIndex( iRAS, cacheIndex );
  Point3<int> index = 
    mDataToIndexCache->Get( cacheIndex[0], cacheIndex[1], cacheIndex[2] );

  oIndex[0] = index.x();
  oIndex[1] = index.y();
  oIndex[2] = index.z();

#endif

  mWorldToIndexTransform.MultiplyVector3( iRAS, oIndex );
}

void
VolumeCollection::RASToMRIIndex ( float iRAS[3], float oIndex[3] ) {

  mWorldToIndexTransform.MultiplyVector3( iRAS, oIndex );
}

void
VolumeCollection::MRIIndexToRAS ( int iIndex[3], float oRAS[3] ) {
  
  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::MRIIndexToRAS ( float iIndex[3], float oRAS[3] ) {
  
  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::RASToDataRAS ( float iRAS[3], float oDataRAS[3] ) {
  
  mDataToWorldTransform->InvMultiplyVector3( iRAS, oDataRAS );
}

bool 
VolumeCollection::IsRASInMRIBounds ( float iRAS[3] ) {

  int index[3];
  RASToMRIIndex( iRAS, index );

  return IsMRIIndexInMRIBounds( index );
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

    if( index[0] >= 0 && index[0] < mMRI->width &&
	index[1] >= 0 && index[1] < mMRI->height &&
	index[2] >= 0 && index[2] < mMRI->depth ) {
      
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
  }
  
  return (float)value;
}

float 
VolumeCollection::GetMRITrilinearValueAtRAS ( float iRAS[3] ) {

  Real value = 0;
  if( NULL != mMRI ) {
    float index[3];
    RASToMRIIndex( iRAS, index );

    if( index[0] >= 0 && index[0] < mMRI->width &&
	index[1] >= 0 && index[1] < mMRI->height &&
	index[2] >= 0 && index[2] < mMRI->depth ) {
      
      MRIsampleVolumeType( mMRI, index[0], index[1], index[2],
			   &value, SAMPLE_TRILINEAR );
    }
  }
  return (float)value;
}

float 
VolumeCollection::GetMRISincValueAtRAS ( float iRAS[3] ) {
  
  Real value = 0;
  if( NULL != mMRI ) {
    float index[3];
    RASToMRIIndex( iRAS, index );

    if( index[0] >= 0 && index[0] < mMRI->width &&
	index[1] >= 0 && index[1] < mMRI->height &&
	index[2] >= 0 && index[2] < mMRI->depth ) {
      
      MRIsampleVolumeType( mMRI, index[0], index[1], index[2],
			   &value, SAMPLE_SINC );
    }
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

  DataChanged();
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

void
VolumeCollection::DoListenToMessage ( string isMessage, void* iData ) {
  
  if( isMessage == "transformChanged" ) {
    CalcWorldToIndexTransform();
  }

  DataCollection::DoListenToMessage( isMessage, iData );
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
  
  Point3<float> centerRAS( iCenterRAS );
  Point3<int> centerIdx;
  RASToMRIIndex( centerRAS.xyz(), centerIdx.xyz() );
  int cXVoxels = (int)ceil( iRadius / GetVoxelXSize() );
  int cYVoxels = (int)ceil( iRadius / GetVoxelYSize() );
  int cZVoxels = (int)ceil( iRadius / GetVoxelZSize() );
  
  if( !ibBrushX ) {
    cXVoxels = 0;
  }
  if( !ibBrushY ) {
    cYVoxels = 0;
  }
  if( !ibBrushZ ) {
    cZVoxels = 0;
  }

  oPoints.push_back( Point3<float>(iCenterRAS) );
  
  for( int nZ = -cZVoxels; nZ <= cZVoxels; nZ ++ ) {
    for( int nY = -cYVoxels; nY <= cYVoxels; nY ++ ) {
      for( int nX = -cXVoxels; nX <= cXVoxels; nX ++ ) {

	Point3<float> ras( iCenterRAS[0] + nX*GetVoxelXSize(),
			   iCenterRAS[1] + nY*GetVoxelYSize(),
			   iCenterRAS[2] + nZ*GetVoxelZSize() );
	
	if( fabs(ras[0] - iCenterRAS[0]) > iRadius ||
	    fabs(ras[1] - iCenterRAS[1]) > iRadius ||
	    fabs(ras[2] - iCenterRAS[2]) > iRadius ) {
	  continue;
	}

	if( IsRASInMRIBounds(ras.xyz()) )
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

  Point3<float> centerRAS( iCenterRAS );
  Point3<int> centerIdx;
  RASToMRIIndex( centerRAS.xyz(), centerIdx.xyz() );
  int cXVoxels = (int)ceil( iRadius / GetVoxelXSize() );
  int cYVoxels = (int)ceil( iRadius / GetVoxelYSize() );
  int cZVoxels = (int)ceil( iRadius / GetVoxelZSize() );
  
  if( !ibBrushX ) {
    cXVoxels = 0;
  }
  if( !ibBrushY ) {
    cYVoxels = 0;
  }
  if( !ibBrushZ ) {
    cZVoxels = 0;
  }

  oPoints.push_back( Point3<float>(iCenterRAS) );
  
  for( int nZ = -cZVoxels; nZ <= cZVoxels; nZ ++ ) {
    for( int nY = -cYVoxels; nY <= cYVoxels; nY ++ ) {
      for( int nX = -cXVoxels; nX <= cXVoxels; nX ++ ) {

	Point3<float> ras( iCenterRAS[0] + nX*GetVoxelXSize(),
			   iCenterRAS[1] + nY*GetVoxelYSize(),
			   iCenterRAS[2] + nZ*GetVoxelZSize() );
	
	float distance = 
	  sqrt( ((ras[0]-iCenterRAS[0]) * (ras[0]-iCenterRAS[0])) + 
		((ras[1]-iCenterRAS[1]) * (ras[1]-iCenterRAS[1])) + 
		((ras[2]-iCenterRAS[2]) * (ras[2]-iCenterRAS[2])) );

	if( distance > iRadius ) {
	  continue;
	}

	if( IsRASInMRIBounds(ras.xyz()) )
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
VolumeCollection::CalcWorldToIndexTransform () {

  // This makes it look like tkmedit when it loads a display
  // transform, is this right???
  //  mWorldToIndexTransform = mDataToIndexTransform;
  //  mWorldToIndexTransform.ApplyTransform( mDataToWorldTransform->Inverse() );
  
  //mWorldToIndexTransform = mDataToWorldTransform->Inverse();
  //mWorldToIndexTransform.ApplyTransform( mDataToIndexTransform );
  
  mWorldToIndexTransform =
    mDataToIndexTransform * mDataToWorldTransform->Inverse();

#if 0
  cerr << "mDataToIndex " << mDataToIndexTransform << endl;
  cerr << "mDataToWorldTransform inv " << mDataToWorldTransform->Inverse() << endl;
  cerr << "mWorldToIndexTransform" << mWorldToIndexTransform << endl;
#endif

  DataChanged();
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

  float voxelSize[3];
  voxelSize[0] = iVolume.GetVoxelXSize();
  voxelSize[1] = iVolume.GetVoxelYSize();
  voxelSize[2] = iVolume.GetVoxelZSize();

  this->DoBegin();

  Volume3<bool>* bVisited =
    new Volume3<bool>( iVolume.mMRI->width, 
		       iVolume.mMRI->height, 
		       iVolume.mMRI->depth, false );

  // Save the initial value.
  float seedValue = iVolume.GetMRINearestValueAtRAS( iRASSeed );

  // Push the seed onto the list. 
  vector<Point3<float> > RASPoints;
  RASPoints.push_back( iRASSeed );
  while( RASPoints.size() > 0 &&
	 !this->DoStopRequested() ) {
    
    Point3<float> ras = RASPoints.back();
    RASPoints.pop_back();


    if( !iVolume.IsRASInMRIBounds( ras.xyz() ) ) {
      continue;
    }

    Point3<int> index;
    iVolume.RASToMRIIndex( ras.xyz(), index.xyz() );
    if( bVisited->Get_Unsafe( index.x(), index.y(), index.z() ) ) {
      continue;
    }
    bVisited->Set_Unsafe( index.x(), index.y(), index.z(), true );

    // Check if this is an edge or an ROI. If so, and our params say
    // not go to here, continue.
    if( iParams.mbStopAtEdges ) {
      if( iVolume.IsRASEdge( ras.xyz() ) ) {
	continue;
      }
    }
    if( iParams.mbStopAtROIs ) {
      if( iVolume.IsOtherRASSelected( ras.xyz(), iVolume.GetSelectedROI() ) ) {
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
      float value = iVolume.GetMRINearestValueAtRAS( ras.xyz() );
      if( fabs( value - seedValue ) > iParams.mFuzziness ) {
	continue;
      }
    }

    // Call the user compare function to give them a chance to bail.
    if( !this->CompareVoxel( ras.xyz() ) ) {
      continue;
    }
    
    // Call the user function.
    this->DoVoxel( ras.xyz() );

    // Add adjacent voxels.
    float beginX = ras.x() - voxelSize[0];
    float endX   = ras.x() + voxelSize[0];
    float beginY = ras.y() - voxelSize[1];
    float endY   = ras.y() + voxelSize[1];
    float beginZ = ras.z() - voxelSize[2];
    float endZ   = ras.z() + voxelSize[2];
    if( !iParams.mb3D && iParams.mbWorkPlaneX ) {
      beginX = endX = ras.x();
    }
    if( !iParams.mb3D && iParams.mbWorkPlaneY ) {
      beginY = endY = ras.y();
    }
    if( !iParams.mb3D && iParams.mbWorkPlaneZ ) {
      beginZ = endZ = ras.z();
    }
    Point3<float> newRAS;
    if( iParams.mbDiagonal ) {
      for( float nZ = beginZ; nZ <= endZ; nZ += voxelSize[2] ) {
	for( float nY = beginY; nY <= endY; nY += voxelSize[1] ) {
	  for( float nX = beginX; nX <= endX; nX += voxelSize[0] ) {
	    newRAS.Set( nX, nY, nZ );
	    RASPoints.push_back( newRAS );
	  }
	}
      }
    } else {
      newRAS.Set( beginX, ras.y(), ras.z() );
      RASPoints.push_back( newRAS );
      newRAS.Set( endX, ras.y(), ras.z() );
      RASPoints.push_back( newRAS );
      newRAS.Set( ras.x(), beginY, ras.z() );
      RASPoints.push_back( newRAS );
      newRAS.Set( ras.x(), endY, ras.z() );
      RASPoints.push_back( newRAS );
      newRAS.Set( ras.x(), ras.y(), beginZ );
      RASPoints.push_back( newRAS );
      newRAS.Set( ras.x(), ras.y(), endZ );
      RASPoints.push_back( newRAS );
    }

  }

  mVolume = NULL;
  mParams = NULL;

  this->DoEnd();
}
