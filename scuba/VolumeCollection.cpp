#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include <vector>
extern "C" {
#include "error.h"
#include "ctrpoints.h"
}
#include "VolumeCollection.h"
#include "DataManager.h"
#include "Point3.h"
#include "Utilities.h"
#include "PathManager.h"
#include "VectorOps.h"

using namespace std;


VolumeCollection::VolumeCollection () :
  DataCollection() {
  mMRI = NULL;
  mMagnitudeMRI = NULL;
  mSelectedVoxels = NULL;
  mVoxelSize[0] = mVoxelSize[1] = mVoxelSize[2] = 0;
  mbUseDataToIndexTransform = true;

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
  commandMgr.AddCommand( *this, "MakeVolumeUsingTemplate", 2, 
			 "collectionID templateCollectionID", 
			 "Makes a volume using an existing volume "
			 "as a template." );
  commandMgr.AddCommand( *this, "SaveVolume", 1, 
			 "collectionID", "Save volume with its file name." );
  commandMgr.AddCommand( *this, "SaveVolumeWithFileName", 2, 
			 "collectionID fileName", "Save volume with "
			 "a given file name." );
  commandMgr.AddCommand( *this, "SetUseVolumeDataToIndexTransform", 2, 
			 "collectionID use", "Use or don't use the volume's "
			 "Data to Index transform (usually RAS transform) "
			 "in displaying data." );
  commandMgr.AddCommand( *this, "GetUseVolumeDataToIndexTransform", 1, 
			 "collectionID", "Returns whether or not a volume "
			 "is using its Data to Index transform "
			 "(usually RAS transform) in displaying data." );
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

DataLocation&
VolumeCollection::MakeLocationFromRAS ( float const iRAS[3] ) {
  
  VolumeLocation* loc = new VolumeLocation( *this, iRAS );
  return *loc;
}

void
VolumeCollection::SetFileName ( string& ifnMRI ) {

  mfnMRI = ifnMRI;
}

void
VolumeCollection::MakeUsingTemplate ( int iCollectionID ) {

  VolumeCollection* vol = NULL;
  try { 
    DataCollection* col = &DataCollection::FindByID( iCollectionID );
    //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
    vol = (VolumeCollection*)col;
  }
  catch (...) {
    throw runtime_error( "Couldn't find template." );
  }

  // Get the mri from the template volume.
  MRI* mri = vol->GetMRI();
  if( NULL == mri ) {
    throw runtime_error( "Couldn't get MRI from template" );
  }

  // Allocate the mri with the size from the template.
  MRI* newMri = MRIallocSequence( mri->width, mri->height, mri->depth,
				  mri->type, mri->nframes );
  if( NULL == newMri ) {
    throw runtime_error( "Couldn't allocate new mri." );
  }
  
  // Copy the header from the template into the new mri.
  MRIcopyHeader( mri, newMri );

  // Save the MRI.
  mMRI = newMri;
  
  // Initialize from it.
  InitializeFromMRI();
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

    InitializeFromMRI();
  }

  return mMRI; 
}

void
VolumeCollection::Save () {

  Save( mfnMRI );
}

void
VolumeCollection::Save ( string ifn ) {

  char* fn = strdup( ifn.c_str() ); 
  int rMRI = MRIwrite( mMRI, fn );
  free( fn );

  if( ERROR_NONE != rMRI ) {
    stringstream ssError;
    ssError << "Couldn't write file " << ifn;
    throw runtime_error( ssError.str() );
  }
}

void
VolumeCollection::InitializeFromMRI () {

  if( NULL == mMRI ) {
    throw runtime_error( "InitializeFromMRI called without an MRI" );
  }

  // Get our surfaceRAS -> index transform.
  MATRIX* voxelFromSurfaceRAS = extract_r_to_i( mMRI );
  if( NULL == voxelFromSurfaceRAS ) {
    throw runtime_error( "Couldn't get voxelFromSurfaceRAS matrix" );
  }

  // Copy it to a Matrix44 and release the MATRIX. Then set our
  // mDataToIndexTransform transform from this matrix. Then calculate
  // the WorldToIndex transform.
  Matrix44 m;
  m.SetMatrix( voxelFromSurfaceRAS );

  MatrixFree( &voxelFromSurfaceRAS );
  
  mDataToIndexTransform.SetMainTransform( m );
  
  CalcWorldToIndexTransform();

  // Update (initialize) our MRI value range.
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
  
  // Init the selection volume.
  InitSelectionVolume();
  
  // Save our voxel sizes.
  mVoxelSize[0] = mMRI->xsize;
  mVoxelSize[1] = mMRI->ysize;
  mVoxelSize[2] = mMRI->zsize;
  
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
VolumeCollection::GetMRINearestValueAtLocation ( VolumeLocation& iLoc ) {

  Real value = 0;
  if( NULL != mMRI ) {

    if( iLoc.mIdx[0] >= 0 && iLoc.mIdx[0] < mMRI->width &&
	iLoc.mIdx[1] >= 0 && iLoc.mIdx[1] < mMRI->height &&
	iLoc.mIdx[2] >= 0 && iLoc.mIdx[2] < mMRI->depth ) {
      
      switch( mMRI->type ) {
      case MRI_UCHAR:
	value = (float)MRIvox(mMRI, iLoc.mIdx[0], iLoc.mIdx[1], iLoc.mIdx[2] );
	break ;
      case MRI_SHORT:
      value = (float)MRISvox(mMRI, iLoc.mIdx[0], iLoc.mIdx[1], iLoc.mIdx[2] );
      break ;
      case MRI_INT:
	value = (float)MRIIvox(mMRI, iLoc.mIdx[0], iLoc.mIdx[1], iLoc.mIdx[2] );
	break ;
      case MRI_FLOAT:
	value = MRIFvox(mMRI, iLoc.mIdx[0], iLoc.mIdx[1], iLoc.mIdx[2] );
	break ;
      default:
	value = 0;
      }
    }
  }
  
  return (float)value;
  
}

bool 
VolumeCollection::IsLocationInBounds ( VolumeLocation& iLoc ) {

  return IsMRIIndexInMRIBounds( iLoc.mIdx );
}

bool 
VolumeCollection::IsLocationSelected ( VolumeLocation& iLoc, int oColor[3] ) {

  // Check the selection volume cache first.
  if( !(mSelectedVoxels->Get_Unsafe( iLoc.mIdx[0], iLoc.mIdx[1], iLoc.mIdx[2] )) )
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
      if( volumeROI->IsVoxelSelected( iLoc.mIdx ) ) {
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


void
VolumeCollection::GetRASRange ( float oRASRange[6] ) {

  if( NULL != mMRI ) {
    oRASRange[0] = mMRI->xstart;
    oRASRange[1] = mMRI->xend;
    oRASRange[2] = mMRI->ystart;
    oRASRange[3] = mMRI->yend;
    oRASRange[4] = mMRI->zstart;
    oRASRange[5] = mMRI->zend;
  }
}

void
VolumeCollection::GetMRIIndexRange ( int oMRIIndexRange[3] ) {

  if( NULL != mMRI ) {
    oMRIIndexRange[0] = mMRI->width;
    oMRIIndexRange[1] = mMRI->height;
    oMRIIndexRange[2] = mMRI->depth;
  }
}

float
VolumeCollection::GetMRIMagnitudeMaxValue () { 
  if( NULL == mMagnitudeMRI ) {
    MakeMagnitudeVolume();
  }
  return mMRIMagMaxValue; 
}


void
VolumeCollection::RASToMRIIndex ( float const iRAS[3], int oIndex[3] ) {
  

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
VolumeCollection::RASToMRIIndex ( float const iRAS[3], float oIndex[3] ) {

  mWorldToIndexTransform.MultiplyVector3( iRAS, oIndex );
}

void
VolumeCollection::MRIIndexToRAS ( int const iIndex[3], float oRAS[3] ) {
  
  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::MRIIndexToRAS ( float const iIndex[3], float oRAS[3] ) {
  
  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::RASToDataRAS ( float const iRAS[3], float oDataRAS[3] ) {
  
  mDataToWorldTransform->InvMultiplyVector3( iRAS, oDataRAS );
}

bool 
VolumeCollection::IsRASInMRIBounds ( float const iRAS[3] ) {

  int index[3];
  RASToMRIIndex( iRAS, index );

  return IsMRIIndexInMRIBounds( index );
}

bool 
VolumeCollection::IsMRIIndexInMRIBounds ( int const iIndex[3] ) {

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
    default:
      break ;
    }
  }

  if( iValue < mMRIMinValue ) {
    mMRIMinValue = iValue;
  }
  if( iValue > mMRIMaxValue ) {
    mMRIMaxValue = iValue;
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
  
  // MakeVolumeUsingTemplate <collectionID> <templateID>
  if( 0 == strcmp( isCommand, "MakeVolumeUsingTemplate" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
      int templateID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad template ID";
	return error;
      }
    
      MakeUsingTemplate( templateID );
    }
  }
  
  // SaveVolume <collectionID>
  if( 0 == strcmp( isCommand, "SaveVolume" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
      Save();
    }
  }
  
  // SaveVolumeWithFileName <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "SaveVolumeWithFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
     
      string fn( iasArgv[2] );
      Save( fn );
    }
  }

  // SetUseVolumeDataToIndexTransform <collectionID> <use>
  if( 0 == strcmp( isCommand, "SetUseVolumeDataToIndexTransform" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      try {
	bool bUse =
	  TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
	SetUseWorldToIndexTransform( bUse );
      }
      catch( runtime_error e ) {
	sResult = "bad use \"" + string(iasArgv[2]) + "\"," + e.what();
	return error;	
      }
    }
  }
  

  // GetUseVolumeDataToIndexTransform <layerID>
  if( 0 == strcmp( isCommand, "GetUseVolumeDataToIndexTransform" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {

      sReturnValues =
	TclCommandManager::ConvertBooleanToReturnValue( mbUseDataToIndexTransform );
      sReturnFormat = "i";
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
VolumeCollection::FindRASPointsInSquare ( float iPointA[3], float iPointB[3],
					  float iPointC[3], float iPointD[3],
					  float iMaxDistance,
					  list<Point3<float> >& oPoints ) {

  Point3<float> squareRAS[4];
  squareRAS[0].Set( iPointA );
  squareRAS[1].Set( iPointB );
  squareRAS[2].Set( iPointC );
  squareRAS[3].Set( iPointD );

  // Find a plane normal from the four points.
  Point3<float> v[2];
  v[0] = squareRAS[1] - squareRAS[0];
  v[1] = squareRAS[2] - squareRAS[0];
  Point3<float> n = VectorOps::Cross( v[0], v[1] );

  // Find the MRI indices. Find the square-to-volume bound.
  Point3<int> squareIdx[4];
  RASToMRIIndex( squareRAS[0].xyz(), squareIdx[0].xyz() );
  RASToMRIIndex( squareRAS[1].xyz(), squareIdx[1].xyz() );
  RASToMRIIndex( squareRAS[2].xyz(), squareIdx[2].xyz() );
  RASToMRIIndex( squareRAS[3].xyz(), squareIdx[3].xyz() );


  Point3<int> volumeBoundIdx[2];
  volumeBoundIdx[0].Set
    ( MIN(MIN(MIN(squareIdx[0][0],squareIdx[1][0]),squareIdx[2][0]),
	  squareIdx[3][0]),
      MIN(MIN(MIN(squareIdx[0][1],squareIdx[1][1]),squareIdx[2][1]),
	  squareIdx[3][1]),
      MIN(MIN(MIN(squareIdx[0][2],squareIdx[1][2]),squareIdx[2][2]),
	  squareIdx[3][2]));
  volumeBoundIdx[1].Set
    ( MAX(MAX(MAX(squareIdx[0][0],squareIdx[1][0]),squareIdx[2][0]),
	  squareIdx[3][0]),
      MAX(MAX(MAX(squareIdx[0][1],squareIdx[1][1]),squareIdx[2][1]),
	  squareIdx[3][1]),
      MAX(MAX(MAX(squareIdx[0][2],squareIdx[1][2]),squareIdx[2][2]),
	  squareIdx[3][2]));

  // If there's only one voxel in the bounds, this is a special case
  // where no voxels's edges will cross the square. So we'll just add
  // this one and return.
  if( volumeBoundIdx[0].x() == volumeBoundIdx[1].x() &&
      volumeBoundIdx[0].y() == volumeBoundIdx[1].y() &&
      volumeBoundIdx[0].z() == volumeBoundIdx[1].z() ) {
    Point3<float> centerRAS;
    MRIIndexToRAS( volumeBoundIdx[0].xyz(), centerRAS.xyz() );
    oPoints.push_back( centerRAS );
    return;
  }

  // For each voxel in the cuboid...
  for( int nZ = volumeBoundIdx[0].z(); nZ <= volumeBoundIdx[1].z(); nZ++ ) {
    for( int nY = volumeBoundIdx[0].y(); nY <= volumeBoundIdx[1].y(); nY++ ) {
      for( int nX = volumeBoundIdx[0].x(); nX <= volumeBoundIdx[1].x(); nX++ ){


	// Create RAS versions of our corners.
	Point3<int> voxelIdx[8];
	voxelIdx[0].Set( nX  , nY  , nZ   );
	voxelIdx[1].Set( nX+1, nY  , nZ   );
	voxelIdx[2].Set( nX  , nY+1, nZ   );
	voxelIdx[3].Set( nX+1, nY+1, nZ   );
	voxelIdx[4].Set( nX  , nY  , nZ+1 );
	voxelIdx[5].Set( nX+1, nY  , nZ+1 );
	voxelIdx[6].Set( nX  , nY+1, nZ+1 );
	voxelIdx[7].Set( nX+1, nY+1, nZ+1 );
	Point3<float> voxelRAS[8];
	for( int nCorner = 0; nCorner < 8; nCorner++ ) {
	  MRIIndexToRAS( voxelIdx[nCorner].xyz(), voxelRAS[nCorner].xyz() );
	}

	Point3<float> segmentRAS[12][2];
	int anSegments[12][2] = { {0, 1}, {4, 5}, {6, 7}, {2, 3},
				  {0, 2}, {1, 3}, {4, 6}, {5, 7},
				  {0, 4}, {1, 5}, {2, 6}, {3, 7} };
	for( int nSegment = 0; nSegment < 12; nSegment++ ) {
	  segmentRAS[nSegment][0].Set( voxelRAS[anSegments[nSegment][0]] );
	  segmentRAS[nSegment][1].Set( voxelRAS[anSegments[nSegment][1]] );
	}

	// Intersect these segments with the plane. If any of them hit...
	for( int nSegment = 0; nSegment < 12; nSegment++ ) {

	  Point3<float> intersectionRAS;
	  VectorOps::IntersectionResult rInt =
	    VectorOps::SegmentIntersectsPlane
	    ( segmentRAS[nSegment][0], segmentRAS[nSegment][1],
	      squareRAS[0], n, intersectionRAS );
	  
	  if( VectorOps::intersect == rInt ) {

	    // Calculate the anglesum of the intersection point with the
	    // four plane corner points. If it is 2pi, this point is
	    // inside the poly.
	    double angleSum = 0;
	    for( int nVector = 0; nVector < 4; nVector++ ) {
	      Point3<float> v1 = squareRAS[nVector]       - intersectionRAS;
	      Point3<float> v2 = squareRAS[(nVector+1)%4] - intersectionRAS;

	      float v1Length = VectorOps::Length(v1);
	      float v2Length = VectorOps::Length(v2);

	      if( fabs(v1Length * v2Length) <= (float)0.0001 ) {
		angleSum = 2*M_PI;
		break;
	      }
	      
	      double rads = VectorOps::RadsBetweenVectors( v1, v2 );
	      angleSum += rads;
	    }
	    
	    if( fabs(angleSum - 2.0*M_PI) <= (float)0.0001 ) {
	      oPoints.push_back( voxelRAS[0] );
	      break;
	    }

	  } // if intersect
	} // for segment
      } // for x
    } // for y
  } // for z
}

void
VolumeCollection::FindRASPointsInCircle ( float iPointA[3], float iPointB[3],
					  float iPointC[3], float iPointD[3],
					  float iMaxDistance,
					  float iCenter[3], float iRadius,
					  list<Point3<float> >& oPoints ) {

  list<Point3<float> > squarePoints;
  FindRASPointsInSquare( iPointA, iPointB, iPointC, iPointD, 
			 iMaxDistance, squarePoints );

  Point3<float> center( iCenter );
  list<Point3<float> >::iterator tPoints;
  for( tPoints = squarePoints.begin(); tPoints != squarePoints.end();
       ++tPoints ) {

    Point3<float> point = *tPoints;
    if( VectorOps::Distance( point, center ) <= iRadius ) {
      oPoints.push_back( point );
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
VolumeCollection::SetDataToWorldTransform ( int iTransformID ) {

  DataCollection::SetDataToWorldTransform( iTransformID );
  CalcWorldToIndexTransform();
}


Matrix44&
VolumeCollection::GetWorldToIndexTransform () {

  return mWorldToIndexTransform.GetMainMatrix(); 
}


void
VolumeCollection::SetUseWorldToIndexTransform ( bool ibUse ) {

  mbUseDataToIndexTransform = ibUse;
  CalcWorldToIndexTransform();
}

void
VolumeCollection::CalcWorldToIndexTransform () {

  if( mbUseDataToIndexTransform ) {

    // Just mult our transforms together.
    Transform44 worldToData = mDataToWorldTransform->Inverse();
    Transform44 tmp = mDataToIndexTransform * worldToData;
    mWorldToIndexTransform = tmp;

  } else {

    Transform44 center;
    center.SetMainTransform( 1, 0, 0, mMRI->width/2,
			     0, 1, 0, mMRI->height/2,
			     0, 0, 1, mMRI->depth/2,
			     0, 0, 0, 1 );
    
    Transform44 worldToData = mDataToWorldTransform->Inverse();
    Transform44 tmp = center * worldToData;
    mWorldToIndexTransform = tmp;
      

  }

#if 0
  cerr << "mDataToIndex " << mDataToIndexTransform << endl;
  cerr << "mDataToWorldTransform inv " << mDataToWorldTransform->Inverse() << endl;
  cerr << "mWorldToIndexTransform" << mWorldToIndexTransform << endl;
#endif

  DataChanged();
}

void
VolumeCollection::ImportControlPoints ( string ifnControlPoints,
					list<Point3<float> >& oControlPoints ){

  int cControlPoints;
  int bUseRealRAS;
  char* fnControlPoints = strdup( ifnControlPoints.c_str() );
  MPoint* aControlPoints = 
    MRIreadControlPoints( fnControlPoints, &cControlPoints, &bUseRealRAS );
  free( fnControlPoints );
  if( NULL == aControlPoints ) {
    throw runtime_error( "Couldn't read control points file." );
  }

  for( int nControlPoint = 0; nControlPoint < cControlPoints; nControlPoint++){

    Point3<float> newControlPoint;
    if( !bUseRealRAS ) {

      // We're getting surface RAS points. Need to convert to normal
      // RAS points.
      Real surfaceRAS[3];
      Real ras[3];
      surfaceRAS[0] = aControlPoints[nControlPoint].x;
      surfaceRAS[1] = aControlPoints[nControlPoint].y;
      surfaceRAS[2] = aControlPoints[nControlPoint].z;
      MRIsurfaceRASToRAS( mMRI, surfaceRAS[0], surfaceRAS[1], surfaceRAS[2],
			  &ras[0], &ras[1], &ras[2] );
      newControlPoint[0] = ras[0];
      newControlPoint[1] = ras[1];
      newControlPoint[2] = ras[2];

    } else {

      // Already in real ras.
      newControlPoint[0] = aControlPoints[nControlPoint].x;
      newControlPoint[1] = aControlPoints[nControlPoint].y;
      newControlPoint[2] = aControlPoints[nControlPoint].z;
    }

    oControlPoints.push_back( newControlPoint );
  }

  free( aControlPoints );
}

void
VolumeCollection::ExportControlPoints ( string ifnControlPoints,
					list<Point3<float> >& iControlPoints) {


  int cControlPoints = iControlPoints.size();
  MPoint* aControlPoints = (MPoint*) calloc( sizeof(MPoint), cControlPoints );
  if( NULL == aControlPoints ) {
    throw runtime_error( "Couldn't allocate control point storage." );
  }

  int nControlPoint = 0;
  list<Point3<float> >::iterator tControlPoint;
  for( tControlPoint = iControlPoints.begin(); 
       tControlPoint != iControlPoints.end();
       ++tControlPoint ) {

    aControlPoints[nControlPoint].x = (*tControlPoint)[0];
    aControlPoints[nControlPoint].y = (*tControlPoint)[1];
    aControlPoints[nControlPoint].z = (*tControlPoint)[2];
    nControlPoint++;
  }

  char* fnControlPoints = strdup( ifnControlPoints.c_str() );
  int rWrite = MRIwriteControlPoints( aControlPoints, cControlPoints, 
				      true, fnControlPoints );
  free( fnControlPoints );

  if( NO_ERROR != rWrite ) {
    throw runtime_error( "Couldn't write control point file." );
  }
}

void
VolumeCollection::MakeHistogram ( list<Point3<float> >& iRASPoints, 
				  int icBins,
				  float& oMinBinValue, float& oBinIncrement,
				  map<int,int>& oBinCounts ) {

  // If no points, return.
  if( iRASPoints.size() == 0 ) {
    return;
  }

  // Find values for all the RAS points. Save the highs and lows.
  float low = 999999;
  float high = -999999;
  list<Point3<float> >::iterator tPoint;
  list<float> values;
  for( tPoint = iRASPoints.begin(); tPoint != iRASPoints.end(); ++tPoint ) {
    Point3<float> point = *tPoint;
    float value = GetMRINearestValueAtRAS( point.xyz() );
    if( value < low ) { low = value; }
    if( value > high ) { high = value; }
    values.push_back( value );
  }

  float binIncrement = (high - low) / (float)icBins;

  // Now go through the values and calculate a bin for them. Increment
  // the count in the proper bin.
  list<float>::iterator tValue;
  for( tValue = values.begin(); tValue != values.end(); ++tValue ) {
    float value = *tValue;
    int nBin = (int)floor( (value - low) / (float)binIncrement );
    if( nBin == icBins && value == high ) { // if value == high, bin will be
      nBin = icBins-1;		            // icBins, should be icBins-1
    }
    oBinCounts[nBin]++;
  }

  oMinBinValue = low;
  oBinIncrement = binIncrement;
}




VolumeCollectionFlooder::VolumeCollectionFlooder () {
  mVolume = NULL;
  mParams = NULL;
}

VolumeCollectionFlooder::~VolumeCollectionFlooder () {
}

VolumeCollectionFlooder::Params::Params () {
  mbStopAtPaths = true;
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

  // Get the source volume.
  VolumeCollection* sourceVol = NULL;
  try { 
    DataCollection* col = 
      &DataCollection::FindByID( iParams.mSourceCollection );
    //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
    sourceVol = (VolumeCollection*)col;
  }
  catch (...) {
    throw runtime_error( "Couldn't find source volume." );
  }
  bool bDifferentSource = sourceVol->GetID() != iVolume.GetID();

  // Save the initial value.
  float seedValue = sourceVol->GetMRINearestValueAtRAS( iRASSeed );

  // Push the seed onto the list. 
  vector<CheckPair> checkPairs;
  Point3<float> seed( iRASSeed );
  CheckPair seedPair( seed, seed );
  checkPairs.push_back( seedPair );
  while( checkPairs.size() > 0 &&
	 !this->DoStopRequested() ) {
    
    CheckPair checkPair = checkPairs.back();
    checkPairs.pop_back();

    Point3<float> ras = checkPair.mCheckRAS;
    Point3<float> sourceRAS = checkPair.mSourceRAS;

    // Check the bound of this volume and the source one.
    if( !iVolume.IsRASInMRIBounds( ras.xyz() ) ) { 
     continue;
    }
    if( bDifferentSource && !sourceVol->IsRASInMRIBounds( ras.xyz() ) ) { 
     continue;
    }

    Point3<int> index;
    iVolume.RASToMRIIndex( ras.xyz(), index.xyz() );
    if( bVisited->Get_Unsafe( index.x(), index.y(), index.z() ) ) {
      continue;
    }
    bVisited->Set_Unsafe( index.x(), index.y(), index.z(), true );

    // Check if this is an path or an ROI. If so, and our params say
    // not go to here, continue.
    if( iParams.mbStopAtPaths ) {
      
      // Create a line from our current point to this check
      // point. Then see if the segment intersects the path. If so,
      // don't proceed.
      
      bool bCross = false;
      Point3<float> x;
      list<Path<float>*>::iterator tPath;
      PathManager& pathMgr = PathManager::GetManager();
      list<Path<float>*>& pathList = pathMgr.GetPathList();
      for( tPath = pathList.begin(); tPath != pathList.end() && !bCross; ++tPath ) {
	Path<float>* path = *tPath;
	if( path->GetNumVertices() > 0 ) {

	  int cVertices = path->GetNumVertices();
	  for( int nCurVertex = 1; nCurVertex < cVertices && !bCross; nCurVertex++ ) {
	    
	    int nBackVertex = nCurVertex - 1;
	    
	    Point3<float>& curVertex  = path->GetVertexAtIndex( nCurVertex );
	    Point3<float>& backVertex = path->GetVertexAtIndex( nBackVertex );

	    Point3<float> segVector = curVertex - backVertex;
	    Point3<float> viewNormal( iParams.mViewNormal );

	    Point3<float> normal = VectorOps::Cross( segVector, viewNormal );

	    VectorOps::IntersectionResult rInt =
	      VectorOps::SegmentIntersectsPlane( sourceRAS, ras,
						 curVertex, normal,
						 x );
	    if( VectorOps::intersect == rInt ) {

	      // Make sure x is in the cuboid formed by curVertex and
	      // backVertex.
	      if( !(x[0] < MIN(curVertex[0],backVertex[0]) ||
		    x[0] > MAX(curVertex[0],backVertex[0]) ||
		    x[1] < MIN(curVertex[1],backVertex[1]) ||
		    x[1] > MAX(curVertex[1],backVertex[1]) ||
		    x[2] < MIN(curVertex[2],backVertex[2]) ||
		    x[2] > MAX(curVertex[2],backVertex[2])) ) {
		bCross = true;
	      }
	    }
	  }
	}
      }

      // If we crossed a path, continue.
      if( bCross )
	continue;
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


    // Check only zero.
    if( iParams.mbOnlyZero ) {
      float value = iVolume.GetMRINearestValueAtRAS( ras.xyz() );
      if( value != 0 ) {
	continue;
      }
    }

    // Check fuzziness.
    if( iParams.mFuzziness > 0 ) {
      float value = sourceVol->GetMRINearestValueAtRAS( ras.xyz() );
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
    CheckPair newPair;
    newPair.mSourceRAS = ras;
    if( iParams.mbDiagonal ) {
      for( float nZ = beginZ; nZ <= endZ; nZ += voxelSize[2] ) {
	for( float nY = beginY; nY <= endY; nY += voxelSize[1] ) {
	  for( float nX = beginX; nX <= endX; nX += voxelSize[0] ) {
	    newRAS.Set( nX, nY, nZ );
	    newPair.mCheckRAS = newRAS;
	    checkPairs.push_back( newPair );
	  }
	}
      }
    } else {
      newRAS.Set( beginX, ras.y(), ras.z() );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

      newRAS.Set( endX, ras.y(), ras.z() );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

      newRAS.Set( ras.x(), beginY, ras.z() );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

      newRAS.Set( ras.x(), endY, ras.z() );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

      newRAS.Set( ras.x(), ras.y(), beginZ );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

      newRAS.Set( ras.x(), ras.y(), endZ );
      newPair.mCheckRAS = newRAS;
      checkPairs.push_back( newPair );

    }

  }

  mVolume = NULL;
  mParams = NULL;

  this->DoEnd();
}

VolumeLocation::VolumeLocation ( VolumeCollection& iVolume,
				 float const iRAS[3] )
  : DataLocation( iRAS ), mVolume( iVolume ) {
  mVolume.RASToMRIIndex( iRAS, mIdx );
}
