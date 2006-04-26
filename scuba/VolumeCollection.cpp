#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
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
  mbAutosave = false;
  mbAutosaveDirty = false;
  mfnAutosave = "";
  mVoxelSize[0] = mVoxelSize[1] = mVoxelSize[2] = 0;
  mbUseDataToIndexTransform = true;
  mbBoundsCacheDirty = true;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollectionFileName", 2, 
			 "collectionID fileName", 
			 "Sets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "LoadVolumeFromFileName", 1, "collectionID", 
			 "Loads the volume from the file name.");
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
  commandMgr.AddCommand( *this, "SetVolumeAutosaveOn", 2, "collectionID on",
			 "Set whether or not autosave is on for this "
			 "volume." );
  commandMgr.AddCommand( *this, "GetVolumeAutosaveOn", 1, "collectionID",
			 "Returns whether or not autosave is on for this "
			 "volume." );
  commandMgr.AddCommand( *this, "GetRASCoordsFromVolumeSurfaceRAS", 4, 
			 "collectionID x y z", "Returns a list of RAS coords "
			 "converted from the input surface RAS coords. This "
			 "is for converting RAS points acquired from a "
			 "surface that is associated with a volume and "
			 "didn't generate coordinates with CRAS info." );
  commandMgr.AddCommand( *this, "GetVolumeSurfaceRASCoordsFromRAS", 4, 
			 "collectionID x y z", "Returns a list of surface "
			 "RAS coords converted from the input RAS coords. This"
			 "is for converting RAS points acquired from a "
			 "surface that is associated with a volume and "
			 "didn't generate coordinates with CRAS info." );
  commandMgr.AddCommand( *this, "GetVolumeAverageValueInROI", 2,
			 "collectionID, roiID", "Returns the average value "
			 "of the voxels in an ROI in a volume." );
  commandMgr.AddCommand( *this, "GetVolumeStandardDeviationInROI", 2,
			 "collectionID, roiID", "Returns the standard "
			 "deviation of the voxels in an ROI in a volume." );
}

VolumeCollection::~VolumeCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();
  try { 
    mriLoader.ReleaseData( &mMRI );
    SendBroadcast( "DataDeleted", NULL );
  } 
  catch( exception& e ) {
    cerr << "Error releasing data: " << e.what() << endl;
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

DataLocation&
VolumeCollection::MakeLocationFromIndex ( int const iIndex[3] ) {
  
  VolumeLocation* loc = new VolumeLocation( *this, iIndex );
  return *loc;
}

void
VolumeCollection::SetFileName ( string& ifnMRI ) {

  mfnMRI = ifnMRI;
  mfnAutosave = MakeAutosaveFileName( mfnMRI );
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

  // Set a temporary filename.
  string fn = "New_Volume.mgh";
  SetFileName( fn );
}

void
VolumeCollection::LoadVolume () {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();

  // If we already have data...
  if( NULL != mMRI ) {

    // Try to load this and see what we get. If it's the same as what
    // we already have, we're fine. If not, keep this one and release
    // the one we have.
    MRI* newMRI = NULL; try { newMRI = mriLoader.GetData( mfnMRI ); }
    catch( exception& e ) { throw logic_error( "Couldn't load MRI" );
    }

    if( newMRI == mMRI ) {
      return;
    }

    // Release old data.
    try { 
      mriLoader.ReleaseData( &mMRI );
    } 
    catch(...) {
      cerr << "Couldn't release data"  << endl;
    }

    // Save new data.
    mMRI = newMRI;
    InitializeFromMRI();
    DataChanged();

  } else {

    // Don't already have it, so get it.
    try { 
      mMRI = mriLoader.GetData( mfnMRI );
    }
    catch( exception& e ) {
      throw logic_error( "Couldn't load MRI" );
    }

    if( msLabel == "" ) {
      SetLabel( mfnMRI );
    }

    InitializeFromMRI();
  }
}

MRI*
VolumeCollection::GetMRI() { 

  // If we don't already have one, load it.
  if( NULL == mMRI ) {
    LoadVolume();
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
VolumeCollection::GetDataRASBounds ( float oRASBounds[6] ) {

  if( NULL != mMRI ) {
    
    if( mbBoundsCacheDirty ) {
      
      mRASBounds[0] = mRASBounds[2] = mRASBounds[4] = 999999;
      mRASBounds[1] = mRASBounds[3] = mRASBounds[5] = -999999;

      int index[3];
      float RAS[3];
      for( index[2] = 0; index[2] < mMRI->depth; index[2]++ ) {
	for( index[1] = 0; index[1] < mMRI->height; index[1]++ ) {
	  for( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {
	    
	    MRIIndexToRAS( index, RAS );

	    if( RAS[0] < mRASBounds[0] ) mRASBounds[0] = RAS[0];
	    if( RAS[0] > mRASBounds[1] ) mRASBounds[1] = RAS[0];
	    if( RAS[1] < mRASBounds[2] ) mRASBounds[2] = RAS[1];
	    if( RAS[1] > mRASBounds[3] ) mRASBounds[3] = RAS[1];
	    if( RAS[2] < mRASBounds[4] ) mRASBounds[4] = RAS[2];
	    if( RAS[2] > mRASBounds[5] ) mRASBounds[5] = RAS[2];
	  }
	}
      }

      mbBoundsCacheDirty = false;
    }
    
    oRASBounds[0] = mRASBounds[0];
    oRASBounds[1] = mRASBounds[1];
    oRASBounds[2] = mRASBounds[2];
    oRASBounds[3] = mRASBounds[3];
    oRASBounds[4] = mRASBounds[4];
    oRASBounds[5] = mRASBounds[5];
    
  } else {
    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] = 
      oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
  }
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


bool 
VolumeCollection::IsInBounds ( VolumeLocation& iLoc ) {

  if( NULL != mMRI ) {
      return ( iLoc.mIdxi[0] >= 0 && iLoc.mIdxi[0] < mMRI->width &&
	       iLoc.mIdxi[1] >= 0 && iLoc.mIdxi[1] < mMRI->height &&
	       iLoc.mIdxi[2] >= 0 && iLoc.mIdxi[2] < mMRI->depth );
  } else {
    return false;
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
  oIndex[0] += 0.5;
  oIndex[1] += 0.5;
  oIndex[2] += 0.5;
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

void
VolumeCollection::DataRASToRAS ( float const iDataRAS[3], float oRAS[3] ) {
  
  mDataToWorldTransform->MultiplyVector3( iDataRAS, oRAS );
}

void
VolumeCollection::SurfaceRASToRAS ( float const iSurfaceRAS[3],
				    float oRAS[3] ) {

  if( NULL != mMRI ) {
    Real surfaceRAS[3], RAS[3];
    surfaceRAS[0] = iSurfaceRAS[0];
    surfaceRAS[1] = iSurfaceRAS[1];
    surfaceRAS[2] = iSurfaceRAS[2];
    MRIsurfaceRASToRAS( mMRI, surfaceRAS[0], surfaceRAS[1], surfaceRAS[2],
		       &RAS[0], &RAS[1], &RAS[2] );
    oRAS[0] = RAS[0];
    oRAS[1] = RAS[1];
    oRAS[2] = RAS[2];
  } else {
    throw runtime_error( "Cannot transform from surface RAS because no "
			 "MRI is loaded." );
  }
}

void
VolumeCollection::RASToSurfaceRAS ( float const iRAS[3],
				    float oSurfaceRAS[3] ) {

  if( NULL != mMRI ) {
    Real surfaceRAS[3], RAS[3];
    RAS[0] = iRAS[0];
    RAS[1] = iRAS[1];
    RAS[2] = iRAS[2];
    MRIRASToSurfaceRAS( mMRI, RAS[0], RAS[1], RAS[2],
			&surfaceRAS[0], &surfaceRAS[1], &surfaceRAS[2] );
    oSurfaceRAS[0] = surfaceRAS[0];
    oSurfaceRAS[1] = surfaceRAS[1];
    oSurfaceRAS[2] = surfaceRAS[2];
  } else {
    throw runtime_error( "Cannot transform from RAS to surface RAS because no"
			 "MRI is loaded." );
  }
}

float 
VolumeCollection::GetMRINearestValue ( VolumeLocation& iLoc ) {

  Real value = 0;
  if( NULL != mMRI ) {

    if( iLoc.mIdxi[0] >= 0 && iLoc.mIdxi[0] < mMRI->width &&
	iLoc.mIdxi[1] >= 0 && iLoc.mIdxi[1] < mMRI->height &&
	iLoc.mIdxi[2] >= 0 && iLoc.mIdxi[2] < mMRI->depth ) {
      
      switch( mMRI->type ) {
      case MRI_UCHAR:
	value = (float)MRIvox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2] );
	break ;
      case MRI_SHORT:
	value =(float)MRISvox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2] );
	break ;
      case MRI_INT:
	value =(float)MRIIvox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2] );
	break ;
      case MRI_FLOAT:
	value = MRIFvox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2] );
	break ;
      default:
	value = 0;
      }
    }
  }
  
  return (float)value;
  
}

float
VolumeCollection::GetMRINearestValueAtIndexUnsafe ( int iIndex[3] ) {

  Real value = 0;
  
  switch( mMRI->type ) {
  case MRI_UCHAR:
    value = (float)MRIvox(mMRI, iIndex[0], iIndex[1], iIndex[2] );
    break ;
  case MRI_SHORT:
    value = (float)MRISvox(mMRI, iIndex[0], iIndex[1], iIndex[2] );
    break ;
  case MRI_INT:
    value = (float)MRIIvox(mMRI, iIndex[0], iIndex[1], iIndex[2] );
    break ;
  case MRI_FLOAT:
    value = MRIFvox(mMRI, iIndex[0], iIndex[1], iIndex[2] );
    break ;
  default:
    value = 0;
  }
  
  return (float)value;
}

float 
VolumeCollection::GetMRITrilinearValue ( VolumeLocation& iLoc ) {

  Real value = 0;
  if( NULL != mMRI ) {

    if( iLoc.mIdxi[0] >= 0 && iLoc.mIdxi[0] < mMRI->width &&
	iLoc.mIdxi[1] >= 0 && iLoc.mIdxi[1] < mMRI->height &&
	iLoc.mIdxi[2] >= 0 && iLoc.mIdxi[2] < mMRI->depth ) {
      
      MRIsampleVolumeType( mMRI, iLoc.mIdxf[0], iLoc.mIdxf[1], iLoc.mIdxf[2],
			   &value, SAMPLE_TRILINEAR );
    }
  }
  return (float)value;
}

float 
VolumeCollection::GetMRISincValue ( VolumeLocation& iLoc ) {
  
  Real value = 0;
  if( NULL != mMRI ) {

    if( iLoc.mIdxi[0] >= 0 && iLoc.mIdxi[0] < mMRI->width &&
	iLoc.mIdxi[1] >= 0 && iLoc.mIdxi[1] < mMRI->height &&
	iLoc.mIdxi[2] >= 0 && iLoc.mIdxi[2] < mMRI->depth ) {
      
      MRIsampleVolumeType( mMRI, iLoc.mIdxf[0], iLoc.mIdxf[1], iLoc.mIdxf[2],
			   &value, SAMPLE_SINC );
    }
  }
  return (float)value;
}

void
VolumeCollection::SetMRIValue ( VolumeLocation& iLoc,
					  float iValue ) {

  if( NULL != mMRI ) {
    switch( mMRI->type ) {
    case MRI_UCHAR:
      MRIvox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] ) =
	(BUFTYPE) iValue;
      break ;
    case MRI_SHORT:
      MRISvox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] ) = 
	(short) iValue;
      break ;
    case MRI_FLOAT:
      MRIFvox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] ) = 
	(float) iValue;
      break ;
    case MRI_LONG:
      MRILvox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] ) = 
	(long) iValue;
      break ;
    case MRI_INT:
      MRIIvox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] ) = 
	(int) iValue;
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
VolumeCollection::GetMRIMagnitudeValue ( VolumeLocation& iLoc ) {

  Real value = 0;

  // If we don't have the magnitude volume, calculate it.
  if( NULL == mMagnitudeMRI ) {
    MakeMagnitudeVolume();
  }

  // Get the value.
  if( NULL != mMagnitudeMRI ) {
    value = MRIFvox( mMagnitudeMRI, iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] );
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
  
  // LoadVolumeFromFileName <collectionID>
  if( 0 == strcmp( isCommand, "LoadVolumeFromFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      LoadVolume();
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
      catch( runtime_error& e ) {
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

  // SetVolumeAutosaveOn <collectionID> <on>
  if( 0 == strcmp( isCommand, "SetVolumeAutosaveOn" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      try {
	bool bOn =
	  TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
	mbAutosave = bOn;
      }
      catch( runtime_error& e ) {
	sResult = "bad on \"" + string(iasArgv[2]) + "\"," + e.what();
	return error;	
      }
    }
  }
  
  // GetVolumeAutosaveOn <collectionID>
  if( 0 == strcmp( isCommand, "GetVolumeAutosaveOn" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {

      sReturnValues =
	TclCommandManager::ConvertBooleanToReturnValue( mbAutosave );
      sReturnFormat = "i";
    }
  }

  // GetRASCoordsFromVolumeSurfaceRAS <colllectionID> <x> <y> <z>
  if( 0 == strcmp( isCommand, "GetRASCoordsFromVolumeSurfaceRAS" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }
    
    if( mID == collectionID ) {

      float surfaceRAS[3];
      surfaceRAS[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
      surfaceRAS[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
      surfaceRAS[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );

      float ras[3];
      SurfaceRASToRAS( surfaceRAS, ras );

      stringstream ssReturnValues;
      ssReturnValues << ras[0] << " " << ras[1] << " " << ras[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }
  
  // GetVolumeSurfaceRASCoordsFromRAS <colllectionID> <x> <y> <z>
  if( 0 == strcmp( isCommand, "GetVolumeSurfaceRASCoordsFromRAS" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }
    
    if( mID == collectionID ) {

      float ras[3];
      ras[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
      ras[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
      ras[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );

      float surfaceRAS[3];
      RASToSurfaceRAS( ras, surfaceRAS );

      stringstream ssReturnValues;
      ssReturnValues << surfaceRAS[0] << " " << surfaceRAS[1] << " " << surfaceRAS[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }

  // GetVolumeAverageValueInROI <colllectionID> <roiID>
  if( 0 == strcmp( isCommand, "GetVolumeAverageValueInROI" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }
    
    if( mID == collectionID ) {

      int roiID;
      try {
	roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      }
      catch( runtime_error& e ) {
	sResult = string("bad roiID: ") + e.what();
	return error;
      }
    
      map<int,ScubaROI*>::iterator tIDROI;
      tIDROI = mROIMap.find( roiID );
      if( tIDROI != mROIMap.end() ) {
	ScubaROIVolume* roi = (ScubaROIVolume*)(*tIDROI).second;

	float averageValue;
	averageValue = GetAverageValue( *roi );
	stringstream ssReturnValues;
	ssReturnValues << averageValue;
	sReturnValues = ssReturnValues.str();
	sReturnFormat = "f";

      } else {
	sResult = string("ROI doesn't belong to that collection.");
	return error;
      }	

      return ok;
    }
  }

  // GetVolumeStandardDeviationInROI <colllectionID> <roiID>
  if( 0 == strcmp( isCommand, "GetVolumeStandardDeviationInROI" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    }
    catch( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }
    
    if( mID == collectionID ) {

      int roiID;
      try {
	roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      }
      catch( runtime_error& e ) {
	sResult = string("bad roiID: ") + e.what();
	return error;
      }
    
      map<int,ScubaROI*>::iterator tIDROI;
      tIDROI = mROIMap.find( roiID );
      if( tIDROI != mROIMap.end() ) {
	ScubaROIVolume* roi = (ScubaROIVolume*)(*tIDROI).second;

	float averageValue = GetAverageValue( *roi );
	float stdDev = GetStandardDeviation( *roi, averageValue );
	stringstream ssReturnValues;
	ssReturnValues << stdDev;
	sReturnValues = ssReturnValues.str();
	sReturnFormat = "f";

      } else {
	sResult = string("ROI doesn't belong to that collection.");
	return error;
      }	

      return ok;
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
  
  int color[] = { 0, 0, 255 };
  roi->SetColor( color );

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
VolumeCollection::Select ( VolumeLocation& iLoc ) {

  if( mSelectedROIID >= 0 ) {

    ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;

    // Selectg in the ROI.
    volumeROI->SelectVoxel( iLoc.mIdxi );

    // Also mark this in the selection voxel.
    mSelectedVoxels->Set_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2],
				 true );
  }
}

void 
VolumeCollection::Unselect ( VolumeLocation& iLoc ) {
  
  if( mSelectedROIID >= 0 ) {

    ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
    //    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    ScubaROIVolume* volumeROI = (ScubaROIVolume*)roi;

    // Unselect in the ROI.
    volumeROI->UnselectVoxel( iLoc.mIdxi );

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

      if( volumeROI->IsVoxelSelected( iLoc.mIdxi ) ) {
	bSelected = true;
	break;
      }
    }
    if( !bSelected ) {
      mSelectedVoxels->Set_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2],
				   false );
    }
  }
}

bool 
VolumeCollection::IsSelected ( VolumeLocation& iLoc, int oColor[3] ) {

  // Check the selection volume cache first.
  if( !(mSelectedVoxels->Get_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] )) )
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
      if( volumeROI->IsVoxelSelected( iLoc.mIdxi ) ) {
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
  catch( runtime_error& e ) {
    cerr << "Error in IsSelected(): " << e.what() << endl;
    return false;
  }
  catch(...) {
    return false;
  }
  
}

bool 
VolumeCollection::IsSelected ( VolumeLocation& iLoc ) {

  // Check the selection volume cache.
  return 
    mSelectedVoxels->Get_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] );
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

#define PRINTOUT 0

void
VolumeCollection::FindRASPointsInSquare ( float iCenter[3],
					  float iPointA[3], float iPointB[3],
					  float iPointC[3], float iPointD[3],
					  float,
					  list<Point3<float> >& oPoints ) {

  Point3<float> planeRAS, squareRAS[4];
  squareRAS[0].Set( iPointA );
  squareRAS[1].Set( iPointB );
  squareRAS[2].Set( iPointC );
  squareRAS[3].Set( iPointD );
  planeRAS.Set( iCenter );

  // Find a plane normal from the four points.
  Point3<float> vRAS[2];
  vRAS[0] = squareRAS[1] - squareRAS[0];
  vRAS[1] = squareRAS[2] - squareRAS[0];
  Point3<float> planeNRAS = VectorOps::Cross( vRAS[0], vRAS[1] );

  // Get our plane and normal in floating index coords.
  Point3<float> planeIdx, planeNIdx;
  RASToMRIIndex( planeRAS.xyz(), planeIdx.xyz() );

  // Only use the rotation part of the RAS->MRI index here for the
  // normal.
  Matrix44 r = mWorldToIndexTransform.GetMainMatrix().ExtractRotation();
  r.MultiplyVector3( planeNRAS.xyz(), planeNIdx.xyz() );

  // Find the MRI indices. Find the square-to-volume bound.
  Point3<float> squareIdxf[4];
  RASToMRIIndex( squareRAS[0].xyz(), squareIdxf[0].xyz() );
  RASToMRIIndex( squareRAS[1].xyz(), squareIdxf[1].xyz() );
  RASToMRIIndex( squareRAS[2].xyz(), squareIdxf[2].xyz() );
  RASToMRIIndex( squareRAS[3].xyz(), squareIdxf[3].xyz() );
  Point3<int> squareIdx[4];
  RASToMRIIndex( squareRAS[0].xyz(), squareIdx[0].xyz() );
  RASToMRIIndex( squareRAS[1].xyz(), squareIdx[1].xyz() );
  RASToMRIIndex( squareRAS[2].xyz(), squareIdx[2].xyz() );
  RASToMRIIndex( squareRAS[3].xyz(), squareIdx[3].xyz() );

#if PRINTOUT
  cerr << "Square is " << squareRAS[0] << " " << squareRAS[1] << " " 
       << squareRAS[2] << " " << squareRAS[3] << endl;
  cerr << "Square is " << squareIdx[0] << " " << squareIdx[1] << " " 
       << squareIdx[2] << " " << squareIdx[3] << endl;
  cerr << "Plane RAS is " << planeRAS << ", N is " << planeNRAS << endl;
  cerr << "Plane idx is " << planeIdx << ", N is " << planeNIdx << endl;
#endif

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

  // For each voxel in the cuboid...
  Point3<float> intersectionIdx;
  for( int nZ = volumeBoundIdx[0].z(); nZ <= volumeBoundIdx[1].z(); nZ++ ) {
    for( int nY = volumeBoundIdx[0].y(); nY <= volumeBoundIdx[1].y(); nY++ ) {
      for( int nX = volumeBoundIdx[0].x(); nX <= volumeBoundIdx[1].x(); nX++ ){
	
	Point3<int> curMRIIndex( nX, nY, nZ );
	VectorOps::IntersectionResult rInt =
	  VoxelIntersectsPlane( curMRIIndex, planeIdx, planeNIdx, 
				intersectionIdx );

	if( VectorOps::intersect == rInt ) {
#if PRINTOUT
	  cerr << "\t\thit";
#endif
	  
	  // Calculate the anglesum of the intersection point with the
	  // four plane corner points. If it is 2pi, this point is
	  // inside the poly.
	  double angleSum = 0;
	  for( int nVector = 0; nVector < 4; nVector++ ) {
	    Point3<float> v1 = squareIdxf[nVector]       - intersectionIdx;
	    Point3<float> v2 = squareIdxf[(nVector+1)%4] - intersectionIdx;
	    
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

	    // Add the RAS.
	    Point3<float> curRAS;
	    MRIIndexToRAS( curMRIIndex.xyz(), curRAS.xyz() );
	    
	    oPoints.push_back( curRAS );
#if PRINTOUT
	    Point3<float> curVoxelF;
	    Point3<int> curVoxelI;
	    RASToMRIIndex( curRAS.xyz(), curVoxelI.xyz() );
	    RASToMRIIndex( curRAS.xyz(), curVoxelF.xyz() );
	    cerr << "\n\t\tadded index " << curVoxelI << " " << curVoxelF << " RAS " << curRAS << endl;
#endif
	  
#if PRINTOUT
	  } else {
	    cerr << " but not in poly, angleSum = " << angleSum << endl;
#endif
	  }

	} // if intersect
#if PRINTOUT
	else {
	  cerr << "   didn't hit" << endl;
	}
      cerr << "DONE x = " << nX << " y = " << nY << " z = " << nZ << endl;
#endif
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

  // Get a list of RAS voxels in the square.
  list<Point3<float> > squarePoints;
  FindRASPointsInSquare( iCenter, iPointA, iPointB, iPointC, iPointD, 
			 iMaxDistance, squarePoints );

  // For each one of those, check if it's within the circle.
  Point3<float> center( iCenter );
  list<Point3<float> >::iterator tPoints;
  for( tPoints = squarePoints.begin(); tPoints != squarePoints.end();
       ++tPoints ) {

    // Get the point. This is actually the corner of the voxel in
    // RAS. We want to take that, convert to index, add 0.5 to it, and
    // convert back to RAS to get the RAS of the center of the voxel.
    Point3<float> pointRAS = *tPoints;
    Point3<float> pointIdx;
    Point3<float> centerRAS;

    RASToMRIIndex( pointRAS.xyz(), pointIdx.xyz() );
    pointIdx[0] += 0.5; pointIdx[1] += 0.5; pointIdx[2] += 0.5;
    MRIIndexToRAS( pointIdx.xyz(), centerRAS.xyz() );

    // If the center voxel is within the radius, add the original
    // corner RAS point.
    if( VectorOps::Distance( centerRAS, center ) <= iRadius ) {
      oPoints.push_back( pointRAS );
    }
  }
}

void 
VolumeCollection::FindRASPointsOnSegment ( float iPointA[3], float iPointB[3],
					 std::list<Point3<float> >& oPoints ) {

  // Convert the end points to MRI index, in float and int.
  Point3<float> segIdxf[2];
  RASToMRIIndex( iPointA, segIdxf[0].xyz() );
  RASToMRIIndex( iPointB, segIdxf[1].xyz() );
  Point3<int> segIdx[2];
  RASToMRIIndex( iPointA, segIdx[0].xyz() );
  RASToMRIIndex( iPointB, segIdx[1].xyz() );

  // Find the bounds in indices of the cuboid made by the line.
  Point3<int> volumeBoundIdx[2];
  volumeBoundIdx[0].Set( MIN(segIdx[0][0],segIdx[1][0]),
			 MIN(segIdx[0][1],segIdx[1][1]),
			 MIN(segIdx[0][2],segIdx[1][2]) );
  volumeBoundIdx[1].Set( MAX(segIdx[0][0],segIdx[1][0]),
			 MAX(segIdx[0][1],segIdx[1][1]),
			 MAX(segIdx[0][2],segIdx[1][2]) );

  // For each voxel in the cuboid...
  Point3<float> intersectionIdx;
  for( int nZ = volumeBoundIdx[0].z(); nZ <= volumeBoundIdx[1].z(); nZ++ ) {
    for( int nY = volumeBoundIdx[0].y(); nY <= volumeBoundIdx[1].y(); nY++ ) {
      for( int nX = volumeBoundIdx[0].x(); nX <= volumeBoundIdx[1].x(); nX++ ){
	
	// If it intersects the segment...
	Point3<int> curMRIIndex( nX, nY, nZ );
	VectorOps::IntersectionResult rInt =
	  VoxelIntersectsSegment( curMRIIndex, segIdxf[0], segIdxf[1], 
				  intersectionIdx );
	
	// Add it to the list to return.
	if( VectorOps::intersect == rInt ) {

	  Point3<float> curRAS;
	  MRIIndexToRAS( curMRIIndex.xyz(), curRAS.xyz() );
	  oPoints.push_back( curRAS );
	}
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
	    VolumeLocation& loc = (VolumeLocation&) MakeLocationFromRAS( ras );

	    float surfaceRAS[3];
	    RASToSurfaceRAS( ras, surfaceRAS );

	    label->lv[nPoint].x = surfaceRAS[0];
	    label->lv[nPoint].y = surfaceRAS[1];
	    label->lv[nPoint].z = surfaceRAS[2];
	    label->lv[nPoint].stat = GetMRINearestValue( loc );
	    label->lv[nPoint].vno = -1;
	    label->lv[nPoint].deleted = false;
	    
	    nPoint++;

	    delete &loc;
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

    float surfaceRAS[3];
    surfaceRAS[0] = label->lv[nPoint].x;
    surfaceRAS[1] = label->lv[nPoint].y;
    surfaceRAS[2] = label->lv[nPoint].z;

    float ras[3];
    SurfaceRASToRAS( surfaceRAS, ras );

    int index[3];
    RASToMRIIndex( ras, index );

    volumeROI->SelectVoxel( index );

    // Also mark this in the selection voxel.
    mSelectedVoxels->Set_Unsafe( index[0], index[1], index[2], true );
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

  MRIcopyHeader( mMRI, segVolume );

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

  // Don't set if we're already using this one.
  if( NULL != mDataToWorldTransform &&
      iTransformID == mDataToWorldTransform->GetID() )
    return;

  DataCollection::SetDataToWorldTransform( iTransformID );
  CalcWorldToIndexTransform();
}


Matrix44&
VolumeCollection::GetWorldToIndexTransform () {

  return mWorldToIndexTransform.GetMainMatrix(); 
}


void
VolumeCollection::SetUseWorldToIndexTransform ( bool ibUse ) {

  // Don't change if it's the current setting.
  if( ibUse == mbUseDataToIndexTransform )
    return;

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
    center.SetMainTransform( 1, 0, 0, (double)mMRI->width/2.0,
			     0, 1, 0, (double)mMRI->height/2.0,
			     0, 0, 1, (double)mMRI->depth/2.0,
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

  // Set all bins to zero.
  for( int nBin = 0; nBin < icBins; nBin++ ) {
    oBinCounts[nBin] = 0;
  }

  // Find values for all the RAS points. Save the highs and lows.
  float low = 999999;
  float high = -999999;
  list<Point3<float> >::iterator tPoint;
  list<float> values;
  for( tPoint = iRASPoints.begin(); tPoint != iRASPoints.end(); ++tPoint ) {
    Point3<float> point = *tPoint;
    VolumeLocation& loc = (VolumeLocation&) MakeLocationFromRAS( point.xyz());
    float value = GetMRINearestValue( loc );
    if( value < low ) { low = value; }
    if( value > high ) { high = value; }
    values.push_back( value );
    delete &loc;
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

void
VolumeCollection::MakeHistogram ( int icBins,
				  float iMinThresh, float iMaxThresh,
				  float& oMinBinValue, float& oBinIncrement,
				  map<int,int>& oBinCounts ) {

  if( iMinThresh < mMRIMinValue || iMinThresh > mMRIMaxValue ||
      iMaxThresh < mMRIMinValue || iMaxThresh > mMRIMaxValue ) {
    throw runtime_error( "Thresh values out of bounds." );
  }

  if( iMinThresh > iMaxThresh ) {
    throw runtime_error( "Min thresh is greater than max thresh." );
  }

  if( icBins <= 0 ) {
    throw runtime_error( "Number of bins is less than zero." );
  }

  float binIncrement = (mMRIMaxValue - mMRIMinValue) / (float)icBins;

  // Set all bins to zero.
  for( int nBin = 0; nBin < icBins; nBin++ ) {
    oBinCounts[nBin] = 0;
  }

  int index[3];
  for( index[2] = 0; index[2] < mMRI->depth; index[2]++ ) {
    for( index[1] = 0; index[1] < mMRI->height; index[1]++ ) {
      for( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {
	
	float value = GetMRINearestValueAtIndexUnsafe( index );
	if( value < iMinThresh || value > iMaxThresh )
	  continue;

	int nBin = (int)floor( (value - mMRIMinValue) / (float)binIncrement );
	if( nBin == icBins && value == iMaxThresh ) { 
	  nBin = icBins-1;
	}
	oBinCounts[nBin]++;

      }
    }
  }

  oMinBinValue = mMRIMinValue;
  oBinIncrement = binIncrement;
}

void
VolumeCollection::DataChanged () {

  mbAutosaveDirty = true;
  DataCollection::DataChanged();
}

string
VolumeCollection::MakeAutosaveFileName ( string& ifn ) {

  // Generate an autosave name.
  string fnAutosave = ifn;
  if ( fnAutosave == "" ) {
    fnAutosave = "newvolume";
  }
  string::size_type c = fnAutosave.find( '/' );
  while( c != string::npos ) {
    fnAutosave[c] = '.';
    c = fnAutosave.find( '/' );
  }
  fnAutosave = "/tmp/" + fnAutosave + ".mgz";

  return fnAutosave;
}

void
VolumeCollection::DeleteAutosave () {

  char* fnAutosave = strdup( mfnAutosave.c_str() );

  struct stat info;
  int rStat = stat( fnAutosave, &info );

  if( 0 == rStat ) {
    if( S_ISREG(info.st_mode) ) {
      remove( fnAutosave );
    }
  }

  free( fnAutosave );
}

void
VolumeCollection::AutosaveIfDirty () {

  if( mbAutosave && mbAutosaveDirty ) {
    Save( mfnAutosave );
    mbAutosaveDirty = false;
  }
}

void
VolumeCollection::PrintVoxelCornerCoords ( ostream& iStream,
					   Point3<int>& iMRIIndex ) {

  // Create RAS versions of our corners.
  Point3<int> voxelIdx[8];
  voxelIdx[0].Set( iMRIIndex.x()    , iMRIIndex.y()    , iMRIIndex.z()   );
  voxelIdx[1].Set( iMRIIndex.x()+1, iMRIIndex.y()    , iMRIIndex.z()   );
  voxelIdx[2].Set( iMRIIndex.x()    , iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[3].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[4].Set( iMRIIndex.x()    , iMRIIndex.y()    , iMRIIndex.z()+1 );
  voxelIdx[5].Set( iMRIIndex.x()+1, iMRIIndex.y()    , iMRIIndex.z()+1 );
  voxelIdx[6].Set( iMRIIndex.x()    , iMRIIndex.y()+1, iMRIIndex.z()+1 );
  voxelIdx[7].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()+1 );
  Point3<float> voxelRAS[8];
  for( int nCorner = 0; nCorner < 8; nCorner++ ) {
    MRIIndexToRAS( voxelIdx[nCorner].xyz(), voxelRAS[nCorner].xyz() );
    iStream << nCorner << ": " << voxelRAS[nCorner] << endl;
  }
}

void
VolumeCollection::GetVoxelsWithValue ( float iValue,
				       list<VolumeLocation>& olLocations ) {

  for( int nZ = 0; nZ < mMRI->depth; nZ++ ) {
    for( int nY = 0; nY < mMRI->height; nY++ ) {
      for( int nX = 0; nX < mMRI->width; nX++ ) {
	
	float value;
	switch( mMRI->type ) {
	case MRI_UCHAR:
	  value = (float)MRIvox(mMRI, nX, nY, nZ );
	  break;
	case MRI_SHORT:
	  value = (float)MRISvox( mMRI, nX, nY, nZ );
	  break;
	case MRI_INT:
	  value = (float)MRIIvox( mMRI, nX, nY, nZ );
	  break;
	case MRI_FLOAT:
	  value = MRIFvox( mMRI, nX, nY, nZ );
	  break;
	default:
	  value = 0;
	}
	
	if( fabs( value - iValue ) < EPSILON ) {
	  int index[3] = { nX, nY, nZ };
	  VolumeLocation& loc = 
	    (VolumeLocation&) MakeLocationFromIndex( index );
	  olLocations.push_back( loc );
	}
      }
    }
  }
}

float 
VolumeCollection::GetRASVolumeOfNVoxels ( int icVoxels ) {

  return GetVoxelXSize() * GetVoxelYSize() * GetVoxelZSize() * (float)icVoxels;
}

float
VolumeCollection::GetPreferredValueIncrement () {

  // Look at the range. If it's > 100, inc is 1, 10-100, inc is .1,
  // 1-10, inc is .01, etc.
  float range = mMRIMaxValue - mMRIMinValue;
  float inc = 1;
       if( range >= 1000000 )        { inc = 1000; }
  else if( range >=  100000 )        { inc =  100; }
  else if( range >=   10000 )        { inc =   10; }
  else if( range >=    1000 )        { inc =    1; }
  else if( range >=      10 )        { inc =    0.1; }
  else if( range >=       1 )        { inc =    0.01; }
  else if( range >=       0.1 )      { inc =    0.001; }
  else if( range >=       0.01 )     { inc =    0.0001; }
  else if( range >=       0.001 )    { inc =    0.00001; }
  else if( range >=       0.0001 )   { inc =    0.000001; }
  else if( range >=       0.00001 )  { inc =    0.0000001; }
  else if( range >=       0.000001 ) { inc =    0.00000001; }
       else                          { inc =    0.000000001; }

       return inc;
}

float 
VolumeCollection::GetAverageValue ( list<VolumeLocation>& ilLocations ) {
  
  if( ilLocations.size() < 1 ) {
    throw runtime_error( "No voxels" );
  }

  float sumValue = 0;
  list<VolumeLocation>::iterator tLocation;
  for( tLocation = ilLocations.begin(); tLocation != ilLocations.end();
       ++tLocation ) {

    VolumeLocation loc = *tLocation;
    sumValue += GetMRINearestValue( loc );
  }

  return sumValue / (float)ilLocations.size();
}

float
VolumeCollection::GetAverageValue ( ScubaROIVolume& iROI ) {
  
  list<VolumeLocation> lLocs;
  int voxel[3];
  for( voxel[2] = 0; voxel[2] < mMRI->depth; voxel[2]++ ) {
    for( voxel[1] = 0; voxel[1] < mMRI->height; voxel[1]++ ) {
      for( voxel[0] = 0; voxel[0] < mMRI->width; voxel[0]++ ) {
	
	if( iROI.IsVoxelSelected( voxel ) ) {

	  VolumeLocation& loc = 
	    (VolumeLocation&) MakeLocationFromIndex( voxel );
	  lLocs.push_back( loc );
	  delete &loc;
	}
      }
    }
  }

  return GetAverageValue( lLocs );
}

float 
VolumeCollection::GetStandardDeviation ( list<VolumeLocation>& ilLocations,
					 float iMean) {

  if( ilLocations.size() < 1 ) {
    throw runtime_error( "No voxels" );
  }

  float sumDeviation = 0;
  list<VolumeLocation>::iterator tLocation;
  for( tLocation = ilLocations.begin(); tLocation != ilLocations.end();
       ++tLocation ) {

    VolumeLocation loc = *tLocation;
    float value = GetMRINearestValue( loc );
    float deviation = (value - iMean) * (value - iMean);
    sumDeviation += deviation;
  }

  return sqrt( sumDeviation / (float)ilLocations.size() );
}

float
VolumeCollection::GetStandardDeviation ( ScubaROIVolume& iROI, float iMean ) {
  
  list<VolumeLocation> lLocs;
  int voxel[3];
  for( voxel[2] = 0; voxel[2] < mMRI->depth; voxel[2]++ ) {
    for( voxel[1] = 0; voxel[1] < mMRI->height; voxel[1]++ ) {
      for( voxel[0] = 0; voxel[0] < mMRI->width; voxel[0]++ ) {
	
	if( iROI.IsVoxelSelected( voxel ) ) {
	  
	  VolumeLocation& loc = 
	    (VolumeLocation&) MakeLocationFromIndex( voxel );
	  lLocs.push_back( loc );
	  delete &loc;
	}
      }
    }
  }
  
  return GetStandardDeviation( lLocs, iMean );
}

VectorOps::IntersectionResult 
VolumeCollection::VoxelIntersectsPlane ( Point3<int>& iMRIIndex, 
					 Point3<float>& iPlaneIdx, 
					 Point3<float>& iPlaneIdxNormal,
					 Point3<float>& oIntersectionIdx ){

  // Create float idx versions of our corners.
  Point3<float> voxelIdx[8];
  voxelIdx[0].Set( iMRIIndex.x()  , iMRIIndex.y()  , iMRIIndex.z()   );
  voxelIdx[1].Set( iMRIIndex.x()+1, iMRIIndex.y()  , iMRIIndex.z()   );
  voxelIdx[2].Set( iMRIIndex.x()  , iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[3].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[4].Set( iMRIIndex.x()  , iMRIIndex.y()  , iMRIIndex.z()+1 );
  voxelIdx[5].Set( iMRIIndex.x()+1, iMRIIndex.y()  , iMRIIndex.z()+1 );
  voxelIdx[6].Set( iMRIIndex.x()  , iMRIIndex.y()+1, iMRIIndex.z()+1 );
  voxelIdx[7].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()+1 );
  
#if PRINTOUT
  cerr << "Looking at idx " << voxelIdx[0] << endl;
#endif
  // Make segments for each edge.
  Point3<float> segmentIdx[12][2];
  int anSegments[12][2] = { {0, 1}, {4, 5}, {6, 7}, {2, 3},
			    {0, 2}, {1, 3}, {4, 6}, {5, 7},
			    {0, 4}, {1, 5}, {2, 6}, {3, 7} };
  for( int nSegment = 0; nSegment < 12; nSegment++ ) {
    segmentIdx[nSegment][0].Set( voxelIdx[anSegments[nSegment][0]] );
    segmentIdx[nSegment][1].Set( voxelIdx[anSegments[nSegment][1]] );
  }
  
  // Intersect these segments with the plane. If any of them hit...
  for( int nSegment = 0; nSegment < 12; nSegment++ ) {
    
    Point3<float> intersectionIdx;
    VectorOps::IntersectionResult rInt =
      VectorOps::SegmentIntersectsPlane
      ( segmentIdx[nSegment][0], segmentIdx[nSegment][1],
	iPlaneIdx, iPlaneIdxNormal, intersectionIdx );
    
#if PRINTOUT
    cerr << "\tTesting " << segmentIdx[nSegment][0] << " to " 
	       <<  segmentIdx[nSegment][1] << "..." << endl;
#endif

    if( VectorOps::intersect == rInt ) {

      oIntersectionIdx = intersectionIdx;
      return rInt;
    }
  }

  return VectorOps::dontIntersect;
}

VectorOps::IntersectionResult 
VolumeCollection::VoxelIntersectsSegment( Point3<int>& iMRIIndex, 
					  Point3<float>& iSegIdxA, 
					  Point3<float>& iSegIdxB, 
					  Point3<float>& oIntersectionIdx ) {

#if PRINTOUT
  cerr << "VoxelIntersectsSegment idx " << iMRIIndex
       << " seg " << iSegIdxA << ", " << iSegIdxB << endl;
#endif
  
  // Create float idx versions of our corners.
  Point3<float> voxelIdx[8];
  voxelIdx[0].Set( iMRIIndex.x()  , iMRIIndex.y()  , iMRIIndex.z()   );
  voxelIdx[1].Set( iMRIIndex.x()+1, iMRIIndex.y()  , iMRIIndex.z()   );
  voxelIdx[2].Set( iMRIIndex.x()  , iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[3].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()   );
  voxelIdx[4].Set( iMRIIndex.x()  , iMRIIndex.y()  , iMRIIndex.z()+1 );
  voxelIdx[5].Set( iMRIIndex.x()+1, iMRIIndex.y()  , iMRIIndex.z()+1 );
  voxelIdx[6].Set( iMRIIndex.x()  , iMRIIndex.y()+1, iMRIIndex.z()+1 );
  voxelIdx[7].Set( iMRIIndex.x()+1, iMRIIndex.y()+1, iMRIIndex.z()+1 );

  // Make plane corners for each face.
  Point3<float> planeIdx[6][4]; // 0-5 face, 0-3 corner
  int anFaces[6][4] = { {0, 1, 3, 2}, {1, 5, 7, 3},
			{0, 4, 5, 1}, {4, 5, 7, 6},
			{0, 4, 6, 2}, {2, 6, 7, 3} };

  for( int nFace = 0; nFace < 6; nFace++ ) {
    planeIdx[nFace][0].Set( voxelIdx[anFaces[nFace][0]] );
    planeIdx[nFace][1].Set( voxelIdx[anFaces[nFace][1]] );
    planeIdx[nFace][2].Set( voxelIdx[anFaces[nFace][2]] );
    planeIdx[nFace][3].Set( voxelIdx[anFaces[nFace][3]] );
  }

  Point3<float> planeNIdx[6];
  planeNIdx[0].Set( 0, 0, 1 ); planeNIdx[3].Set( 0, 0, 1 );
  planeNIdx[1].Set( 1, 0, 0 ); planeNIdx[4].Set( 1, 0, 0 );
  planeNIdx[2].Set( 0, 1, 0 ); planeNIdx[5].Set( 0, 1, 0 );

  for( int nFace = 0; nFace < 6; nFace++ ) {

#if PRINTOUT
    cerr << "\tTesting face p " << planeIdx[nFace][0] 
	 << " N " << planeNIdx[nFace] << "... ";
#endif

    Point3<float> intersectionIdx;
    VectorOps::IntersectionResult rInt =
      VectorOps::SegmentIntersectsPlane
      ( iSegIdxA, iSegIdxB, planeIdx[nFace][0], planeNIdx[nFace], 
	intersectionIdx );
 
    if( VectorOps::intersect == rInt ) {

#if PRINTOUT
      cerr << " hit " << endl
	   << "\t\tIntersection " << intersectionIdx << endl;
#endif
      
      // Calculate the anglesum of the intersection point with the
      // four plane corner points. If it is 2pi, this point is
      // inside the poly.
      double angleSum = 0;
      for( int nCorner = 0; nCorner < 4; nCorner++ ) {
	Point3<float> v1 = planeIdx[nFace][nCorner]       - intersectionIdx;
	Point3<float> v2 = planeIdx[nFace][(nCorner+1)%4] - intersectionIdx;
#if PRINTOUT
	cerr << "\t\tCorner " << planeIdx[nFace][nCorner] << endl;
#endif
	
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
#if PRINTOUT
	cerr << "\t\tIn face" << endl;
#endif
	oIntersectionIdx = intersectionIdx;
	return VectorOps::intersect;
      }
#if PRINTOUT
      else {
	cerr << "\t\tNot in face" << endl;
      }
#endif
    } 
#if PRINTOUT
    else {
      cerr << " miss " << endl;
    }
#endif
  }
#if PRINTOUT
  cerr << "\tDidn't hit. " << endl;
#endif
  
  return VectorOps::dontIntersect;
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
  mViewNormal[0] = mViewNormal[1] = mViewNormal[2] = 0;
  mFuzziness    = 1;
  mMaxDistance  = 0;
  mbDiagonal    = false;
  mbOnlyZero    = false;
  mFuzziness    = seed;
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
VolumeCollectionFlooder::DoVoxel ( float[3] ) {

}

bool 
VolumeCollectionFlooder::CompareVoxel ( float[3] ) {
  return true;
}

void
VolumeCollectionFlooder::Flood ( VolumeCollection& iVolume, 
				 float iRASSeed[3], Params& iParams ) {

  mVolume = &iVolume;
  mParams = &iParams;

  // Get the seed voxel.
  Point3<int> seedVoxel;
  iVolume.RASToMRIIndex( iRASSeed, seedVoxel.xyz() );

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

  // Init a visited volume.
  Volume3<bool>* bVisited =
    new Volume3<bool>( iVolume.mMRI->width, 
		       iVolume.mMRI->height, 
		       iVolume.mMRI->depth, false );

  // Save the initial value.
  VolumeLocation& seedLoc =
    (VolumeLocation&) sourceVol->MakeLocationFromRAS( iRASSeed );
  float seedValue = sourceVol->GetMRINearestValue( seedLoc );

  // Create locations in the source volume space. These will be used
  // in the main loop.
  VolumeLocation& sourceLoc =
    (VolumeLocation&) sourceVol->MakeLocationFromRAS( iRASSeed );
   VolumeLocation& sourceFromLoc =
     (VolumeLocation&) sourceVol->MakeLocationFromRAS( iRASSeed );

#if PRINTOUT
  cerr << "\n\nSeed " << Point3<int>(seedLoc.Index())
       << ", " << Point3<float>(seedLoc.RAS())
       << endl;
#endif

  // If we're not working in 3D, we need to find plane info here. The
  // plane point is just the seed point. The normal is just in
  // relation to the view. Convert those to index coords.
  Point3<float> planeRAS, planeNRAS;
  Point3<float> planeIdx, planeNIdx;
  Point3<float> intersectionIdx;
  if( !iParams.mb3D ) {
    planeRAS.Set( iRASSeed );
    if( iParams.mbWorkPlaneX ) {
      planeNRAS.Set( 1, 0, 0 ); 
    }
    if( iParams.mbWorkPlaneY ) {
      planeNRAS.Set( 0, 1, 0 );
    }
    if( iParams.mbWorkPlaneZ ) {
      planeNRAS.Set( 0, 0, 1 );
    }

    iVolume.RASToMRIIndex( planeRAS.xyz(), planeIdx.xyz() );

    // Only use the rotation part of the RAS->MRI index here for the
    // normal.
    Matrix44 r = 
      iVolume.GetWorldToIndexTransform().ExtractRotation();
    r.MultiplyVector3( planeNRAS.xyz(), planeNIdx.xyz() );
  }

  // Start it up.
  this->DoBegin();

  // Push the seed onto the list. 
  vector<CheckPair> checkPairs;
  CheckPair seedPair( seedVoxel, seedVoxel );
  checkPairs.push_back( seedPair );
  while( checkPairs.size() > 0 &&
	 !this->DoStopRequested() ) {

    // Get this voxel, the from voxel, and RAS versions.
    CheckPair checkPair = checkPairs.back();
    checkPairs.pop_back();

    Point3<int> index = checkPair.mToIndex;
    VolumeLocation& loc = 
      (VolumeLocation&) iVolume.MakeLocationFromIndex( index.xyz() );

    Point3<int> fromIndex = checkPair.mFromIndex;
    VolumeLocation& fromLoc =
      (VolumeLocation&) iVolume.MakeLocationFromIndex( fromIndex.xyz() );
    
    Point3<float> ras;
    iVolume.MRIIndexToRAS( index.xyz(), ras.xyz() );
    Point3<float> fromRAS;
    iVolume.MRIIndexToRAS( fromIndex.xyz(), fromRAS.xyz() );

    // If we have a different source, we need to get a location in its
    // space from the RAS here. If we have the same source, this will
    // just be the same as loc and fromLoc.
    sourceLoc.SetFromRAS( ras.xyz() );
    sourceFromLoc.SetFromRAS( fromRAS.xyz() );
      
    // Check the bound of this volume and the source one.
    if( !iVolume.IsInBounds( loc ) ) { 
      delete &loc;
      delete &fromLoc;
      continue;
    }
    if( bDifferentSource && !sourceVol->IsInBounds( sourceLoc ) ) { 
      delete &loc;
      delete &fromLoc;
      continue;
    }

    if( bVisited->Get_Unsafe( index.x(), index.y(), index.z() ) ) {
      delete &loc;
      delete &fromLoc;
      continue;
    }
    bVisited->Set_Unsafe( index.x(), index.y(), index.z(), true );


    if( iParams.mbStopAtROIs ) {
      if( iVolume.IsOtherRASSelected( ras.xyz(), iVolume.GetSelectedROI() ) ) {
	delete &loc;
	delete &fromLoc;
	continue;
      }
    }
    
    // Check max distance.
    if( iParams.mMaxDistance > 0 ) {
      float distance = sqrt( ((ras[0]-iRASSeed[0]) * (ras[0]-iRASSeed[0])) + 
			     ((ras[1]-iRASSeed[1]) * (ras[1]-iRASSeed[1])) + 
			     ((ras[2]-iRASSeed[2]) * (ras[2]-iRASSeed[2])) );
      if( distance > iParams.mMaxDistance ) {
	delete &loc;
	delete &fromLoc;
	continue;
      }
    }


    // Check only zero.
    if( iParams.mbOnlyZero ) {
      float value = iVolume.GetMRINearestValue( loc );
      if( value != 0 ) {
	delete &loc;
	delete &fromLoc;
	continue;
      }
    }

    // Check fuzziness.
    if( iParams.mFuzziness > 0 ) {
      float value = sourceVol->GetMRINearestValue( sourceLoc );
      switch( iParams.mFuzzinessType ) {
      case Params::seed:
	if( fabs( value - seedValue ) > iParams.mFuzziness ) {
	  delete &loc;
	  delete &fromLoc;
	  continue;
	}
	break;
      case Params::gradient: {
	float fromValue = 
	  sourceVol->GetMRINearestValue( sourceFromLoc );
	if( fabs( value - fromValue ) > iParams.mFuzziness ) {
	  delete &loc;
	  delete &fromLoc;
	  continue;
	}
      } break;
      }
    }

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
	    
	    // Get the two vertices.
	    Point3<float>& curVertex  = path->GetVertexAtIndex( nCurVertex );
	    Point3<float>& backVertex = path->GetVertexAtIndex( nBackVertex );

	    Point3<float> viewNormal( iParams.mViewNormal );

	    // Now convert everything to index coords.
	    Point3<int> curVertexIdx, backVertexIdx, viewIdx,
	      fromIdx, curIdx;
	    iVolume.RASToMRIIndex( curVertex.xyz(), curVertexIdx.xyz() );
	    iVolume.RASToMRIIndex( backVertex.xyz(), backVertexIdx.xyz() );
	    iVolume.RASToMRIIndex( viewNormal.xyz(), viewIdx.xyz() );
	    iVolume.RASToMRIIndex( fromRAS.xyz(), fromIdx.xyz() );
	    iVolume.RASToMRIIndex( ras.xyz(), curIdx.xyz() );

	    // And then get float versions of those index coords
	    // because that's what we do the math in.
	    Point3<float> curVertexIdxf, backVertexIdxf, segVertexIdxf,
	      viewIdxf, normalIdxf, fromIdxf, curIdxf;
	    curVertexIdxf.Set(curVertexIdx[0],curVertexIdx[1],curVertexIdx[2]);
	    backVertexIdxf.Set(backVertexIdx[0],backVertexIdx[1],backVertexIdx[2]);
	    fromIdxf.Set( fromIdx[0], fromIdx[1], fromIdx[2] );
	    curIdxf.Set( curIdx[0], curIdx[1], curIdx[2] );
	    viewIdxf.Set( viewIdx[0], viewIdx[1], viewIdx[2]);

	    // Calculate a vector between the two verts. Calc a plane
	    // normal by cross that vector and the view normal.
	    segVertexIdxf = curVertexIdxf - backVertexIdxf;
	    normalIdxf = VectorOps::Cross( segVertexIdxf, viewIdxf );
	    normalIdxf = VectorOps::Normalize( normalIdxf );

	    VectorOps::IntersectionResult rInt =
	      VectorOps::SegmentIntersectsPlane( fromIdxf, curIdxf,
						 curVertexIdxf, normalIdxf,
						 x );

#if 0
	    cerr << "SegmentIntersectsPlane results: "  << endl
		 << "  from\t" <<fromRAS << "\t" << fromIdx << endl
		 << "  current\t" << ras << "\t" << curIdx << endl
		 << "  curVertex " << nCurVertex << "\t" << curVertex << "\t" << curVertexIdx << endl
		 << "  backVertex " << nBackVertex << "\t" << backVertex << "\t" << backVertexIdx << endl
		 << "  normal\t" << normal << "\t" << normalIdxf << endl
		 << "  int\t" 
		 << VectorOps::IntersectionResultToString(rInt) << "\t" 
		 << VectorOps::IntersectionResultToString(rInt2) << endl;
#endif	    
	    if( VectorOps::intersect == rInt ) {

	      // Make sure x is in the cuboid formed by curVertex and
	      // backVertex.
	      if( !(x[0] < MIN(curVertexIdx[0],backVertexIdx[0]) ||
		    x[0] > MAX(curVertexIdx[0],backVertexIdx[0]) ||
		    x[1] < MIN(curVertexIdx[1],backVertexIdx[1]) ||
		    x[1] > MAX(curVertexIdx[1],backVertexIdx[1]) ||
		    x[2] < MIN(curVertexIdx[2],backVertexIdx[2]) ||
		    x[2] > MAX(curVertexIdx[2],backVertexIdx[2])) ) {
		bCross = true;
	      }
	    }
	  }
	}
      }

      // If we crossed a path, continue.
      if( bCross ) {
	delete &loc;
	delete &fromLoc;
	continue;
      }
    }

    // Call the user compare function to give them a chance to bail.
    if( !this->CompareVoxel( ras.xyz() ) ) {
      delete &loc;
      delete &fromLoc;
      continue;
    }
    
    // Call the user function.
    this->DoVoxel( ras.xyz() );

    // Add adjacent voxels. Look at all adjacent voxels. If we're only
    // working on a single plane, make sure the voxel is in the same
    // plane.
    int beginX = index.x() - 1;    int endX   = index.x() + 1;
    int beginY = index.y() - 1;    int endY   = index.y() + 1;
    int beginZ = index.z() - 1;    int endZ   = index.z() + 1;
#if PRINTOUT
    cerr << "\tbounds: " << beginX << ", " << endX << "  "
	 << beginY << ", " << endY << "  "
	 << beginZ << ", " << endZ
	 << endl;
#endif

    Point3<int> newIndex;
    CheckPair newPair;
    newPair.mFromIndex = index;
    //    if( iParams.mbDiagonal ) {
      for( int nZ = beginZ; nZ <= endZ; nZ += 1 ) {
	for( int nY = beginY; nY <= endY; nY += 1 ) {
	  for( int nX = beginX; nX <= endX; nX += 1 ) {

	    newIndex.Set( nX, nY, nZ );

	    // If not 3D, check if voxel is in same plane.
	    if( !iParams.mb3D ) {
	      VectorOps::IntersectionResult rInt =
		iVolume.VoxelIntersectsPlane( newIndex, planeIdx, planeNIdx,
					      intersectionIdx );
	      if( VectorOps::intersect != rInt ) {
		continue;
	      }
	    }

	    newPair.mToIndex = newIndex;
	    checkPairs.push_back( newPair );
	  }
	}
      }
      //    }
    delete &loc;
    delete &fromLoc;
  }

  this->DoEnd();

  delete &seedLoc;
  delete &sourceLoc;
  delete &sourceFromLoc;
  delete bVisited;

  mVolume = NULL;
  mParams = NULL;
}

VolumeLocation::VolumeLocation ( VolumeCollection& iVolume,
				 float const iRAS[3] )
  : DataLocation( iRAS ), mVolume( &iVolume ) {

  mVolume->RASToMRIIndex( iRAS, mIdxf );
  mVolume->RASToMRIIndex( iRAS, mIdxi );
#if 0
  mIdxi[0] = (int) mIdxf[0];
  mIdxi[1] = (int) mIdxf[1];
  mIdxi[2] = (int) mIdxf[2];
#endif
}

VolumeLocation::VolumeLocation ( VolumeCollection& iVolume,
				 int const iIndex[3] )
  : DataLocation(), mVolume( &iVolume ) {

  mVolume->MRIIndexToRAS( iIndex, mRAS );
  mIdxi[0] = iIndex[0];
  mIdxi[1] = iIndex[1];
  mIdxi[2] = iIndex[2];
  mIdxf[0] = iIndex[0];
  mIdxf[1] = iIndex[1];
  mIdxf[2] = iIndex[2];
}

VolumeLocation::VolumeLocation ( const VolumeLocation& iLoc )
  : DataLocation(iLoc), mVolume(iLoc.GetVolume()) {
  
  mIdxi[0] = iLoc.Index(0);
  mIdxi[1] = iLoc.Index(1);
  mIdxi[2] = iLoc.Index(2);
  mIdxf[0] = iLoc.IndexF(0);
  mIdxf[1] = iLoc.IndexF(1);
  mIdxf[2] = iLoc.IndexF(2);
}

void
VolumeLocation::SetFromRAS ( float const iRAS[3] ) {

  mRAS[0] = iRAS[0];
  mRAS[1] = iRAS[1];
  mRAS[2] = iRAS[2];
  mVolume->RASToMRIIndex( iRAS, mIdxf );
  mVolume->RASToMRIIndex( iRAS, mIdxi );
#if 0
  mIdxi[0] = (int) mIdxf[0];
  mIdxi[1] = (int) mIdxf[1];
  mIdxi[2] = (int) mIdxf[2];
#endif
}
