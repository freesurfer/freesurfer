/**
 * @file  VolumeCollection.cpp
 * @brief MRI volume access and utilities
 *
 * Provides access to MRI volume data and volume ROIs. Also contains
 * functions for determining of RAS points are in a square or circular
 * area, used by tools. Can also generate histogram data given a set
 * of points. Can generate average value and standard deviation of
 * points. Can perform value-based fills that set ranges of volumes
 * to new values.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.122 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <limits>
#include "string_fixed.h"
#include <cassert>
#include <errno.h>
#include <stdexcept>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
extern "C" {
#include <assert.h>
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

// Used in RASToTal to convert MNI Talairach coords to "real" Tal
// coords.
Matrix44 VolumeCollection::mMatrixMNITalToTalGtz( 0.99,       0,      0, 0,
						  0.00,  0.9688,  0.046, 0,
						  0.00, -0.0485, 0.9189, 0,
						  0.00,       0,      0, 1 );
Matrix44 VolumeCollection::mMatrixMNITalToTalLtz( 0.99,       0,      0, 0,
						  0.00,  0.9688,  0.042, 0,
						  0.00, -0.0485,  0.839, 0,
						  0.00,       0,      0, 1 );

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
  mcFrames = 0;
  mbInterpretFramesAsTimePoints = false;
  mcConditions = 0;
  mcTimePoints = 0;
  mcPreStimTimePoints = 0;
  mTimeResolution = 0;
  mbDataContainsErrorValues = false;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetVolumeCollectionFileName", 2,
                         "collectionID fileName",
                         "Sets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "LoadVolumeFromFileName", 1, "collectionID",
                         "Loads the volume from the file name.");
  commandMgr.AddCommand( *this, "GetVolumeCollectionFileName", 1,
                         "collectionID",
                         "Gets the file name for a given volume collection.");
  commandMgr.AddCommand( *this, "GetVolumeNumberOfFrames", 1,
                         "collectionID",
                         "Gets the number of frames in a volume collection.");
  commandMgr.AddCommand( *this, "WriteVolumeROIToLabel", 4,
                         "collectionID roiID useRealRAS fileName",
                         "Writes an ROI to a label file. useRealRAS is a "
                         "flag specifying whether to consider the RAS "
                         "coords as 'real' RAS or to convert them from "
                         "TkRegRAS space." );
  commandMgr.AddCommand( *this, "NewVolumeROIFromLabel", 3,
                         "collectionID useRealRAS fileName",
                         "Creates an ROI from a label file and returns the "
                         "ID of the new ROI. useRealRAS is a flag specifying "
                         "whether to consider the RAS coords as 'real' RAS "
                         "or to convert them from TkRegRAS space." );
  commandMgr.AddCommand( *this, "WriteVolumeROIsToSegmentation", 2,
                         "collectionID fileName",
                         "Writes a series of structure ROIs to a "
                         "segmentation volume." );
  commandMgr.AddCommand( *this, "MakeVolumeUsingTemplate", 3,
                         "collectionID templateCollectionID type",
                         "Makes a volume using an existing volume "
                         "as a template. type should be 0 for uchar, "
			 "1 for int, 2 for long, 3 for float, and 4 "
			 "for short, or -1 to copy from the template "
			 "volume." );
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
  commandMgr.AddCommand( *this, "GetVolumeHistogram", 3, "volID roiID numBins",
			 "Returns a histogram of the volume. Returns the "
			 "format: minBinValue binIncrement{binCount0 "
			 "binCount1 .. binCountN}" );
  commandMgr.AddCommand( *this, "SetVolumeAutosaveOn", 2, "collectionID on",
                         "Set whether or not autosave is on for this "
                         "volume." );
  commandMgr.AddCommand( *this, "GetVolumeAutosaveOn", 1, "collectionID",
                         "Returns whether or not autosave is on for this "
                         "volume." );
  commandMgr.AddCommand( *this, "GetRASCoordsFromVolumeTkRegRAS", 4,
                         "collectionID x y z", "Returns a list of RAS coords "
                         "converted from the input surface RAS coords. This "
                         "is for converting RAS points acquired from a "
                         "surface that is associated with a volume and "
                         "didn't generate coordinates with CRAS info." );
  commandMgr.AddCommand( *this, "GetVolumeTkRegRASCoordsFromRAS", 4,
                         "collectionID x y z", "Returns a list of surface "
                         "RAS coords converted from the input RAS coords. This"
                         "is for converting RAS points acquired from a "
                         "surface that is associated with a volume and "
                         "didn't generate coordinates with CRAS info." );
  commandMgr.AddCommand( *this, "GetVolumeAverageValueInROI", 2,
                         "collectionID roiID", "Returns the average value "
                         "of the voxels in an ROI in a volume." );
  commandMgr.AddCommand( *this, "GetVolumeStandardDeviationInROI", 2,
                         "collectionID roiID", "Returns the standard "
                         "deviation of the voxels in an ROI in a volume." );
  commandMgr.AddCommand( *this, "GetVolumeInterpretFramesAsTime", 1,
                         "collectionID", "Returns whether or not " );
  commandMgr.AddCommand( *this, "GetVolumeNumberOfConditions", 1,
                         "collectionID", "Returns the number of conditions "
                         "in the volume, if available.." );
  commandMgr.AddCommand( *this, "GetVolumeNumberOfTimePoints", 1,
                         "collectionID", "Returns the number of time points "
                         "in the volume, if available.." );
  commandMgr.AddCommand( *this, "GetVolumeTimeResolution", 1,
                         "collectionID", "Returns the time resolution of "
                         "time points in the volume, if available.." );
  commandMgr.AddCommand( *this, "GetVolumeNumberOfPreStimTimePoints", 1,
                         "collectionID", "Returns the number of time points "
                         "before the stimulus, if available." );
  commandMgr.AddCommand( *this, "DoVolumeValueRangeFill", 4,
                         "volID sourceVolID roiID valueRanges",
                         "Performs a series of value range fills using the "
                         "value ranges in sourceVolID to set values in "
                         "this volume. valueRanges should be a list of "
			 "triples: the begining of the range, the end of "
			 "the range, and the new value." );
  
}

VolumeCollection::~VolumeCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();
  try {
    mriLoader.ReleaseData( &mMRI );
    SendBroadcast( "DataDeleted", NULL );
  } catch ( exception& e ) {
    cerr << "Error releasing data: " << e.what() << endl;
  } catch (...) {
    cerr << "Couldn't release data"  << endl;
  }
}

VolumeLocation
VolumeCollection::MakeVolumeLocationFromRAS ( float const iRAS[3] ) const {

  return VolumeLocation( *this, iRAS );
}

VolumeLocation
VolumeCollection::MakeVolumeLocationFromIndex ( int const iIndex[3] ) const {

  return VolumeLocation( *this, iIndex );
}

VolumeLocation
VolumeCollection::MakeVolumeLocationFromRAS ( float const iRAS[3],
					      int iFrame ) const {

  return VolumeLocation( *this, iRAS, iFrame );
}

VolumeLocation
VolumeCollection::MakeVolumeLocationFromIndex ( int const iIndex[3], 
						int iFrame ) const {

  return VolumeLocation( *this, iIndex, iFrame );
}

void
VolumeCollection::SetFileName ( string const& ifnMRI ) {

  // Set our file name and make the autosave name.
  mfnMRI = ifnMRI;
  mfnAutosave = MakeAutosaveFileName( mfnMRI );
}

string const&
VolumeCollection::GetFileName () const {

  return mfnMRI;
}

void
VolumeCollection::MakeUsingTemplate ( int iCollectionID, int iType ) {

  if( iType < -1 || iType > MRI_SHORT ) {
    stringstream ssError;
    ssError << "iType " << iType << " was invalid.";
    throw runtime_error( ssError.str() );
  }

  // Try to find our template volume.
  VolumeCollection* vol = NULL;
  try {
    DataCollection* col = &DataCollection::FindByID( iCollectionID );
    vol = dynamic_cast<VolumeCollection*>(col);
  } catch (...) {
    throw runtime_error( "Couldn't find template." );
  }

  // Get the mri from the template volume.
  MRI* mri = const_cast<MRI*>(vol->GetMRI());
  if ( NULL == mri ) {
    throw runtime_error( "Couldn't get MRI from template" );
  }

  // If they passed -1, copy the type.
  int type = iType;
  if( iType == -1 )
    type = mri->type;

  // Allocate the mri with the size from the template.
  MRI* newMri = MRIallocSequence( mri->width, mri->height, mri->depth,
                                  type, mri->nframes );
  if ( NULL == newMri ) {
    throw runtime_error( "Couldn't allocate new mri." );
  }

  // Copy the header from the template into the new mri.
  MRIcopyHeader( mri, newMri );

  // Save the MRI.
  mMRI = newMri;

  // Initialize from it.
  InitializeFromMRI();

  // Set a temporary filename.
  SetFileName( "New_Volume.mgh" );
}

void
VolumeCollection::LoadVolume () {

  DataManager dataMgr = DataManager::GetManager();
  MRILoader mriLoader = dataMgr.GetMRILoader();

  // If we already have data...
  if ( NULL != mMRI ) {

    // Try to load this and see what we get. If it's the same as what
    // we already have, we're fine. If not, keep this one and release
    // the one we have.
    MRI* newMRI = NULL;
    try {
      newMRI = mriLoader.GetData( mfnMRI );
    } catch ( exception& e ) {
      throw logic_error( "Couldn't load MRI" );
    }

    if ( newMRI == mMRI ) {
      return;
    }

    // Release old data.
    try {
      mriLoader.ReleaseData( &mMRI );
    } catch (...) {
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
    } catch ( exception& e ) {
      throw logic_error( "Couldn't load MRI" );
    }

    if ( msLabel == "" ) {
      SetLabel( mfnMRI );
    }

    InitializeFromMRI();
  }
}

MRI const*
VolumeCollection::GetMRI() const {

  // If we don't already have one, load it.
  if ( NULL == mMRI ) {

    // Make a mutable copy of ourselves so we can call the non-const
    // function LoadVolume.
    VolumeCollection* mutableThis = const_cast<VolumeCollection*>(this);
    mutableThis->LoadVolume();
  }

  return mMRI;
}

void
VolumeCollection::Save () const {

  Save( mfnMRI );
}

void
VolumeCollection::Save ( string const& ifn ) const {

  // Make a local copy of the filename and call MRIwrite.
  char* fn = strdup( ifn.c_str() );
  int rMRI = MRIwrite( mMRI, fn );
  free( fn );

  // If they gave us an error code (who knows what it could mean),
  // throw an error.
  if ( ERROR_NONE != rMRI ) {
    stringstream ssError;
    ssError << "Couldn't write file " << ifn;
    throw runtime_error( ssError.str() );
  }
}

void
VolumeCollection::InitializeFromMRI () {

  if ( NULL == mMRI ) {
    throw runtime_error( "InitializeFromMRI called without an MRI" );
  }

  // Get our ScannerRAS -> index transform.
  MATRIX* voxelFromScannerRAS = extract_r_to_i( mMRI );
  if ( NULL == voxelFromScannerRAS ) {
    throw runtime_error( "Couldn't get voxelFromScannerRAS matrix" );
  }

  if (getenv("FS_DEBUG_XFORMS"))
    {
      printf("voxelFromScannerRAS:\n") ;
      MatrixPrint(stdout, voxelFromScannerRAS) ;
    }
  // Copy it to a Matrix44 and release the MATRIX. Then set our
  // mDataToIndexTransform transform from this matrix. Then calculate
  // the WorldToIndex transform.
  Matrix44 m;
  m.SetMatrix( voxelFromScannerRAS );

  MatrixFree( &voxelFromScannerRAS );

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
  for ( tIDROI = mROIMap.begin();
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

  // Save number of frames.
  mcFrames = mMRI->nframes;

  // Look for time data.
  TryReadingTimeMetadata();
}

void
VolumeCollection::TryReadingTimeMetadata () {

  // We'll only have this data if we actually loaded from a bfloat or
  // bshort set of files, in which case our file name will be in the
  // format /path/to/data/stem_nnn.b{float,short} where stem is a
  // string and nnn is a three digit number. The header file will be
  // /path/to/data/stem.dat. So we'll extract everything up to the
  // last underscore, append a .dat, and look for that file.
  string::size_type nPos = mfnMRI.rfind( "_", string::npos );
  if ( nPos == string::npos ) {
    // Underscore not found in file name, oh well.
    mbInterpretFramesAsTimePoints = false;
    mbDataContainsErrorValues = false;
    return;
  }

  // This file has an underscore. Construct the metadata file name.
  string fnMetadata = mfnMRI.substr( 0, nPos ) + string(".dat");
  ifstream fMetadata( fnMetadata.c_str(), ios::in );
  if ( !fMetadata || fMetadata.bad() ) {
    // File not found, bail.
    mbInterpretFramesAsTimePoints = false;
    mbDataContainsErrorValues = false;
    return;
  }

  // Just go through the file and look for keywords. Parse out the
  // values when we find them.
  string sKeyword;
  int cConditions = -1;
  int cTimePoints = -1;
  float timeResolution;
  float preStimSeconds;
  CovarianceTable covTable;
  while ( !fMetadata.eof() ) {

    fMetadata >> sKeyword;

    if ( sKeyword == "TER" ) {

      fMetadata >> timeResolution;
      if ( fMetadata.fail() ) {
        cerr << "Error in time metadata file: TER is formatted badly" << endl;
        fMetadata.clear();
      }

    } else if ( sKeyword == "TPreStim" ) {

      fMetadata >> preStimSeconds;
      if ( fMetadata.fail() ) {
        cerr << "Error in time metadata file: TPreStim is "
        << "formatted badly" << endl;
        fMetadata.clear();
      }

    } else if ( sKeyword == "nCond" ) {

      fMetadata >> cConditions;
      if ( fMetadata.fail() ) {
        cerr << "Error in time metadata file: nCond is "
        << "formatted badly" << endl;
        fMetadata.clear();
      }

    } else if ( sKeyword == "Nh" ) {

      fMetadata >> cTimePoints;
      if ( fMetadata.fail() ) {
        cerr << "Error in time metadata file: NH is "
        << "formatted badly" << endl;
        fMetadata.clear();
      }

    } else if ( sKeyword == "hCovMtx" ) {

      // By now, the nCond and Nh should have already been stated, so
      // we have our cTimePoints and cConditions. If not, there's a
      // problem with the file.
      if ( cConditions == -1 || cTimePoints == -1 ) {
        throw runtime_error( "Error in time metadata file: hCovMtx "
                             "defined before Nh or nCond" );
      }

      int cRows = cTimePoints*(cConditions-1);
      int cCols = cTimePoints*(cConditions-1);
      for ( int nRow = 0; nRow < cRows; nRow++ ) {
        for ( int nCol = 0; nCol < cCols; nCol++ ) {
          fMetadata >> covTable[nRow][nCol];
          if ( fMetadata.fail() ) {
            cerr << "Error in time metadata file: hCovMtx is "
            << "formatted badly" << endl;
            fMetadata.clear();
          }
        }
      }
    }
  }

  // Save all the data we read.
  mbInterpretFramesAsTimePoints = true;
  mcConditions        = cConditions;
  mcTimePoints        = cTimePoints;
  mTimeResolution     = timeResolution;
  mcPreStimTimePoints = (int) floor(preStimSeconds * timeResolution);

  mbDataContainsErrorValues = true;
  mCovarianceTable          = covTable;
}


float 
VolumeCollection::GetMRIMinValue () const {
  return mMRIMinValue;
}

float 
VolumeCollection::GetMRIMaxValue () const {
  return mMRIMaxValue;
}

float
VolumeCollection::GetMRIMagnitudeMinValue () const {

  if ( NULL == mMagnitudeMRI ) {

    // Make a mutable copy of ourselves so we can call
    // MakeMagnitudeVolume.
    VolumeCollection* mutableThis = const_cast<VolumeCollection*>(this);
    mutableThis->MakeMagnitudeVolume();
  }

  return mMRIMagMinValue;
}

float
VolumeCollection::GetMRIMagnitudeMaxValue () const {

  if ( NULL == mMagnitudeMRI ) {

    // Make a mutable copy of ourselves so we can call
    // MakeMagnitudeVolume.
    VolumeCollection* mutableThis = const_cast<VolumeCollection*>(this);
    mutableThis->MakeMagnitudeVolume();
  }

  return mMRIMagMaxValue;
}

float 
VolumeCollection::GetVoxelXSize () const {
  return mVoxelSize[0];
}

float 
VolumeCollection::GetVoxelYSize () const {
    return mVoxelSize[1];
}

float 
VolumeCollection::GetVoxelZSize () const {
  return mVoxelSize[2];
}

int
VolumeCollection::GetDataType () const {

  if( NULL != mMRI ) {
    return mMRI->type;
  } 

  return -1;
}

int
VolumeCollection::GetNumberOfFrames () const {
  return mcFrames;
}

bool 
VolumeCollection::InterpretFramesAsTimePoints () const {
  return mbInterpretFramesAsTimePoints;
}

int 
VolumeCollection::GetNumberOfConditions () const {
  return mcConditions;
}

int 
VolumeCollection::GetNumberOfTimePoints () const {
  return mcTimePoints;
}

int
VolumeCollection::ConvertConditionAndTimePointToFrame ( int iCondition,
							int iTimePoint ) const {

  if ( !mbInterpretFramesAsTimePoints )
    throw runtime_error( "Volume doesn't use time points." );

  if ( iCondition < 0 || iCondition >= mcConditions )
    throw runtime_error( "Invalid condition." );
  if ( iTimePoint < 0 || iTimePoint >= mcTimePoints )
    throw runtime_error( "Invalid time point." );

  // If we have error data, we have to be cognizant of the way the
  // data is arranged.
  //
  //    WITHOUT ERROR DATA

  //    condition
  //      time point
  //        frame of displayable data

  //                        time
  //    frame   condition   point   type
  //      0         0         0     data
  //      1         0         1     data
  //      2         0         2     data
  //      3         1         0     data
  //      4         1         1     data
  //      5         1         2     data

  //    WITH ERROR DATA

  //    condition 0
  //      time point
  //        null values (ignored)
  //        sigma values
  //    condition 1..n
  //      time point
  //        frame of displayable data
  //        std dev (ignored)

  //                        time
  //    frame   condition   point   type
  //      0         0               null (ignored)
  //      1         0               null (ignored)
  //      2         0               null (ignored)
  //      3         0         0     sigma
  //      4         0         1     sigma
  //      5         0         2     sigma
  //      6         1         0     data
  //      7         1         1     data
  //      8         1         2     data
  //      9         1               std dev (ignored)
  //      10        1               std dev (ignored)
  //      11        1               std dev (ignored)
  //      12        2         0     data
  //      13        2         1     data
  //      14        2         2     data
  //      15        2               std dev (ignored)
  //      16        2               std dev (ignored)
  //      17        2               std dev (ignored)

  if ( mbDataContainsErrorValues ) {
    return (iCondition * 2 * mcTimePoints) + iTimePoint;
  } else {
    return (iCondition * mcTimePoints) + iTimePoint;
  }
}

int
VolumeCollection::ExtractConditionFromFrame ( int iFrame ) const {

  if ( !mbInterpretFramesAsTimePoints )
    throw runtime_error( "Volume doesn't use conditions." );

  if ( iFrame < 0 || iFrame >= mcFrames )
    throw runtime_error( "Invalid frame." );

  if ( mbDataContainsErrorValues ) {
    // Check the time point. This will throw an error if it's an
    // invalid frame. We do this because it's easier to check for an
    // invalid point with the mod calc rather than the div calc.
    ExtractTimePointFromFrame( iFrame );

    int condition = (int) floor( (float)iFrame / (float)(2*mcTimePoints) );
    return condition;

  } else {
    return (int) floor( (float)iFrame / (float)mcTimePoints );
  }
}

int
VolumeCollection::ExtractTimePointFromFrame ( int iFrame ) const {

  if ( !mbInterpretFramesAsTimePoints )
    throw runtime_error( "Volume doesn't use time points." );

  if ( iFrame < 0 || iFrame >= mcFrames )
    throw runtime_error( "Invalid frame." );

  if ( mbDataContainsErrorValues ) {
    int timePoint = iFrame % (2 * mcTimePoints);

    // If this happens, we're in the std dev part of the data, which
    // is not really valid.
    if ( timePoint > mcTimePoints )
      throw runtime_error( "Invalid frame; in error data" );

    return timePoint;
  } else {
    return iFrame % mcTimePoints;
  }
}

void
VolumeCollection::GetDataRASBounds ( float oRASBounds[6] ) const {

  if ( NULL != mMRI ) {

    // If our cache is dirty...
    if ( mbBoundsCacheDirty ) {

      // Recalc our bounds.
      mRASBounds[0] = mRASBounds[2] = mRASBounds[4] = 
	numeric_limits<float>::max();
      mRASBounds[1] = mRASBounds[3] = mRASBounds[5] = 
	numeric_limits<float>::min();

      int index[3];
      int cornerFactor[3];
      float RAS[3];
      for( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ ) {
        index[2] = cornerFactor[2] * mMRI->width;
        for( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ ) {
          index[1] = cornerFactor[1] * mMRI->width;
          for( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ ) {
            index[0] = cornerFactor[0] * mMRI->width;

            MRIIndexToRAS( index, RAS );

            if ( RAS[0] < mRASBounds[0] ) mRASBounds[0] = RAS[0];
            if ( RAS[0] > mRASBounds[1] ) mRASBounds[1] = RAS[0];
            if ( RAS[1] < mRASBounds[2] ) mRASBounds[2] = RAS[1];
            if ( RAS[1] > mRASBounds[3] ) mRASBounds[3] = RAS[1];
            if ( RAS[2] < mRASBounds[4] ) mRASBounds[4] = RAS[2];
            if ( RAS[2] > mRASBounds[5] ) mRASBounds[5] = RAS[2];
          }
        }
      }

      mbBoundsCacheDirty = false;
    }

    // Set the outgoing bounds.
    oRASBounds[0] = mRASBounds[0];
    oRASBounds[1] = mRASBounds[1];
    oRASBounds[2] = mRASBounds[2];
    oRASBounds[3] = mRASBounds[3];
    oRASBounds[4] = mRASBounds[4];
    oRASBounds[5] = mRASBounds[5];

  } else {

    // No MRI, not much we can do.
    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
      oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
  }
}

void
VolumeCollection::GetMRIIndexRange ( int oMRIIndexRange[3] ) const {

  if ( NULL != mMRI ) {
    oMRIIndexRange[0] = mMRI->width;
    oMRIIndexRange[1] = mMRI->height;
    oMRIIndexRange[2] = mMRI->depth;

  } else {
    
    oMRIIndexRange[0] = oMRIIndexRange[1] = oMRIIndexRange[2] = 0;
  }
}

bool
VolumeCollection::IsInBounds ( VolumeLocation const& iLoc ) const {

  if ( NULL != mMRI ) {
    return ( iLoc.mIdxi[0] >= 0 && iLoc.mIdxi[0] < mMRI->width &&
             iLoc.mIdxi[1] >= 0 && iLoc.mIdxi[1] < mMRI->height &&
             iLoc.mIdxi[2] >= 0 && iLoc.mIdxi[2] < mMRI->depth &&
             iLoc.mFrame >= 0 && iLoc.mFrame < mMRI->nframes );
  } else {
    return false;
  }
}

bool
VolumeCollection::IsMRIIdxInBounds ( int const iMRIIdx[3] ) const {

  if ( NULL != mMRI ) {
    return ( iMRIIdx[0] >= 0 && iMRIIdx[0] < mMRI->width &&
             iMRIIdx[1] >= 0 && iMRIIdx[1] < mMRI->height &&
             iMRIIdx[2] >= 0 && iMRIIdx[2] < mMRI->depth );
  } else {
    return false;
  }
}

bool
VolumeCollection::IsMRIIdxInBounds ( float const iMRIIdx[3] ) const {

  if ( NULL != mMRI ) {
    return ( iMRIIdx[0] >= 0 && iMRIIdx[0] < mMRI->width &&
             iMRIIdx[1] >= 0 && iMRIIdx[1] < mMRI->height &&
             iMRIIdx[2] >= 0 && iMRIIdx[2] < mMRI->depth );
  } else {
    return false;
  }
}

bool
VolumeCollection::IsFrameInBounds ( int const iFrame ) const {

  if ( NULL != mMRI ) {
    return ( iFrame >= 0 && iFrame < mMRI->nframes );
  } else {
    return false;
  }
}

void
VolumeCollection::RASToMRIIndex ( float const iRAS[3], int oIndex[3] ) const {


  mWorldToIndexTransform.MultiplyVector3( iRAS, oIndex );
}

void
VolumeCollection::RASToMRIIndex ( float const iRAS[3], float oIndex[3] ) const {

  mWorldToIndexTransform.MultiplyVector3( iRAS, oIndex );
  oIndex[0] += 0.5;
  oIndex[1] += 0.5;
  oIndex[2] += 0.5;
}

void
VolumeCollection::MRIIndexToRAS ( int const iIndex[3], float oRAS[3] ) const {

  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::MRIIndexToRAS ( float const iIndex[3], float oRAS[3] ) const {

  mWorldToIndexTransform.InvMultiplyVector3( iIndex, oRAS );
}

void
VolumeCollection::RASToDataRAS ( float const iRAS[3], float oDataRAS[3] ) const {

  mDataToWorldTransform->InvMultiplyVector3( iRAS, oDataRAS );
}

void
VolumeCollection::DataRASToRAS ( float const iDataRAS[3], float oRAS[3] ) const {

  mDataToWorldTransform->MultiplyVector3( iDataRAS, oRAS );
}

bool
VolumeCollection::IsTalTransformPresent() const {

  if( mMRI )
    return (mMRI->linear_transform != NULL);
  else
    return false;
}

void
VolumeCollection::RASToTal ( float const iRAS[3], float oTal[3] ) const {

  // Only do this if we have a transform.
  if( mMRI && mMRI->linear_transform ) {

    // First calculate the MNI Tal coords using the transform in the
    // volume.
    Real mniX, mniY, mniZ;
    transform_point( mMRI->linear_transform,
		     iRAS[0], iRAS[1], iRAS[2],
		     &mniX, &mniY, &mniZ );

    // Now convert to "real" Tal coords. We have two matrices to use:
    // one if the MNI Tal Z coord is <= 0, and one if it's > 0.
    float mni[3];
    mni[0] = mniX;
    mni[1] = mniY;
    mni[2] = mniZ;

    if( mniZ <= 0 )
      mMatrixMNITalToTalLtz.MultiplyVector3( mni, oTal );
    else
      mMatrixMNITalToTalGtz.MultiplyVector3( mni, oTal );
  }
}

void
VolumeCollection::TkRegRASToRAS ( float const iTkRegRAS[3],
                                  float oRAS[3] ) const {

  if ( NULL != mMRI ) {
    Real TkRegRAS[3], RAS[3];
    TkRegRAS[0] = iTkRegRAS[0];
    TkRegRAS[1] = iTkRegRAS[1];
    TkRegRAS[2] = iTkRegRAS[2];
    MRIsurfaceRASToRAS( mMRI, TkRegRAS[0], TkRegRAS[1], TkRegRAS[2],
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
VolumeCollection::RASToTkRegRAS ( float const iRAS[3],
                                  float oTkRegRAS[3] ) const {

  if ( NULL != mMRI ) {
    Real TkRegRAS[3], RAS[3];
    RAS[0] = iRAS[0];
    RAS[1] = iRAS[1];
    RAS[2] = iRAS[2];
    MRIRASToSurfaceRAS( mMRI, RAS[0], RAS[1], RAS[2],
                        &TkRegRAS[0], &TkRegRAS[1], &TkRegRAS[2] );
    oTkRegRAS[0] = TkRegRAS[0];
    oTkRegRAS[1] = TkRegRAS[1];
    oTkRegRAS[2] = TkRegRAS[2];
  } else {
    throw runtime_error( "Cannot transform from RAS to surface RAS because no"
                         "MRI is loaded." );
  }
}

float
VolumeCollection::GetMRINearestValue ( VolumeLocation const& iLoc ) const {

  Real value = 0;
  if ( NULL != mMRI ) {

    if ( IsInBounds( iLoc ) ) {

      switch ( mMRI->type ) {
      case MRI_UCHAR:
        value = (float)
                MRIseq_vox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2],iLoc.mFrame);
        break ;
      case MRI_SHORT:
        value =(float)
               MRISseq_vox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2],iLoc.mFrame);
        break ;
      case MRI_INT:
        value =(float)
               MRIIseq_vox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2],iLoc.mFrame);
        break ;
      case MRI_FLOAT:
        value =
          MRIFseq_vox(mMRI,iLoc.mIdxi[0],iLoc.mIdxi[1],iLoc.mIdxi[2],iLoc.mFrame);
        break ;
      default:
        value = 0;
      }
    }
  }

  return (float)value;

}

float
VolumeCollection::GetMRINearestValueAtIndexUnsafe ( int const iIndex[3],
						    int inFrame ) const {

  Real value = 0;

  switch ( mMRI->type ) {
  case MRI_UCHAR:
    value = (float)MRIseq_vox(mMRI, iIndex[0], iIndex[1], iIndex[2], inFrame );
    break ;
  case MRI_SHORT:
    value = (float)MRISseq_vox(mMRI, iIndex[0], iIndex[1], iIndex[2], inFrame);
    break ;
  case MRI_INT:
    value = (float)MRIIseq_vox(mMRI, iIndex[0], iIndex[1], iIndex[2], inFrame);
    break ;
  case MRI_FLOAT:
    value = MRIFseq_vox(mMRI, iIndex[0], iIndex[1], iIndex[2], inFrame );
    break ;
  default:
    value = 0;
  }

  return (float)value;
}

float
VolumeCollection::GetMRITrilinearValue ( VolumeLocation const& iLoc ) const {

  Real value = 0;
  if ( NULL != mMRI ) {

    if ( IsInBounds( iLoc ) ) {

      MRIsampleVolumeFrameType( mMRI,
                                iLoc.mIdxf[0], iLoc.mIdxf[1], iLoc.mIdxf[2],
                                iLoc.mFrame, SAMPLE_TRILINEAR, &value );
    }
  }
  return (float)value;
}

float
VolumeCollection::GetMRISincValue ( VolumeLocation const& iLoc ) const {

  Real value = 0;
  if ( NULL != mMRI ) {

    if ( IsInBounds( iLoc ) ) {

      /* Unfortunately, MRIsampleVolumeFrameType doesn't support sinc,
      so we just use the single frame sampler here. */
      MRIsampleVolumeType( mMRI, iLoc.mIdxf[0], iLoc.mIdxf[1], iLoc.mIdxf[2],
                           &value , SAMPLE_SINC );
    }
  }
  return (float)value;
}

void
VolumeCollection::SetMRIValue ( VolumeLocation const& iLoc,
                                float iValue ) {

  if ( NULL != mMRI ) {

    if ( IsInBounds( iLoc ) ) {

      switch ( mMRI->type ) {
      case MRI_UCHAR:
        MRIseq_vox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                    iLoc.mIdxi[2],iLoc.mFrame ) = (BUFTYPE) iValue;
        break ;
      case MRI_SHORT:
        MRISseq_vox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                     iLoc.mIdxi[2],iLoc.mFrame ) = (short) iValue;
        break ;
      case MRI_FLOAT:
        MRIFseq_vox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                     iLoc.mIdxi[2],iLoc.mFrame ) = (float) iValue;
        break ;
      case MRI_LONG:
        MRILseq_vox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                     iLoc.mIdxi[2],iLoc.mFrame ) = (long) iValue;
        break ;
      case MRI_INT:
        MRIIseq_vox( mMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                     iLoc.mIdxi[2],iLoc.mFrame ) = (int) iValue;
        break ;
      default:
        break ;
      }

      if ( iValue < mMRIMinValue ) {
        mMRIMinValue = iValue;
      }
      if ( iValue > mMRIMaxValue ) {
        mMRIMaxValue = iValue;
      }

      DataChanged();
    }
  }
}

float
VolumeCollection::GetMRIMagnitudeValue ( VolumeLocation const& iLoc ) const {

  Real value = 0;

  // If we don't have the magnitude volume, calculate it.
  if ( NULL == mMagnitudeMRI ) {

    // Make a mutable copy of ourselves so we can call
    // MakeMagnitudeVolume.
    VolumeCollection* mutableThis = const_cast<VolumeCollection*>(this);
    mutableThis->MakeMagnitudeVolume();
  }

  // Get the value.
  if ( NULL != mMagnitudeMRI ) {
    value = MRIFseq_vox( mMagnitudeMRI, iLoc.mIdxi[0], iLoc.mIdxi[1],
                         iLoc.mIdxi[2], iLoc.mFrame );
  }
  return (float)value;
}

TclCommandListener::TclCommandResult
VolumeCollection::DoListenToTclCommand ( char* isCommand,
    int iArgc, char** iasArgv ) {

  // SetVolumeCollectionFileName <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "SetVolumeCollectionFileName" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      string fnVolume = iasArgv[2];
      SetFileName( fnVolume );
    }
  }

  // LoadVolumeFromFileName <collectionID>
  if ( 0 == strcmp( isCommand, "LoadVolumeFromFileName" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      LoadVolume();
    }
  }

  // GetVolumeCollectionFileName <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeCollectionFileName" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      sReturnFormat = "s";
      sReturnValues = mfnMRI;
    }
  }

  // GetVolumeNumberOfFrames <colllectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeNumberOfFrames" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetNumberOfFrames();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";

      return ok;
    }
  }

  // WriteVolumeROIToLabel <collectionID> <roiID> <useRealRAS> <fileName>
  if ( 0 == strcmp( isCommand, "WriteVolumeROIToLabel" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      int roiID;
      try {
	roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
	sResult = string("bad roi ID: ") + e.what();
	return error;
      }

      bool bUseRealRAS = false;
      try {
        bUseRealRAS =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = "bad useRealRAS \"" + string(iasArgv[3]) + "\"," + e.what();
        return error;
      }

      try {
        WriteROIToLabel( roiID, bUseRealRAS, string(iasArgv[4]) );
      } catch (...) {
        sResult = "That ROI doesn't belong to this collection";
        return error;
      }
    }
  }

  // NewVolumeROIFromLabel <collectionID> <useRealRAS> <fileName>
  if ( 0 == strcmp( isCommand, "NewVolumeROIFromLabel" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      bool bUseRealRAS = false;
      try {
        bUseRealRAS =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = "bad useRealRAS \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }

      int roiID = NewROIFromLabel( bUseRealRAS, string(iasArgv[3]) );
      stringstream ssReturnValues;
      ssReturnValues << roiID;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // WriteVolumeROIsToSegmentation <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "WriteVolumeROIsToSegmentation" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      WriteROIsToSegmentation( string(iasArgv[2]) );
    }
  }

  // MakeVolumeUsingTemplate <collectionID> <templateID> <type>
  if ( 0 == strcmp( isCommand, "MakeVolumeUsingTemplate" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      int templateID;
      try {
	templateID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
	sResult = string("bad template ID: ") + e.what();
	return error;
      }

      int type;
      try {
	type = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
	sResult = string("bad type: ") + e.what();
	return error;
      }
    
      MakeUsingTemplate( templateID, type );
    }
  }

  // SaveVolume <collectionID>
  if ( 0 == strcmp( isCommand, "SaveVolume" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      Save();
    }
  }

  // SaveVolumeWithFileName <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "SaveVolumeWithFileName" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      string fn( iasArgv[2] );
      Save( fn );
    }
  }

  // SetUseVolumeDataToIndexTransform <collectionID> <use>
  if ( 0 == strcmp( isCommand, "SetUseVolumeDataToIndexTransform" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      try {
        bool bUse =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetUseWorldToIndexTransform( bUse );
      } catch ( runtime_error& e ) {
        sResult = "bad use \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }


  // GetUseVolumeDataToIndexTransform <layerID>
  if ( 0 == strcmp( isCommand, "GetUseVolumeDataToIndexTransform" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mbUseDataToIndexTransform );
      sReturnFormat = "i";
    }
  }

  // GetVolumeHistogram <volID> <roiID> <numBins>
  if( 0 == strcmp( isCommand, "GetVolumeHistogram" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      // Find an ROI if they want.
      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      ScubaROIVolume* volROI = NULL;
      if ( roiID != -1 ) {
        try {
          ScubaROI* roi = &ScubaROI::FindByID( roiID );
          volROI = dynamic_cast<ScubaROIVolume*>(roi);
        } catch (...) {
          throw runtime_error( "Couldn't find ROI or ROI "
                               "wasn't a volume ROI." );
        }
      }

      int cBins;
      try {
        cBins = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad numBins: ") + e.what();
        return error;
      }

      float minBinValue, binIncrement;
      map<int,int> binCounts;
      try {
	MakeHistogram( volROI, cBins, mMRIMinValue, mMRIMaxValue,
		       minBinValue, binIncrement, binCounts );
      } catch ( exception& e ) {
	throw runtime_error( string("Error making the histogram: ") + e.what() );
      }

      // Build the output.
      stringstream ssValues;
      stringstream ssFormat;

      ssFormat << "Lff";
      ssValues << minBinValue << " " << binIncrement << " ";

      ssFormat << "L";
      for ( int nBin = 0; nBin < (int)binCounts.size(); nBin++ ) {
        ssFormat << "i";
        ssValues << binCounts[nBin] << " ";
      }
      ssFormat << "ll";

      sReturnFormat = ssFormat.str();
      sReturnValues = ssValues.str();

    }
  }

  // SetVolumeAutosaveOn <collectionID> <on>
  if ( 0 == strcmp( isCommand, "SetVolumeAutosaveOn" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      try {
        bool bOn =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        mbAutosave = bOn;
      } catch ( runtime_error& e ) {
        sResult = "bad on \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // GetVolumeAutosaveOn <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeAutosaveOn" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mbAutosave );
      sReturnFormat = "i";
    }
  }

  // GetRASCoordsFromVolumeTkRegRAS <colllectionID> <x> <y> <z>
  if ( 0 == strcmp( isCommand, "GetRASCoordsFromVolumeTkRegRAS" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      float TkRegRAS[3];
      try {
	TkRegRAS[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
	TkRegRAS[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
	TkRegRAS[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );
      } catch ( runtime_error& e ) {
	sResult = string("bad x y z: ") + e.what();
	return error;
      }

      float ras[3];
      TkRegRASToRAS( TkRegRAS, ras );

      stringstream ssReturnValues;
      ssReturnValues << ras[0] << " " << ras[1] << " " << ras[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }

  // GetVolumeTkRegRASCoordsFromRAS <colllectionID> <x> <y> <z>
  if ( 0 == strcmp( isCommand, "GetVolumeTkRegRASCoordsFromRAS" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      float ras[3];
      try {
	ras[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
	ras[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
	ras[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );
      } catch ( runtime_error& e ) {
	sResult = string("bad x y z: ") + e.what();
	return error;
      }

      float TkRegRAS[3];
      RASToTkRegRAS( ras, TkRegRAS );

      stringstream ssReturnValues;
      ssReturnValues << TkRegRAS[0] << " " << TkRegRAS[1] << " " << TkRegRAS[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }

  // GetVolumeAverageValueInROI <colllectionID> <roiID>
  if ( 0 == strcmp( isCommand, "GetVolumeAverageValueInROI" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      map<int,ScubaROI*>::iterator tIDROI;
      tIDROI = mROIMap.find( roiID );
      if ( tIDROI != mROIMap.end() ) {
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
  if ( 0 == strcmp( isCommand, "GetVolumeStandardDeviationInROI" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      map<int,ScubaROI*>::iterator tIDROI;
      tIDROI = mROIMap.find( roiID );
      if ( tIDROI != mROIMap.end() ) {
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

  // GetVolumeInterpretFramesAsTime <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeInterpretFramesAsTime" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mbInterpretFramesAsTimePoints  );
      sReturnFormat = "i";
    }
  }

  // GetVolumeNumberOfConditions <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeNumberOfConditions" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << mcConditions;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // GetVolumeNumberOfTimePoints <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeNumberOfTimePoints" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << mcTimePoints;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // GetVolumeTimeResolution <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeTimeResolution" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << mTimeResolution;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "f";
    }
  }

  // GetVolumeNumberOfPreStimTimePoints <collectionID>
  if ( 0 == strcmp( isCommand, "GetVolumeNumberOfPreStimTimePoints" ) ) {
    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << mcPreStimTimePoints;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // DoVolumeValueRangeFill <volID> <sourceVolID> <roiID> <valueRanges>
  if ( 0 == strcmp( isCommand, "DoVolumeValueRangeFill" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collectionID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      int sourceVolID;
      try {
        sourceVolID = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad sourceVolID: ") + e.what();
        return error;
      }

      // Find a volume collection.
      VolumeCollection* sourceVol = NULL;
      try {
        DataCollection* col = &DataCollection::FindByID( sourceVolID );
        sourceVol = dynamic_cast<VolumeCollection*>(col);
      } catch (...) {
        throw runtime_error( "Couldn't find source volume or data collection "
                             "wasn't a volume." );
      }

      // Find an ROI if they want.
      int roiID;
      try {
        roiID = TclCommandManager::ConvertArgumentToInt( iasArgv[3] );
      } catch ( runtime_error& e ) {
        sResult = string("bad roiID: ") + e.what();
        return error;
      }

      ScubaROIVolume* volROI = NULL;
      if ( roiID != -1 ) {
        try {
          ScubaROI* roi = &ScubaROI::FindByID( roiID );
          volROI = (ScubaROIVolume*)roi;
        } catch (...) {
          throw runtime_error( "Couldn't find ROI or ROI "
			       "wasn't a volume ROI." );
        }
      }

      // Parse out the range elements.
      vector<ValueRangeFillElement> lElements;
      try {
	stringstream ssElements( iasArgv[4] );
	while( !ssElements.eof() ) {
	  float beginValue, endValue, newValue;
	  ssElements >> beginValue;
	  ssElements >> endValue;
	  ssElements >> newValue;
	  ValueRangeFillElement element( beginValue, endValue, newValue );
	  lElements.push_back( element );
	}
      } catch (...) {
        throw runtime_error( "Invalid value range list." );
      }

      DoValueRangeFill( *sourceVol, volROI, lElements );
    }
  }

  return DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
VolumeCollection::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "transformChanged" ) {
    CalcWorldToIndexTransform();
  }

  DataCollection::DoListenToMessage( isMessage, iData );
}

void
VolumeCollection::Select ( VolumeLocation const& iLoc ) {

  // If we have an ROI selected...
  if ( mSelectedROIID >= 0 ) {

    try {

      // Try to find the ROI.
      ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
      ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
      
      // Selectg in the ROI.
      volumeROI->SelectVoxel( iLoc.mIdxi );
      
      // Also mark this in the selection volume.
      mSelectedVoxels->
	Set_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2], true );
    } catch ( runtime_error& e ) {
      throw runtime_error ( "Couldn't select coxel, selected ROI "
			    "could not be found" );
    }
  }
}

void
VolumeCollection::Unselect ( VolumeLocation const& iLoc ) {

  if ( mSelectedROIID >= 0 ) {
    
    try {
      
      // Try to find the ROI.
      ScubaROI* roi = &ScubaROI::FindByID( mSelectedROIID );
      ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
      
      // Unselect in the ROI.
      volumeROI->UnselectVoxel( iLoc.mIdxi );
      
      // If there are no more ROIs with this voxel selected, unselect it
      // in the selection volume.
      bool bSelected = false;
      map<int,ScubaROI*>::const_iterator tIDROI;
      for ( tIDROI = mROIMap.begin();
	    tIDROI != mROIMap.end(); ++tIDROI ) {
	int roiID = (*tIDROI).first;

	// Get this ROI.
	ScubaROI* roi = &ScubaROI::FindByID( roiID );
	ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);

	// If it's selected, break.
	if ( volumeROI->IsVoxelSelected( iLoc.mIdxi ) ) {
	  bSelected = true;
	  break;
	}
      }
      if ( !bSelected ) {
	mSelectedVoxels->
	  Set_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2], false );
      }
      
    } catch ( runtime_error& e ) {
      throw runtime_error ( "Couldn't unselect voxel, selected ROI "
			    "could not be found" );
    }
  }
}

bool
VolumeCollection::IsSelected ( VolumeLocation const& iLoc, int oColor[3] ) const {

  // Check the selection volume cache first.
  if ( !(mSelectedVoxels->
	 Get_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] )) )
    return false;

  try {

    // Something is selected.
    bool bSelected = false;
    bool bFirstColor = true;

    // Go through our ROIs...
    map<int,ScubaROI*>::const_iterator tIDROI;
    for ( tIDROI = mROIMap.begin();
          tIDROI != mROIMap.end(); ++tIDROI ) {
      int roiID = (*tIDROI).first;

      // Get the ROI.
      ScubaROI* roi = &ScubaROI::FindByID( roiID );
      ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);

      // If this point is selected, get the color and blend it with
      // previous colors we've gotten.
      if ( volumeROI->IsVoxelSelected( iLoc.mIdxi ) ) {
        bSelected = true;
        int color[3];
        volumeROI->GetDrawColor( color );
        if ( bFirstColor ) {
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

  } catch ( runtime_error& e ) {
    cerr << "Error in IsSelected(): " << e.what() << endl;
    return false;
  } catch (...) {
    return false;
  }

}

bool
VolumeCollection::IsSelected ( VolumeLocation const& iLoc ) const {

  // Check the selection volume cache.
  return
    mSelectedVoxels->Get_Unsafe( iLoc.mIdxi[0], iLoc.mIdxi[1], iLoc.mIdxi[2] );
}

bool
VolumeCollection::IsOtherRASSelected ( float const iRAS[3],
				       int iThisROIID ) const {

  // Check the selection volume cache first.
  int index[3];
  RASToMRIIndex( iRAS, index );
  if ( !(mSelectedVoxels->Get_Unsafe( index[0], index[1], index[2] )) )
    return false;

  // Look at all our ROIs and if they aren't the one we are checking
  // against, check if it's selected. Break if we find one.
  bool bSelected = false;
  map<int,ScubaROI*>::const_iterator tIDROI;
  for ( tIDROI = mROIMap.begin();
        tIDROI != mROIMap.end(); ++tIDROI ) {
    int roiID = (*tIDROI).first;

    ScubaROI* roi = &ScubaROI::FindByID( roiID );
    if ( roiID == iThisROIID ) {
      continue;
    }

    ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
    if ( volumeROI->IsVoxelSelected( index ) ) {
      bSelected = true;
      break;
    }
  }

  return bSelected;
}

void
VolumeCollection::WriteROIToLabel ( int iROIID, bool ibUseRealRAS,
                                    string const& ifnLabel ) const {

  // Look for our ROI.
  map<int,ScubaROI*>::const_iterator tIDROI;
  tIDROI = mROIMap.find( iROIID );
  if ( tIDROI == mROIMap.end() ) 
    throw runtime_error( "ROI doesn't belong to this collection" );
  
  // Get our ROI.
  ScubaROIVolume const* roi = (ScubaROIVolume const*)(*tIDROI).second;
  
  // Get the bounds.
  int bounds[3];
  roi->GetROIBounds( bounds );

  // If there are selected voxels, bail.
  int cSelectedVoxels = roi->NumSelectedVoxels();
  if ( 0 == cSelectedVoxels ) {
    throw runtime_error( "No selected voxels." );
  }

  // Allocate a label for us.
  char* fnLabel = strdup( ifnLabel.c_str() );
  LABEL* label = LabelAlloc( cSelectedVoxels, NULL, fnLabel );
  if ( NULL == label ) {
    throw runtime_error( "Couldn't allocate label" );
  }

  // Set the number of points.
  label->n_points = cSelectedVoxels;
  
  // Go through the selected voxels...
  int nPoint = 0;
  list<Point3<int> > lSelected = roi->GetSelectedVoxelList();
  list<Point3<int> >::const_iterator tSelected;
  for ( tSelected = lSelected.begin();
	tSelected != lSelected.end(); ++tSelected ) {
    
    Point3<int> voxel = *tSelected;
    
    // Convert the voxel to RAS.
    float ras[3];
    MRIIndexToRAS( voxel.xyz(), ras );
    VolumeLocation loc( MakeVolumeLocationFromRAS( ras ) );
    
    // Convert these points to TkRegRAS space if needed.
    if ( !ibUseRealRAS ) {
      RASToTkRegRAS( ras, ras );
    }
    
    // Write the point to the label.
    label->lv[nPoint].x = ras[0];
    label->lv[nPoint].y = ras[1];
    label->lv[nPoint].z = ras[2];
    label->lv[nPoint].stat = GetMRINearestValue( loc );
    label->lv[nPoint].vno = -1;
    label->lv[nPoint].deleted = false;
    
    nPoint++;
  }
  
  // Write the label.
  int error = LabelWrite( label, fnLabel );
  if ( NO_ERROR != error ) {
    throw runtime_error( "Couldn't write label" );
  }
  
  free( fnLabel );
}

int
VolumeCollection::NewROIFromLabel ( bool ibUseRealRAS, 
				    string const& ifnLabel ) {

  // Read the label.
  char* fnLabel = strdup( ifnLabel.c_str() );
  LABEL* label = LabelRead( NULL, fnLabel );
  free( fnLabel );
  if ( NULL == label ) {
    throw runtime_error( "Couldn't read label" );
  }

  // Make a new ROI for us.
  ScubaROIVolume* volumeROI = NULL;
  try {
    int roiID = NewROI();
    ScubaROI* roi = &ScubaROI::FindByID( roiID );
    volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
  } catch (...) {
    throw runtime_error( "Couldn't make ROI" );
  }

  // For each point in the label...
  for ( int nPoint = 0; nPoint < label->n_points; nPoint++ ) {

    // Read the RAS point from the label.
    float ras[3];
    ras[0] = label->lv[nPoint].x;
    ras[1] = label->lv[nPoint].y;
    ras[2] = label->lv[nPoint].z;

    // Convert these points from TkRegRAS if needed.
    if ( !ibUseRealRAS )
      TkRegRASToRAS( ras, ras );

    // Convert to index.
    int index[3];
    RASToMRIIndex( ras, index );

    // Make sure it's in the bounds of the ROI (otherwise the label is
    // giving us points that are out of our own bounds).
    if ( !IsMRIIdxInBounds( index ) )
      continue;

    // Select the voxel.
    volumeROI->SelectVoxel( index );

    // Also mark this in the selection voxel in our cache.
    mSelectedVoxels->Set_Unsafe( index[0], index[1], index[2], true );
  }

  // Get rid of the label.
  LabelFree( &label );

  // Return thenew ID we got.
  return volumeROI->GetID();
}

void
VolumeCollection::WriteROIsToSegmentation ( string const& ifnVolume ) const {

  // Create a volume of the same size as our own.
  MRI* segVolume = MRIallocSequence( mMRI->width, mMRI->height, mMRI->depth,
                                     MRI_UCHAR, mMRI->nframes );
  if ( NULL == segVolume ) {
    throw runtime_error( "Couldn't create seg volume" );
  }

  // Copy the header from our own volume.
  MRIcopyHeader( mMRI, segVolume );

  // Go through the volume...
  int index[3];
  for ( index[2] = 0; index[2] < mMRI->depth; index[2]++ ) {
    for ( index[1] = 0; index[1] < mMRI->height; index[1]++ ) {
      for ( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {

        // For each of our ROIs, if one is selected here and if it's a
        // structure ROI, set the value of the seg volume to the
        // structure index.
        map<int,ScubaROI*>::const_iterator tIDROI;
        for ( tIDROI = mROIMap.begin();
              tIDROI != mROIMap.end(); ++tIDROI ) {
          int roiID = (*tIDROI).first;

          ScubaROI* roi = &ScubaROI::FindByID( roiID );
          ScubaROIVolume* volumeROI = dynamic_cast<ScubaROIVolume*>(roi);
          if ( volumeROI->GetType() == ScubaROI::Structure &&
               volumeROI->IsVoxelSelected( index ) )
            MRIvox( segVolume, index[0], index[1], index[2] ) =
              (BUFTYPE) volumeROI->GetStructure();
        }
      }
    }
  }

  // Write the volume.
  char* fnVolume = strdup( ifnVolume.c_str() );
  int error = MRIwrite( segVolume, fnVolume );
  free( fnVolume );

  // Free the volume.
  MRIfree( &segVolume );

  // Throw an error if we couldn't write.
  if ( NO_ERROR != error ) {
    throw runtime_error( "Couldn't write segmentation." );
  }
}


void
VolumeCollection::SetDataToWorldTransform ( int iTransformID ) {

  // Call the superclass function, then calc our composed transform.
  DataCollection::SetDataToWorldTransform( iTransformID );
  CalcWorldToIndexTransform();
}


Matrix44 const&
VolumeCollection::GetWorldToIndexTransform () const {

  return mWorldToIndexTransform.GetMainMatrix();
}


void
VolumeCollection::SetUseWorldToIndexTransform ( bool ibUse ) {

  // Don't change if it's the current setting.
  if ( ibUse == mbUseDataToIndexTransform )
    return;

  // Set and recalc our composed transform.
  mbUseDataToIndexTransform = ibUse;
  CalcWorldToIndexTransform();
}

bool 
VolumeCollection::GetUseWorldToIndexTransform () const {

  return mbUseDataToIndexTransform;
}

#define PRINTOUT 0

void
VolumeCollection::FindRASPointsInSquare ( float const iCenter[3],
					  float const iPointA[3],
					  float const iPointB[3],
					  float const iPointC[3],
					  float const iPointD[3],
					  float,
				       list<Point3<float> >& oPoints ) const {

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
  for ( int nZ = volumeBoundIdx[0].z(); nZ <= volumeBoundIdx[1].z(); nZ++ ) {
    for ( int nY = volumeBoundIdx[0].y(); nY <= volumeBoundIdx[1].y(); nY++ ) {
      for ( int nX = volumeBoundIdx[0].x(); nX <= volumeBoundIdx[1].x(); nX++ ) {

        Point3<int> curMRIIndex( nX, nY, nZ );
        VectorOps::IntersectionResult rInt =
          VoxelIntersectsPlane( curMRIIndex, planeIdx, planeNIdx,
                                intersectionIdx );

        if ( VectorOps::intersect == rInt ) {
#if PRINTOUT
          cerr << "\t\thit";
#endif

          // Calculate the anglesum of the intersection point with the
          // four plane corner points. If it is 2pi, this point is
          // inside the poly.
          double angleSum = 0;
          for ( int nVector = 0; nVector < 4; nVector++ ) {
            Point3<float> v1 = squareIdxf[nVector]       - intersectionIdx;
            Point3<float> v2 = squareIdxf[(nVector+1)%4] - intersectionIdx;

            float v1Length = VectorOps::Length(v1);
            float v2Length = VectorOps::Length(v2);

            if ( fabs(v1Length * v2Length) <= (float)0.0001 ) {
              angleSum = 2*M_PI;
              break;
            }

            double rads = VectorOps::RadsBetweenVectors( v1, v2 );
            angleSum += rads;
          }

          if ( fabs(angleSum - 2.0*M_PI) <= (float)0.0001 ) {

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
VolumeCollection::FindRASPointsInCircle ( float const iPointA[3], 
					  float const iPointB[3],
					  float const iPointC[3],
					  float const iPointD[3],
					  float iMaxDistance,
					  float const iCenter[3],
					  float iRadius,
				      list<Point3<float> >& oPoints ) const {

  // Get a list of RAS voxels in the square.
  list<Point3<float> > squarePoints;
  FindRASPointsInSquare( iCenter, iPointA, iPointB, iPointC, iPointD,
                         iMaxDistance, squarePoints );

  // For each one of those, check if it's within the circle.
  Point3<float> center( iCenter );
  list<Point3<float> >::iterator tPoints;
  for ( tPoints = squarePoints.begin(); tPoints != squarePoints.end();
        ++tPoints ) {

    Point3<float> pointRAS = *tPoints;

    // If it's within the radius, add the original corner RAS point.
    if ( VectorOps::Distance( pointRAS, center ) < iRadius ) {
      oPoints.push_back( pointRAS );
    }
  }
}

void
VolumeCollection::FindRASPointsOnSegment ( float const iPointA[3],
					   float const iPointB[3],
			       std::list<Point3<float> >& oPoints ) const {

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
  for ( int nZ = volumeBoundIdx[0].z(); nZ <= volumeBoundIdx[1].z(); nZ++ ) {
    for ( int nY = volumeBoundIdx[0].y(); nY <= volumeBoundIdx[1].y(); nY++ ) {
      for ( int nX = volumeBoundIdx[0].x(); nX <= volumeBoundIdx[1].x(); nX++ ) {

        // If it intersects the segment...
        Point3<int> curMRIIndex( nX, nY, nZ );
        VectorOps::IntersectionResult rInt =
          VoxelIntersectsSegment( curMRIIndex, segIdxf[0], segIdxf[1],
                                  intersectionIdx );

        // Add it to the list to return.
        if ( VectorOps::intersect == rInt ) {

          Point3<float> curRAS;
          MRIIndexToRAS( curMRIIndex.xyz(), curRAS.xyz() );
          oPoints.push_back( curRAS );
        }
      }
    }
  }
}

void
VolumeCollection::ImportControlPoints ( string const& ifnControlPoints,
                              list<Point3<float> >& oControlPoints ) const {

  // Read a control points file.
  int cControlPoints;
  int bUseRealRAS;
  char* fnControlPoints = strdup( ifnControlPoints.c_str() );
  MPoint* aControlPoints =
    MRIreadControlPoints( fnControlPoints, &cControlPoints, &bUseRealRAS );
  free( fnControlPoints );
  if ( NULL == aControlPoints ) {
    throw runtime_error( "Couldn't read control points file." );
  }

  // For each one we got...
  for ( int nControlPoint = 0; 
	nControlPoint < cControlPoints; nControlPoint++) {

    Point3<float> newControlPoint;
    if ( !bUseRealRAS ) {

      // We're getting surface RAS points. Need to convert to normal
      // RAS points.
      Real TkRegRAS[3];
      Real ras[3];
      TkRegRAS[0] = aControlPoints[nControlPoint].x;
      TkRegRAS[1] = aControlPoints[nControlPoint].y;
      TkRegRAS[2] = aControlPoints[nControlPoint].z;
      MRIsurfaceRASToRAS( mMRI, TkRegRAS[0], TkRegRAS[1], TkRegRAS[2],
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

    // Push the point onto the list of control points.
    oControlPoints.push_back( newControlPoint );
  }

  free( aControlPoints );
}

void
VolumeCollection::ExportControlPoints ( string const& ifnControlPoints,
                           list<Point3<float> > const& iControlPoints) const {

  // Allocate some control points.
  int cControlPoints = iControlPoints.size();
  MPoint* aControlPoints = (MPoint*) calloc( sizeof(MPoint), cControlPoints );
  if ( NULL == aControlPoints ) {
    throw runtime_error( "Couldn't allocate control point storage." );
  }

  // Fill them out with the points we got.
  int nControlPoint = 0;
  list<Point3<float> >::const_iterator tControlPoint;
  for ( tControlPoint = iControlPoints.begin();
        tControlPoint != iControlPoints.end();
        ++tControlPoint ) {

    aControlPoints[nControlPoint].x = (*tControlPoint)[0];
    aControlPoints[nControlPoint].y = (*tControlPoint)[1];
    aControlPoints[nControlPoint].z = (*tControlPoint)[2];
    nControlPoint++;
  }

  // WRite the file.
  char* fnControlPoints = strdup( ifnControlPoints.c_str() );
  int rWrite = MRIwriteControlPoints( aControlPoints, cControlPoints,
                                      true, fnControlPoints );
  free( fnControlPoints );

  // Throw an error if the write failed.
  if ( NO_ERROR != rWrite ) {
    throw runtime_error( "Couldn't write control point file." );
  }
}

void
VolumeCollection::MakeHistogram ( list<Point3<float> > const& iRASPoints,
				  ScubaROIVolume const* iROI, 
                                  int icBins,
                                  float& oMinBinValue, float& oBinIncrement,
                                  map<int,int>& oBinCounts ) const {

  // If no points, return.
  if ( iRASPoints.size() == 0 ) {
    return;
  }

  // Find values for all the RAS points. Save the highs and lows.
  float low = numeric_limits<float>::max();
  float high = numeric_limits<float>::min();
  list<Point3<float> >::const_iterator tPoint;
  list<float> values;
  for ( tPoint = iRASPoints.begin(); tPoint != iRASPoints.end(); ++tPoint ) {
    Point3<float> point = *tPoint;
    VolumeLocation loc( MakeVolumeLocationFromRAS( point.xyz()) );
    if( iROI && !iROI->IsVoxelSelected( loc.Index() ) ) {
      continue;
    }
    float value = GetMRINearestValue( loc );
    if ( value < low ) {
      low = value;
    }
    if ( value > high ) {
      high = value;
    }
    values.push_back( value );
  }

  // No points in the ROI.
  if( values.size() == 0 ) {
    oMinBinValue = 0;
    oBinIncrement = 0;
    return;
  }

  // Set all bins to zero.
  for ( int nBin = 0; nBin < icBins; nBin++ ) {
    oBinCounts[nBin] = 0;
  }

  // Increment our high a little so that we can fit the high value
  // into the last bin nicely. Use our increment for this.
  high += GetPreferredValueIncrement();

  float binIncrement = (high - low) / (float)icBins;
  
  // Now go through the values and calculate a bin for them. Increment
  // the count in the proper bin.
  list<float>::iterator tValue;
  for ( tValue = values.begin(); tValue != values.end(); ++tValue ) {
    float value = *tValue;
    int nBin = (int)floor( (value - low) / (float)binIncrement );
    assert( nBin >= 0 );
    assert( nBin < icBins );
    oBinCounts[nBin]++;
  }

  oMinBinValue = low;
  oBinIncrement = binIncrement;
}

void
VolumeCollection::MakeHistogram ( ScubaROIVolume const* iROI, 
				  int icBins,
                                  float iMinThresh, float iMaxThresh,
                                  float& oMinBinValue, float& oBinIncrement,
                                  map<int,int>& oBinCounts ) const {

  if ( iMinThresh < mMRIMinValue || iMinThresh > mMRIMaxValue ||
       iMaxThresh < mMRIMinValue || iMaxThresh > mMRIMaxValue ) {
    throw runtime_error( "Thresh values out of bounds." );
  }

  if ( iMinThresh > iMaxThresh ) {
    throw runtime_error( "Min thresh is greater than max thresh." );
  }

  if ( icBins <= 0 ) {
    throw runtime_error( "Number of bins is less than zero." );
  }

  if( iROI && iROI->NumSelectedVoxels() == 0 ) {
    throw runtime_error( "Empty ROI." );
  }

  // If we have an ROI, we have to find the min/max in it. Otherwise
  // we can just use our whole range min/max.
  float low = numeric_limits<float>::max();
  float high = numeric_limits<float>::min();
  if( iROI ) {

    int index[3] = { 0,0,0 };
    for ( index[2] = 0; index[2] < mMRI->depth; index[2]++ )
      for ( index[1] = 0; index[1] < mMRI->height; index[1]++ )
	for ( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {
	  if( iROI && !iROI->IsVoxelSelected( index ) )
	    continue;
	  float value = GetMRINearestValueAtIndexUnsafe( index, 0 );
	  if ( value < low )
	    low = value;
	  if ( value > high )
	    high = value;
	}
    
  } else {
    low = mMRIMinValue;
    high = mMRIMaxValue;
  }

  float binIncrement = (high - low) / (float)icBins;

  // Set all bins to zero.
  for ( int nBin = 0; nBin < icBins; nBin++ ) {
    oBinCounts[nBin] = 0;
  }

  int index[3] = { 0,0,0 };
  for ( index[2] = 0; index[2] < mMRI->depth; index[2]++ ) {
    for ( index[1] = 0; index[1] < mMRI->height; index[1]++ ) {
      for ( index[0] = 0; index[0] < mMRI->width; index[0]++ ) {

	if( iROI && !iROI->IsVoxelSelected( index ) )
	    continue;

        float value = GetMRINearestValueAtIndexUnsafe( index, 0 );
        if ( value < iMinThresh || value > iMaxThresh )
          continue;

        int nBin = (int)floor( (value - mMRIMinValue) / (float)binIncrement );
        if ( nBin == icBins && value == iMaxThresh ) {
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
VolumeCollection::SetAutoSaveOn ( bool ibAutosave ) {
  mbAutosave = ibAutosave;
}

bool 
VolumeCollection::GetAutoSaveOn () const {
  return mbAutosave;
}

bool 
VolumeCollection::IsAutosaveDirty () const {
  return mbAutosaveDirty;
}

string
VolumeCollection::MakeAutosaveFileName ( string const& ifn ) const {

  // Generate an autosave name.
  string fnAutosave = ifn;
  if ( fnAutosave == "" ) {
    fnAutosave = "newvolume";
  }
  string::size_type c = fnAutosave.find( '/' );
  while ( c != string::npos ) {
    fnAutosave[c] = '.';
    c = fnAutosave.find( '/' );
  }
  fnAutosave = "/tmp/" + fnAutosave + ".mgz";

  return fnAutosave;
}

void
VolumeCollection::DeleteAutosave () {

  // Copy our autosave file name.
  char* fnAutosave = strdup( mfnAutosave.c_str() );

  // If it exists, delete the file.
  struct stat info;
  int rStat = stat( fnAutosave, &info );
  if ( 0 == rStat ) {
    if ( S_ISREG(info.st_mode) ) {
      remove( fnAutosave );
    }
  }

  free( fnAutosave );
}

void
VolumeCollection::AutosaveIfDirty () const {

  // If we're dirty and on, save a copy of ourselves in the autosave
  // location.
  if ( mbAutosave && mbAutosaveDirty ) {
    Save( mfnAutosave );
    mbAutosaveDirty = false;
  }
}

VectorOps::IntersectionResult
VolumeCollection::VoxelIntersectsPlane ( Point3<int> const& iMRIIndex,
					 Point3<float> const& iPlaneIdx,
					 Point3<float> const& iPlaneIdxNormal,
				   Point3<float>& oIntersectionIdx ) const {


  // Call the big function and just return the first result if
  // present.
  int cIntersections;
  VectorOps::IntersectionResult results[12];
  Point3<float> intersections[12];
  if ( VoxelIntersectsPlane( iMRIIndex, iPlaneIdx, iPlaneIdxNormal,
                             cIntersections, results, intersections ) ) {

    if ( cIntersections > 0 ) {
      oIntersectionIdx = intersections[0];
      return results[0];
    }
  }

  return VectorOps::dontIntersect;
}

bool
VolumeCollection::VoxelIntersectsPlane ( Point3<int> const& iMRIIndex,
					 Point3<float> const& iPlaneIdx,
					 Point3<float> const& iPlaneIdxNormal,
					 int& ocIntersections,
			    VectorOps::IntersectionResult orIntersection[12],
			       Point3<float> oIntersectionIdx[12] ) const {
  
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
  for ( int nSegment = 0; nSegment < 12; nSegment++ ) {
    segmentIdx[nSegment][0].Set( voxelIdx[anSegments[nSegment][0]] );
    segmentIdx[nSegment][1].Set( voxelIdx[anSegments[nSegment][1]] );
  }

  // Intersect these segments with the plane. If any of them hit...
  int cIntersections = 0;
  for ( int nSegment = 0; nSegment < 12; nSegment++ ) {

    Point3<float> intersectionIdx;
    VectorOps::IntersectionResult rInt =
      VectorOps::SegmentIntersectsPlane
      ( segmentIdx[nSegment][0], segmentIdx[nSegment][1],
        iPlaneIdx, iPlaneIdxNormal, intersectionIdx );

#if PRINTOUT
    cerr << "\tTesting " << segmentIdx[nSegment][0] << " to "
    <<  segmentIdx[nSegment][1] << "..." << endl;
#endif

    // If an intersection, add the results to the output arrays.
    if ( VectorOps::intersect == rInt ) {

      orIntersection[cIntersections] = rInt;
      oIntersectionIdx[cIntersections] = intersectionIdx;
      cIntersections++;
    }
  }

  // Return the number of intersections.
  ocIntersections = cIntersections;

  // Return true if there was an intersection.
  return (cIntersections > 0);
}

VectorOps::IntersectionResult
VolumeCollection::VoxelIntersectsSegment( Point3<int> const& iMRIIndex,
					  Point3<float> const& iSegIdxA,
					  Point3<float> const& iSegIdxB,
				    Point3<float>& oIntersectionIdx ) const {
  
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

  for ( int nFace = 0; nFace < 6; nFace++ ) {
    planeIdx[nFace][0].Set( voxelIdx[anFaces[nFace][0]] );
    planeIdx[nFace][1].Set( voxelIdx[anFaces[nFace][1]] );
    planeIdx[nFace][2].Set( voxelIdx[anFaces[nFace][2]] );
    planeIdx[nFace][3].Set( voxelIdx[anFaces[nFace][3]] );
  }

  Point3<float> planeNIdx[6];
  planeNIdx[0].Set( 0, 0, 1 );
  planeNIdx[3].Set( 0, 0, 1 );
  planeNIdx[1].Set( 1, 0, 0 );
  planeNIdx[4].Set( 1, 0, 0 );
  planeNIdx[2].Set( 0, 1, 0 );
  planeNIdx[5].Set( 0, 1, 0 );

  for ( int nFace = 0; nFace < 6; nFace++ ) {

#if PRINTOUT
    cerr << "\tTesting face p " << planeIdx[nFace][0]
    << " N " << planeNIdx[nFace] << "... ";
#endif

    Point3<float> intersectionIdx;
    VectorOps::IntersectionResult rInt =
      VectorOps::SegmentIntersectsPlane
      ( iSegIdxA, iSegIdxB, planeIdx[nFace][0], planeNIdx[nFace],
        intersectionIdx );

    if ( VectorOps::intersect == rInt ) {

#if PRINTOUT
      cerr << " hit " << endl
      << "\t\tIntersection " << intersectionIdx << endl;
#endif

      // Calculate the anglesum of the intersection point with the
      // four plane corner points. If it is 2pi, this point is
      // inside the poly.
      double angleSum = 0;
      for ( int nCorner = 0; nCorner < 4; nCorner++ ) {
        Point3<float> v1 = planeIdx[nFace][nCorner]       - intersectionIdx;
        Point3<float> v2 = planeIdx[nFace][(nCorner+1)%4] - intersectionIdx;
#if PRINTOUT
        cerr << "\t\tCorner " << planeIdx[nFace][nCorner] << endl;
#endif

        float v1Length = VectorOps::Length(v1);
        float v2Length = VectorOps::Length(v2);

        if ( fabs(v1Length * v2Length) <= (float)0.0001 ) {
          angleSum = 2*M_PI;
          break;
        }

        double rads = VectorOps::RadsBetweenVectors( v1, v2 );
        angleSum += rads;
      }

      if ( fabs(angleSum - 2.0*M_PI) <= (float)0.0001 ) {
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

void
VolumeCollection::GetVoxelsWithValue ( float iValue,
                                  list<VolumeLocation>& olLocations ) const {

  // Just go through the volume and look for this value, and insert
  // locations into the outgoing list.
  for ( int nZ = 0; nZ < mMRI->depth; nZ++ ) {
    for ( int nY = 0; nY < mMRI->height; nY++ ) {
      for ( int nX = 0; nX < mMRI->width; nX++ ) {

	int index[3] = { nX, nY, nZ };
        float value = GetMRINearestValueAtIndexUnsafe( index, 0 );
        if ( fabs( value - iValue ) < numeric_limits<float>::epsilon() ) {
          VolumeLocation loc( MakeVolumeLocationFromIndex( index ) );
          olLocations.push_back( loc );
        }
      }
    }
  }
}

float
VolumeCollection::GetRASVolumeOfNVoxels ( int icVoxels ) const {

  return GetVoxelXSize() * GetVoxelYSize() * GetVoxelZSize() * (float)icVoxels;
}

float
VolumeCollection::GetPreferredValueIncrement () const {

  // Look at the range. If it's > 100, inc is 1, 10-100, inc is .1,
  // 1-10, inc is .01, etc.
  float range = mMRIMaxValue - mMRIMinValue;
  float inc = 1;
  if ( range >= 1000000 )        {
    inc = 1000;
  } else if ( range >=  100000 )        {
    inc =  100;
  } else if ( range >=   10000 )        {
    inc =   10;
  } else if ( range >=    1000 )        {
    inc =    1;
  } else if ( range >=      10 )        {
    inc =    0.1;
  } else if ( range >=       1 )        {
    inc =    0.01;
  } else if ( range >=       0.1 )      {
    inc =    0.001;
  } else if ( range >=       0.01 )     {
    inc =    0.0001;
  } else if ( range >=       0.001 )    {
    inc =    0.00001;
  } else if ( range >=       0.0001 )   {
    inc =    0.000001;
  } else if ( range >=       0.00001 )  {
    inc =    0.0000001;
  } else if ( range >=       0.000001 ) {
    inc =    0.00000001;
  } else                          {
    inc =    0.000000001;
  }

  return inc;
}

float
VolumeCollection::GetAverageValue ( list<VolumeLocation> const& ilLocations ) const {

  if ( ilLocations.size() < 1 ) {
    throw runtime_error( "No voxels" );
  }

  float sumValue = 0;
  list<VolumeLocation>::const_iterator tLocation;
  for ( tLocation = ilLocations.begin(); tLocation != ilLocations.end();
        ++tLocation ) {

    VolumeLocation const& loc = *tLocation;
    sumValue += GetMRINearestValue( loc );
  }

  return sumValue / (float)ilLocations.size();
}

float
VolumeCollection::GetAverageValue ( ScubaROIVolume const& iROI ) const {

  list<VolumeLocation> lLocs;

  list<Point3<int> > lSelected = iROI.GetSelectedVoxelList();
  list<Point3<int> >::const_iterator tSelected;
  for ( tSelected = lSelected.begin();
        tSelected != lSelected.end(); ++tSelected ) {

    Point3<int> voxel = *tSelected;
    VolumeLocation loc( MakeVolumeLocationFromIndex( voxel.xyz() ) );
    lLocs.push_back( loc );
  }

  return GetAverageValue( lLocs );
}

float
VolumeCollection::GetStandardDeviation ( list<VolumeLocation> const& ilLocations,
					 float iMean ) const {

  if ( ilLocations.size() < 1 ) {
    throw runtime_error( "No voxels" );
  }

  float sumDeviation = 0;
  list<VolumeLocation>::const_iterator tLocation;
  for ( tLocation = ilLocations.begin(); tLocation != ilLocations.end();
        ++tLocation ) {

    VolumeLocation const& loc = *tLocation;
    float value = GetMRINearestValue( loc );
    float deviation = (value - iMean) * (value - iMean);
    sumDeviation += deviation;
  }

  return sqrt( sumDeviation / (float)ilLocations.size() );
}

float
VolumeCollection::GetStandardDeviation ( ScubaROIVolume const& iROI, 
					 float iMean ) const {

  list<VolumeLocation> lLocs;
  list<Point3<int> > lSelected = iROI.GetSelectedVoxelList();
  list<Point3<int> >::iterator tSelected;
  for ( tSelected = lSelected.begin();
        tSelected != lSelected.end(); ++tSelected ) {

    Point3<int> voxel = *tSelected;
    VolumeLocation loc( MakeVolumeLocationFromIndex( voxel.xyz() ) );
    lLocs.push_back( loc );
  }

  return GetStandardDeviation( lLocs, iMean );
}

void
VolumeCollection::PrintVoxelCornerCoords ( ostream& iStream,
				       Point3<int> const& iMRIIndex ) const {

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
  for ( int nCorner = 0; nCorner < 8; nCorner++ ) {
    MRIIndexToRAS( voxelIdx[nCorner].xyz(), voxelRAS[nCorner].xyz() );
    iStream << nCorner << ": " << voxelRAS[nCorner] << endl;
  }
}

void
VolumeCollection::DoValueRangeFill ( VolumeCollection const& iSourceVol,
				     ScubaROIVolume const* iROI,
				     std::vector<ValueRangeFillElement> const& lElements) {

  // Do the same walk through of voxels that we do in the get histogram.
  int range[3]; range[0]=0; range[1]=0; range[2]=0;
  GetMRIIndexRange( range );

  Point3<int> index( 0, 0, 0 );
  VolumeLocation destLoc( MakeVolumeLocationFromIndex( index.xyz() ) );
  VolumeLocation srcLoc( iSourceVol.MakeVolumeLocationFromIndex( index.xyz() ) );
  //Point3<float> RAS;

  for( index[2] = 0; index[2] < range[2]; index[2]++ )
    for( index[1] = 0; index[1] < range[1]; index[1]++ )
      for( index[0] = 0; index[0] < range[0]; index[0]++ ) {
	
	destLoc.SetFromIndex( index.xyz() );

	// If they gave us an ROI to use, make sure this is in the
        // ROI.
        if ( NULL != iROI &&
	     !iROI->IsVoxelSelected( index.xyz() ) )
	  continue;

	if( !iSourceVol.IsInBounds( destLoc ) )
	    continue;

	// Get the corresponding value from the source vol.
	srcLoc.SetFromRAS( destLoc.RAS() );
        float value = iSourceVol.GetMRINearestValue( srcLoc );
	
        // We need to see if we should fill this voxel. Go through our
        // list and see if it falls into a range. If so, edit the
        // value in the dest.
        vector<ValueRangeFillElement>::const_iterator tElement;
        for ( tElement = lElements.begin();
	      tElement != lElements.end(); ++tElement ) {
          if ( value > tElement->mBegin && value < tElement->mEnd ) {
            SetMRIValue( destLoc, tElement->mValue );
            break;
          }
        }
      }

}

void
VolumeCollection::UpdateMRIValueRange () {

  if ( NULL != mMRI ) {

    // Get the range in every frame and find the most min and max.
    mMRIMinValue = 99999;
    mMRIMaxValue = -99999;
    float minValue, maxValue;
    for ( int nFrame = 0; nFrame < mMRI->nframes; nFrame++ ) {
      MRIvalRangeFrame( mMRI, &minValue, &maxValue, nFrame );
      if ( minValue < mMRIMinValue )
        mMRIMinValue = minValue;
      if ( maxValue > mMRIMaxValue )
        mMRIMaxValue = maxValue;
    }
  }
}

void
VolumeCollection::CalcWorldToIndexTransform () {

  if ( mbUseDataToIndexTransform ) {

    // Just mult our transforms together.
    Transform44 worldToData = mDataToWorldTransform->Inverse();
    mWorldToIndexTransform = mDataToIndexTransform * worldToData;

  } else {

    // If they don't want to use a transform, just use an identity
    // with enough translation to get us into the center of the world.
    Transform44 center;
    center.SetMainTransform( 1, 0, 0, (double)mMRI->width/2.0,
                             0, 1, 0, (double)mMRI->height/2.0,
                             0, 0, 1, (double)mMRI->depth/2.0,
                             0, 0, 0, 1 );

    Transform44 worldToData = mDataToWorldTransform->Inverse();
    mWorldToIndexTransform = center * worldToData;
  }

  // Notify of changes.
  DataChanged();
}

void
VolumeCollection::MakeMagnitudeVolume () {

  if ( NULL != mMRI ) {
    mMagnitudeMRI =
      MRIallocSequence( mMRI->width, mMRI->height, mMRI->depth,
                        MRI_FLOAT, mMRI->nframes );
    MRI* gradMRI = MRIsobel( mMRI, NULL, mMagnitudeMRI );
    MRIfree( &gradMRI );

    MRIvalRange( mMagnitudeMRI, &mMRIMagMinValue, &mMRIMagMaxValue );
  }
}

ScubaROI*
VolumeCollection::DoNewROI () {

  // Create our specific ROI.
  ScubaROIVolume* roi = new ScubaROIVolume();

  // Set its bounds.
  if ( NULL != mMRI ) {
    int bounds[3];
    bounds[0] = mMRI->width;
    bounds[1] = mMRI->height;
    bounds[2] = mMRI->depth;

    roi->SetROIBounds( bounds );
  }

  // Give it a default color.
  int color[] = { 0, 0, 255 };
  roi->SetFreeColor( color );

  return roi;
}


void
VolumeCollection::InitSelectionVolume () {

  if ( NULL != mMRI ) {
    
    delete mSelectedVoxels;

    // Create our selection volume with the same dimensions as our
    // volume.
    mSelectedVoxels =
      new Volume3<bool>( mMRI->width, mMRI->height, mMRI->depth, false );
  }
}

void
VolumeCollection::DataChanged () {

  mbAutosaveDirty = true;
  DataCollection::DataChanged();
}

VolumeCollectionFlooder::VolumeCollectionFlooder () {
  mVolume = NULL;
  mParams = NULL;
}

VolumeCollectionFlooder::~VolumeCollectionFlooder () {}

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
VolumeCollectionFlooder::DoBegin () {}

void
VolumeCollectionFlooder::DoEnd () {}

bool
VolumeCollectionFlooder::DoStopRequested () {
  return false;
}

void
VolumeCollectionFlooder::DoVoxel ( float[3], int iFrame ) {}

bool
VolumeCollectionFlooder::CompareVoxel ( float[3], int iFrame ) {
  return true;
}

void
VolumeCollectionFlooder::Flood ( VolumeCollection& iVolume,
                                 float iRASSeed[3], int iFrame,
                                 Params& iParams ) {

  mVolume = &iVolume;
  mParams = &iParams;

  // Get the seed voxel.
  Point3<int> seedVoxel;
  iVolume.RASToMRIIndex( iRASSeed, seedVoxel.xyz() );

  // Get the source volume. See if it's different from the one we get
  // as input.
  VolumeCollection* sourceVol = NULL;
  try {
    DataCollection* col =
      &DataCollection::FindByID( iParams.mSourceCollection );
    //    VolumeCollection* vol = dynamic_cast<VolumeCollection*>(col);
    sourceVol = (VolumeCollection*)col;
  } catch (...) {
    throw runtime_error( "Couldn't find source volume." );
  }
  bool bDifferentSource = sourceVol->GetID() != iVolume.GetID();

  // Init a visited volume.
  Volume3<bool>* bVisited =
    new Volume3<bool>( iVolume.mMRI->width,
                       iVolume.mMRI->height,
                       iVolume.mMRI->depth, false );

  // Save the initial value.
  VolumeLocation seedLoc( sourceVol->MakeVolumeLocationFromRAS( iRASSeed ) );
  seedLoc.SetFrame( iFrame );
  float seedValue = sourceVol->GetMRINearestValue( seedLoc );

  // Create locations in the source volume space. These will be used
  // in the main loop.
  VolumeLocation sourceLoc( sourceVol->MakeVolumeLocationFromRAS( iRASSeed ) );
  sourceLoc.SetFrame( iFrame );
  VolumeLocation sourceFromLoc( sourceVol->MakeVolumeLocationFromRAS( iRASSeed ) );
  sourceFromLoc.SetFrame( iFrame );

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
  if ( !iParams.mb3D ) {
    planeRAS.Set( iRASSeed );
    if ( iParams.mbWorkPlaneX ) {
      planeNRAS.Set( 1, 0, 0 );
    }
    if ( iParams.mbWorkPlaneY ) {
      planeNRAS.Set( 0, 1, 0 );
    }
    if ( iParams.mbWorkPlaneZ ) {
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
  while ( checkPairs.size() > 0 &&
          !this->DoStopRequested() ) {

    // Get this voxel, the from voxel, and RAS versions.
    CheckPair checkPair = checkPairs.back();
    checkPairs.pop_back();

    Point3<int> index = checkPair.mToIndex;
    VolumeLocation loc( iVolume.MakeVolumeLocationFromIndex( index.xyz() ) );
    loc.SetFrame( iFrame );

    Point3<int> fromIndex = checkPair.mFromIndex;
    VolumeLocation fromLoc( iVolume.MakeVolumeLocationFromIndex( fromIndex.xyz() ) );
    fromLoc.SetFrame( iFrame );

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
    if ( !iVolume.IsInBounds( loc ) ) {
      continue;
    }
    if ( bDifferentSource && !sourceVol->IsInBounds( sourceLoc ) ) {
      continue;
    }

    // See if we've visited here before.
    if ( bVisited->Get_Unsafe( index.x(), index.y(), index.z() ) ) {
      continue;
    }
    bVisited->Set_Unsafe( index.x(), index.y(), index.z(), true );


    if ( iParams.mbStopAtROIs ) {
      if ( iVolume.IsOtherRASSelected( ras.xyz(), iVolume.GetSelectedROI() ) ) {
        continue;
      }
    }

    // Check max distance.
    if ( iParams.mMaxDistance > 0 ) {
      float distance = sqrt( ((ras[0]-iRASSeed[0]) * (ras[0]-iRASSeed[0])) +
                             ((ras[1]-iRASSeed[1]) * (ras[1]-iRASSeed[1])) +
                             ((ras[2]-iRASSeed[2]) * (ras[2]-iRASSeed[2])) );
      if ( distance > iParams.mMaxDistance ) {
        continue;
      }
    }

    // Check only zero.
    if ( iParams.mbOnlyZero ) {
      float value = iVolume.GetMRINearestValue( loc );
      if ( value != 0 ) {
        continue;
      }
    }

    // Check fuzziness.
    if ( iParams.mFuzziness > 0 ) {
      float value = sourceVol->GetMRINearestValue( sourceLoc );
      switch ( iParams.mFuzzinessType ) {
      case Params::seed:
        if ( fabs( value - seedValue ) > iParams.mFuzziness ) {
          continue;
        }
        break;
      case Params::gradient: {
          float fromValue =
            sourceVol->GetMRINearestValue( sourceFromLoc );
          if ( fabs( value - fromValue ) > iParams.mFuzziness ) {
            continue;
          }
        }
        break;
      }
    }

    // Check if this is an path or an ROI. If so, and our params say
    // not go to here, continue.
    if ( iParams.mbStopAtPaths ) {

      // Create a line from our current point to this check
      // point. Then see if the segment intersects the path. If so,
      // don't proceed.
      bool bCross = false;
      Point3<float> x;
      // For each path...
      vector<Path<float>*>::const_iterator tPath;
      PathManager& pathMgr = PathManager::GetManager();
      vector<Path<float>*> const& pathList = pathMgr.GetPathList();
      for ( tPath = pathList.begin();
            tPath != pathList.end() && !bCross; ++tPath ) {
        Path<float> const* path = *tPath;
        if ( path->GetNumVertices() > 0 ) {

          // For each segment in the path...
          int cVertices = path->GetNumVertices();
          for ( int nCurVertex = 1;
                nCurVertex < cVertices && !bCross; nCurVertex++ ) {

            int nBackVertex = nCurVertex - 1;

            // Get the end points of this segment.
            Point3<float> curVertex = path->GetVertexAtIndex( nCurVertex );
            Point3<float> backVertex = path->GetVertexAtIndex( nBackVertex );

            // Get the view normal.
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
            viewIdxf.Set( viewIdx[0], viewIdx[1], viewIdx[2]);
            fromIdxf.Set( fromIdx[0], fromIdx[1], fromIdx[2] );
            curIdxf.Set( curIdx[0], curIdx[1], curIdx[2] );

            // Calculate a vector between the two path segment
            // points. To create a plane on the segment that's
            // perpendicular to the screen, calc a plane normal by
            // crossing the segment vector and the view normal.
            segVertexIdxf = curVertexIdxf - backVertexIdxf;
            normalIdxf = VectorOps::Cross( segVertexIdxf, viewIdxf );
            normalIdxf = VectorOps::Normalize( normalIdxf );

            // Intersect the segment between the points in the volume
            // we're checking and the plane that's on the path
            // segment.
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
            << "  normal\t" << normalIdxf << endl
            << "  int\t"
            << VectorOps::IntersectionResultToString(rInt) <<  endl;
#endif
            if ( VectorOps::intersect == rInt ) {

              // Make sure x is in the cuboid formed by curVertex and
              // backVertex. We extend the cuboid by 0.9 in all
              // directions to compensate for volumes that are rotated
              // oddly. This number seems to work best.
              if ( !((float)x[0] <
                     (float)MIN(curVertexIdx[0],backVertexIdx[0]) - 0.9 ||
                     (float)x[0] >
                     (float)MAX(curVertexIdx[0],backVertexIdx[0]) + 0.9 ||
                     (float)x[1] <
                     (float)MIN(curVertexIdx[1],backVertexIdx[1]) - 0.9 ||
                     (float)x[1] >
                     (float)MAX(curVertexIdx[1],backVertexIdx[1]) + 0.9 ||
                     (float)x[2] <
                     (float)MIN(curVertexIdx[2],backVertexIdx[2]) - 0.9 ||
                     (float)x[2] >
                     (float)MAX(curVertexIdx[2],backVertexIdx[2]) + 0.9 )) {

                bCross = true;
              }
            }
          }
        }
      }

      // If we crossed a path, continue.
      if ( bCross ) {
        continue;
      }
    }

    // Call the user compare function to give them a chance to bail.
    if ( !this->CompareVoxel( ras.xyz(), iFrame ) ) {
      continue;
    }

    // Call the user function.
    this->DoVoxel( ras.xyz(), iFrame );

    // Add adjacent voxels. Look at all adjacent voxels. If we're only
    // working on a single plane, make sure the voxel is in the same
    // plane.
    int beginX = index.x() - 1;
    int endX   = index.x() + 1;
    int beginY = index.y() - 1;
    int endY   = index.y() + 1;
    int beginZ = index.z() - 1;
    int endZ   = index.z() + 1;
#if PRINTOUT
    cerr << "\tbounds: " << beginX << ", " << endX << "  "
    << beginY << ", " << endY << "  "
    << beginZ << ", " << endZ
    << endl;
#endif

    Point3<int> newIndex;
    CheckPair newPair;
    newPair.mFromIndex = index;
    for ( int nZ = beginZ; nZ <= endZ; nZ += 1 ) {
      for ( int nY = beginY; nY <= endY; nY += 1 ) {
        for ( int nX = beginX; nX <= endX; nX += 1 ) {

          if ( !iParams.mbDiagonal ) {
            if ( !( (nX == index.x() && nY == index.y()) ||
                    (nX == index.x() && nZ == index.z()) ||
                    (nY == index.y() && nZ == index.z())) ) {
              continue;
            }
          }

          newIndex.Set( nX, nY, nZ );

          // If not 3D, check if voxel is in same plane.
          if ( !iParams.mb3D ) {
            VectorOps::IntersectionResult rInt =
              iVolume.VoxelIntersectsPlane( newIndex, planeIdx, planeNIdx,
                                            intersectionIdx );
            if ( VectorOps::intersect != rInt ) {
              continue;
            }
          }

#if 0
          Point3<float> intersectionRAS;
          iVolume.MRIIndexToRAS( intersectionIdx.xyz(),
                                 intersectionRAS.xyz() );

          cerr << "Added " << newIndex << " from " << index << endl;
          cerr << "\tplaneRAS " << planeRAS
          << "\tplaneIdx " << planeIdx
          << endl
          << "\tplaneNRAS " << planeNRAS
          << "\tplaneNIdx " << planeNIdx
          << endl
          << "\tintersectionIdx " << intersectionIdx
          << "\tintersectionRAS " << intersectionRAS
          << endl;
#endif

          newPair.mToIndex = newIndex;
          checkPairs.push_back( newPair );
        }
      }
    }

  }

  this->DoEnd();

  delete bVisited;

  mVolume = NULL;
  mParams = NULL;
}

VolumeLocation::VolumeLocation ( VolumeCollection const& iVolume,
                                 float const iRAS[3] ) :
  DataLocation( iRAS ), 
  mVolume( &iVolume ),
  mFrame( 0 ) {

  mVolume->RASToMRIIndex( iRAS, mIdxf );
  mVolume->RASToMRIIndex( iRAS, mIdxi );
}

VolumeLocation::VolumeLocation ( VolumeCollection const& iVolume,
                                 int const iIndex[3] ) :
  DataLocation(), 
  mVolume( &iVolume ),
  mFrame( 0 ) {

  mVolume->MRIIndexToRAS( iIndex, mRAS );
  mIdxi[0] = iIndex[0];
  mIdxi[1] = iIndex[1];
  mIdxi[2] = iIndex[2];
  mIdxf[0] = iIndex[0];
  mIdxf[1] = iIndex[1];
  mIdxf[2] = iIndex[2];
}

VolumeLocation::VolumeLocation ( VolumeCollection const& iVolume,
                                 float const iRAS[3],
				 int iFrame ) :
  DataLocation( iRAS ),
  mVolume( &iVolume ),
  mFrame( iFrame ) {

  mVolume->RASToMRIIndex( iRAS, mIdxf );
  mVolume->RASToMRIIndex( iRAS, mIdxi );
}

VolumeLocation::VolumeLocation ( VolumeCollection const& iVolume,
                                 int const iIndex[3],
				 int iFrame ) :
  DataLocation(),
  mVolume( &iVolume ),
  mFrame( iFrame ) {

  mVolume->MRIIndexToRAS( iIndex, mRAS );
  mIdxi[0] = iIndex[0];
  mIdxi[1] = iIndex[1];
  mIdxi[2] = iIndex[2];
  mIdxf[0] = iIndex[0];
  mIdxf[1] = iIndex[1];
  mIdxf[2] = iIndex[2];
}

VolumeLocation::VolumeLocation ( VolumeLocation const& iLoc ) :
  DataLocation( iLoc ),
  mVolume( &iLoc.GetVolume() ), 
  mFrame( 0 ) {

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

void
VolumeLocation::SetFromIndex ( int const iIdx[3] ) {

  mIdxi[0] = iIdx[0];
  mIdxi[1] = iIdx[1];
  mIdxi[2] = iIdx[2];
  mIdxf[0] = (float) iIdx[0];
  mIdxf[1] = (float) iIdx[1];
  mIdxf[2] = (float) iIdx[2];

  mVolume->MRIIndexToRAS( iIdx, mRAS );
}

void
VolumeLocation::SetFrame ( int iFrame ) {
  mFrame = iFrame;
}

VolumeLocation&
VolumeLocation::operator= ( VolumeLocation const& iLocation ) {

  // Check for self-assignment.
  if( this == &iLocation )
    return *this;

  // Set our values.
  mVolume = &iLocation.GetVolume();
  SetFromRAS( iLocation.RAS() );
  SetFrame( iLocation.GetFrame() );

  // Return ourself.
  return *this;
}
