#include <stdio.h>

#include "mriVolume.h"
#include "error.h"
#include "xUtilities.h"
#include "mri_conform.h"
#include "mriTransform.h"

//#define LINEAR_CORONAL_RAS_TO_CORONAL_RAS       21
// should be in transform.h if they aren't already

char Volm_ksaErrorStrings[Volm_knNumErrorCodes][Volm_knErrStringLen] = {
  "No error.",
  "Invalid signature.",
  "Invalid parameter.",
  "Invalid index.",
  "Memory allocation failed.",
  "Couldn't read volume.",
  "Couldn't copy volume.",
  "Couldn't read transform.",
  "Couldn't copy transform.",
  "Couldn't normalize volume.",
  "Couldn't export volume to COR format.",
  "MRI volume not present.",
  "Scanner transform not present.",
  "Index to RAS transform not present."
};

Volm_tErr Volm_New ( mriVolumeRef* opVolume ) {
  
  Volm_tErr    eResult   = Volm_tErr_NoErr;
  mriVolumeRef this      = NULL;
  MATRIX*      tmpMatrix = NULL;
  
  DebugEnterFunction( ("Volm_New( opVolume=%p )", opVolume) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (opVolume != NULL),
		     eResult, Volm_tErr_InvalidParamater );
  
  /* allocate us. */
  DebugNote( ("Allocating storage for this") );
  this = (mriVolumeRef) malloc( sizeof( mriVolume ));
  DebugAssertThrowX( (NULL != this), eResult, Volm_tErr_AllocationFailed );
  
  /* init members to defaults */
  DebugNote( ("Setting members to default values") );
  this->mSignature             = Volm_kSignature;
  this->mnDimensionX           = 0;
  this->mnDimensionY           = 0;
  this->mnDimensionZ           = 0;
  this->mnDimensionFrame       = 1;
  this->mpMriValues           = NULL;
  this->mpSnapshot             = NULL;
  this->mpMaxValues            = NULL;
  this->mIdxToRASTransform     = NULL;
  this->mDisplayTransform      = NULL;
  this->mScannerTransform      = NULL;
  bzero( this->msSubjectName, sizeof( this->msSubjectName ) );
  bzero( this->msVolumeName, sizeof( this->msVolumeName ) );
  bzero( this->msOriginalPath, sizeof( this->msOriginalPath ) );
  this->mfColorBrightness      = Volm_kfDefaultBrightness;
  this->mfColorContrast        = Volm_kfDefaultContrast;
  bzero( this->maColorTable, sizeof( this->maColorTable ) );

  /* make our color table with default values */
  eResult = Volm_MakeColorTable( this );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );
  
  /* make our tal transforms */
  Trns_New( &this->mMNITalLtzToRealTalTransform );
  tmpMatrix = MatrixIdentity( 4, NULL );
  stuff_four_by_four( tmpMatrix,
		      0.99,       0,     0, 0,
		      0.00,  0.9688, 0.042, 0,
		      0.00, -0.0485, 0.839, 0,
		      0.00,       0,     0, 1 );
  Trns_CopyARAStoBRAS( this->mMNITalLtzToRealTalTransform, tmpMatrix );
  MatrixIdentity( 4, tmpMatrix );
  Trns_CopyAtoRAS( this->mMNITalLtzToRealTalTransform, tmpMatrix );
  Trns_CopyBtoRAS( this->mMNITalLtzToRealTalTransform, tmpMatrix );
  
  Trns_New( &this->mMNITalGtzToRealTalTransform );
  MatrixIdentity( 4, tmpMatrix );
  stuff_four_by_four( tmpMatrix,
		      0.99,       0,      0, 0,
		      0.00,  0.9688,  0.046, 0,
		      0.00, -0.0485, 0.9189, 0,
		      0.00,       0,      0, 1 );
  Trns_CopyARAStoBRAS( this->mMNITalGtzToRealTalTransform, tmpMatrix );
  MatrixIdentity( 4, tmpMatrix );
  Trns_CopyAtoRAS( this->mMNITalGtzToRealTalTransform, tmpMatrix );
  Trns_CopyBtoRAS( this->mMNITalGtzToRealTalTransform, tmpMatrix );

  /* Make our temp vars. */
  this->mpTmpScreenIdx = VectorAlloc( 4, MATRIX_REAL );
  this->mpTmpMRIIdx    = VectorAlloc( 4, MATRIX_REAL );
  
  /* return us. */
  *opVolume = this;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  
  if( NULL != this )
    Volm_Delete( &this );
  
  EndDebugCatch;
  
  DebugExitFunction;
  
  if( NULL != tmpMatrix )
    MatrixFree( &tmpMatrix );
  
  return eResult;
}

Volm_tErr Volm_Delete ( mriVolumeRef* iopVolume ) {
  
  Volm_tErr    eResult = Volm_tErr_NoErr;
  mriVolumeRef this    = NULL;
  
  DebugEnterFunction( ("Volm_Delete( iopVolume=%p )", iopVolume) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iopVolume), 
		     eResult, Volm_tErr_InvalidParamater );
  
  /* get this and verify it */
  this = *iopVolume;
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* free member data */
  if( NULL != this->mpMriValues ) {
    MRIfree( &(this->mpMriValues) );
  }
  if( NULL != this->mpSnapshot ) {
    free( this->mpSnapshot );
  }
  if( NULL != this->mpMaxValues ) {
    free( this->mpMaxValues );
  }
  if( NULL != this->mDisplayTransform ) {
    Trns_Delete( &(this->mDisplayTransform) );
  }
  if( NULL != this->mMNITalLtzToRealTalTransform ) {
    Trns_Delete( &(this->mMNITalLtzToRealTalTransform) );
  }
  if( NULL != this->mMNITalGtzToRealTalTransform ) {
    Trns_Delete( &(this->mMNITalGtzToRealTalTransform) );
  }
  if( NULL != this->mScannerTransform) {
    Trns_Delete( &(this->mScannerTransform) );
  }
  if( NULL != this->mpTmpScreenIdx ) {
    MatrixFree( &this->mpTmpScreenIdx );
  }
  if( NULL != this->mpTmpMRIIdx ) {
    MatrixFree( &this->mpTmpMRIIdx );
  }
  
  /* free us */
  free( this );
  
  /* return null */
  *iopVolume = NULL;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_DeepClone  ( mriVolumeRef  this, 
			    mriVolumeRef* opVolume ) {
  
  Volm_tErr    eResult    = Volm_tErr_NoErr;
  Trns_tErr    eTransform = Trns_tErr_NoErr;
  mriVolumeRef clone      = NULL;
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != opVolume), eResult, Volm_tErr_InvalidParamater );
  
  /* allocate a clone */
  DebugNote( ("Creating clone") );
  eResult = Volm_New( &clone );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );
  
  /* set dimension */
  DebugNote( ("Setting dimension") );
  clone->mnDimensionX = this->mnDimensionX;
  clone->mnDimensionY = this->mnDimensionY;
  clone->mnDimensionZ = this->mnDimensionZ;
  clone->mnDimensionFrame = this->mnDimensionFrame;
  
  /* copy mri volumes */
  DebugNote( ("Cloning normalized volume") );
  clone->mpMriValues = MRIcopy( this->mpMriValues, NULL );
  DebugAssertThrowX( (NULL != clone->mpMriValues),
		     eResult, Volm_tErr_CouldntCopyVolume );
  
  /* copy the transforms */
  if( NULL != this->mIdxToRASTransform ) {
    DebugNote( ("Copying index to RAS transform") );
    eTransform = Trns_DeepClone( this->mIdxToRASTransform, 
				 &(clone->mIdxToRASTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
		       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if( NULL != this->mDisplayTransform ) {
    DebugNote( ("Copying display transform") );
    eTransform = Trns_DeepClone( this->mDisplayTransform, 
				 &(clone->mDisplayTransform) );
    DebugAssertThrowX( (NULL != clone->mDisplayTransform),
		       eResult, Volm_tErr_CouldntCopyTransform );
    memcpy( clone->mDisplayTransform, this->mDisplayTransform, 
	    sizeof(this->mDisplayTransform) );
  }
  if( NULL != this->mMNITalLtzToRealTalTransform ) {
    DebugNote( ("Copying real tal ltz transform") );
    eTransform = Trns_DeepClone( this->mMNITalLtzToRealTalTransform, 
				 &(clone->mMNITalLtzToRealTalTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
		       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if( NULL != this->mMNITalGtzToRealTalTransform ) {
    DebugNote( ("Copying real tal gtz transform") );
    eTransform = Trns_DeepClone( this->mMNITalGtzToRealTalTransform, 
				 &(clone->mMNITalGtzToRealTalTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
		       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if( NULL != this->mScannerTransform ) {
    DebugNote( ("Copying scanner transform") );
    eTransform = Trns_DeepClone( this->mScannerTransform, 
				 &(clone->mScannerTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
		       eResult, Volm_tErr_CouldntCopyTransform );
  }
  
  /* copy the string */
  xUtil_strncpy( clone->msSubjectName, this->msSubjectName, 
		 sizeof(clone->msSubjectName) );
  xUtil_strncpy( clone->msVolumeName, this->msVolumeName, 
		 sizeof(clone->msVolumeName) );
  xUtil_strncpy( clone->msOriginalPath, this->msOriginalPath, 
		 sizeof(clone->msOriginalPath) );
  
  /* copy color table information */
  clone->mfColorBrightness = this->mfColorBrightness;
  clone->mfColorContrast   = this->mfColorContrast;
  
  /* allocate and copy color table */
  memcpy( clone->maColorTable, this->maColorTable, 
	  sizeof( this->maColorTable ));
  
  /* return the clone */
  *opVolume = clone;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  
  if( NULL != clone )
    Volm_Delete( &clone );
  
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_CreateFromVolume ( mriVolumeRef this,
				  mriVolumeRef iVolume ) {
  
  Volm_tErr eResult          = Volm_tErr_NoErr;
  int       nDimensionX      = 0;
  int       nDimensionY      = 0;
  int       nDimensionZ      = 0;
  int       nType            = 0;
  int       nDimensionFrames = 0;
  MRI*      mriVolume        = NULL;
  
  DebugEnterFunction( ("Volm_CreateFromVolume( this=%p, iVolume=%p )",
		       this, iVolume ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iVolume), eResult, Volm_tErr_InvalidParamater );
  
  /* Get the necessary parameters we need to create a new MRI from the
     volume we got passed in. */
  DebugNote( ("Getting volume dimensions") );
  eResult = Volm_GetDimensions( iVolume, 
				&nDimensionX, &nDimensionY, &nDimensionZ );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  DebugNote( ("Getting number of frames") );
  eResult = Volm_GetNumberOfFrames( iVolume, &nDimensionFrames );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  DebugNote( ("Getting volume type") );
  eResult = Volm_GetType( iVolume, &nType );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Create a new MRI of that size. */
  DebugNote( ("Creating the mri with dimensions %d,%d,%d type %d nframes %d",
	      nDimensionX, nDimensionY, nDimensionZ, nType, nDimensionFrames));
  mriVolume = MRIallocSequence( nDimensionX, nDimensionY, nDimensionZ,
				nType, nDimensionFrames );
  DebugAssertThrowX( (NULL != mriVolume), 
		     eResult, Volm_tErr_CouldntReadVolume );

  /* Set from the new MRI. */
  DebugNote( ("Setting from MRI") );
  eResult = Volm_SetFromMRI_( this, mriVolume );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Set the subject name from the original and set the volume name
     and path to blanks. */
  DebugNote( ("Copying original subject name") );
  Volm_CopySubjectName( iVolume, this->msSubjectName,
		 sizeof(this->msSubjectName) );
  DebugNote( ("Setting volume name to empty") );
  xUtil_strncpy( this->msVolumeName, "", sizeof(this->msVolumeName) );
  DebugNote( ("Setting path to empty") );
  xUtil_strncpy( this->msOriginalPath, "", sizeof(this->msOriginalPath));

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ImportData ( mriVolumeRef this,
			    char*        isSource ) {
  
  Volm_tErr eResult           = Volm_tErr_NoErr;
  MRI*      mriVolume        = NULL;
  
  DebugEnterFunction( ("Volm_ImportData( this=%p, isSource=%s )",
		       this, isSource ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isSource), eResult, Volm_tErr_InvalidParamater );
  
  /* attempt to read the volume in */
  DebugNote( ("Importing volume with MRIread") );
  mriVolume = MRIread( isSource );
  DebugAssertThrowX( (NULL != mriVolume), 
		     eResult, Volm_tErr_CouldntReadVolume );

  DebugNote( ("Setting from MRI") );
  eResult = Volm_SetFromMRI_( this, mriVolume );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* save the volume source */
  DebugNote( ("Copying original path") );
  xUtil_strncpy( this->msOriginalPath, isSource, 
		 sizeof(this->msOriginalPath) );

  /* get the subject and volume name */
  eResult = Volm_ExtractAndSetSubjectName( this, isSource );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );
  eResult = Volm_ExtractAndSetVolumeName( this, isSource );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );
  

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


Volm_tErr Volm_SetFromMRI_ ( mriVolumeRef this,
			     MRI*        iMRI ) {
  
  Volm_tErr eResult           = Volm_tErr_NoErr;
  MATRIX*   identity          = NULL;
  MATRIX*   scannerTransform  = NULL;
  MATRIX*   idxToRASTransform = NULL;
  MATRIX    *m_resample ;
  
  DebugEnterFunction( ("Volm_SetFromMRI_( this=%p, iMRI=%p )", this, iMRI ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iMRI), eResult, Volm_tErr_InvalidParamater );
  
  /* Get the reange. */
  MRIvalRange(iMRI, &this->min_val, &this->max_val) ;
  
  /* find the resampling matrix */
  DebugNote( ("Computing norming matrix ") );
  m_resample = MRIgetConformMatrix( iMRI );
  DebugAssertThrowX( (NULL != m_resample), 
                     eResult, Volm_tErr_CouldntNormalizeVolume );
//#define _VID_DEBUG  
#ifdef _VID_DEBUG
  fprintf(stderr, "m_resample = \n");
  MatrixPrint(stderr, m_resample);
#endif
  
  /* grab the end-start values as our dimension. also nframes. */
  /* 
  this->mnDimensionX = abs(iMRI->xend - iMRI->xstart);
  this->mnDimensionY = abs(iMRI->yend - iMRI->ystart);
  this->mnDimensionZ = abs(iMRI->zend - iMRI->zstart);
  */
  this->mnDimensionX = iMRI->width;
  this->mnDimensionY = iMRI->height;
  this->mnDimensionZ = iMRI->depth;
  this->mnDimensionFrame = iMRI->nframes;

  this->m_resample = m_resample ;
  this->m_resample_orig = MatrixCopy(m_resample, NULL) ;
  DebugNote( ("Calculating inverse of resample matrix") );
  this->m_resample_inv = MatrixInverse( m_resample, NULL );

  /* set the volumes in this */
  this->mpMriValues = iMRI;
  
  /* grab the scanner transform and copy it into our transform */
  if( NULL != this->mScannerTransform ) {
    DebugNote( ("Deleting existing scanner transform") );
    Trns_Delete( &(this->mScannerTransform) );
  }
  DebugNote( ("Creating scanner transform") );
  Trns_New( &this->mScannerTransform );
  DebugNote( ("Getting scanner transform matrix") );
  scannerTransform = MRIgetVoxelToRasXform( iMRI );
  DebugAssertThrowX( (NULL != scannerTransform),
		     eResult, Volm_tErr_AllocationFailed );
  DebugNote( ("Copying scanner transform matrix into transform") );
  Trns_CopyARAStoBRAS( this->mScannerTransform, scannerTransform );
  identity = MatrixIdentity( 4, NULL );
  DebugNote( ("Copying identity matrix into scanner transform") );
  Trns_CopyAtoRAS( this->mScannerTransform, identity );
  Trns_CopyBtoRAS( this->mScannerTransform, identity );
  
  /* do the same for the idx -> ras transform. note that this is the
     same as the scanner transform...? */
  if( NULL != this->mIdxToRASTransform ) {
    DebugNote( ("Deleting existing idx to ras transform") );
    Trns_Delete( &(this->mIdxToRASTransform) );
  }
  
  /*E* dangerous to check slice_direction - for .mgh at least this
    shoule break unless this changes to if(ras_good_flag and not that
    kind of coronal) else check slice_direction - actually, what is
    that check for? */

  //  if (iMRI->slice_direction != MRI_CORONAL)
  if ( 
      (iMRI->ras_good_flag &&
       ( iMRI->x_r != -1.0 
	 || iMRI->x_a != 0.0
	 || iMRI->x_s != 0.0
	 || iMRI->y_r != 0.0
	 || iMRI->y_a != 0.0
	 || iMRI->y_s != -1.0
	 || iMRI->z_r != 0.0
	 || iMRI->z_a != 1.0
	 || iMRI->z_s != 0.0
	 // || iMRI->c_r != 0.0
	 // || iMRI->c_a != 0.0
	 // || iMRI->c_s != 0.0
	 )
       )
      || (!iMRI->ras_good_flag && iMRI->slice_direction != MRI_CORONAL)
      )
    {
      DebugNote( ("Creating idx to ras transform") );
      Trns_New( &this->mIdxToRASTransform );
      DebugNote( ("Getting idx to ras matrix") );
      idxToRASTransform = MRIgetVoxelToRasXform( iMRI );
      // that includes voxelsize
      DebugAssertThrowX( (NULL != idxToRASTransform),
                         eResult, Volm_tErr_AllocationFailed );
      DebugNote( ("Copying idx to ras transform matrix into transform") );
      Trns_CopyAtoRAS( this->mIdxToRASTransform, idxToRASTransform );
      DebugNote( ("Copying identity matrix into idx to ras transform") );
      Trns_CopyARAStoBRAS( this->mIdxToRASTransform, identity ); /* no display xform */
      Trns_CopyBtoRAS( this->mIdxToRASTransform, identity );
    }
  else 
    {
      idxToRASTransform = MatrixAlloc( 4, 4, MATRIX_REAL );
      MatrixClear( idxToRASTransform );
      *MATRIX_RELT(idxToRASTransform,1,1) = -1.0;
      *MATRIX_RELT(idxToRASTransform,2,3) = 1.0;
      *MATRIX_RELT(idxToRASTransform,3,2) = -1.0;
      *MATRIX_RELT(idxToRASTransform,1,4) = 128;
      *MATRIX_RELT(idxToRASTransform,2,4) = -128;
      *MATRIX_RELT(idxToRASTransform,3,4) = 128;
      *MATRIX_RELT(idxToRASTransform,4,4) = 1.0;
      
      DebugNote( ("Creating idx to ras transform") );
      Trns_New( &this->mIdxToRASTransform );
      DebugNote( ("Copying idx to ras transform matrix into transform") );
      Trns_CopyARAStoBRAS( this->mIdxToRASTransform, idxToRASTransform );
      DebugNote( ("Copying identity matrix into idx to ras transform") );
      Trns_CopyAtoRAS( this->mIdxToRASTransform, identity );
      Trns_CopyBtoRAS( this->mIdxToRASTransform, identity );
    }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  if( NULL != scannerTransform )
    MatrixFree( &scannerTransform );
  if( NULL != idxToRASTransform )
    MatrixFree( &idxToRASTransform );
  if( NULL != identity )
    MatrixFree( &identity );
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ExportNormToCOR ( mriVolumeRef this,
				 char*        isPath ) {
  
  Volm_tErr eResult              = Volm_tErr_NoErr;
  char      sPath[mri_knPathLen] = "";
  int       eMRI                 = NO_ERROR;
  
  DebugEnterFunction( ("Volm_ExportNormToCOR( this=%p, isPath=%s )",
		       this, isPath) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* if out isPath is null, use the original path. */
  if( NULL == isPath ) {
    xUtil_strncpy( sPath, this->msOriginalPath, sizeof( sPath ) );
  } else {
    xUtil_strncpy( sPath, isPath, sizeof( sPath ) );
  }
  
  /* write the volume */
  DebugNote( ("Writing volume with MRIwriteType") );
  eMRI = MRIwriteType( this->mpMriValues, sPath,MRI_CORONAL_SLICE_DIRECTORY );
  DebugAssertThrowX( (NO_ERROR == eMRI), 
		     eResult, Volm_tErr_CouldntExportVolumeToCOR );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

//#define _VLDT_DEBUG

Volm_tErr Volm_LoadDisplayTransform ( mriVolumeRef this,
				      char*        isFileName ) {
  
  Volm_tErr eResult   = Volm_tErr_NoErr;
  Trns_tErr eTransform = Trns_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_LoadDisplayTransform( this=%p, isFileName=%s )", 
		       this, isFileName) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isFileName),
		     eResult, Volm_tErr_InvalidParamater );
  
  /* try to read a transform and make sure we got it */
  DebugNote( ("Creating transform from lta transform file %s", isFileName) );
  eTransform = Trns_NewFromLTA( &(this->mDisplayTransform), isFileName );
  //E/ Trns_NewFromLTA copies the matrix right out of
  //LTATransform->xforms[0].m_L and notes type
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
		     eResult, Volm_tErr_CouldntReadTransform );
  
  /* convert (VOX->VOX xform to) RAS->RAS transform to screen->voxel coords */
  switch (this->mDisplayTransform->type)
    {
    case LINEAR_CORONAL_RAS_TO_CORONAL_RAS:
      {
	extern MATRIX *gm_screen2ras ;
        MATRIX *m_coronalras2coronalras, *m_ras2ras, *m_ras2ras_inverse, *m_tmp, *m_ras2vox,
	  *m_thisorientation2coronal, *m_coronal2thisorientation;

        Trns_GetARAStoBRAS(this->mDisplayTransform, &m_coronalras2coronalras) ;
	m_thisorientation2coronal = MRIgetConformMatrix(this->mpMriValues);
	//E/ = MRIgetResampleMatrix(mri, templ) = src_inv*m_templ
	m_coronal2thisorientation = MatrixInverse(m_thisorientation2coronal, NULL);

	//E/ want m_ras2ras = m_coronal2thisorientation_src * m_coronalras_t2coronalras_s

	//E/ NOT * m_thisorientation2coronal_targ - leave that end
	//coronal - we don't have any noncoronal targ orientation in
	//mind.

	m_ras2ras = MatrixMultiply(m_coronal2thisorientation, m_coronalras2coronalras, NULL); //m_ras_c2ras_s
	
        m_ras2vox = MRIgetRasToVoxelXform(this->mpMriValues);
        m_ras2ras_inverse = MatrixInverse(m_ras2ras, NULL) ; // ras_s2ras_c
        m_tmp = MatrixMultiply(m_ras2ras_inverse, gm_screen2ras, NULL) ;
        MatrixMultiply(m_ras2vox, m_tmp, this->m_resample) ;
	// m_resample_t2s = m_ras2vox_s * m_ras_t2ras_s * gm_screen2ras_t
        DebugNote( ("Calculating inverse of resample matrix") );
        MatrixInverse( this->m_resample, this->m_resample_inv );
        printf("resample matrix is:\n") ;
        MatrixPrint(stdout, this->m_resample) ;
        MatrixFree(&m_ras2vox) ; MatrixFree(&m_tmp);MatrixFree(&m_ras2ras_inverse);
        break ;
      }
    case LINEAR_VOX_TO_VOX:
      {
        MATRIX *m_ras2ras, *m_vox2vox ;

        /* assume it's in coronal->coronal coordinates */
        Trns_GetARAStoBRAS(this->mDisplayTransform, &m_vox2vox) ;
        m_ras2ras = MRIvoxelXformToRasXform(this->mpMriValues, this->mpMriValues, m_vox2vox, NULL);

        Trns_CopyARAStoBRAS(this->mDisplayTransform, m_ras2ras) ;
        MatrixFree(&m_vox2vox) ; MatrixFree(&m_ras2ras) ;
#ifdef _VLDT_DEBUG
	fprintf(stderr, "Volm_LoadDisplayTransform: case LINEAR_VOX_TO_VOX\n");
#endif
      }
      /*  NOBREAK!!  break ;*/

    case LINEAR_RAS_TO_RAS:
      //E/ could probably fix it for this case
      //    case LINEAR_RAS_TO_CORONAL_RAS:
      //E/ this case was probably a mistake
      {
	extern MATRIX *gm_screen2ras ;
	MATRIX *m_ras2ras, *m_ras2ras_inverse, *m_tmp, *m_ras2vox ;
	
	if (this->mDisplayTransform->type == LINEAR_VOX_TO_VOX)
	  fprintf(stderr, "Don't really know what to do with a LTA of type LINEAR_VOX_TO_VOX - we'll pretend it's LINEAR_VOX_TO_CONFORM_VOX and see what happens.\n");
	//	if (this->mDisplayTransform->type == LINEAR_RAS_TO_CORONAL_RAS)
	//	  fprintf(stderr, "Don't really know what to do with a LTA of type LINEAR_RAS_TO_CORONAL_RAS - we'll pretend it's LINEAR_CORONAL_RAS_TO_CORONAL_RAS and see what happens.\n");
	if (this->mDisplayTransform->type == LINEAR_RAS_TO_RAS)
	  fprintf(stderr, "Don't really know what to do with a LTA of type LINEAR_RAS_TO_RAS - we'll pretend it's LINEAR_CORONAL_RAS_TO_CORONAL_RAS and see what happens.\n");

	Trns_GetARAStoBRAS(this->mDisplayTransform, &m_ras2ras) ;
	//E/ m_ras2ras is either from v2v case above or what mri_rigid_register wrote out, i.e. ras_c2ras_s

	m_ras2vox = MRIgetRasToVoxelXform(this->mpMriValues) ;
	m_ras2ras_inverse = MatrixInverse(m_ras2ras, NULL) ;
	m_tmp = MatrixMultiply(m_ras2ras_inverse, gm_screen2ras, NULL) ;
	MatrixMultiply(m_ras2vox, m_tmp, this->m_resample) ;
	//E/ m_resample = m_ras2vox_s * m_ras_c2ras_s * gm_screen2ras_c
	//E/ the value of gm_screen2ras is hardcoded in tkmedit.c = "W"
	DebugNote( ("Calculating inverse of resample matrix") );
	MatrixInverse( this->m_resample, this->m_resample_inv );

	//E/ now m_resample_inv = gm_ras_c2screen * m_ras_s2ras_c * m_vox2ras_s

	// Trns_GetARAStoBRAS() delivers whatever you wrote out to the
	// lta file, in e.g. mri_rigid_register, which is ras_s2ras_c.
	// Or ras_c2ras_s.  Uhh..

#ifdef _VLDT_DEBUG
	fprintf(stderr, "Volm_LoadDisplayTransform: case LINEAR_RAS_TO_RAS or fell through\n");
	fprintf(stderr, "m_ras2ras (= m_ras_s2ras_c, I think) = \n");
	MatrixPrint(stderr, m_ras2ras) ;
	fprintf(stderr, "m_ras2vox (= m_ras2vox_s) = \n");
	MatrixPrint(stderr, m_ras2vox) ;
	fprintf(stderr, "m_ras2ras_inverse (= m_ras_c2ras_s, I think) = \n");
	MatrixPrint(stderr, m_ras2ras_inverse ) ;
	fprintf(stderr, "gm_screen2ras(_c) = \n");
	MatrixPrint(stderr, gm_screen2ras) ;
	fprintf(stderr, "m_tmp = m_ras2ras_inverse * gm_screen2ras = \n");
	MatrixPrint(stderr, m_tmp) ;
	fprintf(stderr, "this->m_resample = m_ras2vox * m_ras2ras_inverse * gm_screen2ras = m_ras2vox(_s) * m_ras_c2ras_s * gm_screen2ras(_c) = ");
#endif
	
	printf("resample matrix is:\n") ;
	MatrixPrint(stdout, this->m_resample) ;
	MatrixFree(&m_ras2vox) ; MatrixFree(&m_tmp);MatrixFree(&m_ras2ras_inverse);
	break ;
      }
    default:   /* don't know what to do yet */
      fprintf(stderr, "LTA type isn't LINEAR_VOX_TO_VOX, LINEAR_RAS_TO_RAS, nor LINEAR_CORONAL_RAS_TO_CORONAL_RAS or the new ones - don't know what to do.\n");
      break ;
    }
  
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


/*E*
MATRIX *Volm_VoxelXformToCoronalRasXform(MRI *mri_src, MATRIX *m_vox_c2vox_s, MATRIX *m_ras_c2ras_s)
{
  MATRIX   *V, *W, *m_tmp ;

  //E/ need m_ras_c2ras_s = vox_c2ras_c * m_vox_c2vox_s * ras_s2vox_s
  //E/ = W_s * m_vox_c2vox_s * V_c
  V = MatrixAlloc(4, 4, MATRIX_REAL) ;  // world to voxel transform
  W = MatrixAlloc(4, 4, MATRIX_REAL) ;  // voxel to world transform
 
  //E/ 1mmIsotropicCoronalRas to vox
  *MATRIX_RELT(V, 1, 1) = -1 ; *MATRIX_RELT(V, 1, 4) = 128 ;
  *MATRIX_RELT(V, 2, 3) = -1 ; *MATRIX_RELT(V, 2, 4) = 128 ;
  *MATRIX_RELT(V, 3, 2) = 1 ;  *MATRIX_RELT(V, 3, 4) = 128 ;
  *MATRIX_RELT(V, 4, 4) = 1 ;

  //E/ vox to appropriate Ras
  W = MRIgetVoxelToRasXform(mri_src);

  m_tmp = MatrixMultiply(m_vox_c2vox_s, V, NULL) ;
  m_ras_c2ras_s = MatrixMultiply(W, m_tmp, NULL) ;

  //#define _VVXTCRXDEBUG
#ifdef _VVXTCRXDEBUG
  fprintf(stderr, "Volm_VoxelXformToCoronalRasXform: m_vox_c2vox_s=\n");
  MatrixPrint(stderr, m_vox_c2vox_s);
  fprintf(stderr, "Volm_VoxelXformToCoronalRasXform: V=\n");
  MatrixPrint(stderr, V);
  fprintf(stderr, "Volm_VoxelXformToCoronalRasXform: W=\n");
  MatrixPrint(stderr, W);
  fprintf(stderr, "Volm_VoxelXformToCoronalRasXform: m_ras_c2ras_s=\n");
  MatrixPrint(stderr, m_ras_c2ras_s);
#endif

  MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
  return(m_ras_c2ras_s) ;
}
*/

#if 0
//E/ this case was worthless
MATRIX *Volm_V2CVXtoR2CRX(MRI *mri_src, MATRIX *m_vox_s2vox_c, MATRIX *m_ras_s2ras_c)
{
  MATRIX   *V, *W, *m_tmp ;

  //E/ need m_ras_s2ras_c = vox_c2ras_c * m_vox_s2vox_c * ras_s2vox_s
  //E/ = W_c * m_vox_t2vox_s * V_s
  V = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* world to voxel transform */
  W = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* voxel to world transform */
 
  //E/ vox to 1mmIsotropicCoronalRas
    *MATRIX_RELT(W, 1, 1) = -1 ; *MATRIX_RELT(W, 1, 4) = 128 ;
    *MATRIX_RELT(W, 2, 3) = 1 ; *MATRIX_RELT(W, 2, 4) = -128 ;
    *MATRIX_RELT(W, 3, 2) = -1 ;  *MATRIX_RELT(W, 3, 4) = 128 ;
    *MATRIX_RELT(W, 4, 4) = 1 ;

  //E/ ras to appropriate vox
  V = MRIgetRasToVoxelXform(mri_src);

  m_tmp = MatrixMultiply(m_vox_s2vox_c, V, NULL) ;
  m_ras_s2ras_c = MatrixMultiply(W, m_tmp, NULL) ;

  //#define _VVXTCRXDEBUG
#ifdef _VVXTCRXDEBUG
  fprintf(stderr, "Volm_V2CVXtoR2CRX: m_vox_s2vox_c=\n");
  MatrixPrint(stderr, m_vox_s2vox_c);
  fprintf(stderr, "Volm_V2CVXtoR2CRX: V_s=\n");
  MatrixPrint(stderr, V);
  fprintf(stderr, "Volm_V2CVXtoR2CRX: W_c=\n");
  MatrixPrint(stderr, W);
  fprintf(stderr, "Volm_V2CVXtoR2CRX: m_ras_s2ras_c=\n");
  MatrixPrint(stderr, m_ras_s2ras_c);
#endif

  MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
  return(m_ras_s2ras_c) ;
}
#endif


#if 0
//E/ KLUDGE ALERT, getting rid of this!
/*E* See note above, where this is called.  MRIgetRasToVoxelXform
  which is extract_r_to_i was giving an extra displacement by c_ras.
  Other things call it, where that's probably the _right thing.  So
  this is a rewritten version just for _this case to do. */
MATRIX *Volm_MRIgetRasToVoxelXform(MRI *mri)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras ;
  
  m_voxel_to_ras = Volm_MRIgetVoxelToRasXform(mri) ;
  m_ras_to_voxel = MatrixInverse(m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ;
  return(m_ras_to_voxel) ;
} /* end Volm_MRIgetRasToVoxelXform() */

MATRIX *Volm_MRIgetVoxelToRasXform(MRI *mri)
{
  MATRIX *m;
  float m11, m12, m13, m14;
  float m21, m22, m23, m24;
  float m31, m32, m33, m34;
  float ci, cj, ck;
  
#ifdef _VLDT_DEBUG
  /*E*/ fprintf(stderr, "Volm_MRIgetVoxelToRasXform(): mri = \n");
  /*E*/ MRIdump(mri, stderr);
#endif
  
  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if(m == NULL)
    {
      ErrorReturn(NULL, (ERROR_BADPARM, "Volm_MRIgetVoxelToRasXform(): error allocating matrix"));
    }
  
  if(mri->ras_good_flag)
    {
#ifdef _VLDT_DEBUG
      /*E*/ fprintf(stderr, "Volm_MRIgetVoxelToRasXform(): in ras_good_flag\n");
#endif
      
      m11 = mri->xsize * mri->x_r;  m12 = mri->ysize * mri->y_r;  m13 = mri->zsize * mri->z_r;
      m21 = mri->xsize * mri->x_a;  m22 = mri->ysize * mri->y_a;  m23 = mri->zsize * mri->z_a;
      m31 = mri->xsize * mri->x_s;  m32 = mri->ysize * mri->y_s;  m33 = mri->zsize * mri->z_s;
      
      ci = mri->width;
      cj = mri->height;
      ck = mri->depth;
      
      m14 = - (m11 * ci + m12 * cj + m13 * ck)/2.0;
      m24 = - (m21 * ci + m22 * cj + m23 * ck)/2.0;
      m34 = - (m31 * ci + m32 * cj + m33 * ck)/2.0;
      
      /*E* mri.c's extract_i_to_r() had these lines instead:
	ci = (mri->width - 1.0) / 2.0;
	cj = (mri->height - 1.0) / 2.0;
	ck = (mri->depth - 1.0) / 2.0;
	
	m14 = mri->c_r - (m11 * ci + m12 * cj + m13 * ck);
	m24 = mri->c_a - (m21 * ci + m22 * cj + m23 * ck);
	m34 = mri->c_s - (m31 * ci + m32 * cj + m33 * ck);
	which don't seem to be the right thing, here. */
      
    }
  else if(mri->slice_direction == MRI_CORONAL)
    {
      m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         m14 =  mri->xsize * mri->width / 2.0;
      m21 =  0.0;         m22 =  0.0;         m23 = mri->zsize;  m24 = -mri->zsize * mri->depth / 2.0;
      m31 =  0.0;         m32 = -mri->ysize;  m33 = 0.0;         m34 =  mri->ysize * mri->height / 2.0;
    }
  else if(mri->slice_direction == MRI_SAGITTAL)
    {
      m11 =  0.0;         m12 =  0.0;         m13 = mri->zsize;  m14 = -mri->zsize * mri->depth / 2.0;
      m21 =  mri->xsize;  m22 =  0.0;         m23 = 0.0;         m24 = -mri->xsize * mri->width / 2.0;
      m31 =  0.0;         m32 = -mri->ysize;  m33 = 0.0;         m34 =  mri->ysize * mri->height / 2.0;
    }
  else if(mri->slice_direction == MRI_HORIZONTAL)
    {
      m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         m14 =  mri->xsize * mri->width / 2.0;
      m21 =  0.0;         m22 = -mri->ysize;  m23 = 0.0;         m24 =  mri->ysize * mri->height / 2.0;
      m31 =  0.0;         m32 =  0.0;         m33 = mri->zsize;  m34 = -mri->zsize * mri->depth / 2.0;
    }
  else
    {
      m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         
      m14 =  mri->xsize * mri->width / 2.0;
      m21 =  0.0;         m22 = -mri->ysize;  m23 = 0.0;         
      m24 =  mri->ysize * mri->height / 2.0;
      m31 =  0.0;         m32 =  0.0;         m33 = mri->zsize;  
      m34 = -mri->zsize * mri->depth / 2.0;
    }
  
  stuff_four_by_four(m, m11, m12, m13, m14, 
		     m21, m22, m23, m24, 
		     m31, m32, m33, m34, 
		     0.0, 0.0, 0.0, 1.0);
#ifdef _VLDT_DEBUG
  /*E*/  fprintf(stderr, "Volm_MRIgetVoxelToRasXform(): m = \n");
  /*E*/  MatrixPrint(stderr, m) ;
#endif
  
  return(m);
  
} /* end Volm_MRIgetVoxelToRasXform() */
#endif



Volm_tErr Volm_UnloadDisplayTransform ( mriVolumeRef this ) {
  
  Volm_tErr eResult   = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_UnloadDisplayTransform( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* if we already have a transform, delete it */
  if( NULL != this->mDisplayTransform )
    Trns_Delete( &(this->mDisplayTransform) );
  
  MatrixCopy(this->m_resample_orig, this->m_resample) ;
  DebugNote( ("Calculating inverse of resample matrix") );
  MatrixInverse( this->m_resample, this->m_resample_inv );
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}



void Volm_GetColorAtIdx ( mriVolumeRef this,
			  xVoxelRef    iIdx,
			  xColor3fRef  oColor ) {
  
  float  value = 0;
  int    colorIdx = 0;

# if 0 /* ??????? */
  /* transform idx to display transform */
  if( NULL != this->mDisplayTransform && 0) {
    Volm_ApplyDisplayTransform_( this, iIdx, &disp );
    if( Volm_VerifyIdx_( this, &disp ) == Volm_tErr_NoErr ) 
      {
	Volm_GetSincNormValueAtIdx_( this, &disp, &rValue );
	value = (Volm_tValue)rValue;
      }
  }
#endif
  
  /* Get the value and normalize it to 0-255. */
  Volm_GetValueAtIdx_( this, iIdx, &value );
  colorIdx = (int) (255.0 * (value - this->min_val) / 
		    (this->max_val - this->min_val)) ;

  *oColor = this->maColorTable[colorIdx];
}

void Volm_GetColorAtXYSlice ( mriVolumeRef     this,
			      mri_tOrientation iOrientation,
			      xPoint2nRef      iPoint,
			      int              inSlice,
			      xColor3fRef      oColor ) {
  
  xVoxel   idx;

  Volm_ConvertXYSliceToIdx_( iOrientation, iPoint, inSlice, &idx );
  Volm_GetColorAtIdx( this, &idx, oColor );  
}

void Volm_GetMaxColorAtIdx ( mriVolumeRef     this,
			     xVoxelRef        iIdx,
			     mri_tOrientation iOrientation,
			     xColor3fRef      oColor ) {
  xPoint2n    point;
  int         nSlice;
  
  /* if we haven't find the max values yet, find them. */
  if( NULL == this->mpMaxValues ) {
    Volm_FindMaxValues( this );
  }
  
  /* Convert to an xy slice and get the color. */
  Volm_ConvertIdxToXYSlice_( iIdx, iOrientation, &point, &nSlice );
  Volm_GetMaxColorAtXYSlice( this, iOrientation, &point, nSlice, oColor );
}

void Volm_GetMaxColorAtXYSlice ( mriVolumeRef     this,
				 mri_tOrientation iOrientation,
				 xPoint2nRef      iPoint,
				 int              inSlice,
				 xColor3fRef      oColor ) {
  
  float value = 0;
  int   colorIdx = 0;
  
  /* if we haven't find the max values yet, find them. */
  if( NULL == this->mpMaxValues ) {
    Volm_FindMaxValues( this );
  }
  
  /* Get the value and normalize it to 0-255. Return this entry in the
     color table. */
  Volm_GetMaxValueAtXYSlice_( this, iOrientation, iPoint, &value );
  colorIdx = (int) (255.0 * (value - this->min_val) / 
		    (this->max_val - this->min_val)) ;
  *oColor = this->maColorTable[colorIdx];
}

Volm_tErr Volm_GetDimensions ( mriVolumeRef this,
			       int*         onDimensionX, int*  onDimensionY, int* onDimensionZ ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetDimensions( this=%p, onDimension=%p,%p,%p )",
		       this, onDimensionX, onDimensionY, onDimensionZ ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (onDimensionX != NULL) && (onDimensionY != NULL) &&
                     (onDimensionZ != NULL),
                     eResult, Volm_tErr_InvalidParamater );
  
  /* return the dimension */
  DebugNote( ("Returning the dimension") );
  *onDimensionX = this->mnDimensionX;
  *onDimensionY = this->mnDimensionY;
  *onDimensionZ = this->mnDimensionZ;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetNumberOfFrames ( mriVolumeRef this, 
				   int* onDimensionFrames ) {

  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetNumberOfFrames( this=%p, "
		       "onDimensionFrames=%p )", this, onDimensionFrames ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (onDimensionFrames != NULL),
                     eResult, Volm_tErr_InvalidParamater );
  
  /* return the dimension */
  DebugNote( ("Returning the dimension") );
  *onDimensionFrames = this->mnDimensionFrame;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetType ( mriVolumeRef this, int* onType ) {

  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetType( this=%p, onType=%p )", this, onType ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (onType != NULL), eResult, Volm_tErr_InvalidParamater );
  
  /* Make sure we have a volume */
  DebugAssertThrowX( ( NULL != this->mpMriValues ), 
		     eResult, Volm_tErr_MRIVolumeNotPresent );

  /* return the dimension */
  DebugNote( ("Returning the type") );
  *onType = this->mpMriValues->type;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetValueAtIdx ( mriVolumeRef this,
			       xVoxelRef    iIdx,
			       float*       oValue ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetValueAtIdx( this=%p, iIdx=%p, "
		       "oValue=%p )", this, iIdx, oValue ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oValue != NULL),
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* return the value */
  DebugNote( ("Fetching the value at (%d,%d,%d) and returning it",
	      xVoxl_ExpandInt( iIdx )));
  Volm_GetValueAtIdx_( this, iIdx, oValue );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetValueAtIdxUnsafe ( mriVolumeRef this,
				     xVoxelRef    iIdx,
				     float*       oValue ) {
  
  Volm_GetValueAtIdx_( this, iIdx, oValue );

  return Volm_tErr_NoErr;
}

Volm_tErr Volm_SetValueAtIdx ( mriVolumeRef this,
			       xVoxelRef    iIdx,
			       float        iValue ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_SetValueAtIdx( this=%p, iIdx=%p, "
		       "iValue=%d )", this, iIdx, (int)iValue ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL), eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* set the value */
  DebugNote( ("Setting the norm value") );
  Volm_SetValueAtIdx_( this, iIdx, iValue );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetValueAtIdxFrame  ( mriVolumeRef this,
				     xVoxelRef    iIdx,
				     int          iFrame,
				     float*       oValue ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetValueAtIdxFrame( this=%p, iIdx=%p, iFrame=%d,"
		       "oValue=%p )", this, iIdx, iFrame, oValue ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oValue != NULL),
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  eResult = Volm_VerifyFrame_( this, iFrame );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
    
  /* return the value */
  DebugNote( ("Fetching the value at (%d,%d,%d) frame %d and returning it",
	      xVoxl_ExpandInt( iIdx ), iFrame));
  Volm_GetValueAtIdxFrame_( this, iIdx, iFrame, oValue );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetValueAtIdxFrameUnsafe ( mriVolumeRef this,
					  xVoxelRef    iIdx,
					  int          iFrame,
					  float*       oValue ) {
  
  Volm_GetValueAtIdxFrame_( this, iIdx, iFrame, oValue );

  return Volm_tErr_NoErr;
}

Volm_tErr Volm_SetValueAtIdxFrame ( mriVolumeRef this,
				    xVoxelRef    iIdx,
				    int          iFrame,
				    float        iValue ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_SetValueAtIdx( this=%p, iIdx=%p, "
		       "iValue=%d )", this, iIdx, (int)iValue ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL), eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  eResult = Volm_VerifyFrame_( this, iFrame );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* set the value */
  DebugNote( ("Setting the value") );
  Volm_SetValueAtIdxFrame_( this, iIdx, iFrame, iValue );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertIdxToRAS ( mriVolumeRef this,
				 xVoxelRef    iIdx,
				 xVoxelRef    oRAS ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;
  Real      rasX    = 0;
  Real      rasY    = 0;
  Real      rasZ    = 0;
  

  DebugEnterFunction( ("Volm_ConvertIdxToRAS( this=%p, iIdx=%p, oRAS=%p )", 
		       this, iIdx, oRAS) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oRAS != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* First convert the incoming screen idx to our local MRI idx. */
  DebugNote( ("Converting scresshen idx to MRI idx") );
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &mriIdx );

  /* convert idx to ras */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to RAS with MRIvoxelToWorld",
	      xVoxl_ExpandFloat( &mriIdx )) );
  MRIvoxelToWorld( this->mpMriValues, xVoxl_ExpandFloat( &mriIdx ),
		   &rasX, &rasY, &rasZ );
  //E/ rewriting this - it checked slice_direction before ras_good_flag
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oRAS, (float)rasX, (float)rasY, (float)rasZ );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertRASToIdx ( mriVolumeRef this,
				 xVoxelRef    iRAS,
				 xVoxelRef    oIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  Real      idxX    = 0;
  Real      idxY    = 0;
  Real      idxZ    = 0;
  xVoxel    mriIdx;
  
  DebugEnterFunction( ("Volm_ConvertRASToIdx( this=%p, iRAS=%p, oIdx=%p )", 
		       this, iRAS, oIdx) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iRAS != NULL && oIdx != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  
  /* convert ras to MRI idx */
  DebugNote( ("Converting RAS (%.2f, %.2f, %.2f) to idx with MRIworldToVoxel",
	      xVoxl_ExpandFloat( iRAS )) );
  MRIworldToVoxel( this->mpMriValues, xVoxl_ExpandFloat( iRAS ),
		   &idxX, &idxY, &idxZ );
  
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( &mriIdx, (float)idxX, (float)idxY, (float)idxZ );

  /* Then convert from local MRI idx to screen idx. */
  DebugNote( ("Converting MRI idx to screen idx") );
  Volm_ConvertMRIIdxToScreenIdx_( this, &mriIdx, oIdx );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertIdxToMNITal ( mriVolumeRef this,
				    xVoxelRef    iIdx,
				    xVoxelRef    oMNITal ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;
  Real      talX    = 0;
  Real      talY    = 0;
  Real      talZ    = 0;
  
  DebugEnterFunction( ("Volm_ConvertIdxToMNITal( this=%p, iIdx=%p, "
		       "oMNITal=%p )", this, iIdx, oMNITal) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
#if 0  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oMNITal != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
#endif
  
  /* First convert the incoming screen idx to our local MRI idx. */
  DebugNote( ("Converting screen idx to MRI idx") );
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &mriIdx );

  /* convert idx to tal */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to mni tal with "
	      "MRIvoxelToTalairachVoxel", xVoxl_ExpandFloat( &mriIdx )) );
  MRIvoxelToTalairach( this->mpMriValues, xVoxl_ExpandFloat( &mriIdx ),
		       &talX, &talY, &talZ );
  //E/ written for coronal only - might be confusing for .mgh's (they
  //look coronal (slice_direction is not set, so sometimes 0) but
  //aren't)
  
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oMNITal, (float)talX, (float)talY, (float)talZ );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertIdxToTal ( mriVolumeRef this,
				 xVoxelRef    iIdx,
				 xVoxelRef    oTal ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;
  Real      talX    = 0;
  Real      talY    = 0;
  Real      talZ    = 0;
  xVoxel    tal;
  
  DebugEnterFunction( ("Volm_ConvertIdxToTal( this=%p, iIdx=%p, "
		       "oTal=%p )", this, iIdx, oTal) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
#if 0  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oTal != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
#endif
  
  /* First convert the incoming screen idx to our local MRI idx. */
  DebugNote( ("Converting screen idx to MRI idx") );
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &mriIdx );

  /* convert idx to tal */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to tal with "
	      "MRIvoxelToTalairach", xVoxl_ExpandFloat( &mriIdx )) );
  MRIvoxelToTalairach( this->mpMriValues, xVoxl_ExpandFloat( &mriIdx ),
		       &talX, &talY, &talZ );
  
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( &tal, (float)talX, (float)talY, (float)talZ );
  
  /* convert to real tal. switch on the z to see which transform to use. */
  if( xVoxl_GetFloatZ( &tal ) > 0 ) {
    DebugNote( ("Converting to real tal with >0 transform") );
    Trns_ConvertAtoB( this->mMNITalGtzToRealTalTransform, &tal, oTal );
  } else {
    DebugNote( ("Converting to real tal with <0 transform") );
    Trns_ConvertAtoB( this->mMNITalLtzToRealTalTransform, &tal, oTal );
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertTalToIdx ( mriVolumeRef this,
				 xVoxelRef    iTal,
				 xVoxelRef    oIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  Real      idxX    = 0;
  Real      idxY    = 0;
  Real      idxZ    = 0;
  xVoxel    mniTal;
  int       eMRI    = 0;
  
  DebugEnterFunction( ("Volm_ConvertTalToIdx( this=%p, iTal=%p, "
		       "oIdx=%p )", this, iTal, oIdx) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iTal != NULL && oIdx != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  
  /* convert to mni tal. switch on the z to see which transform to use. */
  if( xVoxl_GetFloatZ( iTal ) > 0 ) {
    DebugNote( ("Converting to mni tal with >0 transform") );
    Trns_ConvertBtoA( this->mMNITalGtzToRealTalTransform, iTal, &mniTal );
  } else {
    DebugNote( ("Converting to mni tal with <=0 transform") );
    Trns_ConvertBtoA( this->mMNITalLtzToRealTalTransform, iTal, &mniTal );
  }
  
  /* convert tal to idx */
  DebugNote( ("Converting tal (%.2f, %.2f, %.2f) to idx with "
	      "MRItalairachToVoxel", xVoxl_ExpandFloat( &mniTal )) );
  eMRI = MRItalairachToVoxel( this->mpMriValues, xVoxl_ExpandFloat( &mniTal ),
			      &idxX, &idxY, &idxZ );
  DebugAssertThrowX( (eMRI == NO_ERROR),eResult, Volm_tErr_CouldntReadVolume );
  
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oIdx, (float)idxX, (float)idxY, (float)idxZ );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertIdxToScanner ( mriVolumeRef this,
				     xVoxelRef    iIdx,
				     xVoxelRef    oScanner ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;
  
  DebugEnterFunction( ("Volm_ConvertIdxToScanner( this=%p, iIdx=%p, "
		       "oScanner=%p )", this, iIdx, oScanner) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
#if 0  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oScanner != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
#endif
  
  /* make sure we have the scanner transform */
  DebugAssertThrowX( (NULL != this->mScannerTransform),
		     eResult, Volm_tErr_ScannerTransformNotPresent );
  
  /* First convert the incoming screen idx to our local MRI idx. */
  DebugNote( ("Converting screen idx to MRI idx") );
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &mriIdx );

  /* convert idx to scanner */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to scanner", 
	      xVoxl_ExpandFloat( &mriIdx )) );
  Trns_ConvertAtoB( this->mScannerTransform, &mriIdx, oScanner );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ConvertIdxToMRIIdx ( mriVolumeRef this,
				    xVoxelRef    iIdx,
				    xVoxelRef    oMRIIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_ConvertIdxToMRIIdx( this=%p, iIdx=%p, "
		       "oMRIIdx=%p )", this, iIdx, oMRIIdx) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
#if 0  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oMRIIdx != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
#endif
  
  /* Convert the incoming screen idx to our local MRI idx. */
  DebugNote( ("Converting screen idx to MRI idx") );
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, oMRIIdx );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_GetIdxToRASTransform ( mriVolumeRef     this,
				      mriTransformRef* opTransform ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_GetIdxToRASTransform( this=%p, opTransform=%p )",
		       this, opTransform) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != opTransform),
		     eResult, Volm_tErr_InvalidParamater );
  
  /* make sure we have the transform */
  DebugAssertThrowX( (NULL != this->mIdxToRASTransform),
		     eResult, Volm_tErr_IdxToRASTransformNotPresent );
  
  /* return a pointer to it */
  *opTransform = this->mIdxToRASTransform;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_VisitAllVoxels ( mriVolumeRef        this,
				Volm_tVisitFunction iFunc,
				void*               ipData ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  Volm_tVisitCommand eVisit = Volm_tVisitComm_Continue;
  tBoolean bGo = TRUE;
  xVoxel idx;
  float value = 0;
  
  DebugEnterFunction( ("Volm_tVisitFunction( this=%p, iFunc=%p, ipData=%p )", 
		       this, iFunc, ipData) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iFunc), eResult, Volm_tErr_InvalidParamater );
  
  /* go thru the volume. get the value and if it matches our index,
     add it to the selection. */
  xVoxl_Set( &idx, 0, 0, 0 );
  bGo = TRUE;
  while( bGo ) {
    
    Volm_GetValueAtIdx_( this, &idx, &value );
    eVisit = iFunc( &idx, value, ipData );
    
    /* look at the result code they gave us. see if they want us to continue
       or skip a row. if skip, increment the voxel coordinate to the end. the
       next call to Increment will pop it over into to the next 
       row or plane. */
    switch( eVisit ) {
    case Volm_tVisitComm_Continue:
      bGo = xVoxl_IncrementUntilLimits( &idx, this->mnDimensionX-1, 
					this->mnDimensionY-1,
                                        this->mnDimensionZ-1);
      continue;
      break;
    case Volm_tVisitComm_SkipRestOfRow:
      idx.mfX = this->mnDimensionX-1;
      bGo = xVoxl_IncrementUntilLimits( &idx, this->mnDimensionX-1,
                                        this->mnDimensionY-1,
					this->mnDimensionZ-1);
      break;
    case Volm_tVisitComm_SkipRestOfPlane:
      idx.mfX = this->mnDimensionX-1;
      idx.mfY = this->mnDimensionY-1;
      bGo = xVoxl_IncrementUntilLimits( &idx, this->mnDimensionX-1,
                                        this->mnDimensionY-1,
					this->mnDimensionZ-1);
      break;
    case Volm_tVisitComm_Stop:
      bGo = FALSE;
      break;
    default:
      break;
    }
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_FindMaxValues ( mriVolumeRef this ) {
  
  Volm_tErr        eResult     = Volm_tErr_NoErr;
  mri_tOrientation orientation = mri_tOrientation_None;
  xPoint2n         point       = {0, 0};
  int              nSlice      = 0;
  float            value       = 0;
  float            maxValue    = 0;
  int              nSize       = 0;
#ifdef Solaris
  int              i           = 0;
#endif
  
  DebugEnterFunction( ("Volm_FindMaxValues( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* make sure we have our normalized volume */
  DebugAssertThrowX( (NULL != this->mpMriValues), 
		     eResult, Volm_tErr_MRIVolumeNotPresent );
  
  /* we need 3 slices worth of data. if we have some allocated, delete them,
     and reallocate them. */
  if( NULL != this->mpMaxValues )
    free( this->mpMaxValues );
  nSize = 
    MAX(this->mnDimensionX, this->mnDimensionY) * 
    MAX(this->mnDimensionX, this->mnDimensionZ) *
    mri_knNumOrientations ;
  this->mpMaxValues = (float*) malloc( nSize * sizeof(float) );
  DebugAssertThrowX( (NULL != this->mpMaxValues),
		     eResult, Volm_tErr_AllocationFailed );

#ifdef Solaris
  DebugNote( ("Zeroing visited volume iteratively") );
  for( i = 0; i < nSize; i++ )
    this->mpMaxValues[i] = 0;
#else
  DebugNote( ("Zeroing visted volume with bzero") );
  bzero( this->mpMaxValues, nSize );
#endif
  
  /* go through the max value buffer. for each orientation */
  for( orientation = 0;
       orientation < mri_knNumOrientations; 
       orientation++ ) {
    
    /* for every point in... */
    for( point.mnY = 0; point.mnY < this->mnDimensionY; point.mnY++ ) {
      for( point.mnX = 0; point.mnX < this->mnDimensionX; point.mnX++ ) {
        
        /* now go thru every slice in the volume at this point.
           for every slice... */
        for( nSlice = 0; nSlice < this->mnDimensionZ; nSlice++ ) {
          
          /* get the value and max value at this point. if the value
	     is greater, make it the new max value. */
          Volm_GetValueAtXYSlice_( this, orientation, &point, nSlice, &value );
          Volm_GetMaxValueAtXYSlice_( this, orientation, &point, &maxValue );
          if( value > maxValue ) {
            Volm_SetMaxValueAtXYSlice_( this, orientation, &point, value ); 
	  }
        }
      }
    }
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_MakeColorTable ( mriVolumeRef this ) {
  
  Volm_tErr   eResult    = Volm_tErr_NoErr;
  float       thresh     = 0;
  float       squash     = 0;
  float       scale      = 0;
  int         value      = 0;
  float       fValue     = 0;
  float       fComponent = 0;
  
  DebugEnterFunction( ("Volm_MakeColorTable( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* for each value... */
  scale = Volm_knMaxValue;
  thresh = this->mfColorBrightness;
  squash = this->mfColorContrast;
  for( value = 0; value < Volm_knNumValues; value++ ) {
    
    /* calculate the color */
    fValue = (float)value;
    fComponent = 
      (1.0 / (1.0 + exp(-squash * ((fValue/scale) - thresh)))) * scale + 0.5;
    if( fComponent > Volm_knMaxValue )
      fComponent = Volm_knMaxValue;
    
    /* normalize to 0 - 1 */
    fComponent = (float)fComponent / (float)Volm_knMaxValue;
    
    /* set the color */
    this->maColorTable[value].mfRed   = fComponent;
    this->maColorTable[value].mfGreen = fComponent;
    this->maColorTable[value].mfBlue  = fComponent;
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


Volm_tErr Volm_SetBrightnessAndContrast ( mriVolumeRef this,
					  float        ifBrightness,
					  float        ifContrast ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_SetBrightnessAndContrast( this=%p, "
		       "ifBrightness=%.2f, ifContrast=%.2f )", this,
		       ifBrightness, ifContrast) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
  
  /* set the values */
  this->mfColorBrightness = ifBrightness;
  this->mfColorContrast   = ifContrast;
  
  /* recalc the color table */
  eResult = Volm_MakeColorTable( this );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_SaveToSnapshot ( mriVolumeRef this ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  int       nSize   = 0;
  
  DebugEnterFunction( ("Volm_SaveToSnapshot( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* if we have snapshot data, delete it */
  if( NULL != this->mpSnapshot )
    free( this->mpSnapshot );
  
  /* make a buffer the same size as the slices buffer. */
  nSize = this->mpMriValues->width * 
    this->mpMriValues->height * 
    this->mpMriValues->depth ;
  switch (this->mpMriValues->type)
    {
    case MRI_SHORT:
      nSize *= sizeof(short) ;
      break ;
    case MRI_FLOAT:
      nSize *= sizeof(float) ;
      break ;
    case MRI_LONG:
      nSize *= sizeof(long) ;
      break ;
    case MRI_INT:
      nSize *= sizeof(int) ;
      break ;
    default:
      break ;
    }
  this->mpSnapshot = (float*) malloc( nSize );
  DebugAssertThrowX( (NULL != this->mpSnapshot),
		     eResult, Volm_tErr_AllocationFailed );
  
  /* copy it in */
  memcpy( this->mpSnapshot, this->mpMriValues->slices, nSize );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_RestoreFromSnapshot ( mriVolumeRef this ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  int       nSize   = 0;
  
  DebugEnterFunction( ("Volm_RestoreFromSnapshot( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* if we don't have a snapshot, return. */
  if( NULL == this->mpSnapshot )
    DebugGotoCleanup;
  
  /* copy the snapshot data into the norm volume buffer data */
  nSize = this->mpMriValues->width * 
    this->mpMriValues->height * 
    this->mpMriValues->depth ;
  switch (this->mpMriValues->type)
    {
    case MRI_SHORT:
      nSize *= sizeof(short) ;
      break ;
    case MRI_FLOAT:
      nSize *= sizeof(float) ;
      break ;
    case MRI_LONG:
      nSize *= sizeof(long) ;
      break ;
    case MRI_INT:
      nSize *= sizeof(int) ;
      break ;
    default:
      break ;
    }
  memcpy( this->mpMriValues->slices, this->mpSnapshot, nSize );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


Volm_tErr Volm_Rotate ( mriVolumeRef     this,
			mri_tOrientation iAxis,
			float            ifDegrees ) {
  
  Volm_tErr eResult  = Volm_tErr_NoErr;
  MRI*      newNorm  = NULL;
  float     fRadians = 0;
  
  DebugEnterFunction( ("Volm_Rotate( this=%p, iAxis=%d, ifDegrees=%.2f )", 
		       this, iAxis, ifDegrees) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* convert to radians */
  fRadians = (ifDegrees / 180.0) * M_PI;
  
  /* switch on the axis and pass to the approriate mri function */
  switch( iAxis ) {
  case mri_tOrientation_Coronal:
    DebugNote( ("Rotating normalized volume with MRIrotateY") );
    newNorm = MRIrotateY( this->mpMriValues, NULL, fRadians );
    break;
  case mri_tOrientation_Horizontal:
    DebugNote( ("Rotating normalized volume with MRIrotateZ") );
    newNorm = MRIrotateZ( this->mpMriValues, NULL, fRadians );
    break;
  case mri_tOrientation_Sagittal:
    DebugNote( ("Rotating normalized volume with MRIrotateX") );
    newNorm = MRIrotateX( this->mpMriValues, NULL, fRadians );
    break;
  default:
    DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
    break;
  }    
  
  /* reassign to the new volume */
  DebugNote( ("Freeing old normalized volume and assigned rotated one") );
  MRIfree( &this->mpMriValues );
  this->mpMriValues = newNorm;
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_Threshold ( mriVolumeRef this,
			   float        iThreshold,
			   tBoolean     ibAbove,
			   float        iNewValue ) {
  
  Volm_tErr   eResult = Volm_tErr_NoErr;
  float        value   = 0;
  xVoxel      idx;
  
  DebugEnterFunction( ("Volm_Threshold( this=%p, iThreshold=%d, ibAbove=%d, "
		       "iNewValue=%d )", this, (int)iThreshold, (int)ibAbove,
		       (int)iNewValue) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* step through the volume... */
  xVoxl_Set( &idx, 0, 0, 0 );
  DebugNote( ("Thresholding normalized volume") );
  while( xVoxl_IncrementUntilLimits( &idx, this->mnDimensionX-1,
                                     this->mnDimensionY-1, 
				     this->mnDimensionZ-1) ) {
    
    Volm_GetValueAtIdx_( this, &idx, &value );
    
    /* if we're going above and this value is above the thresh, or if we're
       going below and this value is below the thresh...*/
    if( (ibAbove && value > iThreshold) ||
	(!ibAbove && value < iThreshold ) ) {
      
      /* set the value to the new level. */
      Volm_SetValueAtIdx_( this, &idx, iNewValue );
    }
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


Volm_tErr Volm_Flip ( mriVolumeRef     this,
		      mri_tOrientation iAxis ) {
  
  Volm_tErr   eResult       = Volm_tErr_NoErr;
  float       value         = 0;
  float       flippedValue  = 0;
  xVoxel      idx;
  xVoxel      flippedIdx;
  int         nXLimit       = 0;
  int         nYLimit       = 0;
  int         nZLimit       = 0;
  
  DebugEnterFunction( ("Volm_Flip( this=%p, iAxis=%d )", this, (int)iAxis) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* get our limits */
  nXLimit = this->mnDimensionX-1 ;
  nYLimit = this->mnDimensionY-1 ;
  nZLimit = this->mnDimensionZ-1;
  switch( iAxis ) {
  case mri_tOrientation_Coronal:
    nZLimit = (this->mnDimensionZ-1) / 2;
    break;
  case mri_tOrientation_Sagittal:
    nXLimit = (this->mnDimensionX-1) / 2;
    break;
  case mri_tOrientation_Horizontal:
    nYLimit = (this->mnDimensionY-1) / 2;
    break;
  default:
    DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
    break;
  }
  
  /* step through the volume... */
  xVoxl_Set( &idx, 0, 0, 0 );
  DebugNote( ("Flipping normalized volume") );
  while( xVoxl_IncrementUntilLimits( &idx, nXLimit, nYLimit, nZLimit )) {
    
    switch( iAxis ) {
    case mri_tOrientation_Coronal:
      xVoxl_Set( &flippedIdx,
                 xVoxl_GetX( &idx ), 
                 xVoxl_GetY( &idx ), 
                 this->mnDimensionZ-1 - xVoxl_GetZ( &idx ) );
      break;
    case mri_tOrientation_Sagittal:
      xVoxl_Set( &flippedIdx,
                 this->mnDimensionX-1 - xVoxl_GetX( &idx ), 
                 xVoxl_GetY( &idx ), 
                 xVoxl_GetZ( &idx ) );
      break;
    case mri_tOrientation_Horizontal:
      xVoxl_Set( &flippedIdx,
                 xVoxl_GetX( &idx ), 
                 this->mnDimensionY-1 - xVoxl_GetY( &idx ), 
                 xVoxl_GetZ( &idx ) );
      break;
    default:
      DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
      break;
    }
    
    /* swap the values */
    Volm_GetValueAtIdx_( this, &idx, &value );
    Volm_GetValueAtIdx_( this, &flippedIdx, &flippedValue );
    Volm_SetValueAtIdx_( this, &idx, flippedValue );
    Volm_SetValueAtIdx_( this, &flippedIdx, value );
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_CopySubjectName ( mriVolumeRef this,
				 char*        oSubjectName,
				 int          inDestLen ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_CopySubjectName( this=%p, oSubjectName=%s, "
		       "inDestLen=%d )", this, oSubjectName, inDestLen) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oSubjectName != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  
  DebugNote( ("Copying subject name") );
  xUtil_strncpy( oSubjectName, this->msSubjectName, inDestLen );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_CopyVolumeName ( mriVolumeRef this,
				char*        oVolumeName,
				int          inDestLen ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_CopyVolumeName( this=%p, oVolumeName=%s, "
		       "inDestLen=%d )", this, oVolumeName, inDestLen) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oVolumeName != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  
  DebugNote( ("Copying volume name") );
  xUtil_strncpy( oVolumeName, this->msVolumeName, inDestLen );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_CopySourceDir ( mriVolumeRef this,
			       char*        oSourceDir,
			       int          inDestLen ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_CopySourceDir( this=%p, oSourceDir=%s, "
		       "inDestLen=%d )", this, oSourceDir, inDestLen) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oSourceDir != NULL), 
		     eResult, Volm_tErr_InvalidParamater );
  
  DebugNote( ("Copying source dir") );
  xUtil_strncpy( oSourceDir, this->msOriginalPath, inDestLen );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ ( mriVolumeRef this ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_( this=%p )", this) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_SetSubjectName ( mriVolumeRef this,
				char*        isName ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_SetSubjectName( this=%p, isName=%s )", 
		       this, isName) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (isName != NULL), eResult, Volm_tErr_InvalidParamater );
  
  /* copy the name */
  xUtil_strncpy( this->msSubjectName, isName, mri_knSubjectNameLen );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_SetVolumeName ( mriVolumeRef this,
			       char*        isName ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_SetVolumeName( this=%p, isName=%s )", 
		       this, isName) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (isName != NULL), eResult, Volm_tErr_InvalidParamater );
  
  /* copy the name */
  xUtil_strncpy( this->msVolumeName, isName, mri_knSubjectNameLen );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}



Volm_tErr Volm_ExtractAndSetSubjectName ( mriVolumeRef this,
					  char*        isSource ) {
  
  Volm_tErr eResult    = Volm_tErr_NoErr;
  int       nChar      = 0;
  int       nWordChar  = 0;
  char*     sWord      = NULL;
  char      sName[256] = "";
  int       nLastSlash = 0;
  
  DebugEnterFunction( ("Volm_ExtractAndSetSubjectName( this=%p, isSource=%s)",
		       this, isSource) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isSource), eResult, Volm_tErr_InvalidParamater );
  
  /* look for 'subjects' in the title */
  DebugNote( ("Looking for subjects/ in source") );
  sWord = strstr( isSource, "subjects/" );
  if( NULL != sWord ) {
    
    /* we're at the s in subjects now. scoot ahead to the slash. */
    DebugNote( ("Scooting ahead to the slash") );
    nChar = 0;
    while( sWord[nChar] != '/' &&
	   sWord[nChar] != '\0' ) {
      nChar++;
    }
    
    /* if found, use the next part as the name */
    DebugNote( ("Copying in chars from subjects/ to the next /") );
    nWordChar = 0;
    nChar++; /* get past the slash */
    while( sWord[nChar] != '/' &&
	   sWord[nChar] != '\0' ) {
      sName[nWordChar] = sWord[nChar];
      nWordChar++;
      nChar++;
    }
    
  } else {
    
    /* look for /mri in the title */
    DebugNote( ("Looking for /mri in source") );
    sWord = strstr( isSource, "/mri" );
    if( NULL != sWord ) {
      
      /* we're at the slash now. go to the slash before this one. while
	 we're not at the beginning.. */
      DebugNote( ("Going to the slash before /mri") );
      while( sWord != isSource ) {
	sWord -= sizeof( char );
	if( *sWord == '/' )
	  break;
      }
      
      /* inc past the slash and use the next part as the name */
      sWord += sizeof( char );
      nWordChar = 0;
      DebugNote( ("Copying in the part after the slash before /mri to /mri") );
      while( *sWord != '/' &&
	     *sWord != '\0' ) {
	sName[nWordChar] = *sWord;
	nWordChar++;
	sWord += sizeof( char );
      }
      
    } else {
      
      /* else just use the last part */
      nChar = 0;
      DebugNote( ("Scooting to last slash") );
      while( isSource[nChar] != '\0' ) {
	if( isSource[nChar] == '/' )
	  nLastSlash = nChar;
	nChar++;
      }
      
      /* if we got it, use it, else use the whole source. */
      if( isSource[nLastSlash] == '/' ) {
	DebugNote( ("Copying everything from last slash to end of source") );
	xUtil_strncpy( sName, &(isSource[nLastSlash+1]), sizeof(sName) );
      } else {
	DebugNote( ("Copying entire source") );
	xUtil_strncpy( sName, isSource, sizeof(sName) );
      }
    }
  }
  
  /* set the subject name */
  DebugNote( ("Copying result into subject name") );
  xUtil_strncpy( this->msSubjectName, sName, sizeof(this->msSubjectName) );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_ExtractAndSetVolumeName ( mriVolumeRef this,
					 char*        isSource ) {
  
  Volm_tErr eResult    = Volm_tErr_NoErr;
  int       nChar      = 0;
  int       nLastSlash = 0;
  
  DebugEnterFunction( ("Volm_ExtractAndSetVolumeName( this=%p, isSource=%s )",
		       this, isSource ) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isSource), eResult, Volm_tErr_InvalidParamater );
  
  /* look for the last / */
  DebugNote( ("Looking for last slash") );
  while( isSource[nChar] != '\0' ) {
    if( isSource[nChar] == '/' )
      nLastSlash = nChar;
    nChar++;
  }
  
  /* if we got one, everything from there+1 into the subjectname, 
     otherwise just use the whole name */
  if( isSource[nLastSlash] == '/' ) {
    DebugNote( ("Copying in source name from char %d", nLastSlash+1) );
    xUtil_strncpy( this->msVolumeName, &(isSource[nLastSlash+1]),
		   sizeof( this->msVolumeName ) );
  } else {
    DebugNote( ("Copying in whole source name") );
    xUtil_strncpy( this->msVolumeName, isSource,sizeof( this->msVolumeName ) );
  }
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}


/* ======================================== BEGIN FUNCTION VERION OF MACROS */
#ifndef VOLM_USE_MACROS

void Volm_ConvertScreenIdxToMRIIdx_ ( mriVolumeRef this,
				      xVoxelRef    iScreenIdx,
				      xVoxelRef    oMRIIdx ) {
  
  /* Stuff the screen index. */
  DebugNote( ("Stuffing iScreenIdx into vector") );
  *MATRIX_RELT(this->mpTmpScreenIdx,1,1) = xVoxl_GetFloatX(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,2,1) = xVoxl_GetFloatY(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,3,1) = xVoxl_GetFloatZ(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,4,1) = 1.0 ;

  /* Multiple by the resample matrix. */
  DebugNote( ("Multiplying screen index by resample matrix") );
  MatrixMultiply( this->m_resample, this->mpTmpScreenIdx, this->mpTmpMRIIdx );

  /* Stuff the outgoing voxel */
  DebugNote( ("Stuffing result vector into voxel") );
  xVoxl_SetFloat( oMRIIdx, 
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,1),
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,2), 
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,3) );

  DebugNote( ("Checking mri idx") );
  if( (xVoxl_GetX(oMRIIdx) < 0 ||
       xVoxl_GetX(oMRIIdx) >= this->mnDimensionX ||
       xVoxl_GetY(oMRIIdx) < 0 ||
       xVoxl_GetY(oMRIIdx) >= this->mnDimensionY ||
       xVoxl_GetZ(oMRIIdx) < 0 ||
       xVoxl_GetZ(oMRIIdx) >= this->mnDimensionZ) ) {

    xVoxl_Set( oMRIIdx, 0, 0, 0 );
  }
}

void Volm_ConvertMRIIdxToScreenIdx_ ( mriVolumeRef this,
				      xVoxelRef    iMRIIdx,
				      xVoxelRef    oScreenIdx ) {
  
  /* Stuff the mri index. */
  DebugNote( ("Stuffing iMRIIdx into vector") );
  *MATRIX_RELT(this->mpTmpMRIIdx,1,1) = xVoxl_GetFloatX(iMRIIdx);
  *MATRIX_RELT(this->mpTmpMRIIdx,2,1) = xVoxl_GetFloatY(iMRIIdx);
  *MATRIX_RELT(this->mpTmpMRIIdx,3,1) = xVoxl_GetFloatZ(iMRIIdx);
  *MATRIX_RELT(this->mpTmpMRIIdx,4,1) = 1.0 ;

  /* Multiple by the inverse resample matrix. */
  DebugNote( ("Multiplying mri index by inverse resample matrix") );
  MatrixMultiply( this->m_resample_inv, 
		  this->mpTmpMRIIdx, this->mpTmpScreenIdx );

  /* Stuff the outgoing voxel */
  DebugNote( ("Stuffing result vector into voxel") );
  xVoxl_SetFloat( oScreenIdx, 
		  *MATRIX_RELT(this->mpTmpScreenIdx,1,1),
		  *MATRIX_RELT(this->mpTmpScreenIdx,1,2), 
		  *MATRIX_RELT(this->mpTmpScreenIdx,1,3) );

  DebugNote( ("Checking screen idx") );
  if( (xVoxl_GetX(oScreenIdx) < 0 ||
       xVoxl_GetX(oScreenIdx) >= 256 ||
       xVoxl_GetY(oScreenIdx) < 0 ||
       xVoxl_GetY(oScreenIdx) >= 256 ||
       xVoxl_GetZ(oScreenIdx) < 0 ||
       xVoxl_GetZ(oScreenIdx) >= 256) ) {

    xVoxl_Set( oScreenIdx, 0, 0, 0 );
  }
}

void Volm_ConvertIdxToXYSlice_ ( xVoxelRef         iIdx,
				 mri_tOrientation  iOrientation,
				 xPoint2nRef       oPoint,
				 int*              onSlice ) {
  
  switch( iOrientation ) {
  case mri_tOrientation_Coronal:
    oPoint->mnX = xVoxl_GetX( iIdx );
    oPoint->mnY = xVoxl_GetY( iIdx );
    *onSlice    = xVoxl_GetZ( iIdx );
    break;
  case mri_tOrientation_Horizontal:
    oPoint->mnX = xVoxl_GetX( iIdx );
    oPoint->mnY = xVoxl_GetZ( iIdx );
    *onSlice    = xVoxl_GetY( iIdx );
    break;
  case mri_tOrientation_Sagittal:
    oPoint->mnX = xVoxl_GetZ( iIdx );
    oPoint->mnY = xVoxl_GetY( iIdx );
    *onSlice    = xVoxl_GetZ( iIdx );
    break;
  default:
    DebugPrintStack;
    DebugPrint( ("Volm_ConvertIdxToXYSlice_ called with invalid "
		 "orientation %d", iOrientation) );
    oPoint->mnX = oPoint->mnY = *onSlice = 0;
    break;
  }
}

void Volm_ConvertXYSliceToIdx_ ( mri_tOrientation  iOrientation,
				 xPoint2nRef       iPoint,
				 int               inSlice,
				 xVoxelRef         oIdx ) {
  
  switch( iOrientation ) {
  case mri_tOrientation_Coronal:
    xVoxl_Set( oIdx, iPoint->mnX, iPoint->mnY, inSlice );
    break;
  case mri_tOrientation_Horizontal:
    xVoxl_Set( oIdx, iPoint->mnX, inSlice, iPoint->mnY );
    break;
  case mri_tOrientation_Sagittal:
    xVoxl_Set( oIdx, inSlice, iPoint->mnY, iPoint->mnX );
    break;
  default:
    DebugPrintStack;
    DebugPrint( ("Volm_ConvertXYSliceToIdx_ called with invalid "
		 "orientation %d", iOrientation) );
    xVoxl_Set( oIdx, 0, 0, 0 );
    break;
  }
}

void Volm_GetValueAtXYSlice_ ( mriVolumeRef this,
			       mri_tOrientation  iOrientation,
			       xPoint2nRef       iPoint,
			       int               inSlice,
			       float*            oValue) {
  
  Volm_ConvertXYSliceToIdx_( iOrientation, iPoint, inSlice, &this->mTmpVoxel );
  Volm_GetValueAtIdx_( this, &this->mTmpVoxel, oValue );
}

void Volm_GetValueAtIdx_ ( mriVolumeRef this,
			   xVoxelRef    iIdx,
			   float*       oValue) {
  
  int x, y, z;

  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to get the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );
  
  x = xVoxl_GetX(&this->mTmpVoxel);
  y = xVoxl_GetY(&this->mTmpVoxel);
  z = xVoxl_GetZ(&this->mTmpVoxel);

  switch( this->mpMriValues->type ) {
    case MRI_UCHAR:
      *oValue = 
	MRIvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) );
      break;
    case MRI_INT:
      *oValue = 
	MRIIvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		 xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) );
      break;
    case MRI_LONG:
      *oValue = 
	MRILvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		 xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) );
      break;
    case MRI_FLOAT:
      *oValue = 
	MRIFvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		 xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) );
      break;
    case MRI_SHORT:
      *oValue = 
	MRISvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		 xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) );
      break;
    default:
      *oValue = 0;
      break ;
    }
}

void Volm_SetValueAtIdx_ ( mriVolumeRef     this,
			   xVoxelRef        iIdx,
			   float            iValue ) {
  
  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to set the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  switch (this->mpMriValues->type)
    {
    default:
      break ;
    case MRI_UCHAR:
      MRIvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
	      xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) ) = 
	(BUFTYPE) iValue;
      break ;
    case MRI_SHORT:
      MRISvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
	       xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) ) = 
	(short) iValue;
      break ;
    case MRI_FLOAT:
      MRIFvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
	       xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) ) = 
	(float) iValue;
      break ;
    case MRI_LONG:
      MRILvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
	       xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) ) = 
	(long) iValue;
      break ;
    case MRI_INT:
      MRIIvox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
	       xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel) ) = 
	(int) iValue;
      break ;
    }
}

void Volm_GetValueAtIdxFrame_ ( mriVolumeRef this,
				xVoxelRef    iIdx,
				int          iFrame,
				float*       oValue) {
  
  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to get the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  switch( this->mpMriValues->type ) {
    case MRI_UCHAR:
      *oValue = 
	MRIseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		    xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		    iFrame );
      break;
    case MRI_INT:
      *oValue = 
	MRIIseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		     xVoxl_GetY(&this->mTmpVoxel),xVoxl_GetZ(&this->mTmpVoxel),
		     iFrame );
      break;
    case MRI_LONG:
      *oValue = 
	MRILseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		     xVoxl_GetY(&this->mTmpVoxel),xVoxl_GetZ(&this->mTmpVoxel),
		     iFrame );
      break;
    case MRI_FLOAT:
      *oValue = 
	MRIFseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		     xVoxl_GetY(&this->mTmpVoxel),xVoxl_GetZ(&this->mTmpVoxel),
		     iFrame );
      break;
    case MRI_SHORT:
      *oValue = 
	MRISseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		     xVoxl_GetY(&this->mTmpVoxel),xVoxl_GetZ(&this->mTmpVoxel),
		     iFrame );
      break;
    default:
      *oValue = 0;
      break ;
    }
}

void Volm_SetValueAtIdxFrame_ ( mriVolumeRef     this,
				xVoxelRef        iIdx,
				int          iFrame,
				float            iValue ) {
  
  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to set the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  switch (this->mpMriValues->type)
    {
    default:
      break ;
    case MRI_UCHAR:
      MRIseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		  xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		  iFrame ) = 
	(BUFTYPE) iValue;
      break ;
    case MRI_SHORT:
      MRISseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		   xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		   iFrame ) = 
	(short) iValue;
      break ;
    case MRI_FLOAT:
      MRIFseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		   xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		   iFrame ) = 
	(float) iValue;
      break ;
    case MRI_LONG:
      MRILseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		   xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		   iFrame ) = 
	(long) iValue;
      break ;
    case MRI_INT:
      MRIIseq_vox( this->mpMriValues, xVoxl_GetX(&this->mTmpVoxel), 
		   xVoxl_GetY(&this->mTmpVoxel), xVoxl_GetZ(&this->mTmpVoxel),
		   iFrame ) = 
	(int) iValue;
      break ;
    }
}

void Volm_GetSincValueAtIdx_ ( mriVolumeRef     this,
			       xVoxelRef        iIdx,
			       Real*            irValue) {
  /* First convert to MRI index, then use the MRI access function to
     get the sinc value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  MRIsincSampleVolume( this->mpMriValues, 
		       xVoxl_GetFloatX(&this->mTmpVoxel),
		       xVoxl_GetFloatY(&this->mTmpVoxel),
		       xVoxl_GetFloatZ(&this->mTmpVoxel),
		       2, irValue);
}

void Volm_GetMaxValueAtXYSlice_ ( mriVolumeRef     this,
				  mri_tOrientation iOrientation, 
				  xPoint2nRef      iPoint,
				  float*           oValue ) {
  
  *oValue = 
    this->mpMaxValues[Volm_GetMaxValueIndex_(this,iOrientation,iPoint)];
}

void Volm_SetMaxValueAtXYSlice_ ( mriVolumeRef     this,
				  mri_tOrientation iOrientation, 
				  xPoint2nRef      iPoint,
				  float            iValue ) {
  
  this->mpMaxValues[Volm_GetMaxValueIndex_(this,iOrientation,iPoint)] = iValue;
}


int Volm_GetMaxValueIndex_ ( mriVolumeRef     this,
			     mri_tOrientation iOrientation, 
			     xPoint2nRef      iPoint ) {
  
  return (iOrientation * this->mnDimensionZ * this->mnDimensionY) +
    (iPoint->mnY * this->mnDimensionX) +
    iPoint->mnX;
}

void Volm_ApplyDisplayTransform_ ( mriVolumeRef     this,
				   xVoxelRef        iIdx,
				   xVoxelRef        oIdx ) {

  Trns_ConvertBtoA( this->mDisplayTransform, iIdx, oIdx );
}


#endif /* VOLM_USE_MACROS */
/* ======================================== END OF FUNCTION VERION OF MACROS */


Volm_tErr Volm_Verify ( mriVolumeRef this ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  DebugEnterFunction( ("Volm_Verify( this=%p )", this) );
  
  DebugAssertThrowX( (NULL != this), eResult, Volm_tErr_InvalidParamater );
  DebugAssertThrowX( (Volm_kSignature == this->mSignature),
		     eResult, Volm_tErr_InvalidSignature );
  
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  DebugPrintStack;
  EndDebugCatch;
  
  DebugExitFunction;
  
  return eResult;
}

Volm_tErr Volm_VerifyIdx ( mriVolumeRef this,
			   xVoxelRef    iIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  
  eResult = Volm_Verify( this );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult ) );
  
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult ) );
  
  DebugCatch;
  EndDebugCatch;
  
  return eResult;
}

Volm_tErr Volm_VerifyIdx_ ( mriVolumeRef this,
			    xVoxelRef    iIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel mriIdx;

  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &mriIdx );

  DebugAssertThrowX( (xVoxl_GetX(&mriIdx) >= 0 &&
		      xVoxl_GetRoundX(&mriIdx) < this->mnDimensionX &&
		      xVoxl_GetY(&mriIdx) >= 0 &&
		      xVoxl_GetRoundY(&mriIdx) < this->mnDimensionY &&
		      xVoxl_GetZ(&mriIdx) >= 0 &&
		      xVoxl_GetRoundZ(&mriIdx) < this->mnDimensionZ),
		     eResult, Volm_tErr_InvalidIdx );
  
  DebugCatch;
  EndDebugCatch;
  return eResult;
}

Volm_tErr Volm_VerifyIdxInMRIBounds ( mriVolumeRef this,
				      xVoxelRef    iScreenIdx ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;

  eResult = Volm_Verify( this );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult ) );
  
  /* Stuff the screen index. */
  DebugNote( ("Stuffing iScreenIdx into vector") );
  *MATRIX_RELT(this->mpTmpScreenIdx,1,1) = xVoxl_GetFloatX(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,2,1) = xVoxl_GetFloatY(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,3,1) = xVoxl_GetFloatZ(iScreenIdx);
  *MATRIX_RELT(this->mpTmpScreenIdx,4,1) = 1.0 ;

  /* Multiple by the resample matrix. */
  DebugNote( ("Multiplying screen index by resample matrix") );
  MatrixMultiply( this->m_resample, this->mpTmpScreenIdx, this->mpTmpMRIIdx );

  /* Stuff the outgoing voxel */
  DebugNote( ("Stuffing result vector into voxel") );
  xVoxl_SetFloat( &this->mTmpVoxel, 
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,1),
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,2), 
		  *MATRIX_RELT(this->mpTmpMRIIdx,1,3) );

  DebugNote( ("Checking mri idx") );
  DebugAssertThrowX((xVoxl_GetX(&this->mTmpVoxel) >= 0 &&
		     xVoxl_GetRoundX(&this->mTmpVoxel) < this->mnDimensionX &&
		     xVoxl_GetY(&this->mTmpVoxel) >= 0 &&
		     xVoxl_GetRoundY(&this->mTmpVoxel) < this->mnDimensionY &&
		     xVoxl_GetZ(&this->mTmpVoxel) >= 0 &&
		     xVoxl_GetRoundZ(&this->mTmpVoxel) < this->mnDimensionZ),
		    eResult, Volm_tErr_InvalidIdx );
  
  DebugCatch;
  EndDebugCatch;
  
  return eResult;
}

Volm_tErr Volm_VerifyFrame_ ( mriVolumeRef this,
			      int          iFrame ) {
  
  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugAssertThrowX( (iFrame >= 0 && iFrame < this->mnDimensionFrame),
		     eResult, Volm_tErr_InvalidIdx );
  
  DebugCatch;
  EndDebugCatch;
  
  return eResult;
}

char* Volm_GetErrorString ( Volm_tErr ieCode ) {
  
  Volm_tErr eCode = ieCode;
  
  if( ieCode < 0 ||
      ieCode >= Volm_knNumErrorCodes ) {
    eCode = Volm_tErr_Invalid;
  }
  
  return Volm_ksaErrorStrings[eCode];
}

