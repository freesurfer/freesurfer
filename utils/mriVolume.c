#include "mriVolume.h"
#include "error.h"
#include "xUtilities.h"
#include "mri_conform.h"

char Volm_ksaErrorStrings[Volm_knNumErrorCodes][Volm_knErrStringLen] = {
  "No error.",
  "Invalid signature.",
  "Invalid parameter.",
  "Invalid index.",
  "Memory allocation failed.",
  "Couldn't read volume.",
  "Couldn't copy volume.",
  "Couldn't copy transform.",
  "Couldn't normalize volume.",
  "Couldn't export volume to COR format.",
  "Normalized values not present.",
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
  this->mnDimension            = 0;
  this->mpNormValues           = NULL;
  this->mpRawValues            = NULL;
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
  if( NULL != this->mpNormValues ) {
    MRIfree( &(this->mpNormValues) );
  }
  if( NULL != this->mpRawValues ) {
    MRIfree( &(this->mpRawValues) );
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
  clone->mnDimension = this->mnDimension;
  
  /* copy mri volumes */
  DebugNote( ("Cloning normalized volume") );
  clone->mpNormValues = MRIcopy( this->mpNormValues, NULL );
  DebugAssertThrowX( (NULL != clone->mpNormValues),
         eResult, Volm_tErr_CouldntCopyVolume );

  if( NULL != this->mpRawValues ) {
    DebugNote( ("Cloning raw volume" ) );
    clone->mpRawValues = MRIcopy( this->mpRawValues, NULL );
    DebugAssertThrowX( (NULL != clone->mpRawValues),
           eResult, Volm_tErr_CouldntCopyVolume );
  }

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

Volm_tErr Volm_ImportData ( mriVolumeRef this,
          char*        isSource ) {

  Volm_tErr eResult           = Volm_tErr_NoErr;
  MRI*      rawVolume         = NULL;
  MRI*      normVolume        = NULL;
  MRI*      tmpVolume         = NULL;
  MATRIX*   identity          = NULL;
  MATRIX*   scannerTransform  = NULL;
  MATRIX*   idxToRASTransform = NULL;

  DebugEnterFunction( ("Volm_ImportData( this=%p, isSource=%s )",
           this, isSource ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isSource), eResult, Volm_tErr_InvalidParamater );

  /* attempt to read the volume in */
  DebugNote( ("Importing volume with MRIread") );
  normVolume = MRIread( isSource );
  DebugAssertThrowX( (NULL != normVolume), 
         eResult, Volm_tErr_CouldntReadVolume );

  /* if the element type is not uchar, then we need to make a copy and keep
     it around to get the raw values, because we're going to conform this
     one to uchar to get decent color values. */
  if( MRI_UCHAR != normVolume->type ) {
    
    /* make a copy for our raw version */
    DebugNote( ("Copying volume with MRIcopy") );
    rawVolume = MRIcopy( normVolume, NULL );
    DebugAssertThrowX( (NULL != rawVolume), 
         eResult, Volm_tErr_CouldntCopyVolume );
  }

  /* conform the volume */
  DebugNote( ("Normalizing volume with MRIconform") );
  tmpVolume = MRIconform( normVolume );
    DebugAssertThrowX( (NULL != tmpVolume), 
           eResult, Volm_tErr_CouldntNormalizeVolume );
  MRIfree( &normVolume );
  normVolume = tmpVolume;

  /* grab the xsize of the conformed volume as our dimension */
  this->mnDimension = normVolume->xend - normVolume->xstart;

  /* set the volumes in this */
  this->mpNormValues = normVolume;
  this->mpRawValues  = rawVolume;

  /* grab the scanner transform and copy it into our transform */
  if( NULL != this->mScannerTransform ) {
    DebugNote( ("Deleting existing scanner transform") );
    Trns_Delete( &(this->mScannerTransform) );
  }
  DebugNote( ("Creating scanner transform") );
  Trns_New( &this->mScannerTransform );
  DebugNote( ("Getting scanner transform matrix") );
  scannerTransform = MRIgetVoxelToRasXform( normVolume );
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

#if 0
  DebugNote( ("Creating idx to ras transform") );
  Trns_New( &this->mIdxToRASTransform );
  DebugNote( ("Getting idx to ras matrix") );
  idxToRASTransform = MRIgetVoxelToRasXform( normVolume );
  DebugAssertThrowX( (NULL != idxToRASTransform),
         eResult, Volm_tErr_AllocationFailed );
  DebugNote( ("Copying idx to ras transform matrix into transform") );
  Trns_CopyARAStoBRAS( this->mIdxToRASTransform, idxToRASTransform );
  DebugNote( ("Copying identity matrix into idx to ras transform") );
  Trns_CopyAtoRAS( this->mIdxToRASTransform, identity );
  Trns_CopyBtoRAS( this->mIdxToRASTransform, identity );
#else
  idxToRASTransform = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixClear( idxToRASTransform );
  *MATRIX_RELT(idxToRASTransform,1,1) = -1.0;
  *MATRIX_RELT(idxToRASTransform,2,3) = 1.0;
  *MATRIX_RELT(idxToRASTransform,3,2) = -1.0;
  *MATRIX_RELT(idxToRASTransform,1,4) = 128.0;
  *MATRIX_RELT(idxToRASTransform,2,4) = -128.0;
  *MATRIX_RELT(idxToRASTransform,3,4) = 128.0;
  *MATRIX_RELT(idxToRASTransform,4,4) = 1.0;

  DebugNote( ("Creating idx to ras transform") );
  Trns_New( &this->mIdxToRASTransform );
  DebugNote( ("Copying idx to ras transform matrix into transform") );
  Trns_CopyARAStoBRAS( this->mIdxToRASTransform, idxToRASTransform );
  DebugNote( ("Copying identity matrix into idx to ras transform") );
  Trns_CopyAtoRAS( this->mIdxToRASTransform, identity );
  Trns_CopyBtoRAS( this->mIdxToRASTransform, identity );
#endif

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
  eMRI = MRIwriteType( this->mpNormValues, sPath,MRI_CORONAL_SLICE_DIRECTORY );
  DebugAssertThrowX( (NO_ERROR == eMRI), 
         eResult, Volm_tErr_CouldntExportVolumeToCOR );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;

  return eResult;
}

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
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
         eResult, Volm_tErr_CouldntReadTransform );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_UnloadDisplayTransform ( mriVolumeRef this ) {

  Volm_tErr eResult   = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_UnloadDisplayTransform( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* if we already have a transform, delete it */
  if( NULL != this->mDisplayTransform )
    Trns_Delete( &(this->mDisplayTransform) );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;

  return eResult;
}



void Volm_GetColorAtIdx ( mriVolumeRef this,
        xVoxelRef    iIdx,
        xColor3fRef  oColor ) {

  Volm_tValue value = 0;
  xVoxel disp;
  Real rValue = 0;

  /* transform idx to display transform */
  if( NULL != this->mDisplayTransform ) {
    Volm_ApplyDisplayTransform_( this, iIdx, &disp );
    if( Volm_VerifyIdx_( this, &disp ) == Volm_tErr_NoErr ) {
      Volm_GetSincNormValueAtIdx_( this, &disp, &rValue );
      value = (Volm_tValue)rValue;
    }
  } else {
    if( Volm_VerifyIdx_( this, iIdx ) == Volm_tErr_NoErr )
      value = Volm_GetNormValueAtIdx_( this, iIdx );
  }
  *oColor = this->maColorTable[value];
}

void Volm_GetColorAtXYSlice ( mriVolumeRef     this,
            mri_tOrientation iOrientation,
            xPoint2nRef      iPoint,
            int              inSlice,
            xColor3fRef      oColor ) {

  Volm_tValue value;
  xVoxel      idx;
  Volm_ConvertXYSliceToIdx_( iOrientation, iPoint, inSlice, &idx );
  value = Volm_GetNormValueAtIdx_( this, &idx );
  *oColor = this->maColorTable[value];
}

void Volm_GetMaxColorAtIdx ( mriVolumeRef     this,
           xVoxelRef        iIdx,
           mri_tOrientation iOrientation,
           xColor3fRef      oColor ) {
  Volm_tValue value   = 0;
  xPoint2n    point;
  int         nSlice;
  
  /* if we haven't find the max values yet, find them. */
  if( NULL == this->mpMaxValues ) {
    Volm_FindMaxValues( this );
  }

  Volm_ConvertIdxToXYSlice_( iIdx, iOrientation, &point, &nSlice );
  value = Volm_GetMaxValueAtXYSlice_( this, iOrientation, &point );
  *oColor = this->maColorTable[value];
}

void Volm_GetMaxColorAtXYSlice ( mriVolumeRef     this,
         mri_tOrientation iOrientation,
         xPoint2nRef      iPoint,
         int              inSlice,
         xColor3fRef      oColor ) {

  Volm_tValue value;

  /* if we haven't find the max values yet, find them. */
  if( NULL == this->mpMaxValues ) {
    Volm_FindMaxValues( this );
  }

  value = Volm_GetMaxValueAtXYSlice_( this, iOrientation, iPoint );
  *oColor = this->maColorTable[value];
}

Volm_tErr Volm_GetDimension ( mriVolumeRef this,
            int*         onDimension ) {

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetDimension( this=%p, onDimension=%p )",
           this, onDimension ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (onDimension != NULL),
         eResult, Volm_tErr_InvalidParamater );

  /* return the dimension */
  DebugNote( ("Returning the dimension") );
  *onDimension = this->mnDimension;

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
  xVoxel disp;
  Real rValue = 0;

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

  /* return the raw value */
  DebugNote( ("Fetching the value at (%d,%d,%d) and return it",
       xVoxl_ExpandInt( iIdx )));
  /* transform idx to display transform */
  if( NULL != this->mDisplayTransform ) {
    Volm_ApplyDisplayTransform_( this, iIdx, &disp );
    if( Volm_VerifyIdx_( this, &disp ) == Volm_tErr_NoErr ) {
      Volm_GetSincNormValueAtIdx_( this, &disp, &rValue );
    }
    *oValue = (Volm_tValue)rValue;
  } else {
    *oValue = (float) Volm_GetRawValueAtIdx_( this, iIdx );
  }
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetValueAtIdxUnsafe ( mriVolumeRef this,
             xVoxelRef    iIdx,
             float*       oValue ) {

  xVoxel disp;
  Real rValue = 0;

  /* return the raw value. transform idx to display transform */
  if( NULL != this->mDisplayTransform ) {
    Volm_ApplyDisplayTransform_( this, iIdx, &disp );
    if( Volm_VerifyIdx_( this, &disp ) == Volm_tErr_NoErr ) {
      Volm_GetSincNormValueAtIdx_( this, &disp, &rValue );
    }
    *oValue = (Volm_tValue)rValue;
  } else {
    *oValue = (float) Volm_GetRawValueAtIdx_( this, iIdx );
  }

  return Volm_tErr_NoErr;
}

Volm_tErr Volm_SetValueAtIdx ( mriVolumeRef this,
             xVoxelRef    iIdx,
             Volm_tValue  iValue ) {

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

  /* set the norm value */
  DebugNote( ("Setting the norm value") );
  Volm_SetNormValueAtIdx_( this, iIdx, iValue );

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

  /* convert idx to ras */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to RAS with MRIvoxelToWorld",
        xVoxl_ExpandFloat( iIdx )) );
  MRIvoxelToWorld( this->mpNormValues, xVoxl_ExpandFloat( iIdx ),
       &rasX, &rasY, &rasZ );
  
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

  DebugEnterFunction( ("Volm_ConvertRASToIdx( this=%p, iRAS=%p, oIdx=%p )", 
           this, iRAS, oIdx) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iRAS != NULL && oIdx != NULL), 
         eResult, Volm_tErr_InvalidParamater );

  /* convert ras to idx */
  DebugNote( ("Converting RAS (%.2f, %.2f, %.2f) to idx with MRIworldToVoxel",
        xVoxl_ExpandFloat( iRAS )) );
  MRIworldToVoxel( this->mpNormValues, xVoxl_ExpandFloat( iRAS ),
       &idxX, &idxY, &idxZ );
  
  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oIdx, (float)idxX, (float)idxY, (float)idxZ );

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
  Real      talX    = 0;
  Real      talY    = 0;
  Real      talZ    = 0;

  DebugEnterFunction( ("Volm_ConvertIdxToMNITal( this=%p, iIdx=%p, "
           "oMNITal=%p )", this, iIdx, oMNITal) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oMNITal != NULL), 
         eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* convert idx to tal */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to mni tal with "
        "MRIvoxelToTalairachVoxel", xVoxl_ExpandFloat( iIdx )) );
  MRIvoxelToTalairach( this->mpNormValues, xVoxl_ExpandFloat( iIdx ),
           &talX, &talY, &talZ );
  
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
  Real      talX    = 0;
  Real      talY    = 0;
  Real      talZ    = 0;
  xVoxel    tal;

  DebugEnterFunction( ("Volm_ConvertIdxToTal( this=%p, iIdx=%p, "
           "oTal=%p )", this, iIdx, oTal) );
  
  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oTal != NULL), 
         eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* convert idx to tal */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to tal with "
        "MRIvoxelToTalairach", xVoxl_ExpandFloat( iIdx )) );
  MRIvoxelToTalairach( this->mpNormValues, xVoxl_ExpandFloat( iIdx ),
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
  eMRI = MRItalairachToVoxel( this->mpNormValues, xVoxl_ExpandFloat( &mniTal ),
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

  DebugEnterFunction( ("Volm_ConvertIdxToScanner( this=%p, iIdx=%p, "
           "oScanner=%p )", this, iIdx, oScanner) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iIdx != NULL && oScanner != NULL), 
         eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyIdx_( this, iIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* make sure we have the scanner transform */
  DebugAssertThrowX( (NULL != this->mScannerTransform),
         eResult, Volm_tErr_ScannerTransformNotPresent );

  /* convert idx to scanner */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to scanner", 
        xVoxl_ExpandFloat( iIdx )) );
  Trns_ConvertAtoB( this->mScannerTransform, iIdx, oScanner );

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
  Volm_tValue value = 0;

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
    
    value = Volm_GetNormValueAtIdx_( this, &idx );
    eVisit = iFunc( &idx, value, ipData );

    /* look at the result code they gave us. see if they want us to continue
       or skip a row. if skip, increment the voxel coordinate to the end. the
       next call to Increment will pop it over into to the next 
       row or plane. */
    switch( eVisit ) {
    case Volm_tVisitComm_Continue:
      bGo = xVoxl_IncrementUntilLimit( &idx, this->mnDimension-1 );
      continue;
      break;
    case Volm_tVisitComm_SkipRestOfRow:
      idx.mfX = this->mnDimension-1;
      bGo = xVoxl_IncrementUntilLimit( &idx, this->mnDimension-1 );
      break;
    case Volm_tVisitComm_SkipRestOfPlane:
      idx.mfX = this->mnDimension-1;
      idx.mfY = this->mnDimension-1;
      bGo = xVoxl_IncrementUntilLimit( &idx, this->mnDimension-1 );
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
  Volm_tValue      value       = 0;
#ifdef Solaris
  int                  i           = 0;
  int                  nSize       = 0;
#endif

  DebugEnterFunction( ("Volm_FindMaxValues( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  
  /* make sure we have our normalized volume */
  DebugAssertThrowX( (NULL != this->mpNormValues), 
         eResult, Volm_tErr_NormValuesNotPresent );

  /* we need 3 slices worth of data. if we have some allocated, delete them,
     and reallocate them. */
  if( NULL != this->mpMaxValues )
    free( this->mpMaxValues );
  this->mpMaxValues = (Volm_tValueRef) malloc( mri_knNumOrientations * 
                 pow( this->mnDimension, 2) *
                 sizeof(Volm_tValue) );
  DebugAssertThrowX( (NULL != this->mpMaxValues),
         eResult, Volm_tErr_AllocationFailed );
#ifdef Solaris
  DebugNote( ("Zeroing visited volume iteratively") );
  nSize = pow( mri_knNumOrientations * this->mnDimension, 2 );
  for( i = 0; i < nSize; i++ )
    this->mpMaxValues[i] = (Volm_tValue)0;
#else
  DebugNote( ("Zeroing visted volume with bzero") );
  bzero( this->mpMaxValues, mri_knNumOrientations * 
   pow(this->mnDimension,2) * sizeof(Volm_tValue) );
#endif

  /* go through the max value buffer. for each orientation */
  for( orientation = 0;
       orientation < mri_knNumOrientations; 
       orientation++ ) {

    /* for every point in... */
    for( point.mnY = 0; point.mnY < this->mnDimension; point.mnY++ ) {
      for( point.mnX = 0; point.mnX < this->mnDimension; point.mnX++ ) {
    
  /* now go thru every slice in the volume at this point.
     for every slice... */
  for( nSlice = 0; nSlice < this->mnDimension; nSlice++ ) {

    /* get the value at this point */
    Volm_GetNormValueAtXYSlice_( orientation, &point, nSlice, &value );

    if( value > Volm_GetMaxValueAtXYSlice_( this, orientation, &point ))
      Volm_SetMaxValueAtXYSlice_( this, orientation, &point, value ); 
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
  nSize = pow( this->mnDimension, 3 );
  this->mpSnapshot = (BUFTYPE*) malloc( nSize * sizeof(BUFTYPE) );
  DebugAssertThrowX( (NULL != this->mpSnapshot),
         eResult, Volm_tErr_AllocationFailed );

  /* copy it in */
  memcpy( this->mpSnapshot, this->mpNormValues->slices, nSize );

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
  nSize = pow( this->mnDimension, 3 );
  memcpy( this->mpNormValues->slices, this->mpSnapshot, nSize );

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
    newNorm = MRIrotateY( this->mpNormValues, NULL, fRadians );
    break;
  case mri_tOrientation_Horizontal:
    DebugNote( ("Rotating normalized volume with MRIrotateZ") );
    newNorm = MRIrotateZ( this->mpNormValues, NULL, fRadians );
    break;
  case mri_tOrientation_Sagittal:
    DebugNote( ("Rotating normalized volume with MRIrotateX") );
    newNorm = MRIrotateX( this->mpNormValues, NULL, fRadians );
    break;
  default:
    DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
    break;
  }    

  /* reassign to the new volume */
  DebugNote( ("Freeing old normalized volume and assigned rotated one") );
  MRIfree( &this->mpNormValues );
  this->mpNormValues = newNorm;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;
  
  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_Threshold ( mriVolumeRef this,
         Volm_tValue  iThreshold,
         tBoolean     ibAbove,
         Volm_tValue  iNewValue ) {

  Volm_tErr   eResult = Volm_tErr_NoErr;
  Volm_tValue value   = 0;
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
  while( xVoxl_IncrementUntilLimit( &idx, this->mnDimension-1 ) ) {
    
    value = Volm_GetNormValueAtIdx_( this, &idx );
  
    /* if we're going above and this value is above the thresh, or if we're
       going below and this value is below the thresh...*/
    if( (ibAbove && value > iThreshold) ||
  (!ibAbove && value < iThreshold ) ) {
      
      /* set the value to the new level. */
      Volm_SetNormValueAtIdx_( this, &idx, iNewValue );
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
  Volm_tValue value         = 0;
  Volm_tValue flippedValue  = 0;
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
  nXLimit = nYLimit = nZLimit = this->mnDimension-1;
  switch( iAxis ) {
  case mri_tOrientation_Coronal:
    nZLimit = (this->mnDimension-1) / 2;
    break;
  case mri_tOrientation_Sagittal:
    nXLimit = (this->mnDimension-1) / 2;
    break;
  case mri_tOrientation_Horizontal:
    nYLimit = (this->mnDimension-1) / 2;
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
     this->mnDimension-1 - xVoxl_GetZ( &idx ) );
      break;
    case mri_tOrientation_Sagittal:
      xVoxl_Set( &flippedIdx,
     this->mnDimension-1 - xVoxl_GetX( &idx ), 
     xVoxl_GetY( &idx ), 
     xVoxl_GetZ( &idx ) );
      break;
    case mri_tOrientation_Horizontal:
      xVoxl_Set( &flippedIdx,
     xVoxl_GetX( &idx ), 
     this->mnDimension-1 - xVoxl_GetY( &idx ), 
     xVoxl_GetZ( &idx ) );
      break;
    default:
      DebugAssertThrowX( (1), eResult, Volm_tErr_InvalidParamater );
      break;
    }
    
    /* swap the values */
    value = Volm_GetNormValueAtIdx_( this, &idx );
    flippedValue = Volm_GetNormValueAtIdx_( this, &flippedIdx );
    Volm_SetNormValueAtIdx_( this, &idx, flippedValue );
    Volm_SetNormValueAtIdx_( this, &flippedIdx, value );
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

#ifndef VOLM_USE_MACROS

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

Volm_tValue Volm_GetNormValueAtXYSlice_ ( mriVolumeRef this,
            mri_tOrientation  iOrientation,
            xPoint2nRef       iPoint,
            int               inSlice ) {

  xVoxel idx;
  Volm_ConvertXYSliceToIdx_( iOrientation, iPoint, inSlice, &idx );
  return Volm_GetNormValueAtIdx_( this, &idx );
}

Volm_tValue Volm_GetNormValueAtIdx_ ( mriVolumeRef this,
              xVoxelRef    iIdx ) {

  if( Volm_tErr_NoErr != Volm_VerifyIdx_( this, iIdx ) )
    return 0;

  return (Volm_tValue) MRIvox( this->mpNormValues, xVoxl_GetX(iIdx),
             xVoxl_GetY(iIdx), xVoxl_GetZ(iIdx) );
}

void Volm_SetNormValueAtIdx_ ( mriVolumeRef     this,
             xVoxelRef        iIdx,
             Volm_tValue      iValue ) {

  if( Volm_tErr_NoErr != Volm_VerifyIdx_( this, iIdx ) )
    return;

  MRIvox( this->mpNormValues, xVoxl_GetX(iIdx), 
    xVoxl_GetY(iIdx), xVoxl_GetZ(iIdx) ) = (BUFTYPE) iValue;

}

#endif /* VOLM_USE_MACROS */

float Volm_GetRawValueAtIdx_ ( mriVolumeRef     this,
             xVoxelRef        iIdx ) {
  
  xVoxel rasVox;
  MATRIX* rasToIdx = NULL;
  MATRIX* ras = NULL;
  MATRIX* idx = NULL;
  int i[3];

  if( Volm_tErr_NoErr != Volm_VerifyIdx_( this, iIdx ) )
    return 0;

  if( NULL != this->mpRawValues ) {
    
    DebugPrint( ("USING RAW VOLUME\n") );

    /* need to convert from normalized idx to ras space then to the
       raw index. first convert to ras. */
    Volm_ConvertIdxToRAS( this, iIdx, &rasVox );


    rasToIdx = MRIgetRasToVoxelXform( this->mpRawValues );
    ras = MatrixAlloc( 4, 1, MATRIX_REAL );
    idx = MatrixAlloc( 4, 1, MATRIX_REAL );

    *MATRIX_RELT(ras,1,1) = xVoxl_GetFloatX( &rasVox );
    *MATRIX_RELT(ras,2,1) = xVoxl_GetFloatY( &rasVox );
    *MATRIX_RELT(ras,3,1) = xVoxl_GetFloatZ( &rasVox );
    *MATRIX_RELT(ras,4,1) = 1;

    /* then to index */
    MatrixMultiply( rasToIdx, ras, idx );
    
    i[0] = rint( *MATRIX_RELT(idx,1,1) );
    i[1] = rint( *MATRIX_RELT(idx,2,1) );
    i[2] = rint( *MATRIX_RELT(idx,3,1) );

    /* check if in bounds */
    if( i[0] < 0 ||
  i[0] >= this->mpRawValues->width ||
  i[1] < 0 ||
  i[1] >= this->mpRawValues->height ||
  i[2] < 0 ||
  i[2] >= this->mpRawValues->depth )
      return 0;

    switch( this->mpRawValues->type ) {
    case MRI_UCHAR:
      return MRIvox( this->mpRawValues, i[0], i[1], i[2] );
      break;
    case MRI_INT:
      return MRIIvox( this->mpRawValues, i[0], i[1], i[2] );
      break;
    case MRI_LONG:
      return MRILvox( this->mpRawValues, i[0], i[1], i[2] );
      break;
    case MRI_FLOAT:
      return MRIFvox( this->mpRawValues, i[0], i[1], i[2] );
      break;
    case MRI_SHORT:
      return MRISvox( this->mpRawValues, i[0], i[1], i[2] );
      break;
    default:
      return 0;
    }
  } else {
     return MRIvox( this->mpNormValues, xVoxl_GetX(iIdx), 
         xVoxl_GetY(iIdx), xVoxl_GetZ(iIdx) );
  }
}

#ifndef VOLM_USE_MACROS

Volm_tValue Volm_GetMaxValueAtXYSlice_ ( mriVolumeRef     this,
           mri_tOrientation iOrientation, 
           xPoint2nRef      iPoint ) {

  return this->mpMaxValues[Volm_GetMaxValueIndex_(this,iOrientation,iPoint)];
}

void Volm_SetMaxValueAtXYSlice_ ( mriVolumeRef     this,
          mri_tOrientation iOrientation, 
          xPoint2nRef      iPoint,
          Volm_tValue      iValue ) {

  this->mpMaxValues[Volm_GetMaxValueIndex_(this,iOrientation,iPoint)] = iValue;
}


int Volm_GetMaxValueIndex_ ( mriVolumeRef     this,
           mri_tOrientation iOrientation, 
           xPoint2nRef      iPoint ) {

  return (iOrientation * this->mnDimension * this->mnDimension) +
    (iPoint->mnY * this->mnDimension) +
    iPoint->mnX;
}

#endif /* VOLM_USE_MACROS */

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

  DebugAssertThrowX( (xVoxl_GetX(iIdx) >= 0 &&
          xVoxl_GetRoundX(iIdx) < this->mnDimension &&
          xVoxl_GetY(iIdx) >= 0 &&
          xVoxl_GetRoundY(iIdx) < this->mnDimension &&
          xVoxl_GetZ(iIdx) >= 0 &&
          xVoxl_GetRoundZ(iIdx) < this->mnDimension),
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

