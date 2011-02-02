/**
 * @file  mriVolume.c
 * @brief early form of mri volume utils
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.93 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mriVolume.h"
#include "error.h"
#include "xUtilities.h"
#include "mri_conform.h"
#include "mriTransform.h"
#include "xList.h"
#include "macros.h" // FEQUAL
#include "proto.h" // nint
#include "fsenv.h"


// should be in transform.h if they aren't already

char Volm_ksaErrorStrings[Volm_knNumErrorCodes][Volm_knErrStringLen] =
{
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
  "Index to RAS transform not present.",
  "Flood user function can only return Continue or Stop",
  "While trying to find a suitable destination volume for the display transform, couldn't find FREESURFER_HOME environment variable.",
  "While trying to find a suitable destination volume for the display transform, couldn't find FREESURFER_HOME/average/mni305.cor.mgz.",
  "Error reading FREESURFER_HOME/average/mni305.cor.mgz while trying to population destination information in the display transform."
};

Volm_tErr Volm_New ( mriVolumeRef* opVolume )
{

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
  this->mpMaxValuesX           = NULL;
  this->mpMaxValuesY           = NULL;
  this->mpMaxValuesZ           = NULL;
  this->mSampleType            = Volm_tSampleType_Nearest;
  this->mMRIIdxToAnaIdxTransform     = NULL;
  this->mDisplayTransform      = NULL;
  this->mScannerTransform      = NULL;
  memset( this->msSubjectName, 0, sizeof( this->msSubjectName ) );
  memset( this->msVolumeName, 0, sizeof( this->msVolumeName ) );
  memset( this->msOriginalPath, 0, sizeof( this->msOriginalPath ) );
  this->mfColorBrightness      = Volm_kfDefaultBrightness;
  this->mfColorContrast        = Volm_kfDefaultContrast;
  this->mfColorMin             = 0;
  this->mfColorMax             = 0;
  memset( this->mafColorTable, 0, sizeof( this->mafColorTable ) );
  memset( this->manColorTable, 0, sizeof( this->manColorTable ) );

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

  if ( NULL != this )
  {
    Volm_Delete( &this );
  }

  EndDebugCatch;

  DebugExitFunction;

  if ( NULL != tmpMatrix )
  {
    MatrixFree( &tmpMatrix );
  }

  return eResult;
}

Volm_tErr Volm_Delete ( mriVolumeRef* iopVolume )
{

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
  if ( NULL != this->mpMriValues )
  {
    MRIfree( &(this->mpMriValues) );
  }
  if ( NULL != this->mpSnapshot )
  {
    free( this->mpSnapshot );
  }
  if ( NULL != this->mpMaxValuesX )
  {
    free( this->mpMaxValuesX );
  }
  if ( NULL != this->mpMaxValuesY )
  {
    free( this->mpMaxValuesY );
  }
  if ( NULL != this->mpMaxValuesZ )
  {
    free( this->mpMaxValuesZ );
  }
  if ( NULL != this->mDisplayTransform )
  {
    TransformFree( &(this->mDisplayTransform) );
  }
  if ( NULL != this->mMNITalLtzToRealTalTransform )
  {
    Trns_Delete( &(this->mMNITalLtzToRealTalTransform) );
  }
  if ( NULL != this->mMNITalGtzToRealTalTransform )
  {
    Trns_Delete( &(this->mMNITalGtzToRealTalTransform) );
  }
  if ( NULL != this->mScannerTransform)
  {
    Trns_Delete( &(this->mScannerTransform) );
  }
  if ( NULL != this->mpTmpScreenIdx )
  {
    MatrixFree( &this->mpTmpScreenIdx );
  }
  if ( NULL != this->mpTmpMRIIdx )
  {
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
                            mriVolumeRef* opVolume )
{

  Volm_tErr    eResult    = Volm_tErr_NoErr;
  Trns_tErr    eTransform = Trns_tErr_NoErr;
  mriVolumeRef clone      = NULL;

  DebugEnterFunction( ("Volm_DeepClone( this=%p, opVolume=%p )",
                       this, opVolume) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != opVolume), eResult, Volm_tErr_InvalidParamater );

  /* allocate a clone */
  DebugNote( ("Creating clone") );
  eResult = Volm_New( &clone );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );

  /* Copy simple data */
  clone->mnDimensionX      = this->mnDimensionX;
  clone->mnDimensionY      = this->mnDimensionY;
  clone->mnDimensionZ      = this->mnDimensionZ;
  clone->mnDimensionFrame  = this->mnDimensionFrame;
  clone->mSampleType       = this->mSampleType;
  clone->mfColorBrightness = this->mfColorBrightness;
  clone->mfColorContrast   = this->mfColorContrast;
  clone->mfColorMin        = this->mfColorMin;
  clone->mfColorMax        = this->mfColorMax;
  clone->mfMinValue        = this->mfMinValue;
  clone->mfMaxValue        = this->mfMaxValue;

  /* copy mri volumes */
  DebugNote( ("Cloning normalized volume") );
  clone->mpMriValues = MRIcopy( this->mpMriValues, NULL );
  DebugAssertThrowX( (NULL != clone->mpMriValues),
                     eResult, Volm_tErr_CouldntCopyVolume );

  /* copy the transforms */
  if ( NULL != this->mMRIIdxToAnaIdxTransform )
  {
    DebugNote( ("Copying MRI index to ana idx transform") );
    eTransform = Trns_DeepClone( this->mMRIIdxToAnaIdxTransform,
                                 &(clone->mMRIIdxToAnaIdxTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
                       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if ( NULL != this->mDisplayTransform )
  {
    DebugNote( ("Copying display transform") );
    clone->mDisplayTransform =
      TransformCopy( this->mDisplayTransform, NULL );
    DebugAssertThrowX( (NULL != clone->mDisplayTransform),
                       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if ( NULL != this->mMNITalLtzToRealTalTransform )
  {
    DebugNote( ("Copying real tal ltz transform") );
    eTransform = Trns_DeepClone( this->mMNITalLtzToRealTalTransform,
                                 &(clone->mMNITalLtzToRealTalTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
                       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if ( NULL != this->mMNITalGtzToRealTalTransform )
  {
    DebugNote( ("Copying real tal gtz transform") );
    eTransform = Trns_DeepClone( this->mMNITalGtzToRealTalTransform,
                                 &(clone->mMNITalGtzToRealTalTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
                       eResult, Volm_tErr_CouldntCopyTransform );
  }
  if ( NULL != this->mScannerTransform )
  {
    DebugNote( ("Copying scanner transform") );
    eTransform = Trns_DeepClone( this->mScannerTransform,
                                 &(clone->mScannerTransform) );
    DebugAssertThrowX( (Trns_tErr_NoErr == eTransform),
                       eResult, Volm_tErr_CouldntCopyTransform );
  }

  /* copy the strings */
  xUtil_strncpy( clone->msSubjectName, this->msSubjectName,
                 sizeof(clone->msSubjectName) );
  xUtil_strncpy( clone->msVolumeName, this->msVolumeName,
                 sizeof(clone->msVolumeName) );
  xUtil_strncpy( clone->msOriginalPath, this->msOriginalPath,
                 sizeof(clone->msOriginalPath) );

  /* allocate and copy color tables */
  memmove( clone->mafColorTable, this->mafColorTable,
           sizeof( this->manColorTable ));
  memmove( clone->manColorTable, this->manColorTable,
           sizeof( this->mafColorTable ));

  /* Copy the resample information then calc the resample matrices,
     that will point the m_resample{} matrices to the right place in
     the clone. */
  clone->mResampleMethod = this->mResampleMethod;
  clone->mResampleToRAS =
    MatrixCopy( this->mResampleToRAS, NULL );
  clone->mResampleToRASOrig =
    MatrixCopy( this->mResampleToRASOrig, NULL );
  clone->mResampleToRASInverse =
    MatrixCopy( this->mResampleToRASInverse, NULL );
  clone->mResampleToSlice =
    MatrixCopy( this->mResampleToSlice, NULL );
  clone->mResampleToSliceOrig =
    MatrixCopy( this->mResampleToSliceOrig, NULL );
  clone->mResampleToSliceInverse =
    MatrixCopy( this->mResampleToSliceInverse, NULL );

  Volm_SetResampleMethod( clone, clone->mResampleMethod );

  /* return the clone */
  *opVolume = clone;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );

  if ( NULL != clone )
  {
    Volm_Delete( &clone );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_CreateFromVolume ( mriVolumeRef this,
                                  mriVolumeRef iVolume )
{

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
  MRIcopyHeader( iVolume->mpMriValues, mriVolume );
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
                            char*        isSource )
{

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
                             MRI*        iMRI )
{

  Volm_tErr eResult           = Volm_tErr_NoErr;
  MATRIX*   identity          = NULL;
  MATRIX*   scannerTransform  = NULL;
  MATRIX*   m_resample;
  MATRIX*   m_resample_inv ;

  DebugEnterFunction( ("Volm_SetFromMRI_( this=%p, iMRI=%p )", this, iMRI ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iMRI), eResult, Volm_tErr_InvalidParamater );

  /* Clear old volume if we have one. */
  if ( NULL != this->mpMriValues &&
       this->mpMriValues != iMRI )
  {
    MRIfree( &this->mpMriValues );
  }

  /* Get the reange. */
  MRIvalRange(iMRI, &this->mfMinValue, &this->mfMaxValue) ;

  /* Set the color min/max to these values. */
  this->mfColorMin = this->mfMinValue;
  this->mfColorMax = this->mfMaxValue;


  /* Grab the dimensions and number of frames. */
  this->mnDimensionX = iMRI->width;
  this->mnDimensionY = iMRI->height;
  this->mnDimensionZ = iMRI->depth;
  this->mnDimensionFrame = iMRI->nframes;

  /* set the volumes in this */
  this->mpMriValues = iMRI;


  /* Initialize the resample matrix. */
  DebugNote( ("Computing norming matrix ") );
  m_resample = MRIgetConformMatrix( iMRI );
  DebugAssertThrowX( (NULL != m_resample),
                     eResult, Volm_tErr_CouldntNormalizeVolume );
  DebugNote( ("Setting resample RAS matrix") );
  this->mResampleToRAS = m_resample;
  DebugNote( ("Calculating inverse of resample matrix") );
  m_resample_inv = MatrixInverse( m_resample, NULL );
  DebugNote( ("Setting inverse resample RAS matrix") );
  this->mResampleToRASInverse = m_resample_inv;
  DebugNote( ("Copying original resample RAS matrix") );
  this->mResampleToRASOrig = MatrixCopy(m_resample, NULL) ;


  /* Calculate the resample matrices that don't use the RAS
     transform. This just does the pixel size transform. */
  DebugNote( ("Setting resample to slice matrix") );
  this->mResampleToSlice = MatrixIdentity( 4, NULL );
  stuff_four_by_four( this->mResampleToSlice,
                      1.0/this->mpMriValues->xsize, 0, 0, 0,
                      0, 1.0/this->mpMriValues->ysize, 0, 0,
                      0, 0, 1.0/this->mpMriValues->zsize, 0,
                      0, 0, 0, 1 );
  DebugNote( ("Calculating inverse of resample to slice matrix") );
  this->mResampleToSliceInverse =MatrixInverse( this->mResampleToSlice, NULL );
  DebugNote( ("Copying original resample to slice matrix") );
  this->mResampleToSliceOrig = MatrixCopy( this->mResampleToSlice, NULL );


  /* grab the scanner transform and copy it into our transform */
  if ( NULL != this->mScannerTransform )
  {
    DebugNote( ("Deleting existing scanner transform") );
    Trns_Delete( &(this->mScannerTransform) );
  }
  DebugNote( ("Creating scanner transform") );
  Trns_New( &this->mScannerTransform );
  DebugNote( ("Getting scanner transform matrix") );
  scannerTransform = extract_i_to_r( iMRI );
  DebugAssertThrowX( (NULL != scannerTransform),
                     eResult, Volm_tErr_AllocationFailed );
  DebugNote( ("Copying scanner transform matrix into transform") );
  Trns_CopyARAStoBRAS( this->mScannerTransform, scannerTransform );
  identity = MatrixIdentity( 4, NULL );
  DebugNote( ("Copying identity matrix into scanner transform") );
  Trns_CopyAtoRAS( this->mScannerTransform, identity );
  Trns_CopyBtoRAS( this->mScannerTransform, identity );


  /* Set the initial resample method. */
  Volm_SetResampleMethod( this, Volm_tResampleMethod_RAS );


  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  if ( NULL != scannerTransform )
  {
    MatrixFree( &scannerTransform );
  }
  if ( NULL != identity )
  {
    MatrixFree( &identity );
  }

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_CalculateMRIIdxToAnaIdx_ ( mriVolumeRef this )
{

  Volm_tErr eResult           = Volm_tErr_NoErr;
  MATRIX*   identity          = NULL;
  MATRIX*   idxToRASTransform = NULL;
  MATRIX*   BtoRAS            = NULL;

  DebugEnterFunction( ("Volm_CalculateMRIIdxToAnaIdx_( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );

  /* This function creates a matrix that in which A is MRI Idx and B
     is conformed idx or Ana idx. The RAS space is scanner RAS
     space. Then you can use A->RAS and RAS->A to convert between
     scanner RAS and MRI idx, and B->RAS and RAS->B to convert between
     scanner RAS and Ana Idx.

     Furthermore, this transform is never used in mriVolume.c, as we
     now use the MRIworldToVoxel etc functions. But this transform is
     used by tkmedit, esp in the functional transform code.

     Here's the Tosa diagram.


                           extract_i_to_r()
                           = MRIxfmCRS2XYZ()

       "MRI Idx" or (i,j,k) (A)   --------> RAS
                                             |
                                             |    Identity
                                             |
                                             V
      "Ana Idx" or conformed (B)  --------> RAS

                    extract_i_to_r * m_resample
                    = MRIxfmCRS2XYZ() * MRIGetConformMtx()

  */

  /* Delete existing transform. */
  if ( NULL != this->mMRIIdxToAnaIdxTransform )
  {
    DebugNote( ("Deleting existing MRI idx to ana idx transform") );
    Trns_Delete( &(this->mMRIIdxToAnaIdxTransform) );
  }

  /* Build the idxToRAS transform. */
  DebugNote( ("Creating MRI idx to ana idx transform") );
  Trns_New( &this->mMRIIdxToAnaIdxTransform );

  identity = MatrixIdentity( 4, NULL );

  /* AtoRAS is extract_i_to_r. This goes from MRI Idx or (i,j,k) to
     scanner RAS. This includes the voxel size calculation. */
  DebugNote( ("Getting idx to ras matrix") );
  idxToRASTransform = extract_i_to_r( this->mpMriValues );
  DebugAssertThrowX( (NULL != idxToRASTransform),
                     eResult, Volm_tErr_AllocationFailed );
  DebugNote( ("Copying idx to ras transform matrix into AtoRAS") );
  Trns_CopyAtoRAS( this->mMRIIdxToAnaIdxTransform, idxToRASTransform );

  /* ARSToBRAS is identity. */
  DebugNote( ("Copying identity matrix into ARAStoBRAS") );
  Trns_CopyARAStoBRAS( this->mMRIIdxToAnaIdxTransform, identity );

  /* BtoRAS = AtoRAS * m_resample */
  /*        = extract_i_to_r * m_resample */
  /*        = MRIxfmCRS2XYZ() * MRIGetConformMtx() */
  /* This goes from Ana Idx or conformed index to scanner RAS. */
  DebugNote( ("Copying calc'd matrix into BtoRAS") );
  BtoRAS = MatrixMultiply(idxToRASTransform, this->m_resample, NULL);
  Trns_CopyBtoRAS(this->mMRIIdxToAnaIdxTransform, BtoRAS);

  DebugNote( ("Freeing BtoRAS tmp matrix") );
  MatrixFree(&BtoRAS);

  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  if ( NULL != identity )
  {
    MatrixFree( &identity );
  }

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetMaxValueAtMRIIdx_ ( mriVolumeRef     this,
                                      xVoxelRef        iMRIIdx,
                                      mri_tOrientation iOrientation,
                                      float*           oValue )
{

  Volm_tErr eResult              = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetMaxValueAtMRIIdx_( this=%p, iMRIIdx=%p, "
                       "iOrientation=%d, oValue=%p )", this, iMRIIdx,
                       iOrientation, oValue) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  if ( NULL == this->mpMaxValuesX )
  {
    Volm_FindMaxValues( this );
  }

  switch ( iOrientation )
  {
  case mri_tOrientation_Sagittal:
    *oValue = this->mpMaxValuesX[xVoxl_GetY(iMRIIdx)][xVoxl_GetZ(iMRIIdx)];
    break;
  case mri_tOrientation_Horizontal:
    *oValue = this->mpMaxValuesY[xVoxl_GetX(iMRIIdx)][xVoxl_GetZ(iMRIIdx)];
    break;
  case mri_tOrientation_Coronal:
    *oValue = this->mpMaxValuesZ[xVoxl_GetX(iMRIIdx)][xVoxl_GetY(iMRIIdx)];
    break;
  default:
    *oValue = 0;
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_ExportNormToCOR ( mriVolumeRef this,
                                 char*        isPath )
{

  Volm_tErr eResult              = Volm_tErr_NoErr;
  char      sPath[mri_knPathLen] = "";
  int       eMRI                 = NO_ERROR;

  DebugEnterFunction( ("Volm_ExportNormToCOR( this=%p, isPath=%s )",
                       this, isPath) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* if out isPath is null, use the original path. */
  if ( NULL == isPath )
  {
    xUtil_strncpy( sPath, this->msOriginalPath, sizeof( sPath ) );
  }
  else
  {
    xUtil_strncpy( sPath, isPath, sizeof( sPath ) );
  }

  /* write the volume */
  DebugNote( ("Writing volume with MRIwriteType") );
#if 0
  eMRI = MRIwriteType( this->mpMriValues, sPath,MRI_CORONAL_SLICE_DIRECTORY );
#else
  eMRI = MRIwrite( this->mpMriValues, sPath );
#endif
  DebugAssertThrowX( (NO_ERROR == eMRI),
                     eResult, Volm_tErr_CouldntExportVolumeToCOR );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_Save ( mriVolumeRef this,
                      char*        isFileName,
                      tBoolean     ibSaveFileName )
{

  Volm_tErr eResult              = Volm_tErr_NoErr;
  char      sFileName[mri_knPathLen] = "";
  int       eMRI                 = NO_ERROR;

  DebugEnterFunction( ("Volm_Save( this=%p, isFileName=%s, ibSaveFileName )",
                       this, isFileName, ibSaveFileName) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* if out isFileName is null, use the original. */
  if ( NULL == isFileName )
  {
    xUtil_strncpy( sFileName, this->msOriginalPath, sizeof( sFileName ) );
  }
  else
  {
    xUtil_strncpy( sFileName, isFileName, sizeof( sFileName ) );
  }

  /* write the volume */
  DebugNote( ("Writing volume with MRIwrite") );
  eMRI = MRIwrite( this->mpMriValues, sFileName );
  DebugAssertThrowX( (NO_ERROR == eMRI),
                     eResult, Volm_tErr_CouldntExportVolumeToCOR );

  /* Save the filename if they want. */
  if ( NULL != isFileName && ibSaveFileName )
  {
    DebugNote( ("Saving file name") );
    xUtil_strncpy( this->msOriginalPath, isFileName,
                   sizeof( this->msOriginalPath ) );
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_LoadDisplayTransform ( mriVolumeRef this,
                                      char*        isFileName )
{
  Volm_tErr  eResult    = Volm_tErr_NoErr;
  TRANSFORM* transform  = NULL;
  int        eTransform = ERROR_NONE;
  LTA*       lta        = NULL;
  FSENV*     fsenv      = NULL;
  char       sMNI305FileName[1000] = "";
  int        eStat;
  struct stat statMNI305;
  MRI*       MNI305volume = NULL;

  DebugEnterFunction( ("Volm_LoadDisplayTransform( this=%p, isFileName=%s )",
                       this, isFileName) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != isFileName),
                     eResult, Volm_tErr_InvalidParamater );

  /* try to read a transform and make sure we got it */
  DebugNote( ("Reading transform from %s", isFileName) );
  transform = TransformRead( isFileName );
  DebugAssertThrowX( (NULL != transform),
                     eResult, Trns_tErr_LTAImportFailed );

  /* If this is an RAS_TO_RAS, take a look at the valid flags for
     source and dest. If they aren't valid, we have to provide some
     data. */
  if( transform->type == LINEAR_RAS_TO_RAS )
  {

    lta = (LTA *)transform->xform ;

    /* If the source is invalid, we can provide the data from this MRI
       volume. */
    if( !lta->xforms[0].src.valid )
    {

      fprintf( stderr, "INFO: Setting transform source geometry data to MRI volume's geometry\n" );

      getVolGeom( this->mpMriValues, &(lta->xforms[0].src) );
    }

    /* If the dest is invalid, and the file name is talairach.xfm, we
       can provide the volume from
       $FREESURFER_HOME/average/mni305.cor.mgz if it's available. */
    if( !lta->xforms[0].dst.valid &&
        NULL != strstr( isFileName, "talairach.xfm" ) )
    {

      fprintf( stderr, "INFO: Setting transform dest geometry data to mni305's geometry\n" );

      /* Get FREESURFER_HOME. */
      fsenv = FSENVgetenv();
      DebugAssertThrowX( (NULL != fsenv),
                         eResult, Volm_tErr_FsenvNotFoundWhileSearchingForTransformDest );
      DebugAssertThrowX( (NULL != fsenv->FREESURFER_HOME),
                         eResult, Volm_tErr_FsenvNotFoundWhileSearchingForTransformDest );

      /* Build the file name. */
      snprintf( sMNI305FileName, sizeof(sMNI305FileName),
                "%s/average/mni305.cor.mgz", fsenv->FREESURFER_HOME );

      /* Make sure the file exists. */
      eStat = stat( sMNI305FileName, &statMNI305 );
      DebugAssertThrowX( (0 == eStat),
                         eResult, Volm_tErr_MNI305NotFoundWhilePopulatingTramsformDest );

      /* Load it and read the header. */
      MNI305volume = MRIreadHeader( sMNI305FileName, MRI_VOLUME_TYPE_UNKNOWN );
      DebugAssertThrowX( (NULL != MNI305volume),
                         eResult, Volm_tErr_ErrorReadingMNI305WhilePopulatingTramsformDest );

      /* Copy the geometry data. */
      getVolGeom( MNI305volume, &(lta->xforms[0].dst) );
    }
  }

  /* This doesn't invert the transform (as you might think), but
     calculates the inverses that we'll use later. */
  eTransform = TransformInvert( transform, this->mpMriValues );
  if( 0 != eTransform )
  {
    fprintf( stderr, "TransformInvert returned %d\n", eTransform );
  }

  /* Success! Out with the old and in with the new. */
  if( NULL != this->mDisplayTransform )
  {
    TransformFree( &this->mDisplayTransform );
  }
  this->mDisplayTransform = transform;

  //  TransformRas2Vox(transform, this->mpMriValues, this->mpMriValues) ;
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  if( NULL != fsenv )
  {
    FSENVfree( &fsenv );
  }
  if( NULL != MNI305volume )
  {
    MRIfree( &MNI305volume );
  }

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_UnloadDisplayTransform ( mriVolumeRef this )
{

  Volm_tErr eResult   = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_UnloadDisplayTransform( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* if we already have a transform, delete it */
  if ( NULL != this->mDisplayTransform )
  {
    TransformFree( &(this->mDisplayTransform) );
    this->mDisplayTransform = NULL;
  }

  MatrixCopy(this->m_resample_orig, this->m_resample) ;
  DebugNote( ("Calculating inverse of resample matrix") );
  MatrixInverse( this->m_resample, this->m_resample_inv );
  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_CopyGeometryInformation ( mriVolumeRef this,
    VOL_GEOM*    ioVolumeGeometry )
{

  Volm_tErr eResult   = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_CopyGeometryInformation( this=%p, "
                       "ioVolumeGeometry=%p )", this, ioVolumeGeometry) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != ioVolumeGeometry),
                     eResult, Volm_tErr_InvalidParamater );

  /* Copy the information out. */
  getVolGeom( this->mpMriValues, ioVolumeGeometry );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


void Volm_GetIntColorAtIdx ( mriVolumeRef this,
                             xVoxelRef    iIdx,
                             xColor3nRef  oColor )
{

  float  x, y, z;
  float  value = 0;
  int    colorIdx = 0;
  xVoxel disp;

  /* transform idx to display transform */
  if ( NULL != this->mDisplayTransform )
  {

    /* Transform using the inverse of the display transform. */
    TransformSampleInverse( this->mDisplayTransform,
                            xVoxl_ExpandFloat( iIdx ), &x, &y, &z );
    xVoxl_SetFloat( &disp, x, y, z );

    if ( Volm_VerifyIdx_( this, &disp ) == Volm_tErr_NoErr )
    {
      Volm_GetSampledValueAtIdx_( this, &disp, &value );
    }

  }
  else
  {

    /* Get the sampled value. Cap it to the color min and max. If in
       between that, calculate the color index between 0 and 255. */
    Volm_GetSampledValueAtIdx_( this, iIdx, &value );
  }

  if ( value <= this->mfColorMin )
  {
    colorIdx = 0;
  }
  else if ( value >= this->mfColorMax )
  {
    colorIdx = 255;
  }
  else
  {
    colorIdx = (int) (Volm_knNumColorTableEntries *
                      (value - this->mfColorMin) /
                      (this->mfColorMax - this->mfColorMin)) ;
  }

  *oColor = this->manColorTable[colorIdx];
}


void Volm_GetMaxIntColorAtIdx ( mriVolumeRef     this,
                                xVoxelRef        iIdx,
                                mri_tOrientation iOrientation,
                                xColor3nRef      oColor )
{
  xVoxel MRI;
  float  value    = 0;
  int    colorIdx = 0;


  Volm_ConvertIdxToMRIIdx( this, iIdx, &MRI );

  Volm_GetMaxValueAtMRIIdx_( this, &MRI, iOrientation, &value );

  colorIdx = (int) (255.0 * (value - this->mfMinValue) /
                    (this->mfMaxValue - this->mfMinValue)) ;
  *oColor = this->manColorTable[colorIdx];
}

Volm_tErr Volm_GetDimensions ( mriVolumeRef this,
                               int*         onDimensionX,
                               int*         onDimensionY,
                               int*         onDimensionZ )
{

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
                                   int* onDimensionFrames )
{

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

Volm_tErr Volm_GetType ( mriVolumeRef this, int* onType )
{

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

Volm_tErr Volm_GetValueMinMax ( mriVolumeRef this,
                                float*       ofMin,
                                float*       ofMax )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetValueMinMax( this=%p, ofMin=%p, ofMax )",
                       this, ofMin, ofMax ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (ofMin != NULL), eResult, Volm_tErr_InvalidParamater );
  DebugAssertThrowX( (ofMax != NULL), eResult, Volm_tErr_InvalidParamater );

  /* return the min and max */
  DebugNote( ("Returning the min") );
  *ofMin = this->mfMinValue;
  *ofMax = this->mfMaxValue;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetResampleMethod ( mriVolumeRef          this,
                                   Volm_tResampleMethod* oMethod )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetResampleMethod( this=%p, "
                       "oMethod=%p )", this, oMethod ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oMethod != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* return the method */
  DebugNote( ("Returning the method") );
  *oMethod = this->mResampleMethod;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_SetResampleMethod ( mriVolumeRef         this,
                                   Volm_tResampleMethod iMethod )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetResampleMethod( this=%p, iMethod=%d )",
                       this, (int)iMethod ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iMethod >= 0 && iMethod < Volm_knNumResampleMethods),
                     eResult, Volm_tErr_InvalidParamater );

  /* Set the resample method. */
  this->mResampleMethod = iMethod;

  /* Point the matrices correctly. */
  switch ( this->mResampleMethod )
  {
  case Volm_tResampleMethod_RAS:
    this->m_resample      = this->mResampleToRAS;
    this->m_resample_orig = this->mResampleToRASOrig;
    this->m_resample_inv  = this->mResampleToRASInverse;
    break;
  case Volm_tResampleMethod_Slice:
    this->m_resample      = this->mResampleToSlice;
    this->m_resample_orig = this->mResampleToSliceOrig;
    this->m_resample_inv  = this->mResampleToSliceInverse;
    break;
  default:
    break;
  }


  /* Recalc the mri idx to and idx transform. */
  Volm_CalculateMRIIdxToAnaIdx_( this );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetSampleType ( mriVolumeRef      this,
                               Volm_tSampleType* oType )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetSampleType( this=%p, "
                       "oType=%p )", this, oType ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oType != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* return the type */
  DebugNote( ("Returning the type") );
  *oType = this->mSampleType;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_SetSampleType ( mriVolumeRef     this,
                               Volm_tSampleType iType )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetSampleType( this=%p, iType=%d )",
                       this, (int)iType ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iType >= 0 && iType < Volm_knNumSampleTypes),
                     eResult, Volm_tErr_InvalidParamater );

  /* Set the sample type. */
  this->mSampleType = iType;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetValueAtIdx ( mriVolumeRef this,
                               xVoxelRef    iIdx,
                               float*       oValue )
{

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
                                     float*       oValue )
{

  Volm_GetValueAtIdx_( this, iIdx, oValue );

  return Volm_tErr_NoErr;
}

Volm_tErr Volm_SetValueAtIdx ( mriVolumeRef this,
                               xVoxelRef    iIdx,
                               float        iValue )
{

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

  /* update our min max. */
  if ( iValue > this->mfMaxValue )
  {
    this->mfMaxValue = iValue;
  }
  if ( iValue < this->mfMinValue )
  {
    this->mfMinValue = iValue;
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetValueAtIdxFrame  ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     int          iFrame,
                                     float*       oValue )
{

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
    float*       oValue )
{

  Volm_GetValueAtIdxFrame_( this, iIdx, iFrame, oValue );

  return Volm_tErr_NoErr;
}

Volm_tErr Volm_SetValueAtIdxFrame ( mriVolumeRef this,
                                    xVoxelRef    iIdx,
                                    int          iFrame,
                                    float        iValue )
{

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

Volm_tErr Volm_GetValueAtMRIIdx ( mriVolumeRef this,
                                  xVoxelRef    iMRIIdx,
                                  float*       oValue )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetValueAtMRIIdx( this=%p, iMRIIdx=%p, "
                       "oValue=%p )", this, iMRIIdx, oValue ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iMRIIdx != NULL && oValue != NULL),
                     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyMRIIdx_( this, iMRIIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* return the value */
  Volm_GetValueAtMRIIdx_( this, iMRIIdx, oValue );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_SetValueAtMRIIdx ( mriVolumeRef this,
                                  xVoxelRef    iMRIIdx,
                                  float        iValue )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetValueAtIdx( this=%p, iMRIIdx=%p, "
                       "iValue=%d )", this, iMRIIdx, (int)iValue ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iMRIIdx != NULL), eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyMRIIdx_( this, iMRIIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* set the value */
  Volm_SetValueAtMRIIdx_( this, iMRIIdx, iValue );

  /* update our min max. */
  if ( iValue > this->mfMaxValue )
  {
    this->mfMaxValue = iValue;
  }
  if ( iValue < this->mfMinValue )
  {
    this->mfMinValue = iValue;
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_HasTalTransform ( mriVolumeRef this,
                                 tBoolean*    obHasTransform )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_HasTalTransform( this=%p, obHasTransform=%p)",
                       this, obHasTransform) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (obHasTransform != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* Return true if we have a tal transform, false if not. */
  if ( NULL != this->mpMriValues->linear_transform )
  {
    *obHasTransform = TRUE;
  }
  else
  {
    *obHasTransform = FALSE;
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_ConvertIdxToRAS ( mriVolumeRef this,
                                 xVoxelRef    iIdx,
                                 xVoxelRef    oRAS )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;
  double      rasX    = 0;
  double      rasY    = 0;
  double      rasZ    = 0;


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
                                 xVoxelRef    oIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      idxX    = 0;
  double      idxY    = 0;
  double      idxZ    = 0;
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
                                    xVoxelRef    oMNITal )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    scannerRAS;
  double      talX    = 0;
  double      talY    = 0;
  double      talZ    = 0;

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

  /* Convert the index coords to scanner RAS coords. */
  Volm_ConvertIdxToScanner( this, iIdx, &scannerRAS );

  /* Transform them to MNI talaiarch with the transform in the MRI. */
  transform_point( this->mpMriValues->linear_transform,
                   xVoxl_ExpandFloat( &scannerRAS ),
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
                                 xVoxelRef    oTal )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    scannerRAS;
  double      talX    = 0;
  double      talY    = 0;
  double      talZ    = 0;
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

  /* Convert the index coords to scanner RAS coords. */
  Volm_ConvertIdxToScanner( this, iIdx, &scannerRAS );

  /* Transform them to MNI talaiarch with the transform in the MRI. */
  transform_point( this->mpMriValues->linear_transform,
                   xVoxl_ExpandFloat( &scannerRAS ),
                   &talX, &talY, &talZ );

  /* convert to real tal. switch on the z to see which transform to use. */
  xVoxl_SetFloat( &tal, talX, talY, talZ );
  if ( xVoxl_GetFloatZ( &tal ) > 0 )
  {
    DebugNote( ("Converting to real tal with >0 transform") );
    Trns_ConvertAtoB( this->mMNITalGtzToRealTalTransform, &tal, oTal );
  }
  else
  {
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
                                 xVoxelRef    oIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      idxX    = 0;
  double      idxY    = 0;
  double      idxZ    = 0;
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
  if ( xVoxl_GetFloatZ( iTal ) > 0 )
  {
    DebugNote( ("Converting to mni tal with >0 transform") );
    Trns_ConvertBtoA( this->mMNITalGtzToRealTalTransform, iTal, &mniTal );
  }
  else
  {
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
                                     xVoxelRef    oScanner )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  xVoxel    mriIdx;

  DebugEnterFunction( ("Volm_ConvertIdxToScanner( this=%p, iIdx=%p, "
                       "oScanner=%p )", this, iIdx, oScanner) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

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

Volm_tErr Volm_ConvertMRIIdxToScanner ( mriVolumeRef this,
                                        xVoxelRef    iMRIIdx,
                                        xVoxelRef    oScanner )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_ConvertIdxToScanner( this=%p, iMRIIdx=%p, "
                       "oScanner=%p )", this, iMRIIdx, oScanner) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* make sure we have the scanner transform */
  DebugAssertThrowX( (NULL != this->mScannerTransform),
                     eResult, Volm_tErr_ScannerTransformNotPresent );

  /* convert idx to scanner */
  DebugNote( ("Converting idx (%.2f, %.2f, %.2f) to scanner",
              xVoxl_ExpandFloat( iMRIIdx )) );
  Trns_ConvertAtoB( this->mScannerTransform, iMRIIdx, oScanner );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_ConvertIdxToMRIIdx ( mriVolumeRef this,
                                    xVoxelRef    iIdx,
                                    xVoxelRef    oMRIIdx )
{

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

Volm_tErr Volm_ConvertMRIIdxToIdx  ( mriVolumeRef this,
                                     xVoxelRef    iMRIIdx,
                                     xVoxelRef    oIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_ConvertMRIIdxToIdx( this=%p, iMRIIdx=%p, "
                       "oIdx=%p )", this, iMRIIdx, oIdx) );

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

  /* Convert the incoming local MRI to  screen idx. */
  DebugNote( ("Converting MRI idx to screen idx") );
  Volm_ConvertMRIIdxToScreenIdx_( this, iMRIIdx, oIdx );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_ConvertMRIIdxToRAS ( mriVolumeRef this,
                                    xVoxelRef    iMRIIdx,
                                    xVoxelRef    oRAS )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      rasX    = 0;
  double      rasY    = 0;
  double      rasZ    = 0;


  DebugEnterFunction( ("Volm_ConvertIdxToRAS( this=%p, iMRIIdx=%p, oRAS=%p )",
                       this, iMRIIdx, oRAS) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iMRIIdx != NULL && oRAS != NULL),
                     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyMRIIdx_( this, iMRIIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Convert idx to ras */
  DebugNote( ("Converting idx to RAS with MRIvoxelToWorld") );
  MRIvoxelToWorld( this->mpMriValues, xVoxl_ExpandFloat( iMRIIdx ),
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

Volm_tErr Volm_ConvertRASToMRIIdx ( mriVolumeRef this,
                                    xVoxelRef    iRAS,
                                    xVoxelRef    oMRIIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      idxX    = 0;
  double      idxY    = 0;
  double      idxZ    = 0;

  DebugEnterFunction( ("Volm_ConvertRASToIdx( this=%p, iRAS=%p, oMRIIdx=%p )",
                       this, iRAS, oMRIIdx) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iRAS != NULL && oMRIIdx != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* convert ras to MRI idx */
  DebugNote( ("Converting RAS to mri idx with MRIworldToVoxel") );
  MRIworldToVoxel( this->mpMriValues, xVoxl_ExpandFloat( iRAS ),
                   &idxX, &idxY, &idxZ );

  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oMRIIdx, (float)idxX, (float)idxY, (float)idxZ );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_ConvertMRIIdxToSurfaceRAS ( mriVolumeRef this,
    xVoxelRef    iMRIIdx,
    xVoxelRef    oSurfaceRAS )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      rasX    = 0;
  double      rasY    = 0;
  double      rasZ    = 0;


  DebugEnterFunction( ("Volm_ConvertIdxToSurfaceRAS( this=%p, iMRIIdx=%p, "
                       "oRAS=%p )", this, iMRIIdx, oSurfaceRAS) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iMRIIdx != NULL && oSurfaceRAS != NULL),
                     eResult, Volm_tErr_InvalidParamater );
  eResult = Volm_VerifyMRIIdx_( this, iMRIIdx );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Convert idx to ras */
  DebugNote( ("Converting idx to surface RAS with MRIvoxelToSurfaceRAS") );
  MRIvoxelToSurfaceRAS( this->mpMriValues, xVoxl_ExpandFloat( iMRIIdx ),
                        &rasX, &rasY, &rasZ );

  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oSurfaceRAS, (float)rasX, (float)rasY, (float)rasZ );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_ConvertSurfaceRASToMRIIdx ( mriVolumeRef this,
    xVoxelRef    iSurfaceRAS,
    xVoxelRef    oMRIIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  double      idxX    = 0;
  double      idxY    = 0;
  double      idxZ    = 0;

  DebugEnterFunction( ("Volm_ConvertRASToIdx( this=%p, iSurfaceRAS=%p, "
                       "oMRIIdx=%p )", this, iSurfaceRAS, oMRIIdx) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iSurfaceRAS != NULL && oMRIIdx != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* convert ras to MRI idx */
  DebugNote( ("Converting surface RAS to mri idx with MRIsurfaceRASToVoxel") );
  MRIsurfaceRASToVoxel( this->mpMriValues, xVoxl_ExpandFloat( iSurfaceRAS ),
                        &idxX, &idxY, &idxZ );

  /* stuff results */
  DebugNote( ("Stuffing result into xVoxel") );
  xVoxl_SetFloat( oMRIIdx, (float)idxX, (float)idxY, (float)idxZ );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetMRIIdxToAnaIdxTransform ( mriVolumeRef     this,
    mriTransformRef* opTransform )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_GetMRIIdxToAnaIdxTransform( this=%p, "
                       "opTransform=%p )", this, opTransform) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != opTransform),
                     eResult, Volm_tErr_InvalidParamater );

  /* make sure we have the transform */
  DebugAssertThrowX( (NULL != this->mMRIIdxToAnaIdxTransform),
                     eResult, Volm_tErr_IdxToRASTransformNotPresent );

  /* return a pointer to it */
  *opTransform = this->mMRIIdxToAnaIdxTransform;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_Flood ( mriVolumeRef        this,
                       Volm_tFloodParams*  iParams )
{

  Volm_tErr    eResult       = Volm_tErr_NoErr;
  xList_tErr   eList         = xList_tErr_NoErr;
  xListRef     list          = NULL;
  int          nSize         = 0;
  tBoolean*    visited       = NULL;
  xVoxelRef    curVoxel      = NULL;
  void*        pCurVoxel     = NULL;
  xVoxelRef    newVoxel      = NULL;
  Volm_tVisitCommand eVisit  = Volm_tVisitComm_Continue;
  int          nVisitedIndex = 0;
  int          nZ            = 0;
  int          nBeginX       = 0;
  int          nEndX         = 0;
  int          nBeginY       = 0;
  int          nEndY         = 0;
  int          nBeginZ       = 0;
  int          nEndZ         = 0;
  int          nY            = 0;
  int          nX            = 0;
  xVoxel       idx;
  float        fValue        = 0;
  float        fDistance     = 0;
#ifdef Solaris
  int          i             = 0;
#endif

  DebugEnterFunction( ("Volm_Flood( this=%p, iParams=%p )", this, iParams) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iParams->mfSourceValue >= this->mfMinValue &&
                      iParams->mfSourceValue <= this->mfMaxValue),
                     eResult, Volm_tErr_InvalidParamater );
  if ( NULL == iParams->mComparatorFunc )
  {
    DebugAssertThrowX
    ((iParams->mComparatorType>Volm_tValueComparator_Invalid &&
      iParams->mComparatorType < Volm_knNumValueComparators),
     eResult, Volm_tErr_InvalidParamater );
  }
  if ( FALSE == iParams->mb3D )
  {
    DebugAssertThrowX( (iParams->mOrientation > mri_tOrientation_None &&
                        iParams->mOrientation < mri_knNumOrientations),
                       eResult, Volm_tErr_InvalidParamater );
  }
  DebugAssertThrowX( (iParams->mpFunction != NULL),
                     eResult, Volm_tErr_InvalidParamater );

  /* Init a volume to keep track of visited voxels */
  DebugNote( ("Allocating visited volume") );
  nSize = this->mnDimensionX * this->mnDimensionY * this->mnDimensionZ;
  visited = (tBoolean*) malloc( sizeof(tBoolean) * nSize) ;
  DebugAssertThrowX( (NULL != visited), eResult, Volm_tErr_AllocationFailed );

  DebugNote( ("Zeroing visited volume") );
  memset( visited, 0, nSize * sizeof( tBoolean ) );

  /* Initialize the list of voxels. Create a new voxel, copy our
     source into it, and put it in the list. */
  xList_New( &list );
  xVoxl_New( &newVoxel );
  xVoxl_Copy( newVoxel, &iParams->mSourceIdx );
  xList_PushItem( list, newVoxel );

  /* Start going through the list. */
  pCurVoxel = NULL;
  eList = xList_tErr_NoErr;
  while ( xList_tErr_NoErr == eList )
  {

    /* Try to pop a voxel off the list. If we get one... */
    eList = xList_PopItem ( list, &pCurVoxel );
    if ( eList == xList_tErr_NoErr && NULL != pCurVoxel )
    {

      /* To avoid type punning, which is bad. */
      curVoxel = (xVoxelRef) pCurVoxel;

      /* Make sure it's in bounds. */
      if ( Volm_tErr_NoErr != Volm_VerifyMRIIdx_( this, curVoxel ) )
      {
        continue;
      }

      /* If we've already been here, continue. If not, mark our visited
      volume and continue. */
      nVisitedIndex = xVoxl_ExpandToIndex( curVoxel, this->mnDimensionZ,
                                           this->mnDimensionY );
      if ( visited[nVisitedIndex] )
      {
        continue;
      }
      visited[nVisitedIndex] = TRUE;

      /* Get the value and do the comparison. If we have a type and
      not a function, just look at the values. Otherwise call the
      comparator function they gave us. If it's okay, go on, if
      not, continue. */
      Volm_GetValueAtMRIIdx_( this, curVoxel, &fValue );
      if ( NULL != iParams->mComparatorFunc )
      {

        /* Call the comparator function they gave us. If it returns
           false, continue. */
        if ( FALSE == iParams->mComparatorFunc( curVoxel, fValue,
                                                iParams->mComparatorFuncData ))
        {
          continue;
        }

      }
      else
      {

        /* Compare ths source value and the original value. */
        switch ( iParams->mComparatorType )
        {
        case Volm_tValueComparator_LTE:
          if ( !(fValue <= iParams->mfSourceValue + iParams->mfFuzziness) )
          {
            continue;
          }
          break;
        case Volm_tValueComparator_EQ:
          if ( !(iParams->mfSourceValue - iParams->mfFuzziness <= fValue &&
                 fValue <= iParams->mfSourceValue + iParams->mfFuzziness) )
          {
            continue;
          }
          break;
        case Volm_tValueComparator_GTE:
          if ( !(fValue >= iParams->mfSourceValue - iParams->mfFuzziness) )
          {
            continue;
          }
          break;
        default:
          break;
        }
      }

      /* Check distance if >0. if it's over our max distance, exit. */
      if ( iParams->mfMaxDistance > 0 )
      {
        fDistance =
          sqrt( ((xVoxl_GetX(curVoxel) - xVoxl_GetX(&(iParams->mSourceIdx))) *
                 (xVoxl_GetX(curVoxel) - xVoxl_GetX(&(iParams->mSourceIdx)))) +
                ((xVoxl_GetY(curVoxel) - xVoxl_GetY(&(iParams->mSourceIdx))) *
                 (xVoxl_GetY(curVoxel) - xVoxl_GetY(&(iParams->mSourceIdx)))) +
                ((xVoxl_GetZ(curVoxel) - xVoxl_GetZ(&(iParams->mSourceIdx))) *
                 (xVoxl_GetZ(curVoxel) - xVoxl_GetZ(&(iParams->mSourceIdx)))));
        if ( fDistance > iParams->mfMaxDistance )
        {
          continue;
        }
      }

      /* This is good, so call the user function. Handle their return code. */
      eVisit =
        iParams->mpFunction( curVoxel, fValue, iParams->mpFunctionData );
      switch ( eVisit )
      {
      case Volm_tVisitComm_SkipRestOfRow:
      case Volm_tVisitComm_SkipRestOfPlane:
        DebugThrowX( eResult, Volm_tErr_FloodVisitCommandNotSupported );
        break;
      case Volm_tVisitComm_Stop:
        goto cleanup;
        break;
      default:
        break;
      }

      /* Calc our bounds */
      nBeginX  = MAX( xVoxl_GetX(curVoxel) - 1, 0 );
      nEndX    = MIN( xVoxl_GetX(curVoxel) + 1, this->mnDimensionX );
      nBeginY  = MAX( xVoxl_GetY(curVoxel) - 1, 0 );
      nEndY    = MIN( xVoxl_GetY(curVoxel) + 1, this->mnDimensionY );
      nBeginZ  = MAX( xVoxl_GetZ(curVoxel) - 1, 0 );
      nEndZ    = MIN( xVoxl_GetZ(curVoxel) + 1, this->mnDimensionZ );

      /* If we're not in 3d, limit one of them based on the orientation */
      if ( !iParams->mb3D )
      {
        switch ( iParams->mOrientation )
        {
        case mri_tOrientation_Coronal:
          nBeginZ = nEndZ = xVoxl_GetZ(curVoxel);
          break;
        case mri_tOrientation_Horizontal:
          nBeginY = nEndY = xVoxl_GetY(curVoxel);
          break;
        case mri_tOrientation_Sagittal:
          nBeginX = nEndX = xVoxl_GetX(curVoxel);
          break;
        default:
          goto cleanup;
          break;
        }
      }

      /* Check the surrounding voxels. For each one, create a new
      voxel and add it to the list.  */
      for ( nZ = nBeginZ; nZ <= nEndZ; nZ++ )
      {
        for ( nY = nBeginY; nY <= nEndY; nY++ )
        {
          for ( nX = nBeginX; nX <= nEndX; nX++ )
          {
            xVoxl_Set( &idx, nX, nY, nZ );
            if ( Volm_VerifyIdx_( this, &idx ) == Volm_tErr_NoErr )
            {
              xVoxl_New( &newVoxel );
              xVoxl_Set( newVoxel, nX, nY, nZ );
              xList_PushItem( list, newVoxel );
            }
          }
        }
      }
      /* Delete the current voxel. */
      xVoxl_Delete( &curVoxel );
    }
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  if ( NULL != visited )
  {
    DebugNote( ("Freeing visited volume") );
    free( visited );
  }

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_VisitAllVoxels ( mriVolumeRef        this,
                                Volm_tVisitFunction iFunc,
                                void*               ipData )
{

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
  while ( bGo )
  {

    Volm_GetValueAtIdx_( this, &idx, &value );
    eVisit = iFunc( &idx, value, ipData );

    /* look at the result code they gave us. see if they want us to continue
       or skip a row. if skip, increment the voxel coordinate to the end. the
       next call to Increment will pop it over into to the next
       row or plane. */
    switch ( eVisit )
    {
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

Volm_tErr Volm_FindMaxValues ( mriVolumeRef this )
{

  Volm_tErr        eResult     = Volm_tErr_NoErr;
  float            value       = 0;
  int              nX          = 0;
  int              nY          = 0;
  int              nZ          = 0;
  xVoxel           MRI;
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

  /* Allocate our buffers. */
  if ( NULL != this->mpMaxValuesX )
  {
    free( this->mpMaxValuesX );
  }
  if ( NULL != this->mpMaxValuesY )
  {
    free( this->mpMaxValuesY );
  }
  if ( NULL != this->mpMaxValuesZ )
  {
    free( this->mpMaxValuesZ );
  }

  this->mpMaxValuesX = (float**) calloc( this->mnDimensionY, sizeof(float*) );
  DebugAssertThrowX( (NULL != this->mpMaxValuesX),
                     eResult, Volm_tErr_AllocationFailed );
  for ( nY = 0; nY < this->mnDimensionY; nY++ )
  {
    this->mpMaxValuesX[nY] =
      (float*) calloc( this->mnDimensionZ, sizeof(float*) );
    DebugAssertThrowX( (NULL != this->mpMaxValuesX[nY]),
                       eResult, Volm_tErr_AllocationFailed );
  }

  this->mpMaxValuesY = (float**) calloc( this->mnDimensionX, sizeof(float*) );
  DebugAssertThrowX( (NULL != this->mpMaxValuesY),
                     eResult, Volm_tErr_AllocationFailed );
  for ( nX = 0; nX < this->mnDimensionX; nX++ )
  {
    this->mpMaxValuesY[nX] =
      (float*) calloc( this->mnDimensionZ, sizeof(float*) );
    DebugAssertThrowX( (NULL != this->mpMaxValuesY[nX]),
                       eResult, Volm_tErr_AllocationFailed );
  }

  this->mpMaxValuesZ = (float**) calloc( this->mnDimensionX, sizeof(float*) );
  DebugAssertThrowX( (NULL != this->mpMaxValuesZ),
                     eResult, Volm_tErr_AllocationFailed );
  for ( nX = 0; nX < this->mnDimensionX; nX++ )
  {
    this->mpMaxValuesZ[nX] =
      (float*) calloc( this->mnDimensionY, sizeof(float*) );
    DebugAssertThrowX( (NULL != this->mpMaxValuesZ[nX]),
                       eResult, Volm_tErr_AllocationFailed );
  }

  /* Go through the volume and select high values. */
  for ( nZ = 0; nZ < this->mnDimensionZ; nZ++ )
  {
    for ( nY = 0; nY < this->mnDimensionY; nY++ )
    {
      for ( nX = 0; nX < this->mnDimensionX; nX++ )
      {
        xVoxl_Set( &MRI, nX, nY, nZ );
        Volm_GetValueAtMRIIdx_( this, &MRI, &value );
        if ( value > this->mpMaxValuesX[nY][nZ] )
        {
          this->mpMaxValuesX[nY][nZ] = value;
        }
        if ( value > this->mpMaxValuesY[nX][nZ] )
        {
          this->mpMaxValuesY[nX][nZ] = value;
        }
        if ( value > this->mpMaxValuesZ[nX][nY] )
        {
          this->mpMaxValuesZ[nX][nY] = value;
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

Volm_tErr Volm_MakeColorTable ( mriVolumeRef this )
{

  Volm_tErr   eResult    = Volm_tErr_NoErr;
  float       min        = 0;
  float       max        = 0;
  float       thresh     = 0;
  float       squash     = 0;
  int         entry      = 0;
  float       value      = 0;
  float       fComponent = 0;

  DebugEnterFunction( ("Volm_MakeColorTable( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* for each value... */
  min = this->mfColorMin;
  max = this->mfColorMax;
  thresh = this->mfColorBrightness;
  squash = this->mfColorContrast;
  if ( FEQUAL( min, max ) )
  {
    for ( entry = 0; entry < Volm_knNumColorTableEntries; entry++ )
    {
      this->mafColorTable[entry].mfRed   = 0;
      this->mafColorTable[entry].mfGreen = 0;
      this->mafColorTable[entry].mfBlue  = 0;
      this->manColorTable[entry].mnRed   = 0;
      this->manColorTable[entry].mnGreen = 0;
      this->manColorTable[entry].mnBlue  = 0;
    }
  }
  else
  {

    for ( entry = 0; entry < Volm_knNumColorTableEntries; entry++ )
    {

      /* Get the corresponding value for this entry, where entry 0 = min
      and entry Volm_knNumColorTableEntries = max. */
      value = ((entry * (max - min)) / Volm_knNumColorTableEntries) + min;

      /* This is a standard sigma 1 / (1+ exp((x-thresh)*-squash) )
      where x=(value-min)/(max-min) to normalize it to 0-1 within the
      range of min->max. This function will populate the 256
      lookuptable entries with intensity values from min->max so
      that by lowering the range of min->max, you increase the
      granularity of intensities in the visible values. */
      fComponent =
        (1.0 / (1.0 + exp( (((value-min)/(max-min))-thresh) * -squash )));

      /* set the float color */
      this->mafColorTable[entry].mfRed   = fComponent;
      this->mafColorTable[entry].mfGreen = fComponent;
      this->mafColorTable[entry].mfBlue  = fComponent;

      /* normalize to 0 - 255 pixel value */
      fComponent = (float)fComponent * 255.0;

      /* set the integer color */
      this->manColorTable[entry].mnRed   = (int)fComponent;
      this->manColorTable[entry].mnGreen = (int)fComponent;
      this->manColorTable[entry].mnBlue  = (int)fComponent;

    }
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_SetBrightnessAndContrast ( mriVolumeRef this,
    float        ifBrightness,
    float        ifContrast )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetBrightnessAndContrast( this=%p, "
                       "ifBrightness=%.2f, ifContrast=%.2f )", this,
                       ifBrightness, ifContrast) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

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

Volm_tErr Volm_SetColorMinMax ( mriVolumeRef this,
                                float        ifMin,
                                float        ifMax )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetColorMinMax( this=%p, "
                       "ifMin=%.2f, ifMax=%.2f )", this, ifMin, ifMax) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* No longer make the check so that a script can set the values
     outside of the bounds to be consistent over different volumes. */
#if 0
  DebugNote( ("Checking if min is okay") );
  DebugAssertThrowX( (rint(ifMin) >= rint(this->mfMinValue) &&
                      rint(ifMin) <= rint(this->mfMaxValue)),
                     eResult, Volm_tErr_InvalidParamater );
  DebugNote( ("Checking if max is okay") );
  DebugAssertThrowX( (rint(ifMax) >= rint(this->mfMinValue) &&
                      rint(ifMax) <= rint(this->mfMaxValue)),
                     eResult, Volm_tErr_InvalidParamater );
#endif

  /* set the values */
  this->mfColorMin = ifMin;
  this->mfColorMax = ifMax;

  /* recalc the color table */
  eResult = Volm_MakeColorTable( this );
  DebugAssertThrow( (Volm_tErr_NoErr == eResult) );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_GetBrightnessAndContrast ( mriVolumeRef this,
    float*       ofBrightness,
    float*       ofContrast )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_getBrightnessAndContrast( this=%p, "
                       "ofBrightness=%p, ofContrast=%p )", this,
                       ofBrightness, ofContrast) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );
  DebugAssertThrowX( (NULL != ofBrightness && NULL != ofContrast),
                     eResult, Volm_tErr_InvalidParamater );

  /* return the values */
  *ofBrightness = this->mfColorBrightness;
  *ofContrast   = this->mfColorContrast;

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_SaveToSnapshot ( mriVolumeRef this )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  int       nSize   = 0;

  DebugEnterFunction( ("Volm_SaveToSnapshot( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* if we have snapshot data, delete it */
  if ( NULL != this->mpSnapshot )
  {
    free( this->mpSnapshot );
  }

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
  memmove( this->mpSnapshot, this->mpMriValues->slices, nSize );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_RestoreFromSnapshot ( mriVolumeRef this )
{

  Volm_tErr eResult = Volm_tErr_NoErr;
  int       nSize   = 0;

  DebugEnterFunction( ("Volm_RestoreFromSnapshot( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* if we don't have a snapshot, return. */
  if ( NULL == this->mpSnapshot )
  {
    DebugGotoCleanup;
  }

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
  memmove( this->mpMriValues->slices, this->mpSnapshot, nSize );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_Rotate ( mriVolumeRef     this,
                        mri_tOrientation iAxis,
                        float            ifDegrees )
{

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
  switch ( iAxis )
  {
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
                           float        iNewValue )
{

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
  while ( xVoxl_IncrementUntilLimits( &idx, this->mnDimensionX-1,
                                      this->mnDimensionY-1,
                                      this->mnDimensionZ-1) )
  {

    Volm_GetValueAtMRIIdx_( this, &idx, &value );

    /* if we're going above and this value is above the thresh, or if we're
       going below and this value is below the thresh...*/
    if ( (ibAbove && value > iThreshold) ||
         (!ibAbove && value < iThreshold ) )
    {

      /* set the value to the new level. */
      Volm_SetValueAtMRIIdx_( this, &idx, iNewValue );
    }
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


Volm_tErr Volm_Flip ( mriVolumeRef     this,
                      mri_tOrientation iAxis )
{

  Volm_tErr   eResult       = Volm_tErr_NoErr;
  float       value         = 0;
  float       flippedValue  = 0;
  xVoxel      idx;
  xVoxel      flippedIdx;
  int         nXLimit       = 0;
  int         nYLimit       = 0;
  int         nZLimit       = 0;

  flippedIdx.mfX=0;
  flippedIdx.mfY=0;
  flippedIdx.mfZ=0;

  DebugEnterFunction( ("Volm_Flip( this=%p, iAxis=%d )", this, (int)iAxis) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* get our limits */
  nXLimit = this->mnDimensionX-1 ;
  nYLimit = this->mnDimensionY-1 ;
  nZLimit = this->mnDimensionZ-1;
  switch ( iAxis )
  {
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
  while ( xVoxl_IncrementUntilLimits( &idx, nXLimit, nYLimit, nZLimit ))
  {

    switch ( iAxis )
    {
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

Volm_tErr Volm_SetAllValues ( mriVolumeRef this,
                              float        iNewValue )
{

  Volm_tErr   eResult = Volm_tErr_NoErr;

  DebugEnterFunction( ("Volm_SetAllValues( this=%p, iNewValue=%d )",
                       this, (int)iNewValue) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugNote( ("Calling MRIsetValues") );
  MRIsetValues( this->mpMriValues, nint(iNewValue) );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_CopySubjectName ( mriVolumeRef this,
                                 char*        oSubjectName,
                                 int          inDestLen )
{

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
                                int          inDestLen )
{

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
                               int          inDestLen )
{

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

Volm_tErr Volm_ ( mriVolumeRef this )
{

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
                                char*        isName )
{

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
                               char*        isName )
{

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
    char*        isSource )
{

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

  /* look for /mri in the title */
  DebugNote( ("Looking for /mri in source") );
  sWord = strstr( isSource, "/mri" );
  if ( NULL != sWord )
  {

    /* we're at the slash now. go to the slash before this one. while
       we're not at the beginning.. */
    DebugNote( ("Going to the slash before /mri") );
    while ( sWord != isSource )
    {
      sWord -= sizeof( char );
      if ( *sWord == '/' )
      {
        break;
      }
    }

    /* inc past the slash and use the next part as the name */
    sWord += sizeof( char );
    nWordChar = 0;
    DebugNote( ("Copying in the part after the slash before /mri to /mri") );
    while ( *sWord != '/' &&
            *sWord != '\0' )
    {
      sName[nWordChar] = *sWord;
      nWordChar++;
      sWord += sizeof( char );
    }

  }
  else
  {

    /* look for 'subjects' in the title */
    DebugNote( ("Looking for subjects/ in source") );
    sWord = strstr( isSource, "subjects/" );
    if ( NULL != sWord )
    {

      /* we're at the s in subjects now. scoot ahead to the slash. */
      DebugNote( ("Scooting ahead to the slash") );
      nChar = 0;
      while ( sWord[nChar] != '/' &&
              sWord[nChar] != '\0' )
      {
        nChar++;
      }

      /* if found, use the next part as the name */
      DebugNote( ("Copying in chars from subjects/ to the next /") );
      nWordChar = 0;
      nChar++; /* get past the slash */
      while ( sWord[nChar] != '/' &&
              sWord[nChar] != '\0' )
      {
        sName[nWordChar] = sWord[nChar];
        nWordChar++;
        nChar++;
      }

    }
    else
    {

      /* else just use the last part */
      nChar = 0;
      DebugNote( ("Scooting to last slash") );
      while ( isSource[nChar] != '\0' )
      {
        if ( isSource[nChar] == '/' )
        {
          nLastSlash = nChar;
        }
        nChar++;
      }

      /* if we got it, use it, else use the whole source. */
      if ( isSource[nLastSlash] == '/' )
      {
        DebugNote( ("Copying everything from last slash to end of source") );
        xUtil_strncpy( sName, &(isSource[nLastSlash+1]), sizeof(sName) );
      }
      else
      {
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
    char*        isSource )
{

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
  while ( isSource[nChar] != '\0' )
  {
    if ( isSource[nChar] == '/' )
    {
      nLastSlash = nChar;
    }
    nChar++;
  }

  /* if we got one, everything from there+1 into the subjectname,
     otherwise just use the whole name */
  if ( isSource[nLastSlash] == '/' )
  {
    DebugNote( ("Copying in source name from char %d", nLastSlash+1) );
    xUtil_strncpy( this->msVolumeName, &(isSource[nLastSlash+1]),
                   sizeof( this->msVolumeName ) );
  }
  else
  {
    DebugNote( ("Copying in whole source name") );
    xUtil_strncpy( this->msVolumeName, isSource,sizeof( this->msVolumeName ) );
  }

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_SetMinVoxelSizeToOne ( mriVolumeRef this )
{

  Volm_tErr eResult    = Volm_tErr_NoErr;
  float scale, fov_x, fov_y, fov_z  ;

  DebugEnterFunction( ("Volm_SetMinVoxelSizeToOne( this=%p )", this ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Calculate the scale factor and size the xyzsize values. */
  scale = 1.0/MIN(MIN(this->mpMriValues->xsize,
                      this->mpMriValues->ysize),
                  this->mpMriValues->zsize) ;
  printf("scaling voxel sizes up by %2.2f\n", scale) ;
  this->mpMriValues->xsize *= scale ;
  this->mpMriValues->ysize *= scale ;
  this->mpMriValues->zsize *= scale ;

  /* Find the fov for each dimension. */
  fov_x = this->mpMriValues->xsize * this->mpMriValues->width;
  fov_y = this->mpMriValues->ysize * this->mpMriValues->height;
  fov_z = this->mpMriValues->zsize * this->mpMriValues->depth;

  /* Recalc the bounds. */
  this->mpMriValues->xend = fov_x / 2.0;
  this->mpMriValues->xstart = -this->mpMriValues->xend;
  this->mpMriValues->yend = fov_y / 2.0;
  this->mpMriValues->ystart = -this->mpMriValues->yend;
  this->mpMriValues->zend = fov_z / 2.0;
  this->mpMriValues->zstart = -this->mpMriValues->zend;

  /* Find the overall fov. */
  this->mpMriValues->fov = (fov_x > fov_y ?
                            (fov_x > fov_z ? fov_x : fov_z) :
                              (fov_y > fov_z ? fov_y : fov_z) );

  /* Set from the newly fov'd MRI. */
  DebugNote( ("Setting from MRI") );
  eResult = Volm_SetFromMRI_( this, this->mpMriValues );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  DebugCatch;
  DebugCatchError( eResult, Volm_tErr_NoErr, Volm_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

Volm_tErr Volm_Conform ( mriVolumeRef this )
{

  Volm_tErr eResult   = Volm_tErr_NoErr;
  MRI*      conformed = NULL;

  DebugEnterFunction( ("Volm_Conform( this=%p )", this ) );

  DebugNote( ("Verifying volume") );
  eResult = Volm_Verify( this );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

  /* Call the conform function. */
  DebugNote( ("Conforming MRI") );
  conformed = MRIconform( this->mpMriValues );
  DebugAssertThrowX( (NULL != conformed), eResult,
                     Volm_tErr_CouldntReadVolume );

  /* Set from the newly conformed MRI. */
  DebugNote( ("Setting from MRI") );
  eResult = Volm_SetFromMRI_( this, conformed );
  DebugAssertThrow( (eResult == Volm_tErr_NoErr) );

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
                                      xVoxelRef    oMRIIdx )
{

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
  if ( (xVoxl_GetX(oMRIIdx) < 0 ||
        xVoxl_GetX(oMRIIdx) >= this->mnDimensionX ||
        xVoxl_GetY(oMRIIdx) < 0 ||
        xVoxl_GetY(oMRIIdx) >= this->mnDimensionY ||
        xVoxl_GetZ(oMRIIdx) < 0 ||
        xVoxl_GetZ(oMRIIdx) >= this->mnDimensionZ) )
  {

    xVoxl_Set( oMRIIdx, 0, 0, 0 );
  }
}

void Volm_ConvertMRIIdxToScreenIdx_ ( mriVolumeRef this,
                                      xVoxelRef    iMRIIdx,
                                      xVoxelRef    oScreenIdx )
{

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
  if ( (xVoxl_GetX(oScreenIdx) < 0 ||
        xVoxl_GetX(oScreenIdx) >= 256 ||
        xVoxl_GetY(oScreenIdx) < 0 ||
        xVoxl_GetY(oScreenIdx) >= 256 ||
        xVoxl_GetZ(oScreenIdx) < 0 ||
        xVoxl_GetZ(oScreenIdx) >= 256) )
  {

    xVoxl_Set( oScreenIdx, 0, 0, 0 );
  }
}


void Volm_GetValueAtIdx_ ( mriVolumeRef this,
                           xVoxelRef    iIdx,
                           float*       oValue)
{

  int x, y, z;

  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to get the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  x = xVoxl_GetX(&this->mTmpVoxel);
  y = xVoxl_GetY(&this->mTmpVoxel);
  z = xVoxl_GetZ(&this->mTmpVoxel);

  switch ( this->mpMriValues->type )
  {
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

void Volm_GetValueAtMRIIdx_ ( mriVolumeRef this,
                              xVoxelRef    iMRIIdx,
                              float*       oValue)
{

  switch ( this->mpMriValues->type )
  {
  case MRI_UCHAR:
    *oValue =
      MRIvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
              xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) );
    break;
  case MRI_INT:
    *oValue =
      MRIIvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
               xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) );
    break;
  case MRI_LONG:
    *oValue =
      MRILvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
               xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) );
    break;
  case MRI_FLOAT:
    *oValue =
      MRIFvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
               xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) );
    break;
  case MRI_SHORT:
    *oValue =
      MRISvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
               xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) );
    break;
  default:
    *oValue = 0;
    break ;
  }
}


void Volm_GetValueAtIdxFrame_ ( mriVolumeRef this,
                                xVoxelRef    iIdx,
                                int          iFrame,
                                float*       oValue)
{

  /* First convert to MRI index, and then switch on the volume data
     type and use the proper MRIvox access function to get the
     value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  switch ( this->mpMriValues->type )
  {
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

void Volm_GetSampledValueAtIdx_ ( mriVolumeRef      this,
                                  xVoxelRef         iIdx,
                                  float*            oValue)
{

  double value;

  /* First convert to MRI index, then use the MRI access function to
     get the sinc value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  MRIsampleVolumeType( this->mpMriValues,
                       xVoxl_GetFloatX(&this->mTmpVoxel),
                       xVoxl_GetFloatY(&this->mTmpVoxel),
                       xVoxl_GetFloatZ(&this->mTmpVoxel),
                       &value,
                       (int)this->mSampleType );

  *oValue = (float)value;
}



void Volm_GetSampledValueAtIdxFrame_ ( mriVolumeRef      this,
                                       xVoxelRef         iIdx,
                                       int               iFrame,
                                       float*            oValue )
{

  double value;

  /* First convert to MRI index, then use the MRI access function to
     get the sinc value. */
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel );

  MRIsampleVolumeFrameType( this->mpMriValues,
                            xVoxl_GetFloatX(&this->mTmpVoxel),
                            xVoxl_GetFloatY(&this->mTmpVoxel),
                            xVoxl_GetFloatZ(&this->mTmpVoxel),
                            iFrame,
                            (int)this->mSampleType,
                            &value );

  *oValue = (float)value;
}



void Volm_SetValueAtIdx_ ( mriVolumeRef     this,
                           xVoxelRef        iIdx,
                           float            iValue )
{

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

void Volm_SetValueAtMRIIdx_ ( mriVolumeRef     this,
                              xVoxelRef        iMRIIdx,
                              float            iValue )
{

  switch (this->mpMriValues->type)
  {
  default:
    break ;
  case MRI_UCHAR:
    MRIvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
            xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) ) =
              (BUFTYPE) iValue;
    break ;
  case MRI_SHORT:
    MRISvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
             xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) ) =
               (short) iValue;
    break ;
  case MRI_FLOAT:
    MRIFvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
             xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) ) =
               (float) iValue;
    break ;
  case MRI_LONG:
    MRILvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
             xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) ) =
               (long) iValue;
    break ;
  case MRI_INT:
    MRIIvox( this->mpMriValues, xVoxl_GetX(iMRIIdx),
             xVoxl_GetY(iMRIIdx), xVoxl_GetZ(iMRIIdx) ) =
               (int) iValue;
    break ;
  }
}

void Volm_SetValueAtIdxFrame_ ( mriVolumeRef     this,
                                xVoxelRef        iIdx,
                                int          iFrame,
                                float            iValue )
{

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
#endif /* VOLM_USE_MACROS */
/* ======================================== END OF FUNCTION VERION OF MACROS */


Volm_tErr Volm_Verify ( mriVolumeRef this )
{

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
                           xVoxelRef    iIdx )
{

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
                            xVoxelRef    iIdx )
{

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

Volm_tErr Volm_VerifyMRIIdx_ ( mriVolumeRef this,
                               xVoxelRef    iMRIIdx )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugAssertThrowX( (xVoxl_GetX(iMRIIdx) >= 0 &&
                      xVoxl_GetRoundX(iMRIIdx) < this->mnDimensionX &&
                      xVoxl_GetY(iMRIIdx) >= 0 &&
                      xVoxl_GetRoundY(iMRIIdx) < this->mnDimensionY &&
                      xVoxl_GetZ(iMRIIdx) >= 0 &&
                      xVoxl_GetRoundZ(iMRIIdx) < this->mnDimensionZ),
                     eResult, Volm_tErr_InvalidIdx );

  DebugCatch;
  EndDebugCatch;
  return eResult;
}

Volm_tErr Volm_VerifyIdxInMRIBounds ( mriVolumeRef this,
                                      xVoxelRef    iScreenIdx )
{

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
                              int          iFrame )
{

  Volm_tErr eResult = Volm_tErr_NoErr;

  DebugAssertThrowX( (iFrame >= 0 && iFrame < this->mnDimensionFrame),
                     eResult, Volm_tErr_InvalidIdx );

  DebugCatch;
  EndDebugCatch;

  return eResult;
}

char* Volm_GetErrorString ( Volm_tErr ieCode )
{

  Volm_tErr eCode = ieCode;

  if ( ieCode < 0 ||
       ieCode >= Volm_knNumErrorCodes )
  {
    eCode = Volm_tErr_Invalid;
  }

  return Volm_ksaErrorStrings[eCode];
}
