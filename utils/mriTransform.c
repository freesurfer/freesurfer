/**
 * @file  mriTransform.c
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.19 $
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


#include <stdlib.h>
#include <math.h>
#include "mriTransform.h"
#include "macros.h"

char *Trns_ksaErrorStrings [Trns_knNumErrorCodes] =
{

  "No error.",
  "Invalid pointer to object.",
  "Invalid signature.",
  "Invalid parameter.",
  "Allocation failed.",
  "Tranformation matrices not inited.",
  "LTA import (LTAread) failed.",
  "Invalid error code."
};


Trns_tErr Trns_New ( mriTransformRef* opTransform )
{

  Trns_tErr       eResult = Trns_tErr_NoErr;
  mriTransformRef this    = NULL;

  /* allocate us */
  this = (mriTransformRef) malloc( sizeof( mriTransform ) );
  if ( NULL == this )
  {
    eResult = Trns_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  this->mSignature = Trns_kSignature;

  /* set our matrices to null. */
  this->mAtoRAS     = NULL;
  this->mBtoRAS     = NULL;
  this->mARAStoBRAS = NULL;
  this->mRAStoA     = NULL;
  this->mRAStoB     = NULL;
  this->mBRAStoARAS = NULL;
  this->mAtoB       = NULL;
  this->mBtoA       = NULL;

  /* allocate our temp matricies */
  this->mCoord1 = MatrixAlloc( 4, 1, MATRIX_REAL );
  this->mCoord2 = MatrixAlloc( 4, 1, MATRIX_REAL );

  /* return us. */
  *opTransform = this;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_New: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_NewFromLTA ( mriTransformRef* opTransform,
                            char*            isLTAFileName )
{

  Trns_tErr       eResult     = Trns_tErr_NoErr;
  mriTransformRef this        = NULL;
  LTA*            LTATransform = NULL;
  MATRIX*         identity    = NULL;

  /* allocate us */
  this = (mriTransformRef) malloc( sizeof( mriTransform ) );
  if ( NULL == this )
  {
    eResult = Trns_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  this->mSignature = Trns_kSignature;

  /* set our matrices to null. */
  this->mAtoRAS     = NULL;
  this->mBtoRAS     = NULL;
  this->mARAStoBRAS = NULL;
  this->mRAStoA     = NULL;
  this->mRAStoB     = NULL;
  this->mBRAStoARAS = NULL;
  this->mAtoB       = NULL;
  this->mBtoA       = NULL;

  /* try to read a transform and make sure we got it */
  DebugNote( ("Trying to read transform file %s", isLTAFileName) );
  LTATransform = LTAreadEx( isLTAFileName );
  DebugAssertThrowX( (NULL != LTATransform),
                     eResult, Trns_tErr_LTAImportFailed );

  this->type = LTATransform->type ;  /* if RAS will be converted
                                        to voxel later */

  /* copy the matrix out of it */
  Trns_CopyARAStoBRAS( this, LTATransform->xforms[0].m_L );

  /* copy identities for the rest */
  identity = MatrixIdentity( 4, NULL );
  Trns_CopyAtoRAS( this, identity );
  Trns_CopyBtoRAS( this, identity );

  /* allocate our temp matricies */
  this->mCoord1 = MatrixAlloc( 4, 1, MATRIX_REAL );
  this->mCoord2 = MatrixAlloc( 4, 1, MATRIX_REAL );

  /* return us. */
  *opTransform = this;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_New: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  if ( NULL != LTATransform )
  {
    LTAfree( &LTATransform );
  }

  if ( NULL != identity )
  {
    MatrixFree( &identity );
  }

  return eResult;
}

Trns_tErr Trns_Delete ( mriTransformRef* iopTransform )
{

  Trns_tErr       eResult = Trns_tErr_NoErr;
  mriTransformRef this    = NULL;

  /* get us. */
  this = *iopTransform;

  /* verify us. */
  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* delete any matricies we have. */
  if ( NULL != this->mAtoRAS )
  {
    MatrixFree( &(this->mAtoRAS) );
  }
  if ( NULL != this->mBtoRAS )
  {
    MatrixFree( &(this->mBtoRAS) );
  }
  if ( NULL != this->mARAStoBRAS )
  {
    MatrixFree( &(this->mARAStoBRAS) );
  }
  if ( NULL != this->mRAStoA )
  {
    MatrixFree( &(this->mRAStoA) );
  }
  if ( NULL != this->mRAStoB )
  {
    MatrixFree( &(this->mRAStoB) );
  }
  if ( NULL != this->mBRAStoARAS )
  {
    MatrixFree( &(this->mBRAStoARAS) );
  }
  if ( NULL != this->mAtoB )
  {
    MatrixFree( &(this->mAtoB) );
  }
  if ( NULL != this->mBtoA )
  {
    MatrixFree( &(this->mBtoA) );
  }
  MatrixFree( &(this->mCoord1) );
  MatrixFree( &(this->mCoord2) );

  /* trash our sig */
  this->mSignature = 0x1;

  /* delete us */
  free( this );
  this = NULL;

  /* return null. */
  *iopTransform = NULL;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_Delete: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_DeepClone ( mriTransformRef  this,
                           mriTransformRef* opTransform )
{

  Trns_tErr       eResult = Trns_tErr_NoErr;
  mriTransformRef clone   = NULL;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* allocate the clone */
  eResult = Trns_New( &clone );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* copy any of our matrices that are defined */
  if ( NULL != this->mAtoRAS )
  {
    Trns_CopyAtoRAS( clone, this->mAtoRAS );
  }
  if ( NULL != this->mBtoRAS )
  {
    Trns_CopyBtoRAS( clone, this->mBtoRAS );
  }
  if ( NULL != this->mARAStoBRAS )
  {
    Trns_CopyARAStoBRAS( clone, this->mARAStoBRAS );
  }

  /* return the clone. */
  *opTransform = clone;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_DeepClone: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyAtoRAS ( mriTransformRef this,
                            MATRIX*         iAtoRAS )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* create a copy of this matrix */
  this->mAtoRAS = MatrixCopy( iAtoRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if ( eResult != Trns_tErr_NoErr )
  {
    Trns_Signal( "Trns_CopyAtoRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_CopyAtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyBtoRAS ( mriTransformRef this,
                            MATRIX*         iBtoRAS )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* create a copy of this matrix */
  this->mBtoRAS = MatrixCopy( iBtoRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if ( eResult != Trns_tErr_NoErr )
  {
    Trns_Signal( "Trns_CopyBtoRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_CopyBtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyARAStoBRAS ( mriTransformRef this,
                                MATRIX*         iARAStoBRAS )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* create a copy of this matrix */
  this->mARAStoBRAS = MatrixCopy( iARAStoBRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if ( eResult != Trns_tErr_NoErr )
  {
    Trns_Signal( "Trns_CopyARAStoBRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_CopyARAStoBRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetAtoRAS ( mriTransformRef this,
                           MATRIX**        opMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the matrix */
  *opMatrix = this->mAtoRAS;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetAtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetBtoRAS ( mriTransformRef this,
                           MATRIX**         opMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the matrix */
  *opMatrix = this->mBtoRAS;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetBtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetARAStoBRAS ( mriTransformRef this,
                               MATRIX**        opMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the matrix */
  *opMatrix = this->mARAStoBRAS;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetARAStoBRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}
Trns_tErr Trns_GetAtoB ( mriTransformRef this,
                         MATRIX**        opMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the matrix */
  *opMatrix = this->mAtoB;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetAtoB: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}
Trns_tErr Trns_GetBtoA ( mriTransformRef this,
                         MATRIX**        opMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the matrix */
  *opMatrix = this->mBtoA;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetBtoA: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

#define mp(l,m) fprintf(stderr, "%s\n", l ); MatrixPrint(stderr,m);

Trns_tErr Trns_ApplyTransform ( mriTransformRef this,
                                MATRIX*         iTransform )
{

  Trns_tErr eResult         = Trns_tErr_NoErr;
  MATRIX*   mTranslation    = NULL;
  MATRIX*   mTranslationInv = NULL;
  MATRIX*   mRotation       = NULL;
  MATRIX*   mRotationInv    = NULL;
  MATRIX*   mScale          = NULL;
  MATRIX*   mNewTranslation = NULL;
  MATRIX*   mNewRotation    = NULL;
  MATRIX*   mNewScale       = NULL;
  MATRIX*   mTmp1           = NULL;
  MATRIX*   mTmp2           = NULL;
  MATRIX*   mNew            = NULL;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* init our matricies */
  mTranslation    = MatrixAlloc( 4, 4, MATRIX_REAL );
  mTranslationInv = MatrixAlloc( 4, 4, MATRIX_REAL );
  mRotation       = MatrixAlloc( 4, 4, MATRIX_REAL );
  mRotationInv    = MatrixAlloc( 4, 4, MATRIX_REAL );
  mScale          = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewTranslation = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewRotation    = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNew            = MatrixAlloc( 4, 4, MATRIX_REAL );
  mTmp1           = MatrixAlloc( 4, 4, MATRIX_REAL );
  mTmp2           = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewScale       = MatrixAlloc( 4, 4, MATRIX_REAL );

  //  mp("orig",this->mARAStoBRAS);

  /* get our matrix components */
  Trns_ExtractTranslationMatrix( this->mARAStoBRAS, mTranslation );
  Trns_ExtractRotationMatrix(    this->mARAStoBRAS, mRotation );
  Trns_ExtractScaleMatrix(       this->mARAStoBRAS, mScale );
  //  mp("t",mTranslation);
  //  mp("r",mRotation);
  //  mp("s",mScale);

  /* we find the inverse rotation so we can apply our new transforms in a
     screen-relative space instead of local space. */
  /* find the inverse rotation */
  MatrixInverse( mRotation, mRotationInv );
  MatrixInverse( mTranslation, mTranslationInv );

  /* new_t = r * new_t * inv_r */
  Trns_ExtractTranslationMatrix( iTransform, mNewTranslation );
  MatrixMultiply( mNewTranslation, mRotationInv, mTmp1 );
  MatrixMultiply( mRotation, mTmp1, mNewTranslation );
  //  mp("new_t",mNewTranslation);

  /* new_r = r * new_r * inv_r */
  Trns_ExtractRotationMatrix( iTransform, mNewRotation );
  MatrixMultiply( mNewRotation, mRotationInv, mTmp1 );
  MatrixMultiply( mRotation, mTmp1, mNewRotation );
  //  mp("new_r",mNewRotation);

  Trns_ExtractScaleMatrix( iTransform, mNewScale );
  /*
  MatrixMultiply( mRotation, mNewScale, mTmp1 );
  MatrixMultiply( mTmp1, mRotationInv, mNewScale );
  */
  //  mp("new_s",mNewScale);

  /* compose them back together in the proper order with the new ones */
  MatrixIdentity( 4, mNew );

  /* new_t t new_r r new_s s */
  MatrixMultiply( mNewScale, mScale, mTmp1 );
  MatrixMultiply( mRotation, mTmp1, mTmp2 );
  MatrixMultiply( mNewRotation, mTmp2, mTmp1 );
  MatrixMultiply( mTranslation, mTmp1, mTmp2 );
  MatrixMultiply( mNewTranslation, mTmp2, mNew );

  /* new_t t new_s s new_r r */
  MatrixMultiply( mNewRotation, mRotation, mTmp1 );
  MatrixMultiply( mScale, mTmp1, mTmp2 );
  MatrixMultiply( mNewScale, mTmp2, mTmp1 );
  MatrixMultiply( mTranslation, mTmp1, mTmp2 );
  MatrixMultiply( mNewTranslation, mTmp2, mNew );

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS( this, mNew );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ApplyTransform: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  if ( mTranslation )
  {
    MatrixFree( &mTranslation );
  }
  if ( mTranslationInv )
  {
    MatrixFree( &mTranslationInv );
  }
  if ( mRotation )
  {
    MatrixFree( &mRotation );
  }
  if ( mRotationInv )
  {
    MatrixFree( &mRotationInv );
  }
  if ( mScale )
  {
    MatrixFree( &mScale );
  }
  if ( mNewTranslation )
  {
    MatrixFree( &mNewTranslation );
  }
  if ( mNewRotation )
  {
    MatrixFree( &mNewRotation );
  }
  if ( mNewScale )
  {
    MatrixFree( &mNewScale );
  }
  if ( mNew )
  {
    MatrixFree( &mNew );
  }
  if ( mTmp1 )
  {
    MatrixFree( &mTmp1 );
  }
  if ( mTmp2 )
  {
    MatrixFree( &mTmp2 );
  }

  return eResult;
}

Trns_tErr Trns_ExtractTranslationMatrix ( MATRIX* iTransform,
    MATRIX* oTranslation )
{

  if ( NULL == iTransform
       || NULL == oTranslation )
  {
    return Trns_tErr_InvalidParameter;
  }

  /* identity matrix with just the translation components */
  MatrixIdentity( 4, oTranslation );
  *MATRIX_RELT(oTranslation,1,4) = *MATRIX_RELT(iTransform,1,4);
  *MATRIX_RELT(oTranslation,2,4) = *MATRIX_RELT(iTransform,2,4);
  *MATRIX_RELT(oTranslation,3,4) = *MATRIX_RELT(iTransform,3,4);

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractRotationMatrix    ( MATRIX* iTransform,
    MATRIX* oRotation )
{

  MATRIX* mScale  = NULL;
  MATRIX* mTmp    = NULL;
  int     n       = 0;

  if ( NULL == iTransform
       || NULL == oRotation )
  {
    return Trns_tErr_InvalidParameter;
  }

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy( iTransform, NULL );
  *MATRIX_RELT(mTmp,1,4) = 0;
  *MATRIX_RELT(mTmp,2,4) = 0;
  *MATRIX_RELT(mTmp,3,4) = 0;

  /* we want to cancel out the scale portion of this matrix. so now we'll
     extract the scale matrix from this matrix */
  mScale = MatrixAlloc( 4, 4, MATRIX_REAL );
  Trns_ExtractScaleMatrix( mTmp, mScale );

  /* now modify it so that all the factors are one-overed. */
  for ( n = 1; n <= 3; n ++ )
    if ( !FZERO(*MATRIX_RELT(mScale,n,n)) )
    {
      *MATRIX_RELT(mScale,n,n) = 1.0 / *MATRIX_RELT(mScale,n,n);
    }
    else
    {
      *MATRIX_RELT(mScale,n,n) = 1.0;
    }

  /* rotation is that matrix composed with original matrix sans translation */
  MatrixMultiply( mTmp, mScale, oRotation );

  MatrixFree( &mTmp );
  MatrixFree( &mScale );

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractScaleMatrix ( MATRIX* iTransform,
                                    MATRIX* oScale )
{

  MATRIX* mTmp    = NULL;
  VECTOR* vTmp    = NULL;
  float   fFactor = 0;

  if ( NULL == iTransform
       || NULL == oScale )
  {
    return Trns_tErr_InvalidParameter;
  }

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy( iTransform, mTmp );
  *MATRIX_RELT(mTmp,1,4) = 0;
  *MATRIX_RELT(mTmp,2,4) = 0;
  *MATRIX_RELT(mTmp,3,4) = 0;

  /* find the x component. create a 1,0,0 vector and compose it. */
  vTmp = VectorAlloc(4,MATRIX_REAL);
  VECTOR_ELT(vTmp,1) = 1.0;
  VECTOR_ELT(vTmp,2) = 0;
  VECTOR_ELT(vTmp,3) = 0;
  VECTOR_ELT(vTmp,4) = 0.0;
  MatrixMultiply(mTmp,vTmp,vTmp);

  /* the x scale factor is the lengh of the transformed vector */
  fFactor = VectorLen( vTmp );
  MatrixIdentity( 4, oScale );
  *MATRIX_RELT(oScale,1,1) = fFactor;

  /* do the same for y with a 0,1,0 vector */
  VECTOR_ELT(vTmp,1) = 0;
  VECTOR_ELT(vTmp,2) = 1.0;
  VECTOR_ELT(vTmp,3) = 0;
  VECTOR_ELT(vTmp,4) = 0.0;
  MatrixMultiply(mTmp,vTmp,vTmp);
  fFactor = VectorLen( vTmp );
  *MATRIX_RELT(oScale,2,2) = fFactor;

  /* do the same for z with a 0,0,1 vector */
  VECTOR_ELT(vTmp,1) = 0;
  VECTOR_ELT(vTmp,2) = 0;
  VECTOR_ELT(vTmp,3) = 1.00;
  VECTOR_ELT(vTmp,4) = 0.0;
  MatrixMultiply(mTmp,vTmp,vTmp);
  fFactor = VectorLen( vTmp );
  *MATRIX_RELT(oScale,3,3) = fFactor;

  MatrixFree( &mTmp );
  VectorFree( &vTmp );

  return Trns_tErr_NoErr;
}


Trns_tErr Trns_Translate ( mriTransformRef this,
                           float           ifAmount,
                           tAxis           iAxis )
{


  Trns_tErr eResult    = Trns_tErr_NoErr;
  MATRIX*   mTransform = NULL;
  MATRIX*   mOld       = NULL;
  MATRIX*   mNew       = NULL;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  mTransform = MatrixAlloc( 4, 4, MATRIX_REAL );
  mOld = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNew = MatrixAlloc( 4, 4, MATRIX_REAL );

  /* create the matrix */
  mTransform = MatrixIdentity( 4, NULL );
  switch ( iAxis )
  {
  case tAxis_X:
    *MATRIX_RELT( mTransform, 1, 4 ) = ifAmount;
    break;
  case tAxis_Y:
    *MATRIX_RELT( mTransform, 2, 4 ) = ifAmount;
    break;
  case tAxis_Z:
    *MATRIX_RELT( mTransform, 3, 4 ) = ifAmount;
    break;
  default:
    break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy( this->mARAStoBRAS, mOld );
  MatrixMultiply( mOld, mTransform, mNew );

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS( this, mNew );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_Translate: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  if ( mTransform != NULL )
  {
    MatrixFree( &mTransform );
  }
  if ( mOld != NULL )
  {
    MatrixFree( &mNew );
  }
  if ( mNew != NULL )
  {
    MatrixFree( &mNew );
  }

  return eResult;
}

Trns_tErr Trns_Rotate ( mriTransformRef this,
                        float           ifDegrees,
                        tAxis           iAxis )
{

  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX*   mTransform = NULL;
  MATRIX*   mOld       = NULL;
  MATRIX*   mNew       = NULL;
  float     fRadians  = 0;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  mTransform = MatrixAlloc( 4, 4, MATRIX_REAL );
  mOld = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNew = MatrixAlloc( 4, 4, MATRIX_REAL );

  /* create the matrix */
  fRadians = ifDegrees * M_PI / 180.0;
  switch ( iAxis )
  {
  case tAxis_X:
    mTransform = MatrixAllocRotation( 4, fRadians, X_ROTATION );
    break;
  case tAxis_Y:
    mTransform = MatrixAllocRotation( 4, fRadians, Y_ROTATION );
    break;
  case tAxis_Z:
    mTransform = MatrixAllocRotation( 4, fRadians, Z_ROTATION );
    break;
  default:
    break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy( this->mARAStoBRAS, mOld );
  MatrixMultiply( mOld, mTransform, mNew );

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS( this, mNew );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_Rotate: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_Scale ( mriTransformRef this,
                       float           ifFactor,
                       tAxis           iAxis )
{

  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX*   mTransform = NULL;
  MATRIX*   mOld       = NULL;
  MATRIX*   mNew       = NULL;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  mTransform = MatrixAlloc( 4, 4, MATRIX_REAL );
  mOld = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNew = MatrixAlloc( 4, 4, MATRIX_REAL );

  /* create the matrix */
  mTransform = MatrixIdentity( 4, NULL );
  switch ( iAxis )
  {
  case tAxis_X:
    *MATRIX_RELT( mTransform, 1, 1 ) = ifFactor;
    break;
  case tAxis_Y:
    *MATRIX_RELT( mTransform, 2, 2 ) = ifFactor;
    break;
  case tAxis_Z:
    *MATRIX_RELT( mTransform, 3, 3 ) = ifFactor;
    break;
  default:
    break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy( this->mARAStoBRAS, mOld );
  MatrixMultiply( mOld, mTransform, mNew );

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS( this, mNew );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_Scale: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Trns_tErr Trns_ConvertAtoB ( mriTransformRef this,
                             xVoxelRef       iAVoxel,
                             xVoxelRef       oBVoxel )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mAtoB )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(this->mCoord1,1,1) = xVoxl_GetFloatX( iAVoxel );
  *MATRIX_RELT(this->mCoord1,2,1) = xVoxl_GetFloatY( iAVoxel );
  *MATRIX_RELT(this->mCoord1,3,1) = xVoxl_GetFloatZ( iAVoxel );
  *MATRIX_RELT(this->mCoord1,4,1) = 1;

  /* do the transform */
  MatrixMultiply( this->mAtoB, this->mCoord1, this->mCoord2 );

  /* set the voxel to the matrix */
  xVoxl_SetFloat( oBVoxel,
                  *MATRIX_RELT(this->mCoord2,1,1),
                  *MATRIX_RELT(this->mCoord2,2,1),
                  *MATRIX_RELT(this->mCoord2,3,1) );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertAtoB: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertAtoRAS ( mriTransformRef this,
                               xVoxelRef       iAVoxel,
                               xVoxelRef       oRASVoxel )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mAtoRAS )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(this->mCoord1,1,1) = xVoxl_GetFloatX( iAVoxel );
  *MATRIX_RELT(this->mCoord1,2,1) = xVoxl_GetFloatY( iAVoxel );
  *MATRIX_RELT(this->mCoord1,3,1) = xVoxl_GetFloatZ( iAVoxel );
  *MATRIX_RELT(this->mCoord1,4,1) = 1;

  /* do the transform */
  MatrixMultiply( this->mAtoRAS, this->mCoord1, this->mCoord2 );

  /* set the voxel to the matrix */
  xVoxl_SetFloat( oRASVoxel,
                  *MATRIX_RELT(this->mCoord2,1,1),
                  *MATRIX_RELT(this->mCoord2,2,1),
                  *MATRIX_RELT(this->mCoord2,3,1) );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertAtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoA ( mriTransformRef this,
                             xVoxelRef       iBVoxel,
                             xVoxelRef       oAVoxel )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mBtoA )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(this->mCoord1,1,1) = xVoxl_GetFloatX( iBVoxel );
  *MATRIX_RELT(this->mCoord1,2,1) = xVoxl_GetFloatY( iBVoxel );
  *MATRIX_RELT(this->mCoord1,3,1) = xVoxl_GetFloatZ( iBVoxel );
  *MATRIX_RELT(this->mCoord1,4,1) = 1;

  /* do the transform */
  MatrixMultiply( this->mBtoA, this->mCoord1, this->mCoord2 );

  /* set the voxel to the matrix */
  xVoxl_SetFloat( oAVoxel,
                  *MATRIX_RELT(this->mCoord2,1,1),
                  *MATRIX_RELT(this->mCoord2,2,1),
                  *MATRIX_RELT(this->mCoord2,3,1) );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertBtoA: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoRAS ( mriTransformRef this,
                               xVoxelRef       iBVoxel,
                               xVoxelRef       oRASVoxel )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mBtoRAS )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(this->mCoord1,1,1) = xVoxl_GetFloatX( iBVoxel );
  *MATRIX_RELT(this->mCoord1,2,1) = xVoxl_GetFloatY( iBVoxel );
  *MATRIX_RELT(this->mCoord1,3,1) = xVoxl_GetFloatZ( iBVoxel );
  *MATRIX_RELT(this->mCoord1,4,1) = 1;

  /* do the transform */
  MatrixMultiply( this->mBtoRAS, this->mCoord1, this->mCoord2 );

  /* set the voxel to the matrix */
  xVoxl_SetFloat( oRASVoxel,
                  *MATRIX_RELT(this->mCoord2,1,1),
                  *MATRIX_RELT(this->mCoord2,2,1),
                  *MATRIX_RELT(this->mCoord2,3,1) );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertBtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Trns_tErr Trns_ConvertBRAStoB ( mriTransformRef this,
                                xVoxelRef       iBRASVoxel,
                                xVoxelRef       oBVoxel )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mRAStoB )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(this->mCoord1,1,1) = xVoxl_GetFloatX( iBRASVoxel );
  *MATRIX_RELT(this->mCoord1,2,1) = xVoxl_GetFloatY( iBRASVoxel );
  *MATRIX_RELT(this->mCoord1,3,1) = xVoxl_GetFloatZ( iBRASVoxel );
  *MATRIX_RELT(this->mCoord1,4,1) = 1;

  /* do the transform */
  MatrixMultiply(this->mRAStoB, this->mCoord1, this->mCoord2 );

  /* set the voxel to the matrix */
  xVoxl_SetFloat( oBVoxel,
                  *MATRIX_RELT(this->mCoord2,1,1),
                  *MATRIX_RELT(this->mCoord2,2,1),
                  *MATRIX_RELT(this->mCoord2,3,1) );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertBRAStoB: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


Trns_tErr Trns_ConvertMatrixAtoB ( mriTransformRef this,
                                   MATRIX*         iAMatrix,
                                   MATRIX*         oBMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mAtoB )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mAtoB, iAMatrix, oBMatrix );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertAtoB: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixAtoRAS ( mriTransformRef this,
                                     MATRIX*         iAMatrix,
                                     MATRIX*         oRASMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mAtoRAS )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mAtoRAS, iAMatrix, oRASMatrix );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertAtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoA ( mriTransformRef this,
                                   MATRIX*         iBMatrix,
                                   MATRIX*         oAMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mBtoA )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mBtoA, iBMatrix, oAMatrix );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertBtoA: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoRAS ( mriTransformRef this,
                                     MATRIX*         iBMatrix,
                                     MATRIX*         oRASMatrix )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if ( NULL == this->mBtoRAS )
  {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mBtoRAS, iBMatrix, oRASMatrix );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_ConvertBtoRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

/* don't use the next function when copying mAtoB or mBtoA */
/* they will be replaced by recalculated ones              */
Trns_tErr Trns_CalcMatricies_ ( mriTransformRef this )
{

  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX*   tmp     = NULL;

  /* for each of our first three matricies, if we have them, calc an
     inverse. */
  if ( NULL != this->mAtoRAS )
  {
    this->mRAStoA = MatrixInverse( this->mAtoRAS, NULL );
  }
  if ( NULL != this->mBtoRAS )
  {
    this->mRAStoB = MatrixInverse( this->mBtoRAS, NULL );
  }
  if ( NULL != this->mARAStoBRAS )
  {
    this->mBRAStoARAS = MatrixInverse( this->mARAStoBRAS, NULL );
  }

  /* if we have everything now, calc the composed conversions */
  if ( NULL != this->mAtoRAS
       && NULL != this->mRAStoA
       && NULL != this->mBtoRAS
       && NULL != this->mRAStoB
       && NULL != this->mARAStoBRAS
       && NULL != this->mBRAStoARAS )
  {

    /* a_index -> a_ras -> b_ras -> b_index */
    tmp = MatrixMultiply( this->mARAStoBRAS, this->mAtoRAS, NULL );
    this->mAtoB = MatrixMultiply( this->mRAStoB, tmp, NULL );

    /* b_index -> b_ras -> a_ras -> a_index */
    tmp = MatrixMultiply( this->mBRAStoARAS, this->mBtoRAS, NULL );
    this->mBtoA = MatrixMultiply( this->mRAStoA, tmp, NULL );
  }

  goto cleanup;

  goto error;
error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_CalcMatricies_: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_Verify ( mriTransformRef this )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this )
  {
    eResult = Trns_tErr_InvalidPtr;
    goto cleanup;
  }

  /* check signature */
  if ( Trns_kSignature != this->mSignature )
  {
    eResult = Trns_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;
}

void Trns_Signal ( char* isFuncName, int inLineNum, Trns_tErr ieCode )
{

  DebugPrint( ("Signal in %s, line %d: %d, %s",
               isFuncName, inLineNum, ieCode, Trns_GetErrorString(ieCode) ) );
}

char* Trns_GetErrorString ( Trns_tErr ieCode )
{

  Trns_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= Trns_knNumErrorCodes )
  {
    eCode = Trns_tErr_InvalidErrorCode;
  }

  return Trns_ksaErrorStrings [eCode];
}

void Trns_DebugPrint_ ( mriTransformRef this )
{

  if ( NULL != this->mAtoRAS )
  {
    DebugPrint( ("A to RAS:\n" ) );
    MatrixPrint( stderr, this->mAtoRAS );
  }

  if ( NULL != this->mBtoRAS )
  {
    DebugPrint( ("B to RAS:\n" ) );
    MatrixPrint( stderr, this->mBtoRAS );
  }

  if ( NULL != this->mARAStoBRAS )
  {
    DebugPrint( ("ARAS to BRAS:\n" ) );
    MatrixPrint( stderr, this->mARAStoBRAS );
  }

  if ( NULL != this->mRAStoA )
  {
    DebugPrint( ("RAS to A:\n" ) );
    MatrixPrint( stderr, this->mRAStoA );
  }

  if ( NULL != this->mRAStoB )
  {
    DebugPrint( ("RAS to B:\n" ) );
    MatrixPrint( stderr, this->mRAStoB );
  }

  if ( NULL != this->mBRAStoARAS )
  {
    DebugPrint( ("BRAS to ARAS:\n" ) );
    MatrixPrint( stderr, this->mBRAStoARAS );
  }

  if ( NULL != this->mAtoB )
  {
    DebugPrint( ("A to B:\n" ) );
    MatrixPrint( stderr, this->mAtoB );
  }
}
Trns_tErr Trns_GetType     ( mriTransformRef this, int *ptype)
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* return the type */
  *ptype = this->type;

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_GetType: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyAtoB ( mriTransformRef this,
                          MATRIX*         iAtoB )
{

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if ( Trns_tErr_NoErr != eResult )
  {
    goto error;
  }

  if (NULL != this->mAtoB)
  {
    MatrixFree(&this->mAtoB) ;
  }
  if (NULL != this->mBtoA)
  {
    MatrixFree(&this->mBtoA) ;
  }

  /* create a copy of this matrix */
  this->mAtoB = MatrixCopy( iAtoB, NULL );
  this->mBtoA = MatrixInverse( iAtoB, NULL );

  DebugAssertThrowX( (NULL != this->mBtoA),
                     eResult, Trns_tErr_LTAImportFailed );

  goto cleanup;

error:

  if ( Trns_tErr_NoErr != eResult )
  {
    DebugPrint( ("Error %d in Trns_CopyARAStoBRAS: %s\n",
                 eResult, Trns_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}
