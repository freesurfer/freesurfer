#include <stdlib.h>
#include "mriTransform.h"

char Trns_ksaErrorStrings [Trns_knNumErrorCodes][256] = {

  "No error.",
  "Invalid pointer to object.",
  "Invalid signature.",
  "Invalid parameter.",
  "Allocation failed.",
  "Tranformation matrices not inited.",
  "Invalid error code."
};


Trns_tErr Trns_New ( mriTransformRef* opTransform ) {

  Trns_tErr       eResult = Trns_tErr_NoErr;
  mriTransformRef this    = NULL;

  /* allocate us */
  this = (mriTransformRef) malloc( sizeof( mriTransform ) );
  if( NULL == this ) {
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_New: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_Delete ( mriTransformRef* iopTransform ) {

  Trns_tErr       eResult = Trns_tErr_NoErr;
  mriTransformRef this    = NULL;

  /* get us. */
  this = *iopTransform;

  /* verify us. */
  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult ) 
    goto error;

  /* delete any matricies we have. */
  if( NULL != this->mAtoRAS )
    MatrixFree( &(this->mAtoRAS) );
  if( NULL != this->mBtoRAS )
    MatrixFree( &(this->mBtoRAS) );
  if( NULL != this->mARAStoBRAS )
    MatrixFree( &(this->mARAStoBRAS) );
  if( NULL != this->mRAStoA )
    MatrixFree( &(this->mRAStoA) );
  if( NULL != this->mRAStoB )
    MatrixFree( &(this->mRAStoB) );
  if( NULL != this->mBRAStoARAS )
    MatrixFree( &(this->mBRAStoARAS) );
  if( NULL != this->mAtoB )
    MatrixFree( &(this->mAtoB) );
  if( NULL != this->mBtoA )
    MatrixFree( &(this->mBtoA) );
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_Delete: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_CopyAtoRAS ( mriTransformRef this,
          MATRIX*         iAtoRAS ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* create a copy of this matrix */
  this->mAtoRAS = MatrixCopy( iAtoRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if( eResult != Trns_tErr_NoErr ) {
    Trns_Signal( "Trns_CopyAtoRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_CopyAtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_CopyBtoRAS ( mriTransformRef this,
          MATRIX*         iBtoRAS ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* create a copy of this matrix */
  this->mBtoRAS = MatrixCopy( iBtoRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if( eResult != Trns_tErr_NoErr ) {
    Trns_Signal( "Trns_CopyBtoRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_CopyBtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_CopyARAStoBRAS ( mriTransformRef this,
        MATRIX*         iARAStoBRAS ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* create a copy of this matrix */
  this->mARAStoBRAS = MatrixCopy( iARAStoBRAS, NULL );

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_( this );
  if( eResult != Trns_tErr_NoErr ) {
    Trns_Signal( "Trns_CopyARAStoBRAS", __LINE__, eResult );
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_CopyARAStoBRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_GetAtoRAS ( mriTransformRef this,
         MATRIX**        opMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* return the matrix */
  *opMatrix = this->mAtoRAS;

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_GetAtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_GetBtoRAS ( mriTransformRef this,
         MATRIX**         opMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* return the matrix */
  *opMatrix = this->mBtoRAS;

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_GetBtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_GetARAStoBRAS ( mriTransformRef this,
             MATRIX**        opMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* return the matrix */
  *opMatrix = this->mARAStoBRAS;

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_GetARAStoBRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ApplyTransform ( mriTransformRef this,
        MATRIX*         iTransform ) {

  Trns_tErr eResult         = Trns_tErr_NoErr;
  MATRIX*   mTranslation    = NULL;
  MATRIX*   mRotation       = NULL;
  MATRIX*   mScale          = NULL;
  MATRIX*   mNewTranslation = NULL;
  MATRIX*   mNewRotation    = NULL;
  MATRIX*   mNewScale       = NULL;
  MATRIX*   mNew            = NULL;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* init our matricies */
  mTranslation    = MatrixAlloc( 4, 4, MATRIX_REAL );
  mRotation       = MatrixAlloc( 4, 4, MATRIX_REAL );
  mScale          = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewTranslation = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewRotation    = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNewScale       = MatrixAlloc( 4, 4, MATRIX_REAL );
  mNew            = MatrixAlloc( 4, 4, MATRIX_REAL );

  /* get our matrix components */
  Trns_ExtractTranslationMatrix( this->mARAStoBRAS, mTranslation );
  Trns_ExtractRotationMatrix(    this->mARAStoBRAS, mRotation );
  Trns_ExtractScaleMatrix(       this->mARAStoBRAS, mScale );
  Trns_ExtractTranslationMatrix( iTransform,        mNewTranslation );
  Trns_ExtractRotationMatrix(    iTransform,        mNewRotation );
  Trns_ExtractScaleMatrix(       iTransform,        mNewScale );

  /* compose them back together in the proper order with the
     new transforms */
  MatrixIdentity( 4, mNew );
  MatrixMultiply( mNewTranslation, mNew, mNew );
  MatrixMultiply( mTranslation,    mNew, mNew );
  MatrixMultiply( mNewScale,       mNew, mNew );
  MatrixMultiply( mScale,          mNew, mNew );
  MatrixMultiply( mNewRotation,    mNew, mNew );
  MatrixMultiply( mRotation,       mNew, mNew );

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS( this, mNew );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ApplyTransform: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  if( mTranslation ) 
    MatrixFree( &mTranslation );
  if( mRotation ) 
    MatrixFree( &mRotation );
  if( mScale ) 
    MatrixFree( &mScale );
  if( mNewTranslation ) 
    MatrixFree( &mNewTranslation );
  if( mNewRotation ) 
    MatrixFree( &mNewRotation );
  if( mNewScale ) 
    MatrixFree( &mNewScale );
  if( mNew ) 
    MatrixFree( &mNew );

  return eResult;
}

Trns_tErr Trns_ExtractTranslationMatrix ( MATRIX* iTransform,
            MATRIX* oTranslation ) {

  if( NULL == iTransform
      || NULL == oTranslation ) 
    return Trns_tErr_InvalidParameter;
  
  /* identity matrix with just the translation components */
  MatrixIdentity( 4, oTranslation );
  *MATRIX_RELT(oTranslation,1,4) = *MATRIX_RELT(iTransform,1,4);
  *MATRIX_RELT(oTranslation,2,4) = *MATRIX_RELT(iTransform,2,4);
  *MATRIX_RELT(oTranslation,3,4) = *MATRIX_RELT(iTransform,3,4);

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractRotationMatrix    ( MATRIX* iTransform,
            MATRIX* oRotation ) {

  MATRIX* mTmp    = NULL;
  MATRIX* mTmp2   = NULL;
  VECTOR* vTmp    = NULL;
  float   fFactor = 0;

  if( NULL == iTransform
      || NULL == oRotation ) 
    return Trns_tErr_InvalidParameter;

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy( iTransform, NULL );
  *MATRIX_RELT(mTmp,1,4) = 0;
  *MATRIX_RELT(mTmp,2,4) = 0;
  *MATRIX_RELT(mTmp,3,4) = 0;

  /* create a 1,0,0 vector and compose it. */
  vTmp = VectorAlloc(4,MATRIX_REAL);
  VECTOR_ELT(vTmp,1) = 1.0; VECTOR_ELT(vTmp,2) = 0;
  VECTOR_ELT(vTmp,3) = 0;   VECTOR_ELT(vTmp,4) = 1.0;
  MatrixMultiply(mTmp,vTmp,vTmp);
  
  /* create a matrix with 1/length in the diagonal */
  fFactor = 1.0 / VectorLen( vTmp );
  mTmp2 = MatrixIdentity( 4, mTmp2 );
  *MATRIX_RELT(mTmp2,1,1) = fFactor;
  *MATRIX_RELT(mTmp2,2,2) = fFactor;
  *MATRIX_RELT(mTmp2,3,3) = fFactor;

  /* rotation is that matrix composed with original matrix sans translation */
  MatrixMultiply( mTmp2, mTmp, oRotation );

  MatrixFree( &mTmp );
  MatrixFree( &mTmp2 );
  VectorFree( &vTmp );

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractScaleMatrix       ( MATRIX* iTransform,
            MATRIX* oScale ) {

  MATRIX* mTmp    = NULL;
  VECTOR* vTmp    = NULL;
  float   fFactor = 0;

  if( NULL == iTransform
      || NULL == oScale ) 
    return Trns_tErr_InvalidParameter;

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy( iTransform, mTmp );
  *MATRIX_RELT(mTmp,1,4) = 0;
  *MATRIX_RELT(mTmp,2,4) = 0;
  *MATRIX_RELT(mTmp,3,4) = 0;

  /* create a 1,0,0 vector and compose it. */
  vTmp = VectorAlloc(4,MATRIX_REAL);
  VECTOR_ELT(vTmp,1) = 1.0; VECTOR_ELT(vTmp,2) = 0;
  VECTOR_ELT(vTmp,3) = 0;   VECTOR_ELT(vTmp,4) = 1.0;
  MatrixMultiply(mTmp,vTmp,vTmp);

  /* the scale is an identity matrix with the length in the diagonal */
  fFactor = VectorLen( vTmp );
  MatrixIdentity( 4, oScale );
  *MATRIX_RELT(oScale,1,1) = fFactor;
  *MATRIX_RELT(oScale,2,2) = fFactor;
  *MATRIX_RELT(oScale,3,3) = fFactor;

  MatrixFree( &mTmp );
  VectorFree( &vTmp );

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ConvertAtoB ( mriTransformRef this,
           xVoxelRef       iAVoxel,
           xVoxelRef       oBVoxel ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mAtoB ) {
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertAtoB: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertAtoRAS ( mriTransformRef this,
           xVoxelRef       iAVoxel,
           xVoxelRef       oRASVoxel ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mAtoRAS ) {
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertAtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoA ( mriTransformRef this,
           xVoxelRef       iBVoxel,
           xVoxelRef       oAVoxel ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mBtoA ) {
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertBtoA: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoRAS ( mriTransformRef this,
           xVoxelRef       iBVoxel,
           xVoxelRef       oRASVoxel ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mBtoRAS ) {
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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertBtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}


Trns_tErr Trns_ConvertMatrixAtoB ( mriTransformRef this,
           MATRIX*         iAMatrix,
           MATRIX*         oBMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mAtoB ) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mAtoB, iAMatrix, oBMatrix );

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertAtoB: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixAtoRAS ( mriTransformRef this,
             MATRIX*         iAMatrix,
             MATRIX*         oRASMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mAtoRAS ) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mAtoRAS, iAMatrix, oRASMatrix );

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertAtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoA ( mriTransformRef this,
           MATRIX*         iBMatrix,
           MATRIX*         oAMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mBtoA ) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mBtoA, iBMatrix, oAMatrix );

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertBtoA: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoRAS ( mriTransformRef this,
             MATRIX*         iBMatrix,
             MATRIX*         oRASMatrix ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify( this );
  if( Trns_tErr_NoErr != eResult )
    goto error;

  /* if we don't have our trans matricies, return an error. */
  if( NULL == this->mBtoRAS ) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply( this->mBtoRAS, iBMatrix, oRASMatrix );

  goto cleanup;

 error:

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_ConvertBtoRAS: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_CalcMatricies_ ( mriTransformRef this ) {

  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX*   tmp     = NULL;

  /* for each of our first three matricies, if we have them, calc an
     inverse. */
  if( NULL != this->mAtoRAS )
    this->mRAStoA = MatrixInverse( this->mAtoRAS, NULL );
  if( NULL != this->mBtoRAS )
    this->mRAStoB = MatrixInverse( this->mBtoRAS, NULL );
  if( NULL != this->mARAStoBRAS )
    this->mBRAStoARAS = MatrixInverse( this->mARAStoBRAS, NULL );

  /* if we have everything now, calc the composed conversions */
  if( NULL != this->mAtoRAS
      && NULL != this->mRAStoA
      && NULL != this->mBtoRAS
      && NULL != this->mRAStoB
      && NULL != this->mARAStoBRAS
      && NULL != this->mBRAStoARAS ) {

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

  if( Trns_tErr_NoErr != eResult ) {
    DebugPrint "Error %d in Trns_CalcMatricies_: %s\n",
      eResult, Trns_GetErrorString( eResult ) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

Trns_tErr Trns_Verify ( mriTransformRef this ) {

  Trns_tErr eResult = Trns_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this ) {
    eResult = Trns_tErr_InvalidPtr;
    goto cleanup;
  }
  
  /* check signature */
  if ( Trns_kSignature != this->mSignature ) {
    eResult = Trns_tErr_InvalidSignature;
    goto cleanup;
  }

 cleanup:

  return eResult;
}

void Trns_Signal ( char* isFuncName, int inLineNum, Trns_tErr ieCode ) {

  DebugPrint "Signal in %s, line %d: %d, %s", 
    isFuncName, inLineNum, ieCode, Trns_GetErrorString(ieCode) EndDebugPrint;
}

char* Trns_GetErrorString ( Trns_tErr ieCode ) {

  Trns_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= Trns_knNumErrorCodes ) {
    eCode = Trns_tErr_InvalidErrorCode;
  }

  return Trns_ksaErrorStrings [eCode];
}

void Trns_DebugPrint_ ( mriTransformRef this ) {

  if( NULL != this->mAtoRAS ) {
    DebugPrint "A to RAS:\n" EndDebugPrint;
    MatrixPrint( stderr, this->mAtoRAS );
  }
  
  if( NULL != this->mBtoRAS ) {
    DebugPrint "B to RAS:\n" EndDebugPrint;
    MatrixPrint( stderr, this->mBtoRAS );
  }
  
  if( NULL != this->mARAStoBRAS ) {
    DebugPrint "ARAS to BRAS:\n" EndDebugPrint;
    MatrixPrint( stderr, this->mARAStoBRAS );
  }
  
  if( NULL != this->mAtoB ) {
    DebugPrint "A to B:\n" EndDebugPrint;
    MatrixPrint( stderr, this->mAtoB );
  }
}
