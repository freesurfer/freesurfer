/**
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <math.h>
#include <stdlib.h>

#include <math.h>
#include <stdlib.h>
#include "macros.h"

#include "mriTransform.h"

const char *Trns_ksaErrorStrings[Trns_knNumErrorCodes] = {
    "No error.",
    "Invalid pointer to object.",
    "Invalid signature.",
    "Invalid parameter.",
    "Allocation failed.",
    "Tranformation matrices not inited.",
    "LTA import (LTAread) failed.",
    "Invalid error code."};

Trns_tErr Trns_New(mriTransformRef *opTransform)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  mriTransformRef ref = NULL;

  /* allocate us */
  ref = (mriTransformRef)malloc(sizeof(mriTransform));
  if (NULL == ref) {
    eResult = Trns_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  ref->mSignature = Trns_kSignature;

  /* set our matrices to null. */
  ref->mAtoRAS = NULL;
  ref->mBtoRAS = NULL;
  ref->mARAStoBRAS = NULL;
  ref->mRAStoA = NULL;
  ref->mRAStoB = NULL;
  ref->mBRAStoARAS = NULL;
  ref->mAtoB = NULL;
  ref->mBtoA = NULL;

  /* allocate our temp matricies */
  ref->mCoord1 = MatrixAlloc(4, 1, MATRIX_REAL);
  ref->mCoord2 = MatrixAlloc(4, 1, MATRIX_REAL);

  /* return us. */
  *opTransform = ref;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_New: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_NewFromLTA(mriTransformRef *opTransform, char *isLTAFileName)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  mriTransformRef ref = NULL;
  LTA *LTATransform = NULL;
  MATRIX *identity = NULL;

  /* allocate us */
  ref = (mriTransformRef)malloc(sizeof(mriTransform));
  if (NULL == ref) {
    eResult = Trns_tErr_AllocationFailed;
    goto error;
  }

  /* set signature */
  ref->mSignature = Trns_kSignature;

  /* set our matrices to null. */
  ref->mAtoRAS = NULL;
  ref->mBtoRAS = NULL;
  ref->mARAStoBRAS = NULL;
  ref->mRAStoA = NULL;
  ref->mRAStoB = NULL;
  ref->mBRAStoARAS = NULL;
  ref->mAtoB = NULL;
  ref->mBtoA = NULL;

  /* try to read a transform and make sure we got it */
  DebugNote(("Trying to read transform file %s", isLTAFileName));
  LTATransform = LTAreadEx(isLTAFileName);
  DebugAssertThrowX((NULL != LTATransform), eResult, Trns_tErr_LTAImportFailed);

  ref->type = LTATransform->type; /* if RAS will be converted
                                      to voxel later */

  /* copy the matrix out of it */
  Trns_CopyARAStoBRAS(ref, LTATransform->xforms[0].m_L);

  /* copy identities for the rest */
  identity = MatrixIdentity(4, NULL);
  Trns_CopyAtoRAS(ref, identity);
  Trns_CopyBtoRAS(ref, identity);

  /* allocate our temp matricies */
  ref->mCoord1 = MatrixAlloc(4, 1, MATRIX_REAL);
  ref->mCoord2 = MatrixAlloc(4, 1, MATRIX_REAL);

  /* return us. */
  *opTransform = ref;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_New: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  if (NULL != LTATransform) {
    LTAfree(&LTATransform);
  }

  if (NULL != identity) {
    MatrixFree(&identity);
  }

  return eResult;
}

Trns_tErr Trns_Delete(mriTransformRef *iopTransform)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  mriTransformRef ref = NULL;

  /* get us. */
  ref = *iopTransform;

  /* verify us. */
  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* delete any matricies we have. */
  if (NULL != ref->mAtoRAS) {
    MatrixFree(&(ref->mAtoRAS));
  }
  if (NULL != ref->mBtoRAS) {
    MatrixFree(&(ref->mBtoRAS));
  }
  if (NULL != ref->mARAStoBRAS) {
    MatrixFree(&(ref->mARAStoBRAS));
  }
  if (NULL != ref->mRAStoA) {
    MatrixFree(&(ref->mRAStoA));
  }
  if (NULL != ref->mRAStoB) {
    MatrixFree(&(ref->mRAStoB));
  }
  if (NULL != ref->mBRAStoARAS) {
    MatrixFree(&(ref->mBRAStoARAS));
  }
  if (NULL != ref->mAtoB) {
    MatrixFree(&(ref->mAtoB));
  }
  if (NULL != ref->mBtoA) {
    MatrixFree(&(ref->mBtoA));
  }
  MatrixFree(&(ref->mCoord1));
  MatrixFree(&(ref->mCoord2));

  /* trash our sig */
  ref->mSignature = 0x1;

  /* delete us */
  free(ref);
  ref = NULL;

  /* return null. */
  *iopTransform = NULL;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_Delete: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_DeepClone(mriTransformRef ref, mriTransformRef *opTransform)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  mriTransformRef clone = NULL;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* allocate the clone */
  eResult = Trns_New(&clone);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* copy any of our matrices that are defined */
  if (NULL != ref->mAtoRAS) {
    Trns_CopyAtoRAS(clone, ref->mAtoRAS);
  }
  if (NULL != ref->mBtoRAS) {
    Trns_CopyBtoRAS(clone, ref->mBtoRAS);
  }
  if (NULL != ref->mARAStoBRAS) {
    Trns_CopyARAStoBRAS(clone, ref->mARAStoBRAS);
  }

  /* return the clone. */
  *opTransform = clone;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_DeepClone: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyAtoRAS(mriTransformRef ref, MATRIX *iAtoRAS)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* create a copy of ref matrix */
  ref->mAtoRAS = MatrixCopy(iAtoRAS, NULL);

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_(ref);
  if (eResult != Trns_tErr_NoErr) {
    Trns_Signal("Trns_CopyAtoRAS", __LINE__, eResult);
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_CopyAtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyBtoRAS(mriTransformRef ref, MATRIX *iBtoRAS)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* create a copy of ref matrix */
  ref->mBtoRAS = MatrixCopy(iBtoRAS, NULL);

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_(ref);
  if (eResult != Trns_tErr_NoErr) {
    Trns_Signal("Trns_CopyBtoRAS", __LINE__, eResult);
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_CopyBtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyARAStoBRAS(mriTransformRef ref, MATRIX *iARAStoBRAS)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* create a copy of ref matrix */
  ref->mARAStoBRAS = MatrixCopy(iARAStoBRAS, NULL);

  /* calc the rest of our matricies */
  eResult = Trns_CalcMatricies_(ref);
  if (eResult != Trns_tErr_NoErr) {
    Trns_Signal("Trns_CopyARAStoBRAS", __LINE__, eResult);
    eResult = Trns_tErr_NoErr;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_CopyARAStoBRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetAtoRAS(mriTransformRef ref, MATRIX **opMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the matrix */
  *opMatrix = ref->mAtoRAS;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetAtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetBtoRAS(mriTransformRef ref, MATRIX **opMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the matrix */
  *opMatrix = ref->mBtoRAS;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetBtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_GetARAStoBRAS(mriTransformRef ref, MATRIX **opMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the matrix */
  *opMatrix = ref->mARAStoBRAS;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetARAStoBRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}
Trns_tErr Trns_GetAtoB(mriTransformRef ref, MATRIX **opMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the matrix */
  *opMatrix = ref->mAtoB;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetAtoB: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}
Trns_tErr Trns_GetBtoA(mriTransformRef ref, MATRIX **opMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the matrix */
  *opMatrix = ref->mBtoA;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetBtoA: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

#define mp(l, m)              \
  fprintf(stderr, "%s\n", l); \
  MatrixPrint(stderr, m);

Trns_tErr Trns_ApplyTransform(mriTransformRef ref, MATRIX *iTransform)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX *mTranslation = NULL;
  MATRIX *mTranslationInv = NULL;
  MATRIX *mRotation = NULL;
  MATRIX *mRotationInv = NULL;
  MATRIX *mScale = NULL;
  MATRIX *mNewTranslation = NULL;
  MATRIX *mNewRotation = NULL;
  MATRIX *mNewScale = NULL;
  MATRIX *mTmp1 = NULL;
  MATRIX *mTmp2 = NULL;
  MATRIX *mNew = NULL;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* init our matricies */
  mTranslation = MatrixAlloc(4, 4, MATRIX_REAL);
  mTranslationInv = MatrixAlloc(4, 4, MATRIX_REAL);
  mRotation = MatrixAlloc(4, 4, MATRIX_REAL);
  mRotationInv = MatrixAlloc(4, 4, MATRIX_REAL);
  mScale = MatrixAlloc(4, 4, MATRIX_REAL);
  mNewTranslation = MatrixAlloc(4, 4, MATRIX_REAL);
  mNewRotation = MatrixAlloc(4, 4, MATRIX_REAL);
  mNew = MatrixAlloc(4, 4, MATRIX_REAL);
  mTmp1 = MatrixAlloc(4, 4, MATRIX_REAL);
  mTmp2 = MatrixAlloc(4, 4, MATRIX_REAL);
  mNewScale = MatrixAlloc(4, 4, MATRIX_REAL);

  //  mp("orig",ref->mARAStoBRAS);

  /* get our matrix components */
  Trns_ExtractTranslationMatrix(ref->mARAStoBRAS, mTranslation);
  Trns_ExtractRotationMatrix(ref->mARAStoBRAS, mRotation);
  Trns_ExtractScaleMatrix(ref->mARAStoBRAS, mScale);
  //  mp("t",mTranslation);
  //  mp("r",mRotation);
  //  mp("s",mScale);

  /* we find the inverse rotation so we can apply our new transforms in a
     screen-relative space instead of local space. */
  /* find the inverse rotation */
  MatrixInverse(mRotation, mRotationInv);
  MatrixInverse(mTranslation, mTranslationInv);

  /* new_t = r * new_t * inv_r */
  Trns_ExtractTranslationMatrix(iTransform, mNewTranslation);
  MatrixMultiply(mNewTranslation, mRotationInv, mTmp1);
  MatrixMultiply(mRotation, mTmp1, mNewTranslation);
  //  mp("new_t",mNewTranslation);

  /* new_r = r * new_r * inv_r */
  Trns_ExtractRotationMatrix(iTransform, mNewRotation);
  MatrixMultiply(mNewRotation, mRotationInv, mTmp1);
  MatrixMultiply(mRotation, mTmp1, mNewRotation);
  //  mp("new_r",mNewRotation);

  Trns_ExtractScaleMatrix(iTransform, mNewScale);
  /*
  MatrixMultiply( mRotation, mNewScale, mTmp1 );
  MatrixMultiply( mTmp1, mRotationInv, mNewScale );
  */
  //  mp("new_s",mNewScale);

  /* compose them back together in the proper order with the new ones */
  MatrixIdentity(4, mNew);

  /* new_t t new_r r new_s s */
  MatrixMultiply(mNewScale, mScale, mTmp1);
  MatrixMultiply(mRotation, mTmp1, mTmp2);
  MatrixMultiply(mNewRotation, mTmp2, mTmp1);
  MatrixMultiply(mTranslation, mTmp1, mTmp2);
  MatrixMultiply(mNewTranslation, mTmp2, mNew);

  /* new_t t new_s s new_r r */
  MatrixMultiply(mNewRotation, mRotation, mTmp1);
  MatrixMultiply(mScale, mTmp1, mTmp2);
  MatrixMultiply(mNewScale, mTmp2, mTmp1);
  MatrixMultiply(mTranslation, mTmp1, mTmp2);
  MatrixMultiply(mNewTranslation, mTmp2, mNew);

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS(ref, mNew);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ApplyTransform: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  if (mTranslation) {
    MatrixFree(&mTranslation);
  }
  if (mTranslationInv) {
    MatrixFree(&mTranslationInv);
  }
  if (mRotation) {
    MatrixFree(&mRotation);
  }
  if (mRotationInv) {
    MatrixFree(&mRotationInv);
  }
  if (mScale) {
    MatrixFree(&mScale);
  }
  if (mNewTranslation) {
    MatrixFree(&mNewTranslation);
  }
  if (mNewRotation) {
    MatrixFree(&mNewRotation);
  }
  if (mNewScale) {
    MatrixFree(&mNewScale);
  }
  if (mNew) {
    MatrixFree(&mNew);
  }
  if (mTmp1) {
    MatrixFree(&mTmp1);
  }
  if (mTmp2) {
    MatrixFree(&mTmp2);
  }

  return eResult;
}

Trns_tErr Trns_ExtractTranslationMatrix(MATRIX *iTransform, MATRIX *oTranslation)
{
  if (NULL == iTransform || NULL == oTranslation) {
    return Trns_tErr_InvalidParameter;
  }

  /* identity matrix with just the translation components */
  MatrixIdentity(4, oTranslation);
  *MATRIX_RELT(oTranslation, 1, 4) = *MATRIX_RELT(iTransform, 1, 4);
  *MATRIX_RELT(oTranslation, 2, 4) = *MATRIX_RELT(iTransform, 2, 4);
  *MATRIX_RELT(oTranslation, 3, 4) = *MATRIX_RELT(iTransform, 3, 4);

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractRotationMatrix(MATRIX *iTransform, MATRIX *oRotation)
{
  MATRIX *mScale = NULL;
  MATRIX *mTmp = NULL;
  int n = 0;

  if (NULL == iTransform || NULL == oRotation) {
    return Trns_tErr_InvalidParameter;
  }

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy(iTransform, NULL);
  *MATRIX_RELT(mTmp, 1, 4) = 0;
  *MATRIX_RELT(mTmp, 2, 4) = 0;
  *MATRIX_RELT(mTmp, 3, 4) = 0;

  /* we want to cancel out the scale portion of ref matrix. so now we'll
     extract the scale matrix from ref matrix */
  mScale = MatrixAlloc(4, 4, MATRIX_REAL);
  Trns_ExtractScaleMatrix(mTmp, mScale);

  /* now modify it so that all the factors are one-overed. */
  for (n = 1; n <= 3; n++)
    if (!FZERO(*MATRIX_RELT(mScale, n, n))) {
      *MATRIX_RELT(mScale, n, n) = 1.0 / *MATRIX_RELT(mScale, n, n);
    }
    else {
      *MATRIX_RELT(mScale, n, n) = 1.0;
    }

  /* rotation is that matrix composed with original matrix sans translation */
  MatrixMultiply(mTmp, mScale, oRotation);

  MatrixFree(&mTmp);
  MatrixFree(&mScale);

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_ExtractScaleMatrix(MATRIX *iTransform, MATRIX *oScale)
{
  MATRIX *mTmp = NULL;
  VECTOR *vTmp = NULL;
  float fFactor = 0;

  if (NULL == iTransform || NULL == oScale) {
    return Trns_tErr_InvalidParameter;
  }

  /* create a matrix identical to the transform but without the translation */
  mTmp = MatrixCopy(iTransform, mTmp);
  *MATRIX_RELT(mTmp, 1, 4) = 0;
  *MATRIX_RELT(mTmp, 2, 4) = 0;
  *MATRIX_RELT(mTmp, 3, 4) = 0;

  /* find the x component. create a 1,0,0 vector and compose it. */
  vTmp = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(vTmp, 1) = 1.0;
  VECTOR_ELT(vTmp, 2) = 0;
  VECTOR_ELT(vTmp, 3) = 0;
  VECTOR_ELT(vTmp, 4) = 0.0;
  MatrixMultiply(mTmp, vTmp, vTmp);

  /* the x scale factor is the lengh of the transformed vector */
  fFactor = VectorLen(vTmp);
  MatrixIdentity(4, oScale);
  *MATRIX_RELT(oScale, 1, 1) = fFactor;

  /* do the same for y with a 0,1,0 vector */
  VECTOR_ELT(vTmp, 1) = 0;
  VECTOR_ELT(vTmp, 2) = 1.0;
  VECTOR_ELT(vTmp, 3) = 0;
  VECTOR_ELT(vTmp, 4) = 0.0;
  MatrixMultiply(mTmp, vTmp, vTmp);
  fFactor = VectorLen(vTmp);
  *MATRIX_RELT(oScale, 2, 2) = fFactor;

  /* do the same for z with a 0,0,1 vector */
  VECTOR_ELT(vTmp, 1) = 0;
  VECTOR_ELT(vTmp, 2) = 0;
  VECTOR_ELT(vTmp, 3) = 1.00;
  VECTOR_ELT(vTmp, 4) = 0.0;
  MatrixMultiply(mTmp, vTmp, vTmp);
  fFactor = VectorLen(vTmp);
  *MATRIX_RELT(oScale, 3, 3) = fFactor;

  MatrixFree(&mTmp);
  VectorFree(&vTmp);

  return Trns_tErr_NoErr;
}

Trns_tErr Trns_Translate(mriTransformRef ref, float ifAmount, tAxis iAxis)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX *mTransform = NULL;
  MATRIX *mOld = NULL;
  MATRIX *mNew = NULL;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  mTransform = MatrixAlloc(4, 4, MATRIX_REAL);
  mOld = MatrixAlloc(4, 4, MATRIX_REAL);
  mNew = MatrixAlloc(4, 4, MATRIX_REAL);

  /* create the matrix */
  mTransform = MatrixIdentity(4, NULL);
  switch (iAxis) {
    case tAxis_X:
      *MATRIX_RELT(mTransform, 1, 4) = ifAmount;
      break;
    case tAxis_Y:
      *MATRIX_RELT(mTransform, 2, 4) = ifAmount;
      break;
    case tAxis_Z:
      *MATRIX_RELT(mTransform, 3, 4) = ifAmount;
      break;
    default:
      break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy(ref->mARAStoBRAS, mOld);
  MatrixMultiply(mOld, mTransform, mNew);

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS(ref, mNew);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_Translate: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  if (mTransform != NULL) {
    MatrixFree(&mTransform);
  }
  if (mOld != NULL) {
    MatrixFree(&mNew);
  }
  if (mNew != NULL) {
    MatrixFree(&mNew);
  }

  return eResult;
}

Trns_tErr Trns_Rotate(mriTransformRef ref, float ifDegrees, tAxis iAxis)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX *mTransform = NULL;
  MATRIX *mOld = NULL;
  MATRIX *mNew = NULL;
  float fRadians = 0;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  mTransform = MatrixAlloc(4, 4, MATRIX_REAL);
  mOld = MatrixAlloc(4, 4, MATRIX_REAL);
  mNew = MatrixAlloc(4, 4, MATRIX_REAL);

  /* create the matrix */
  fRadians = ifDegrees * M_PI / 180.0;
  switch (iAxis) {
    case tAxis_X:
      mTransform = MatrixAllocRotation(4, fRadians, X_ROTATION);
      break;
    case tAxis_Y:
      mTransform = MatrixAllocRotation(4, fRadians, Y_ROTATION);
      break;
    case tAxis_Z:
      mTransform = MatrixAllocRotation(4, fRadians, Z_ROTATION);
      break;
    default:
      break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy(ref->mARAStoBRAS, mOld);
  MatrixMultiply(mOld, mTransform, mNew);

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS(ref, mNew);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_Rotate: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_Scale(mriTransformRef ref, float ifFactor, tAxis iAxis)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX *mTransform = NULL;
  MATRIX *mOld = NULL;
  MATRIX *mNew = NULL;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  mTransform = MatrixAlloc(4, 4, MATRIX_REAL);
  mOld = MatrixAlloc(4, 4, MATRIX_REAL);
  mNew = MatrixAlloc(4, 4, MATRIX_REAL);

  /* create the matrix */
  mTransform = MatrixIdentity(4, NULL);
  switch (iAxis) {
    case tAxis_X:
      *MATRIX_RELT(mTransform, 1, 1) = ifFactor;
      break;
    case tAxis_Y:
      *MATRIX_RELT(mTransform, 2, 2) = ifFactor;
      break;
    case tAxis_Z:
      *MATRIX_RELT(mTransform, 3, 3) = ifFactor;
      break;
    default:
      break;
  }

  /* compose the new matrix with the old one */
  MatrixCopy(ref->mARAStoBRAS, mOld);
  MatrixMultiply(mOld, mTransform, mNew);

  /* set the new matrix */
  eResult = Trns_CopyARAStoBRAS(ref, mNew);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_Scale: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertAtoB(mriTransformRef ref, xVoxelRef iAVoxel, xVoxelRef oBVoxel)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mAtoB) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(ref->mCoord1, 1, 1) = xVoxl_GetFloatX(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 2, 1) = xVoxl_GetFloatY(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 3, 1) = xVoxl_GetFloatZ(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 4, 1) = 1;

  /* do the transform */
  MatrixMultiply(ref->mAtoB, ref->mCoord1, ref->mCoord2);

  /* set the voxel to the matrix */
  xVoxl_SetFloat(
      oBVoxel, *MATRIX_RELT(ref->mCoord2, 1, 1), *MATRIX_RELT(ref->mCoord2, 2, 1), *MATRIX_RELT(ref->mCoord2, 3, 1));

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertAtoB: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertAtoRAS(mriTransformRef ref, xVoxelRef iAVoxel, xVoxelRef oRASVoxel)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mAtoRAS) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(ref->mCoord1, 1, 1) = xVoxl_GetFloatX(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 2, 1) = xVoxl_GetFloatY(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 3, 1) = xVoxl_GetFloatZ(iAVoxel);
  *MATRIX_RELT(ref->mCoord1, 4, 1) = 1;

  /* do the transform */
  MatrixMultiply(ref->mAtoRAS, ref->mCoord1, ref->mCoord2);

  /* set the voxel to the matrix */
  xVoxl_SetFloat(oRASVoxel,
                 *MATRIX_RELT(ref->mCoord2, 1, 1),
                 *MATRIX_RELT(ref->mCoord2, 2, 1),
                 *MATRIX_RELT(ref->mCoord2, 3, 1));

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertAtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoA(mriTransformRef ref, xVoxelRef iBVoxel, xVoxelRef oAVoxel)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mBtoA) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(ref->mCoord1, 1, 1) = xVoxl_GetFloatX(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 2, 1) = xVoxl_GetFloatY(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 3, 1) = xVoxl_GetFloatZ(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 4, 1) = 1;

  /* do the transform */
  MatrixMultiply(ref->mBtoA, ref->mCoord1, ref->mCoord2);

  /* set the voxel to the matrix */
  xVoxl_SetFloat(
      oAVoxel, *MATRIX_RELT(ref->mCoord2, 1, 1), *MATRIX_RELT(ref->mCoord2, 2, 1), *MATRIX_RELT(ref->mCoord2, 3, 1));

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertBtoA: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBtoRAS(mriTransformRef ref, xVoxelRef iBVoxel, xVoxelRef oRASVoxel)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mBtoRAS) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(ref->mCoord1, 1, 1) = xVoxl_GetFloatX(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 2, 1) = xVoxl_GetFloatY(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 3, 1) = xVoxl_GetFloatZ(iBVoxel);
  *MATRIX_RELT(ref->mCoord1, 4, 1) = 1;

  /* do the transform */
  MatrixMultiply(ref->mBtoRAS, ref->mCoord1, ref->mCoord2);

  /* set the voxel to the matrix */
  xVoxl_SetFloat(oRASVoxel,
                 *MATRIX_RELT(ref->mCoord2, 1, 1),
                 *MATRIX_RELT(ref->mCoord2, 2, 1),
                 *MATRIX_RELT(ref->mCoord2, 3, 1));

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertBtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertBRAStoB(mriTransformRef ref, xVoxelRef iBRASVoxel, xVoxelRef oBVoxel)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mRAStoB) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* set our coord matrix */
  *MATRIX_RELT(ref->mCoord1, 1, 1) = xVoxl_GetFloatX(iBRASVoxel);
  *MATRIX_RELT(ref->mCoord1, 2, 1) = xVoxl_GetFloatY(iBRASVoxel);
  *MATRIX_RELT(ref->mCoord1, 3, 1) = xVoxl_GetFloatZ(iBRASVoxel);
  *MATRIX_RELT(ref->mCoord1, 4, 1) = 1;

  /* do the transform */
  MatrixMultiply(ref->mRAStoB, ref->mCoord1, ref->mCoord2);

  /* set the voxel to the matrix */
  xVoxl_SetFloat(
      oBVoxel, *MATRIX_RELT(ref->mCoord2, 1, 1), *MATRIX_RELT(ref->mCoord2, 2, 1), *MATRIX_RELT(ref->mCoord2, 3, 1));

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertBRAStoB: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixAtoB(mriTransformRef ref, MATRIX *iAMatrix, MATRIX *oBMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mAtoB) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply(ref->mAtoB, iAMatrix, oBMatrix);

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertAtoB: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixAtoRAS(mriTransformRef ref, MATRIX *iAMatrix, MATRIX *oRASMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mAtoRAS) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply(ref->mAtoRAS, iAMatrix, oRASMatrix);

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertAtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoA(mriTransformRef ref, MATRIX *iBMatrix, MATRIX *oAMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mBtoA) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply(ref->mBtoA, iBMatrix, oAMatrix);

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertBtoA: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_ConvertMatrixBtoRAS(mriTransformRef ref, MATRIX *iBMatrix, MATRIX *oRASMatrix)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* if we don't have our trans matricies, return an error. */
  if (NULL == ref->mBtoRAS) {
    eResult = Trns_tErr_TransformationMatrixNotInited;
    goto error;
  }

  /* do the transform */
  MatrixMultiply(ref->mBtoRAS, iBMatrix, oRASMatrix);

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_ConvertBtoRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

/* don't use the next function when copying mAtoB or mBtoA */
/* they will be replaced by recalculated ones              */
Trns_tErr Trns_CalcMatricies_(mriTransformRef ref)
{
  Trns_tErr eResult = Trns_tErr_NoErr;
  MATRIX *tmp = NULL;

  /* for each of our first three matricies, if we have them, calc an
     inverse. */
  if (NULL != ref->mAtoRAS) {
    ref->mRAStoA = MatrixInverse(ref->mAtoRAS, NULL);
  }
  if (NULL != ref->mBtoRAS) {
    ref->mRAStoB = MatrixInverse(ref->mBtoRAS, NULL);
  }
  if (NULL != ref->mARAStoBRAS) {
    ref->mBRAStoARAS = MatrixInverse(ref->mARAStoBRAS, NULL);
  }

  /* if we have everything now, calc the composed conversions */
  if (NULL != ref->mAtoRAS && NULL != ref->mRAStoA && NULL != ref->mBtoRAS && NULL != ref->mRAStoB &&
      NULL != ref->mARAStoBRAS && NULL != ref->mBRAStoARAS) {
    /* a_index -> a_ras -> b_ras -> b_index */
    tmp = MatrixMultiply(ref->mARAStoBRAS, ref->mAtoRAS, NULL);
    ref->mAtoB = MatrixMultiply(ref->mRAStoB, tmp, NULL);

    /* b_index -> b_ras -> a_ras -> a_index */
    tmp = MatrixMultiply(ref->mBRAStoARAS, ref->mBtoRAS, NULL);
    ref->mBtoA = MatrixMultiply(ref->mRAStoA, tmp, NULL);
  }

  goto cleanup;

  goto error;
error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_CalcMatricies_: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_Verify(mriTransformRef ref)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  /* check for null ptr */
  if (NULL == ref) {
    eResult = Trns_tErr_InvalidPtr;
    goto cleanup;
  }

  /* check signature */
  if (Trns_kSignature != ref->mSignature) {
    eResult = Trns_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;
}

void Trns_Signal(const char *isFuncName, int inLineNum, Trns_tErr ieCode)
{
  DebugPrint(("Signal in %s, line %d: %d, %s", isFuncName, inLineNum, ieCode, Trns_GetErrorString(ieCode)));
}

const char *Trns_GetErrorString(Trns_tErr ieCode)
{
  Trns_tErr eCode = ieCode;

  if (ieCode < 0 || ieCode >= Trns_knNumErrorCodes) {
    eCode = Trns_tErr_InvalidErrorCode;
  }

  return Trns_ksaErrorStrings[eCode];
}

void Trns_DebugPrint_(mriTransformRef ref)
{
  if (NULL != ref->mAtoRAS) {
    DebugPrint(("A to RAS:\n"));
    MatrixPrint(stderr, ref->mAtoRAS);
  }

  if (NULL != ref->mBtoRAS) {
    DebugPrint(("B to RAS:\n"));
    MatrixPrint(stderr, ref->mBtoRAS);
  }

  if (NULL != ref->mARAStoBRAS) {
    DebugPrint(("ARAS to BRAS:\n"));
    MatrixPrint(stderr, ref->mARAStoBRAS);
  }

  if (NULL != ref->mRAStoA) {
    DebugPrint(("RAS to A:\n"));
    MatrixPrint(stderr, ref->mRAStoA);
  }

  if (NULL != ref->mRAStoB) {
    DebugPrint(("RAS to B:\n"));
    MatrixPrint(stderr, ref->mRAStoB);
  }

  if (NULL != ref->mBRAStoARAS) {
    DebugPrint(("BRAS to ARAS:\n"));
    MatrixPrint(stderr, ref->mBRAStoARAS);
  }

  if (NULL != ref->mAtoB) {
    DebugPrint(("A to B:\n"));
    MatrixPrint(stderr, ref->mAtoB);
  }
}
Trns_tErr Trns_GetType(mriTransformRef ref, int *ptype)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  /* return the type */
  *ptype = ref->type;

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_GetType: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}

Trns_tErr Trns_CopyAtoB(mriTransformRef ref, MATRIX *iAtoB)
{
  Trns_tErr eResult = Trns_tErr_NoErr;

  eResult = Trns_Verify(ref);
  if (Trns_tErr_NoErr != eResult) {
    goto error;
  }

  if (NULL != ref->mAtoB) {
    MatrixFree(&ref->mAtoB);
  }
  if (NULL != ref->mBtoA) {
    MatrixFree(&ref->mBtoA);
  }

  /* create a copy of ref matrix */
  ref->mAtoB = MatrixCopy(iAtoB, NULL);
  ref->mBtoA = MatrixInverse(iAtoB, NULL);

  DebugAssertThrowX((NULL != ref->mBtoA), eResult, Trns_tErr_LTAImportFailed);

  goto cleanup;

error:

  if (Trns_tErr_NoErr != eResult) {
    DebugPrint(("Error %d in Trns_CopyARAStoBRAS: %s\n", eResult, Trns_GetErrorString(eResult)));
  }

cleanup:

  return eResult;
}
