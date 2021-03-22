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


#ifndef mriTransform_h
#define mriTransform_h

#include "matrix.h"
#include "mri.h"    /* transform.h requires mri.h */
#include "transform.h"
#include "xDebug.h"
#include "xVoxel.h"

typedef enum
{

  Trns_tErr_NoErr = 0,
  Trns_tErr_InvalidPtr,
  Trns_tErr_InvalidSignature,
  Trns_tErr_InvalidParameter,
  Trns_tErr_AllocationFailed,
  Trns_tErr_TransformationMatrixNotInited,
  Trns_tErr_LTAImportFailed,
  Trns_tErr_InvalidErrorCode,
  Trns_knNumErrorCodes

} Trns_tErr;

typedef struct
{

  long mSignature;

  /* our various matrices. */
  MATRIX* mAtoRAS;     /* these go from a sepcific system to RAS */
  MATRIX* mBtoRAS;
  MATRIX* mARAStoBRAS; /* this goes from RAS in one system to
             RAS in the other */

  MATRIX* mRAStoA;     /* these are inverses of the above three */
  MATRIX* mRAStoB;
  MATRIX* mBRAStoARAS;

  MATRIX* mAtoB;       /* these are compositions of the first three and */
  MATRIX* mBtoA;       /* second three */

  MATRIX* mCoord1;     /* temp matricies for transforms. */
  MATRIX* mCoord2;
  int     type ;       /* RAS or voxel (see transform.h) */

}
mriTransform, *mriTransformRef;

#define Trns_kSignature 0x09876543

Trns_tErr Trns_New        ( mriTransformRef* opTransform );
Trns_tErr Trns_NewFromLTA ( mriTransformRef* opTransform,
                            char*            isLTAFileName );
Trns_tErr Trns_Delete     ( mriTransformRef* iopTransform );
Trns_tErr Trns_DeepClone  ( mriTransformRef  ThisKRT,
                            mriTransformRef* opTransform );

/* these are the minimum that should be set in a normal
   situation. Note that every time you copy one of these, the other
   matrices (the inverses of these tree and the compositions, AtoB and
   BtoA) will be recalculated in Trns_CalcMatricies_. */
Trns_tErr Trns_CopyAtoRAS     ( mriTransformRef ThisKRT,
                                MATRIX*         iAtoRAS );
Trns_tErr Trns_CopyBtoRAS     ( mriTransformRef ThisKRT,
                                MATRIX*         iBtoRAS );
Trns_tErr Trns_CopyARAStoBRAS ( mriTransformRef ThisKRT,
                                MATRIX*         iARAStoBRAS );

/* Alternatively, use ThisKRT if you only have one matrix to use. It will
   set AtoB and calc BtoA and then not touch the rest of the
   matrices. */
Trns_tErr Trns_CopyAtoB     ( mriTransformRef ThisKRT,
                              MATRIX*         iAtoB );

/* access internal matrices */
Trns_tErr Trns_GetAtoRAS     ( mriTransformRef ThisKRT,
                               MATRIX**        opMatrix );
Trns_tErr Trns_GetType     ( mriTransformRef ThisKRT,
                             int *ptype) ;
Trns_tErr Trns_GetBtoRAS     ( mriTransformRef ThisKRT,
                               MATRIX**        opMatrix );
Trns_tErr Trns_GetARAStoBRAS ( mriTransformRef ThisKRT,
                               MATRIX**        opMatrix );
Trns_tErr Trns_GetAtoB ( mriTransformRef ThisKRT,
                         MATRIX**        opMatrix );
Trns_tErr Trns_GetBtoA ( mriTransformRef ThisKRT,
                         MATRIX**        opMatrix );

/* apply a transformation to the ARAStoBRAS matrix */
Trns_tErr Trns_ApplyTransform ( mriTransformRef ThisKRT,
                                MATRIX*         iTransform );

/* extracts the individual portions of a transformation matrix. necessary
   for local transformation. does not allocate new matricies. note that the
   scale extraction assumes all scales have been uniform in the three axes. */
Trns_tErr Trns_ExtractTranslationMatrix ( MATRIX* iTransform,
    MATRIX* oTranslation );
Trns_tErr Trns_ExtractRotationMatrix    ( MATRIX* iTransform,
    MATRIX* oRotation );
Trns_tErr Trns_ExtractScaleMatrix       ( MATRIX* iTransform,
    MATRIX* oScale );

/* */
Trns_tErr Trns_Translate  ( mriTransformRef ThisKRT,
                            float           ifAmount,
                            tAxis           iAxis );
Trns_tErr Trns_Rotate     ( mriTransformRef ThisKRT,
                            float           ifDegrees,
                            tAxis           iAxis );
Trns_tErr Trns_Scale      ( mriTransformRef ThisKRT,
                            float           ifFactor,
                            tAxis           iAxis );

/* converts from a voxel in A space to one in B space, i.e.
   A -> A_RAS, A_RAS -> B_RAS, B_RAS -> B */
Trns_tErr Trns_ConvertAtoB   ( mriTransformRef ThisKRT,
                               xVoxelRef       iAVoxel,
                               xVoxelRef       oBVoxel );
Trns_tErr Trns_ConvertAtoRAS ( mriTransformRef ThisKRT,
                               xVoxelRef       iAVoxel,
                               xVoxelRef       oRASVoxel );
Trns_tErr Trns_ConvertBtoA   ( mriTransformRef ThisKRT,
                               xVoxelRef       iBVoxel,
                               xVoxelRef       oAVoxel );
Trns_tErr Trns_ConvertBtoRAS ( mriTransformRef ThisKRT,
                               xVoxelRef       iBVoxel,
                               xVoxelRef       oRASVoxel );
Trns_tErr Trns_ConvertBRAStoB( mriTransformRef ThisKRT,
                               xVoxelRef       iBRASVoxel,
                               xVoxelRef       oBVoxel );

/* converts matricies between a and b */
Trns_tErr Trns_ConvertMatrixAtoB ( mriTransformRef ThisKRT,
                                   MATRIX*         iAMatrix,
                                   MATRIX*         oBMatrix );
Trns_tErr Trns_ConvertMatrixAtoRAS ( mriTransformRef ThisKRT,
                                     MATRIX*         iAMatrix,
                                     MATRIX*         oRASMatrix );
Trns_tErr Trns_ConvertMatrixBtoA ( mriTransformRef ThisKRT,
                                   MATRIX*         iBMatrix,
                                   MATRIX*         oAMatrix );
Trns_tErr Trns_ConvertMatrixBtoRAS ( mriTransformRef ThisKRT,
                                     MATRIX*         iBMatrix,
                                     MATRIX*         oRASMatrix );

/* internal function that calculates all the matricies */
Trns_tErr Trns_CalcMatricies_ ( mriTransformRef ThisKRT );


/* debugging support */
Trns_tErr Trns_Verify         ( mriTransformRef ThisKRT );
void      Trns_DebugPrint_    ( mriTransformRef ThisKRT );
void      Trns_Signal         ( const char*           inMsg,
                                int             inLineNum,
                                Trns_tErr       ieCode );
const char* Trns_GetErrorString ( Trns_tErr       ieCode );

#endif
