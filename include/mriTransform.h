#ifndef mriTransform_h
#define mriTransform_h

#include "matrix.h"
#include "xDebug.h"
#include "xVoxel.h"

/*
  future:

    Trns_TranslateA( transform, transVector );
    Trns_RotateB( transform, rotationVector, degrees );
    Trns_ScaleA( transform, scaleVector );
*/

typedef enum {

  Trns_tErr_NoErr = 0,
  Trns_tErr_InvalidPtr,
  Trns_tErr_InvalidSignature,
  Trns_tErr_AllocationFailed,
  Trns_tErr_TransformationMatrixNotInited,
  Trns_tErr_InvalidErrorCode,
  Trns_knNumErrorCodes

} Trns_tErr;

typedef struct {

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
  
} mriTransform, *mriTransformRef;

#define Trns_kSignature 0x09876543

Trns_tErr Trns_New    ( mriTransformRef* opTransform );
Trns_tErr Trns_Delete ( mriTransformRef* iopTransform );

/* these are the minimum that should be set in a normal situation */
Trns_tErr Trns_CopyAtoRAS     ( mriTransformRef this,
        MATRIX*         iAtoRAS );
Trns_tErr Trns_CopyBtoRAS     ( mriTransformRef this,
        MATRIX*         iBtoRAS );
Trns_tErr Trns_CopyARAStoBRAS ( mriTransformRef this,
        MATRIX*         iARAStoBRAS );

/* converts from a voxel in A space to one in B space, i.e. 
   A -> A_RAS, A_RAS -> B_RAS, B_RAS -> B */
Trns_tErr Trns_ConvertAtoB ( mriTransformRef this,
           xVoxelRef       iAVoxel,
           xVoxelRef       oBVoxel );
Trns_tErr Trns_ConvertBtoA ( mriTransformRef this,
           xVoxelRef       iBVoxel,
           xVoxelRef       oAVoxel );

/* converts between a and ras */
Trns_tErr Trns_ConvertAtoRAS ( mriTransformRef this,
             xVoxelRef       iAVoxel,
             xVoxelRef       oARASVoxel );
Trns_tErr Trns_ConvertRAStoA ( mriTransformRef this,
             xVoxelRef       iARASVoxel,
             xVoxelRef       oAVoxel );

/* internal function that calculates all the matricies */
Trns_tErr Trns_CalcMatricies_ ( mriTransformRef this );


/* debugging support */
Trns_tErr Trns_Verify ( mriTransformRef this );
void Trns_DebugPrint_ ( mriTransformRef this );
void Trns_Signal ( char* inMsg, int inLineNum, Trns_tErr ieCode );
char* Trns_GetErrorString ( Trns_tErr ieCode );

#endif
