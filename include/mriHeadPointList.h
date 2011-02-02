/**
 * @file  mriHeadPointList.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.6 $
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


#ifndef mriHeadPointList_h
#define mriHeadPointList_h

#include "mriTypes.h"
#include "xVoxel.h"
#include "mriTransform.h"

typedef enum
{

  HPtL_tErr_NoErr = 0,
  HPtL_tErr_InvalidObject,
  HPtL_tErr_InvalidParameter,
  HPtL_tErr_InvalidSignature,
  HPtL_tErr_AllocationFailed,
  HPtL_tErr_ErrorOpeningHeadPointFile,
  HPtL_tErr_ErrorParsingHeadPointFile,
  HPtL_tErr_ErrorOpeningTransformFile,
  HPtL_tErr_ErrorParsingTransformFile,
  HPtL_tErr_ErrorCreatingTransform,
  HPtL_tErr_ErrorAccessingClientTransform,
  HPtL_tErr_ErrorAccessingTransform,
  HPtL_tErr_LastPoint,
  HPtL_tErr_InvalidErrorCode,
  HPtL_knNumErrorCodes

} HPtL_tErr;

typedef enum
{

  HPtL_tIterationPlane_X = 0,
  HPtL_tIterationPlane_Y,
  HPtL_tIterationPlane_Z,
  HPtL_tIterationPlane_All,
  HPtL_knNumIterationPlanes
} HPtL_tIterationPlane;

typedef struct
{

  char   msLabel[256];
  int    mnIndex;
  xVoxel mPoint;           /* original point */
  xVoxel mClientPoint;     /* client space point */

}
HPtL_tHeadPoint, *HPtL_tHeadPointRef;

#define HPtL_kSignature 0x009087fb

typedef struct
{

  long mSignature;

  /* source */
  char msPointFile[256];
  char msTransformFile[256];

  /* points storage */
  HPtL_tHeadPoint* maPoints;
  int              mnNumPoints;

  /* transform. a is client, b is head point space. */
  mriTransformRef mTransform;

  /* iterator state */
  HPtL_tIterationPlane mIterPlane;  /* in client space */
  float                mfIterPlaneNumber;
  float                mfIterPlaneRange;
  int                  mnCurPoint;

}
mriHeadPointList, *mriHeadPointListRef;

HPtL_tErr HPtL_New    ( mriHeadPointListRef* oList,
                        char*                isListName,
                        char*                isTransformName,
                        MATRIX*              iClientTransform );
HPtL_tErr HPtL_Delete ( mriHeadPointListRef* opList );

HPtL_tErr HPtL_ReadHeadListFile_ ( mriHeadPointListRef this,
                                   char*               isListName );
HPtL_tErr HPtL_CreateTransform_  ( mriHeadPointListRef this,
                                   char*               isTransformName,
                                   MATRIX*             iClientTransform );

HPtL_tErr HPtL_ConvertListToClientSpace_ ( mriHeadPointListRef this );

HPtL_tErr HPtL_WriteTransform     ( mriHeadPointListRef this,
                                    char*               isDest );
HPtL_tErr HPtL_WriteHeadPointFile ( mriHeadPointListRef this,
                                    char*               isDest );

HPtL_tErr HPtL_ResetIterator ( mriHeadPointListRef  this,
                               HPtL_tIterationPlane iPlane,
                               float                ifPlaneNumber,
                               float                ifPlaneRange );
HPtL_tErr HPtL_NextPoint     ( mriHeadPointListRef this,
                               HPtL_tHeadPointRef* opPoint );

HPtL_tErr HPtL_FindNearestPoint ( mriHeadPointListRef  this,
                                  HPtL_tIterationPlane iPlane,
                                  float                ifPlaneRange,
                                  xVoxelRef            iWhere,
                                  HPtL_tHeadPointRef*  opPoint );
HPtL_tErr HPtL_FindFlattenedNearestPoint ( mriHeadPointListRef  this,
    HPtL_tIterationPlane iPlane,
    xVoxelRef            iWhere,
    HPtL_tHeadPointRef*  opPoint );

HPtL_tErr HPtL_RestoreTransform ( mriHeadPointListRef this );
HPtL_tErr HPtL_ApplyTransform   ( mriHeadPointListRef this,
                                  MATRIX*             iTransform );
HPtL_tErr HPtL_Translate ( mriHeadPointListRef this,
                           float               ifDistance,
                           tAxis               iAxis );
HPtL_tErr HPtL_Rotate    ( mriHeadPointListRef this,
                           float               ifDegrees,
                           tAxis               iAxis );
HPtL_tErr HPtL_Scale      ( mriHeadPointListRef this,
                            float               ifFactor,
                            tAxis               iAxis );

HPtL_tErr HPtL_AlignPointToClientVoxel ( mriHeadPointListRef this,
    HPtL_tHeadPointRef  iPoint,
    xVoxelRef           iClientVox );

HPtL_tErr HPtL_Verify         ( mriHeadPointListRef this );
char*     HPtL_GetErrorString ( HPtL_tErr           ieCode );

#endif
