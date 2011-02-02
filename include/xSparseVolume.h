/**
 * @file  xSparseVolume.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.4 $
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


#ifndef xSparseVolume_h
#define xSparseVolume_h

#include "xTypes.h"
#include "xVoxel.h"

typedef enum
{
  xSVol_tErr_NoErr = 0,
  xSVol_tErr_InvalidObject,
  xSVol_tErr_InvalidParameter,
  xSVol_tErr_InvalidSignature,
  xSVol_tErr_AllocationFailed,
  xSVol_tErr_ItemNotFound,
  xSVol_tErr_IndexOutOfBounds,
  xSVol_tErr_IndexNotAllocated,
  xSVol_tErr_InvalidErrorCode,
  xSVol_knNumErrorCodes
} xSVol_tErr;

#define xSVol_kSignature 0x098abcdd

typedef struct
{

  tSignature mSignature;

  /* stored mStorage[z][y][x] */
  void**** mStorage;

  /* dimensions of volume */
  int mnXDim;
  int mnYDim;
  int mnZDim;

}
xSparseVolume, *xSparseVolumeRef;

typedef void(*xSVol_tDeleteEntryFuncPtr)( void* );

/* only allocates an array of z values */
xSVol_tErr xSVol_New    ( xSparseVolumeRef* opVolume,
                          int               inXDim,
                          int               inYDim,
                          int               inZDim );
xSVol_tErr xSVol_Delete ( xSparseVolumeRef*         iopVolume,
                          xSVol_tDeleteEntryFuncPtr ipDeleteFunc );

/* checks if the given location exists. if so, returns the value of the
   data stored there (a pointer). if not, returns NULL and
   xSVol_tErr_ItemNotFound. */
xSVol_tErr xSVol_Get ( xSparseVolumeRef this,
                       xVoxelRef        iWhere,
                       void**           oppItem );

/* sets the pointer value of the location. if necessary, allocates additional
   storage in the volume. starts out with only a z array allocated. when it
   gets a set request for 1,2,3, it will allocate a y array at mStorage[3], and
   then an x array at mStorage[3][2]. a subsequent set request for 2,2,3 will
   not allocate any more memory, while one for 2,3,3 will allocate another y
   and x array. */
xSVol_tErr xSVol_Set ( xSparseVolumeRef this,
                       xVoxelRef        iWhere,
                       void*            ipItem );

/* purges all memory. calls the given delete function on each object, if
   given. */
xSVol_tErr xSVol_Purge ( xSparseVolumeRef          this,
                         xSVol_tDeleteEntryFuncPtr ipDeleteFunc );

/* calls given function an all entries. calls the z function every time a new
   z index is entered or left, the y function etc */
xSVol_tErr xSVol_VisitAll ( xSparseVolumeRef this,
                            void(*ipVisitFunction)(void* ipItem),
                            void(*ipXFunction)(int inIndex,tBoolean bEntering),
                            void(*ipYFunction)(int inIndex,tBoolean bEntering),
                            void(*ipZFunction)(int inIndex,tBoolean bEntering));

xSVol_tErr xSVol_VerifyIndex_ ( xSparseVolumeRef this,
                                xVoxelRef        iIndex );

char* xSVol_GetErrorString ( xSVol_tErr ieCode );
xSVol_tErr xSVol_Verify ( xSparseVolumeRef this );


#endif
