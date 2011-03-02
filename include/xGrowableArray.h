/**
 * @file  xGrowableArray.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef xGrowableArray_h
#define xGrowableArray_h

typedef enum
{

  xGArr_tErr_NoErr,
  xGArr_tErr_InvalidObject,
  xGArr_tErr_InvalidSignature,
  xGArr_tErr_AllocationFailed,
  xGArr_tErr_LastItem,
  xGArr_tErr_InvalidErrorCode,
  xGArr_knNumErrorCodes
} xGArr_tErr;

#define xGArr_kSignature 0x8765433

typedef struct
{

  long mSignature;

  /* storage */
  int mnNumItems;
  int mnMaxNumItems;
  int mnItemSizeBytes;
  int mnMaxSizeBytes;
  void* mpData;

  /* iterator */
  int mnNext;

}
xGrowableArray, *xGrowableArrayRef;

xGArr_tErr xGArr_New    ( xGrowableArrayRef* opList,
                          int                inSize,
                          int                inNumItems );
xGArr_tErr xGArr_Delete ( xGrowableArrayRef* iopList );

xGArr_tErr xGArr_Add    ( xGrowableArrayRef this,
                          void*             ipSrc );

xGArr_tErr xGArr_ResetIterator ( xGrowableArrayRef this );
xGArr_tErr xGArr_NextItem      ( xGrowableArrayRef this,
                                 void*             opDest );

xGArr_tErr xGArr_Clear  ( xGrowableArrayRef this );

xGArr_tErr xGArr_Verify ( xGrowableArrayRef this );
char* xGArr_GetErrorString ( xGArr_tErr ieCode );

#endif
