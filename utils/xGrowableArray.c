/**
 * @file  xGrowableArray.c
 * @brief arrays that grow
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.10 $
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


#include <stdlib.h>
#include <string.h>
#include "xGrowableArray.h"
#include "xDebug.h"

char *xGArr_ksaErrorStrings [xGArr_knNumErrorCodes] =
{

  "No error.",
  "Invalid object.",
  "Invalid signature.",
  "Allocation failed.",
  "Last item.",
  "Invalid error code."
};


xGArr_tErr xGArr_New ( xGrowableArrayRef* opList,
                       int                inSize,
                       int                inInitialNumItems )
{

  xGArr_tErr        eResult = xGArr_tErr_NoErr;
  xGrowableArrayRef this    = NULL;

  this = (xGrowableArrayRef) malloc( sizeof( xGrowableArray ) );
  ;
  if ( NULL == this )
  {
    eResult = xGArr_tErr_AllocationFailed;
    goto error;
  }

  /* init sig */
  this->mSignature = xGArr_kSignature;

  /* initial values */
  this->mnNumItems      = 0;
  this->mnMaxNumItems   = inInitialNumItems;
  this->mnItemSizeBytes = inSize;
  this->mnMaxSizeBytes  = inInitialNumItems * inSize;
  this->mnNext          = 0;

  /* allocate initial storage */
  this->mpData = malloc( this->mnMaxSizeBytes );
  if ( NULL == this->mpData )
  {
    eResult = xGArr_tErr_AllocationFailed;
    goto error;
  }

  /* return us */
  *opList = this;

  goto cleanup;

error:

  if ( NULL != this )
  {
    free( this );
  }

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_New: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;
}

xGArr_tErr xGArr_Delete ( xGrowableArrayRef* iopList )
{

  xGArr_tErr        eResult = xGArr_tErr_NoErr;
  xGrowableArrayRef this    = NULL;

  if ( iopList == NULL )
  {
    eResult = xGArr_tErr_InvalidObject;
    goto error;
  }

  this = *iopList;

  eResult = xGArr_Verify( this );
  if ( xGArr_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* trash sig */
  this->mSignature = 0x1;

  /* free us */
  free( this->mpData );
  free( this );

  *iopList = NULL;

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_Delete: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;
}


xGArr_tErr xGArr_Add ( xGrowableArrayRef this,
                       void*      ipSrc )
{

  xGArr_tErr eResult     = xGArr_tErr_NoErr;
  void*      pNewStorage = NULL;

  eResult = xGArr_Verify( this );
  if ( xGArr_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* if our num items is not less than our max items.... */
  if ( !(this->mnNumItems < this->mnMaxNumItems) )
  {

    /* allocate twice the storage */
    pNewStorage = malloc( this->mnMaxSizeBytes * 2 );
    if ( NULL == pNewStorage )
    {
      eResult = xGArr_tErr_AllocationFailed;
      goto error;
    }

    /* copy our current data in */
    memmove( pNewStorage, this->mpData, this->mnMaxSizeBytes );

    /* delete old storage and point to new one */
    free( this->mpData );
    this->mpData = pNewStorage;
    pNewStorage = NULL;

    /* save our new storage size */
    this->mnMaxSizeBytes *= 2;
    this->mnMaxNumItems = this->mnMaxSizeBytes / this->mnItemSizeBytes;
  }

  /* copy data into this location */
  memmove( &((char*)this->mpData)[ this->mnNumItems * this->mnItemSizeBytes ],
           ipSrc, this->mnItemSizeBytes );

  /* inc num items */
  this->mnNumItems ++;

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_Add: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;
}


xGArr_tErr xGArr_ResetIterator ( xGrowableArrayRef this )
{

  xGArr_tErr eResult = xGArr_tErr_NoErr;

  eResult = xGArr_Verify( this );
  if ( xGArr_tErr_NoErr != eResult )
  {
    goto error;
  }

  this->mnNext = 0;

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_ResetIterator: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;
}

xGArr_tErr xGArr_NextItem ( xGrowableArrayRef this,
                            void*             opDest )
{

  xGArr_tErr eResult = xGArr_tErr_NoErr;

  eResult = xGArr_Verify( this );
  if ( xGArr_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* see if we're at the end */
  if ( this->mnNext >= (this->mnNumItems) )
  {
    eResult = xGArr_tErr_LastItem;
    goto cleanup;
  }

  /* return data at this location */
  memmove( opDest,
           &((char*)this->mpData)[ this->mnNext * this->mnItemSizeBytes ],
           this->mnItemSizeBytes );

  /* inc interator */
  this->mnNext++;

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_NextItem: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;
}

xGArr_tErr xGArr_Clear  ( xGrowableArrayRef this )
{

  xGArr_tErr eResult = xGArr_tErr_NoErr;

  eResult = xGArr_Verify( this );
  if ( xGArr_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* free our storage and realloc it */
  free( this->mpData );
  this->mpData = malloc( this->mnMaxSizeBytes );
  if ( NULL == this->mpData )
  {
    eResult = xGArr_tErr_AllocationFailed;
    goto error;
  }

  /* reset counter */
  this->mnNumItems = 0;
  this->mnNext     = 0;

  goto cleanup;

error:

  if ( xGArr_tErr_NoErr != eResult )
  {
    DebugPrint( ( "Error %d in xGArr_Clear: %s\n",
                  eResult, xGArr_GetErrorString( eResult ) ));
  }

cleanup:

  return eResult;


}


xGArr_tErr xGArr_Verify ( xGrowableArrayRef this )
{

  xGArr_tErr eResult = xGArr_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this )
  {
    eResult = xGArr_tErr_InvalidObject;
    goto cleanup;
  }

  /* check signature */
  if ( xGArr_kSignature != this->mSignature )
  {
    eResult = xGArr_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;

}

char* xGArr_GetErrorString ( xGArr_tErr ieCode )
{

  xGArr_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= xGArr_knNumErrorCodes )
  {
    eCode = xGArr_tErr_InvalidErrorCode;
  }

  return xGArr_ksaErrorStrings [eCode];
}
