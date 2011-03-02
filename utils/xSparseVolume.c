/**
 * @file  xSparseVolume.c
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.9 $
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
#include <math.h>
#include <string.h>
#include "xSparseVolume.h"
#include "xDebug.h"

char *xSVol_ksaErrorStrings [xSVol_knNumErrorCodes] =
{

  "No error.",
  "Invalid pointer to object.",
  "Invalid parameter.",
  "Invalid signature.",
  "Memory allocation failed.",
  "Item not found.",
  "Index out of bounds.",
  "Index not allocated.",
  "Invalid error code."
};

xSVol_tErr xSVol_New( xSparseVolumeRef* opVolume,
                      int               inXDim,
                      int               inYDim,
                      int               inZDim )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;
  xSparseVolumeRef this = NULL;

  if ( NULL == opVolume ||
       inXDim < 0 ||
       inYDim < 0 ||
       inZDim < 0 )
  {
    eResult = xSVol_tErr_InvalidParameter;
    goto error;
  }

  /* allocate us */
  this = (xSparseVolumeRef) malloc( sizeof( xSparseVolume ));
  if ( NULL == this )
  {
    eResult = xSVol_tErr_AllocationFailed;
    goto error;
  }

  /* set signature and defaults */
  this->mSignature = xSVol_kSignature;
  this->mnXDim = inXDim;
  this->mnYDim = inYDim;
  this->mnZDim = inZDim;

  /* allocate the zs */
  this->mStorage = (void****) malloc( sizeof(void***) * this->mnZDim );
  if ( NULL == this->mStorage )
  {
    eResult = xSVol_tErr_AllocationFailed;
    goto error;
  }
  memset( this->mStorage, 0, sizeof(void***) * this->mnZDim );

  /* return us */
  *opVolume = this;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_New: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_Delete ( xSparseVolumeRef*         iopVolume,
                          xSVol_tDeleteEntryFuncPtr ipDeleteFunc )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;
  xSparseVolumeRef this = NULL;

  if ( NULL == iopVolume )
  {
    eResult = xSVol_tErr_InvalidParameter;
    goto error;
  }

  this = *iopVolume;
  eResult = xSVol_Verify( this );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* first purge the volume */
  eResult = xSVol_Purge( this, ipDeleteFunc );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* trash signature */
  this->mSignature = 0x1;

  /* delete us */
  free( this );

  /* return null */
  *iopVolume = NULL;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_Delete: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_Get ( xSparseVolumeRef this,
                       xVoxelRef        iWhere,
                       void**           oppItem )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;
  void* pItem = NULL;

  eResult = xSVol_Verify( this );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  if ( NULL == iWhere ||
       NULL == oppItem )
  {
    eResult = xSVol_tErr_InvalidParameter;
    goto error;
  }

  /* look at the index */
  eResult = xSVol_VerifyIndex_( this, iWhere );

  /* if it's not allocated, item is null. */
  if ( xSVol_tErr_IndexNotAllocated == eResult )
  {
    pItem = NULL;
    eResult = xSVol_tErr_NoErr;

  }
  else if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;

  }
  else
  {

    /* find the item. */
    pItem = this->mStorage[ xVoxl_GetZ( iWhere ) ]
            [ xVoxl_GetY( iWhere ) ]
            [ xVoxl_GetX( iWhere ) ];
  }

  /* return the item */
  *oppItem = pItem;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_Get: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_Set ( xSparseVolumeRef this,
                       xVoxelRef        iWhere,
                       void*            ipItem )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;

  eResult = xSVol_Verify( this );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* look at the index */
  eResult = xSVol_VerifyIndex_( this, iWhere );

  /* if it's not allocated, allocate it. */
  if ( xSVol_tErr_IndexNotAllocated == eResult )
  {

    /* allocate a y array */
    if ( NULL == this->mStorage[xVoxl_GetZ(iWhere)] )
    {
      this->mStorage[xVoxl_GetZ(iWhere)] =
        (void***) malloc( sizeof(void**) * this->mnYDim );
      if ( NULL == this->mStorage )
      {
        eResult = xSVol_tErr_AllocationFailed;
        goto error;
      }
      memset( this->mStorage[xVoxl_GetZ(iWhere)], 0,
              sizeof(void*) * this->mnYDim );
    }

    /* allocate an x array */
    if ( NULL == this->mStorage[xVoxl_GetZ(iWhere)][xVoxl_GetY(iWhere)] )
    {
      this->mStorage[xVoxl_GetZ(iWhere)][xVoxl_GetY(iWhere)] =
        (void**) malloc( sizeof(void*) * this->mnXDim );
      if ( NULL == this->mStorage )
      {
        eResult = xSVol_tErr_AllocationFailed;
        goto error;
      }
      memset( this->mStorage[xVoxl_GetZ(iWhere)][xVoxl_GetY(iWhere)], 0,
              sizeof(void**) * this->mnXDim );
    }

    /* clear the error */
    eResult = xSVol_tErr_NoErr;

  }
  else if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* set the item. */
  this->mStorage[ xVoxl_GetZ( iWhere ) ]
  [ xVoxl_GetY( iWhere ) ]
  [ xVoxl_GetX( iWhere ) ] = ipItem;

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_Set: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_Purge ( xSparseVolumeRef          this,
                         xSVol_tDeleteEntryFuncPtr ipDeleteFunc )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;
  int nZ = 0;
  int nY = 0;

  eResult = xSVol_Verify( this );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* visit every node with the delete function */
  eResult = xSVol_VisitAll( this, ipDeleteFunc, NULL, NULL, NULL );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* for every y and z, delete the arrays */
  for ( nZ = 0; nZ < this->mnZDim; nZ++ )
  {
    if ( NULL != this->mStorage[nZ] )
    {
      for ( nY = 0; nY < this->mnYDim; nY++ )
      {
        if ( NULL != this->mStorage[nZ][nY] )
        {
          free( this->mStorage[nZ][nY] );
        }
      }
      free( this->mStorage[nZ] );
    }
  }

  /* zero the z array */
  memset( this->mStorage, 0, sizeof(void***) * this->mnZDim );

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_Purge: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_VisitAll ( xSparseVolumeRef this,
                            void(*ipVisitFunction)(void* ipItem),
                            void(*ipXFunction)(int inIndex,tBoolean bEntering),
                            void(*ipYFunction)(int inIndex,tBoolean bEntering),
                            void(*ipZFunction)(int inIndex,tBoolean bEntering)
                          )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;
  int nX = 0;
  int nY = 0;
  int nZ = 0;

  eResult = xSVol_Verify( this );
  if ( xSVol_tErr_NoErr != eResult )
  {
    goto error;
  }

  for ( nZ = 0; nZ < this->mnZDim; nZ++ )
  {

    if ( NULL != this->mStorage[nZ] )
    {
      if ( NULL != ipZFunction )
      {
        ipZFunction( nZ, TRUE );
      }
      for ( nY = 0; nY < this->mnYDim; nY++ )
      {

        if ( NULL != this->mStorage[nZ][nY] )
        {
          if ( NULL != ipYFunction )
          {
            ipYFunction( nY, TRUE );
          }
          for ( nX = 0; nX < this->mnXDim; nX++ )
          {

            if ( NULL != this->mStorage[nZ][nY][nX] )
            {
              if ( NULL != ipXFunction )
              {
                ipXFunction( nX, TRUE );
              }
              if ( NULL != ipVisitFunction )
              {
                ipVisitFunction( this->mStorage[nZ][nY][nX] );
              }
              if ( NULL != ipXFunction )
              {
                ipXFunction( nX, FALSE );
              }
            }
          }
          if ( NULL != ipYFunction )
          {
            ipYFunction( nY, FALSE );
          }
        }
      }
      if ( NULL != ipZFunction )
      {
        ipZFunction( nZ, FALSE );
      }
    }
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in xSVol_: %s\n",
               eResult, xSVol_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

xSVol_tErr xSVol_VerifyIndex_ ( xSparseVolumeRef this,
                                xVoxelRef        iIndex )
{

  /* check for bounds */
  if ( xVoxl_GetX(iIndex) < 0 ||
       xVoxl_GetX(iIndex) >= this->mnXDim ||
       xVoxl_GetY(iIndex) < 0 ||
       xVoxl_GetY(iIndex) >= this->mnYDim ||
       xVoxl_GetZ(iIndex) < 0 ||
       xVoxl_GetZ(iIndex) >= this->mnZDim )
  {
    return xSVol_tErr_IndexOutOfBounds;
  }

  /* check if it exists */
  if ( NULL != this->mStorage[xVoxl_GetZ(iIndex)] &&
       NULL != this->mStorage[xVoxl_GetZ(iIndex)][xVoxl_GetY(iIndex)] )
  {
    return xSVol_tErr_NoErr;
  }
  else
  {
    return xSVol_tErr_IndexNotAllocated;
  }
}

char* xSVol_GetErrorString ( xSVol_tErr ieCode )
{

  xSVol_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= xSVol_knNumErrorCodes )
  {
    eCode = xSVol_tErr_InvalidErrorCode;
  }

  return xSVol_ksaErrorStrings [eCode];
}

xSVol_tErr xSVol_Verify ( xSparseVolumeRef this )
{

  xSVol_tErr eResult = xSVol_tErr_NoErr;

  /* check for null ptr */
  if ( NULL == this )
  {
    eResult = xSVol_tErr_InvalidObject;
    goto cleanup;
  }

  /* check signature */
  if ( xSVol_kSignature != this->mSignature )
  {
    eResult = xSVol_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;

}
