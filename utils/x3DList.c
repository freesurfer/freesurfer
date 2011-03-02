/**
 * @file  x3DList.c
 * @brief general purpose utils
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
#include <math.h>
#include "x3DList.h"
#include "xList.h"
#include "xDebug.h"

char *x3Lst_ksaErrorString [x3Lst_knNumErrorCodes] =
{
  "No error.",
  "Invalid space ptr.",
  "Invalid signature.",
  "Error accessing the list.",
  "Item not found in space.",
  "Invalid space (probably not initialized).",
  "Invalid location (location is out of bounds of the space).",
  "Invalid plane number, out of range.",
  "Couldn't allocate space, memory full?",
  "Couldn't allocate list, memory full?",
  "Invalid error code."
};

/* init the space */
x3Lst_tErr x3Lst_New ( x3DListRef* op3DList,
                       int         inListSize )
{

  x3DListRef   this    = NULL;
  x3Lst_tErr   eResult = x3Lst_tErr_NoErr;
  x3Lst_tPlane plane   = 0;
  int          list    = 0;
  xList_tErr   eList   = xList_tErr_NoErr;

  /* allocate us */
  this = (x3DListRef) malloc( sizeof(x3DList) );
  if ( NULL == this )
  {
    eResult = x3Lst_tErr_CouldntAllocateSpace;
    goto error;
  }

  /* set the sig */
  this->mSignature = x3Lst_kSignature;

  /* set the plane size */
  this->mnPlaneSize = inListSize;

  /* allocate our lists. */
  for ( plane = 0; plane < x3Lst_knNumPlanes; plane++ )
  {
    this->mPlane[plane] = (xListRef*) calloc( sizeof(xList),
                          this->mnPlaneSize );
    if ( NULL == this->mPlane[plane] )
    {
      eResult = x3Lst_tErr_CouldntAllocateSpace;
      goto error;
    }

    /* init all the lists */
    for ( list = 0; list < this->mnPlaneSize; list++ )
    {
      eList = xList_New( &(this->mPlane[plane][list]) );
      if ( xList_tErr_NoErr != eList )
      {
        eResult = x3Lst_tErr_CouldntAllocateList;
        goto error;
      }
    }
  }

  /* return us */
  *op3DList = this;

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_New: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_Delete ( x3DListRef* iop3DList )
{

  x3DListRef   this    = NULL;
  x3Lst_tErr   eResult = x3Lst_tErr_NoErr;
  x3Lst_tPlane plane   = 0;
  int          list    = 0;
  xList_tErr   eList   = xList_tErr_NoErr;

  if ( NULL == iop3DList )
  {
    eResult = x3Lst_tErr_InvalidPtr;
    goto error;
  }

  this = *iop3DList;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* for each plane.. */
  for ( plane = 0; plane < x3Lst_knNumPlanes; plane++ )
  {

    /* delete all the lists */
    for ( list = 0; list < this->mnPlaneSize; list++ )
    {
      eList = xList_Delete( &(this->mPlane[plane][list]) );
      if ( xList_tErr_NoErr != eList )
      {
        eResult = x3Lst_tErr_CouldntAllocateList;
        goto error;
      }
    }

    /* delete the plane */
    free( this->mPlane[plane] );
  }

  /* trash the sig */
  this->mSignature = 0x1;

  /* delete us */
  free( this );

  /* return null */
  *iop3DList = NULL;

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_GetPlaneSize ( x3DListRef this,
                                int*       onPlaneSize )
{

  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  if ( NULL == onPlaneSize )
  {
    eResult = x3Lst_tErr_InvalidPtr;
    goto error;
  }

  /* Return the value */
  *onPlaneSize = this->mnPlaneSize;

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_GetPlaneSize: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_AddItem ( x3DListRef this,
                           xVoxelRef  iWhere,
                           void*      ipItem )
{


  x3Lst_tErr eResult = x3Lst_tErr_NoErr;
  xList_tErr eList   = xList_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  eResult = x3Lst_VerifyLocation( this, iWhere );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* add it to the every plane at the approriate voxel */
  eList = xList_PushItem( this->mPlane[x3Lst_tPlane_X][xVoxl_GetX( iWhere )],
                          ipItem );
  if ( xList_tErr_NoErr != eList )
  {
    eResult = x3Lst_tErr_ErrorAccessingList;
    goto error;
  }

  eList = xList_PushItem( this->mPlane[x3Lst_tPlane_Y][xVoxl_GetY( iWhere )],
                          ipItem );
  if ( xList_tErr_NoErr != eList )
  {
    eResult = x3Lst_tErr_ErrorAccessingList;
    goto error;
  }

  eList = xList_PushItem( this->mPlane[x3Lst_tPlane_Z][xVoxl_GetZ( iWhere )],
                          ipItem );
  if ( xList_tErr_NoErr != eList )
  {
    eResult = x3Lst_tErr_ErrorAccessingList;
    goto error;
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_AddItem: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_RemoveItem ( x3DListRef this,
                              xVoxelRef  iWhere,
                              void**     iopItemToRemove )
{

  x3Lst_tErr eResult       = x3Lst_tErr_NoErr;
  xList_tErr eList         = xList_tErr_NoErr;
  void*      pItemToRemove = NULL;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  eResult = x3Lst_VerifyLocation( this, iWhere );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* remove it from every plane at the approriate voxel */
  pItemToRemove = *iopItemToRemove;
  eList = xList_RemoveItem( this->mPlane[x3Lst_tPlane_X][xVoxl_GetX(iWhere)],
                            &pItemToRemove );
  if ( xList_tErr_NoErr != eList )
  {
    goto error;
  }

  pItemToRemove = *iopItemToRemove;
  eList = xList_RemoveItem( this->mPlane[x3Lst_tPlane_Y][xVoxl_GetY(iWhere)],
                            &pItemToRemove );
  if ( xList_tErr_NoErr != eList )
  {
    goto error;
  }

  pItemToRemove = *iopItemToRemove;
  eList = xList_RemoveItem( this->mPlane[x3Lst_tPlane_Z][xVoxl_GetZ(iWhere)],
                            &pItemToRemove );
  if ( xList_tErr_NoErr != eList )
  {
    goto error;
  }

  /* return the removed item */
  *iopItemToRemove = pItemToRemove;

  goto cleanup;

error:

  if ( xList_tErr_NoErr != eList )
  {
    if ( xList_tErr_ItemNotInList == eList )
    {
      eResult = x3Lst_tErr_ItemNotInSpace;
    }
    else
    {
      eResult = x3Lst_tErr_ErrorAccessingList;
    }
  }

  if ( x3Lst_tErr_NoErr != eResult
       && x3Lst_tErr_ItemNotInSpace != eResult )
  {
    DebugPrint( ("Error %d in x3Lst_RemoveItem: %s\n",
                 eResult, x3Lst_GetErrorString( eResult ) ) );
  }

cleanup:

  return eResult;
}


x3Lst_tErr x3Lst_Clear ( x3DListRef this )
{


  x3Lst_tErr   eResult = x3Lst_tErr_NoErr;
  xList_tErr   eList   = xList_tErr_NoErr;
  x3Lst_tPlane plane   = 0;
  int          list    = 0;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* for each plane.. */
  for ( plane = 0; plane < x3Lst_knNumPlanes; plane++ )
  {

    /* clear all the lists */
    for ( list = 0; list < this->mnPlaneSize; list++ )
    {
      eList = xList_Clear( this->mPlane[plane][list] );
      if ( xList_tErr_NoErr != eList )
      {
        eResult = x3Lst_tErr_ErrorAccessingList;
        goto error;
      }
    }
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_Clear: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_IsInList ( x3DListRef this,
                            void*      ipData,
                            tBoolean*  outIsInList )
{


  x3Lst_tErr   eResult  = x3Lst_tErr_NoErr;
  xList_tErr   eList   = xList_tErr_NoErr;
  x3Lst_tPlane plane    = 0;
  int          list     = 0;
  tBoolean     isInList = FALSE;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* for each plane.. */
  for ( plane = 0; plane < x3Lst_knNumPlanes; plane++ )
  {

    /* ask if it's in this list. */
    for ( list = 0; list < this->mnPlaneSize; list++ )
    {
      eList = xList_IsInList( this->mPlane[plane][list],
                              ipData,
                              &isInList );
      if ( xList_tErr_NoErr != eList )
      {
        eResult = x3Lst_tErr_ErrorAccessingList;
        goto error;
      }

      if ( isInList )
      {
        *outIsInList = isInList;
      }
    }
  }

  *outIsInList = FALSE;

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_IsInList: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_GetItemsInXPlane ( x3DListRef this,
                                    int        inXPlane,
                                    xListRef*  opList )
{


  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  eResult = x3Lst_GetItemsInPlane_( this, x3Lst_tPlane_X,
                                    inXPlane, opList );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_GetItemsInXPlane: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_GetItemsInYPlane ( x3DListRef this,
                                    int        inYPlane,
                                    xListRef*  opList )
{


  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  eResult = x3Lst_GetItemsInPlane_( this, x3Lst_tPlane_Y,
                                    inYPlane, opList );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_GetItemsInYPlane: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_GetItemsInZPlane ( x3DListRef this,
                                    int        inZPlane,
                                    xListRef*  opList )
{


  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  eResult = x3Lst_GetItemsInPlane_( this, x3Lst_tPlane_Z,
                                    inZPlane, opList );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_GetItemsInZPlane: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_GetItemsInPlane_ ( x3DListRef   this,
                                    x3Lst_tPlane iPlane,
                                    int          inPlaneIndex,
                                    xListRef*    opList )
{


  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  /* make sure the index is in bounds */
  if ( inPlaneIndex < 0 ||
       inPlaneIndex >= this->mnPlaneSize )
  {
    eResult = x3Lst_tErr_InvalidPlaneNumber;
    goto error;
  }

  /* return a ptr to the list */
  *opList = this->mPlane[iPlane][inPlaneIndex];

  goto cleanup;

error:

  DebugPrint( ("Error %d in x3Lst_GetItemsInPlane_: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_SetComparator ( x3DListRef this,
                                 xList_tCompare(*iComparator)(void*,void*) )
{

  x3Lst_tErr   eResult = x3Lst_tErr_NoErr;
  x3Lst_tPlane plane   = 0;
  int          list    = 0;
  xList_tErr   eList   = xList_tErr_NoErr;

  eResult = x3Lst_Verify( this );
  if ( x3Lst_tErr_NoErr != eResult )
  {
    goto error;
  }

  /* for each plane.. */
  for ( plane = 0; plane < x3Lst_knNumPlanes; plane++ )
  {

    /* set comparator in all the lists */
    for ( list = 0; list < this->mnPlaneSize; list++ )
    {
      eList = xList_SetComparator( this->mPlane[plane][list], iComparator );
      if ( xList_tErr_NoErr != eList )
      {
        goto error;
      }
    }
  }

  goto cleanup;

error:

  if ( xList_tErr_NoErr != eList )
  {
    eResult = x3Lst_tErr_ErrorAccessingList;
  }

  DebugPrint( ("Error %d in x3Lst_: %s\n",
               eResult, x3Lst_GetErrorString( eResult ) ) );

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_Verify ( x3DListRef this )
{

  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  if ( NULL == this )
  {
    eResult = x3Lst_tErr_InvalidPtr;
    goto cleanup;
  }

  if ( x3Lst_kSignature != this->mSignature )
  {
    eResult = x3Lst_tErr_InvalidSignature;
    goto cleanup;
  }

cleanup:

  return eResult;
}

x3Lst_tErr x3Lst_VerifyLocation ( x3DListRef this,
                                  xVoxelRef  iWhere )
{

  x3Lst_tErr eResult = x3Lst_tErr_NoErr;

  if ( xVoxl_GetX(iWhere) < 0 || xVoxl_GetX(iWhere) >= this->mnPlaneSize ||
       xVoxl_GetY(iWhere) < 0 || xVoxl_GetY(iWhere) >= this->mnPlaneSize ||
       xVoxl_GetZ(iWhere) < 0 || xVoxl_GetZ(iWhere) >= this->mnPlaneSize )
  {
    eResult = x3Lst_tErr_InvalidLocation;
  }

  return eResult;
}

char* x3Lst_GetErrorString ( x3Lst_tErr ieCode )
{

  x3Lst_tErr eCode = ieCode;

  if ( ieCode < 0
       || ieCode >= x3Lst_knNumErrorCodes )
  {
    eCode = x3Lst_tErr_InvalidErrorCode;
  }

  return x3Lst_ksaErrorString[eCode];
}
