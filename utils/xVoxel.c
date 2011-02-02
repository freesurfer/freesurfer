/**
 * @file  xVoxel.c
 * @brief general-purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:20 $
 *    $Revision: 1.14 $
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


#include <stdlib.h>
#include <math.h>
#include "xDebug.h"
#include "xVoxel.h"


void xVoxl_New ( xVoxelRef* oppVoxel )
{

  xVoxelRef this = (xVoxelRef) calloc (1,sizeof(xVoxel));
  *oppVoxel = this;
}

void xVoxl_Delete ( xVoxelRef* ioppVoxel )
{

  xVoxelRef this = NULL;

  if ( NULL == ioppVoxel )
  {
    return;
  }

  this = *ioppVoxel;

  free ( this );
  *ioppVoxel = NULL;
}

void xVoxl_Copy ( xVoxelRef ipVoxDest, xVoxelRef ipVoxSrc )
{

  xVoxl_SetFloat ( ipVoxDest,
                   xVoxl_GetFloatX ( ipVoxSrc ),
                   xVoxl_GetFloatY ( ipVoxSrc ),
                   xVoxl_GetFloatZ ( ipVoxSrc ) );
}

char xVoxl_IsEqualInt ( xVoxelRef ipVox1, xVoxelRef ipVox2 )
{

  if ( xVoxl_GetX(ipVox1) == xVoxl_GetX(ipVox2)
       && xVoxl_GetY(ipVox1) == xVoxl_GetY(ipVox2)
       && xVoxl_GetZ(ipVox1) == xVoxl_GetZ(ipVox2) )
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

char xVoxl_IsEqualFloat ( xVoxelRef ipVox1, xVoxelRef ipVox2 )
{

  if ( fabs( xVoxl_GetFloatX(ipVox1) - xVoxl_GetFloatX(ipVox2) ) < 0.0001 &&
       fabs( xVoxl_GetFloatY(ipVox1) - xVoxl_GetFloatY(ipVox2) ) < 0.0001 &&
       fabs( xVoxl_GetFloatZ(ipVox1) - xVoxl_GetFloatZ(ipVox2) ) < 0.0001 )
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* declare non-macro versions of this stuff. */
#ifndef XVOXL_USE_MACROS

void xVoxl_Set ( xVoxelRef this, int x, int y, int z )
{

  this->mfX = x;
  this->mfY = y;
  this->mfZ = z;
}


int xVoxl_GetX ( xVoxelRef this )
{

  return (int)(floor( this->mfX + 0.5 ));
}


int xVoxl_GetY ( xVoxelRef this )
{

  return (int)(floor( this->mfY + 0.5 ));
}


int xVoxl_GetZ ( xVoxelRef this )
{

  return (int)(floor( this->mfZ + 0.5 ));
}


int xVoxl_GetRoundX ( xVoxelRef this )
{

  return rint( this->mfX );
}


int xVoxl_GetRoundY ( xVoxelRef this )
{

  return rint( this->mfY );
}


int xVoxl_GetRoundZ ( xVoxelRef this )
{

  return rint( this->mfZ );
}


void xVoxl_SetFloat ( xVoxelRef this, float x, float y, float z )
{

  this->mfX = x;
  this->mfY = y;
  this->mfZ = z;
}

void xVoxl_SetX ( xVoxelRef this, int x )
{

  this->mfX = (float)x;
}

void xVoxl_SetY ( xVoxelRef this, int y )
{

  this->mfY = (float)y;
}

void xVoxl_SetZ ( xVoxelRef this, int z )
{

  this->mfZ = (float)z;
}

void xVoxl_SetFloatX ( xVoxelRef this, float x )
{

  this->mfX = x;
}

void xVoxl_SetFloatY ( xVoxelRef this, float y )
{

  this->mfY = y;
}

void xVoxl_SetFloatZ ( xVoxelRef this, float z )
{

  this->mfZ = z;
}


float xVoxl_GetFloatX ( xVoxelRef this )
{

  return this->mfX;
}


float xVoxl_GetFloatY ( xVoxelRef this )
{

  return this->mfY;
}


float xVoxl_GetFloatZ ( xVoxelRef this )
{

  return this->mfZ;
}

#endif /* XVOXL_USE_MACROS */

tBoolean xVoxl_IncrementWithMinUntilLimit ( xVoxelRef this,
    float inMin, float inLimit )
{

  if ( this->mfX < inLimit )
  {
    this->mfX += 1.0;
    return TRUE;
  }
  else if ( this->mfY < inLimit )
  {
    this->mfX = inMin;
    this->mfY += 1.0;
    return TRUE;
  }
  else if ( this->mfZ < inLimit )
  {
    this->mfX = inMin;
    this->mfY = inMin;
    this->mfZ += 1.0;
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

tBoolean xVoxl_IncrementUntilLimit ( xVoxelRef this, float inLimit )
{

  return xVoxl_IncrementWithMinUntilLimit( this, 0, inLimit );
}

tBoolean xVoxl_IncrementUntilLimits ( xVoxelRef this, float inXLimit,
                                      float inYLimit, float inZLimit )
{

  if ( this->mfX < inXLimit )
  {
    this->mfX += 1.0;
    return TRUE;
  }
  else if ( this->mfY < inYLimit )
  {
    this->mfX = 0;
    this->mfY += 1.0;
    return TRUE;
  }
  else if ( this->mfZ < inZLimit )
  {
    this->mfX = 0;
    this->mfY = 0;
    this->mfZ += 1.0;
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

tBoolean xVoxl_IncrementWithMinsUntilLimits ( xVoxelRef this,
    float     inXMin,
    float     inYMin,
    float     inXLimit,
    float     inYLimit,
    float     inZLimit )
{
  if ( this->mfX < inXLimit )
  {
    this->mfX += 1.0;
    return TRUE;
  }
  else if ( this->mfY < inYLimit )
  {
    this->mfX = inXMin;
    this->mfY += 1.0;
    return TRUE;
  }
  else if ( this->mfZ < inZLimit )
  {
    this->mfX = inXMin;
    this->mfY = inYMin;
    this->mfZ += 1.0;
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

int xVoxl_ExpandToIndex ( xVoxelRef this, int inDimensionX, int inDimensionY )
{

  return ( this->mfZ * inDimensionX * inDimensionY ) +
         ( this->mfY * inDimensionX  ) + this->mfX;
}

void
xVoxl_PrintDebug ( xVoxelRef this )
{

  DebugPrint( ("Voxel: %d, %d, %d\n",
               xVoxl_GetX(this), xVoxl_GetY(this), xVoxl_GetZ(this) ) );
}
