/**
 * @file  xVoxel.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/11/11 18:40:47 $
 *    $Revision: 1.12 $
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


#ifndef xxVoxl_H
#define xxVoxl_H

#if defined(__cplusplus)
extern "C" {
#endif

#include "xTypes.h"

/* Enable this to turn on macros, see details below. */
#define XVOXL_USE_MACROS

typedef struct
{

  float mfX, mfY, mfZ;

}
xVoxel, *xVoxelRef;


void xVoxl_New    ( xVoxelRef* oppVoxel );
void xVoxl_Delete ( xVoxelRef* ioppVoxel );

void xVoxl_Copy    ( xVoxelRef ipVoxDest, xVoxelRef ipVoxSrc );

char xVoxl_IsEqualInt   ( xVoxelRef ipVox1, xVoxelRef ipVox2 );
char xVoxl_IsEqualFloat ( xVoxelRef ipVox1, xVoxelRef ipVox2 );

/* These access functions are implemented as functions and macros.
   The functions are slower but can be turned on for debugging purposes,
   while the macros should be used for final builds. */

#ifndef XVOXL_USE_MACROS

void xVoxl_Set      ( xVoxelRef ThisKRT,   int x,   int y,   int z );
void xVoxl_SetFloat ( xVoxelRef ThisKRT, float x, float y, float z );

void xVoxl_SetX ( xVoxelRef ThisKRT, int x );
void xVoxl_SetY ( xVoxelRef ThisKRT, int y );
void xVoxl_SetZ ( xVoxelRef ThisKRT, int z );

void xVoxl_SetFloatX ( xVoxelRef ThisKRT, float x );
void xVoxl_SetFloatY ( xVoxelRef ThisKRT, float y );
void xVoxl_SetFloatZ ( xVoxelRef ThisKRT, float z );

int xVoxl_GetX ( xVoxelRef ThisKRT );
int xVoxl_GetY ( xVoxelRef ThisKRT );
int xVoxl_GetZ ( xVoxelRef ThisKRT );

int xVoxl_GetRoundX ( xVoxelRef ThisKRT );
int xVoxl_GetRoundY ( xVoxelRef ThisKRT );
int xVoxl_GetRoundZ ( xVoxelRef ThisKRT );

float xVoxl_GetFloatX ( xVoxelRef ThisKRT );
float xVoxl_GetFloatY ( xVoxelRef ThisKRT );
float xVoxl_GetFloatZ ( xVoxelRef ThisKRT );

#else /* macros versions */

/* The do/while garbage is to allow the compiler to parse a semicolon
after the macro call. */

#define xVoxl_Set(this,x,y,z) \
  do { \
    (this)->mfX = (x); \
    (this)->mfY = (y); \
    (this)->mfZ = (z); \
  }while(0)

#define xVoxl_GetX(this) (int)(floor( (this)->mfX + 0.5 ))
#define xVoxl_GetY(this) (int)(floor( (this)->mfY + 0.5 ))
#define xVoxl_GetZ(this) (int)(floor( (this)->mfZ + 0.5 ))

#define xVoxl_GetRoundX(this)   rint( (this)->mfX )
#define xVoxl_GetRoundY(this)   rint( (this)->mfY )
#define xVoxl_GetRoundZ(this)   rint( (this)->mfZ )

#define xVoxl_SetFloat(this,x,y,z) \
  do { \
    (this)->mfX = (x); \
    (this)->mfY = (y); \
    (this)->mfZ = (z); \
  }while(0)

#define xVoxl_SetX(this,x)   (this)->mfX = (x)
#define xVoxl_SetY(this,y)   (this)->mfY = (y)
#define xVoxl_SetZ(this,z)   (this)->mfZ = (z)

#define xVoxl_SetFloatX(this,x)   (this)->mfX = (x)
#define xVoxl_SetFloatY(this,y)   (this)->mfY = (y)
#define xVoxl_SetFloatZ(this,z)   (this)->mfZ = (z)

#define xVoxl_GetFloatX(this)   (this)->mfX
#define xVoxl_GetFloatY(this)   (this)->mfY
#define xVoxl_GetFloatZ(this)   (this)->mfZ

#endif /* end of macro versions */

tBoolean xVoxl_IncrementUntilLimit          ( xVoxelRef ThisKRT,
    float      inLimit );
tBoolean xVoxl_IncrementWithMinUntilLimit   ( xVoxelRef ThisKRT,
    float     inMin,
    float     inLimit );
tBoolean xVoxl_IncrementUntilLimits         ( xVoxelRef ThisKRT,
    float     inXLimit,
    float     inYLimit,
    float     inZLimit );
tBoolean xVoxl_IncrementWithMinsUntilLimits  ( xVoxelRef ThisKRT,
    float     inXMin,
    float     inYMin,
    float     inXLimit,
    float     inYLimit,
    float     inZLimit );

int xVoxl_ExpandToIndex ( xVoxelRef ThisKRT, int inDimensionX, int inDimensionY );

void xVoxl_PrintDebug ( xVoxelRef  ThisKRT );

/* easy way to expand the voxel struct into x/y/z coords */
#define xVoxl_ExpandFloat(v) xVoxl_GetFloatX(v), xVoxl_GetFloatY(v), xVoxl_GetFloatZ(v)
#define xVoxl_ExpandInt(v)   xVoxl_GetX(v),      xVoxl_GetY(v),      xVoxl_GetZ(v)
#define xVoxl_ExpandRint(v)   (int)rint(xVoxl_GetX(v)), (int)rint(xVoxl_GetY(v)), (int)rint(xVoxl_GetZ(v))

#if defined(__cplusplus)
};
#endif


#endif







