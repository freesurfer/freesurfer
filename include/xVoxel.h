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


#ifndef xxVoxl_H
#define xxVoxl_H

#include "xTypes.h"

/* Enable item to turn on macros, see details below. */
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

#define xVoxl_Set(item,x,y,z) \
  do { \
    (item)->mfX = (x); \
    (item)->mfY = (y); \
    (item)->mfZ = (z); \
  }while(0)

#define xVoxl_GetX(item) (int)(floor( (item)->mfX + 0.5 ))
#define xVoxl_GetY(item) (int)(floor( (item)->mfY + 0.5 ))
#define xVoxl_GetZ(item) (int)(floor( (item)->mfZ + 0.5 ))

#define xVoxl_GetRoundX(item)   rint( (item)->mfX )
#define xVoxl_GetRoundY(item)   rint( (item)->mfY )
#define xVoxl_GetRoundZ(item)   rint( (item)->mfZ )

#define xVoxl_SetFloat(item,x,y,z) \
  do { \
    (item)->mfX = (x); \
    (item)->mfY = (y); \
    (item)->mfZ = (z); \
  }while(0)

#define xVoxl_SetX(item,x)   (item)->mfX = (x)
#define xVoxl_SetY(item,y)   (item)->mfY = (y)
#define xVoxl_SetZ(item,z)   (item)->mfZ = (z)

#define xVoxl_SetFloatX(item,x)   (item)->mfX = (x)
#define xVoxl_SetFloatY(item,y)   (item)->mfY = (y)
#define xVoxl_SetFloatZ(item,z)   (item)->mfZ = (z)

#define xVoxl_GetFloatX(item)   (item)->mfX
#define xVoxl_GetFloatY(item)   (item)->mfY
#define xVoxl_GetFloatZ(item)   (item)->mfZ

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

#endif
