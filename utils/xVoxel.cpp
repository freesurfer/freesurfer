/**
 * @brief general-purpose utils
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

#include "xVoxel.h"
#include <math.h>
#include <stdlib.h>
#include "xDebug.h"

void xVoxl_New(xVoxelRef *oppVoxel)
{
  xVoxelRef ref = (xVoxelRef)calloc(1, sizeof(xVoxel));
  *oppVoxel = ref;
}

void xVoxl_Delete(xVoxelRef *ioppVoxel)
{
  xVoxelRef ref = NULL;

  if (NULL == ioppVoxel) {
    return;
  }

  ref = *ioppVoxel;

  free(ref);
  *ioppVoxel = NULL;
}

void xVoxl_Copy(xVoxelRef ipVoxDest, xVoxelRef ipVoxSrc)
{
  xVoxl_SetFloat(ipVoxDest, xVoxl_GetFloatX(ipVoxSrc), xVoxl_GetFloatY(ipVoxSrc), xVoxl_GetFloatZ(ipVoxSrc));
}

char xVoxl_IsEqualInt(xVoxelRef ipVox1, xVoxelRef ipVox2)
{
  if (xVoxl_GetX(ipVox1) == xVoxl_GetX(ipVox2) && xVoxl_GetY(ipVox1) == xVoxl_GetY(ipVox2) &&
      xVoxl_GetZ(ipVox1) == xVoxl_GetZ(ipVox2)) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

char xVoxl_IsEqualFloat(xVoxelRef ipVox1, xVoxelRef ipVox2)
{
  if (fabs(xVoxl_GetFloatX(ipVox1) - xVoxl_GetFloatX(ipVox2)) < 0.0001 &&
      fabs(xVoxl_GetFloatY(ipVox1) - xVoxl_GetFloatY(ipVox2)) < 0.0001 &&
      fabs(xVoxl_GetFloatZ(ipVox1) - xVoxl_GetFloatZ(ipVox2)) < 0.0001) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

/* declare non-macro versions of ref stuff. */
#ifndef XVOXL_USE_MACROS

void xVoxl_Set(xVoxelRef ref, int x, int y, int z)
{
  ref->mfX = x;
  ref->mfY = y;
  ref->mfZ = z;
}

int xVoxl_GetX(xVoxelRef ref) { return (int)(floor(ref->mfX + 0.5)); }

int xVoxl_GetY(xVoxelRef ref) { return (int)(floor(ref->mfY + 0.5)); }

int xVoxl_GetZ(xVoxelRef ref) { return (int)(floor(ref->mfZ + 0.5)); }

int xVoxl_GetRoundX(xVoxelRef ref) { return rint(ref->mfX); }

int xVoxl_GetRoundY(xVoxelRef ref) { return rint(ref->mfY); }

int xVoxl_GetRoundZ(xVoxelRef ref) { return rint(ref->mfZ); }

void xVoxl_SetFloat(xVoxelRef ref, float x, float y, float z)
{
  ref->mfX = x;
  ref->mfY = y;
  ref->mfZ = z;
}

void xVoxl_SetX(xVoxelRef ref, int x) { ref->mfX = (float)x; }

void xVoxl_SetY(xVoxelRef ref, int y) { ref->mfY = (float)y; }

void xVoxl_SetZ(xVoxelRef ref, int z) { ref->mfZ = (float)z; }

void xVoxl_SetFloatX(xVoxelRef ref, float x) { ref->mfX = x; }

void xVoxl_SetFloatY(xVoxelRef ref, float y) { ref->mfY = y; }

void xVoxl_SetFloatZ(xVoxelRef ref, float z) { ref->mfZ = z; }

float xVoxl_GetFloatX(xVoxelRef ref) { return ref->mfX; }

float xVoxl_GetFloatY(xVoxelRef ref) { return ref->mfY; }

float xVoxl_GetFloatZ(xVoxelRef ref) { return ref->mfZ; }

#endif /* XVOXL_USE_MACROS */

tBoolean xVoxl_IncrementWithMinUntilLimit(xVoxelRef ref, float inMin, float inLimit)
{
  if (ref->mfX < inLimit) {
    ref->mfX += 1.0;
    return TRUE;
  }
  else if (ref->mfY < inLimit) {
    ref->mfX = inMin;
    ref->mfY += 1.0;
    return TRUE;
  }
  else if (ref->mfZ < inLimit) {
    ref->mfX = inMin;
    ref->mfY = inMin;
    ref->mfZ += 1.0;
    return TRUE;
  }
  else {
    return FALSE;
  }
}

tBoolean xVoxl_IncrementUntilLimit(xVoxelRef ref, float inLimit)
{
  return xVoxl_IncrementWithMinUntilLimit(ref, 0, inLimit);
}

tBoolean xVoxl_IncrementUntilLimits(xVoxelRef ref, float inXLimit, float inYLimit, float inZLimit)
{
  if (ref->mfX < inXLimit) {
    ref->mfX += 1.0;
    return TRUE;
  }
  else if (ref->mfY < inYLimit) {
    ref->mfX = 0;
    ref->mfY += 1.0;
    return TRUE;
  }
  else if (ref->mfZ < inZLimit) {
    ref->mfX = 0;
    ref->mfY = 0;
    ref->mfZ += 1.0;
    return TRUE;
  }
  else {
    return FALSE;
  }
}

tBoolean xVoxl_IncrementWithMinsUntilLimits(
    xVoxelRef ref, float inXMin, float inYMin, float inXLimit, float inYLimit, float inZLimit)
{
  if (ref->mfX < inXLimit) {
    ref->mfX += 1.0;
    return TRUE;
  }
  else if (ref->mfY < inYLimit) {
    ref->mfX = inXMin;
    ref->mfY += 1.0;
    return TRUE;
  }
  else if (ref->mfZ < inZLimit) {
    ref->mfX = inXMin;
    ref->mfY = inYMin;
    ref->mfZ += 1.0;
    return TRUE;
  }
  else {
    return FALSE;
  }
}

int xVoxl_ExpandToIndex(xVoxelRef ref, int inDimensionX, int inDimensionY)
{
  return (ref->mfZ * inDimensionX * inDimensionY) + (ref->mfY * inDimensionX) + ref->mfX;
}

void xVoxl_PrintDebug(xVoxelRef ref)
{
  DebugPrint(("Voxel: %d, %d, %d\n", xVoxl_GetX(ref), xVoxl_GetY(ref), xVoxl_GetZ(ref)));
}
