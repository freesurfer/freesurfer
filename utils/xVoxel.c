#include <stdlib.h>
#include <math.h>
#include "xDebug.h"
#include "xVoxel.h"


void xVoxl_New ( xVoxelRef* oppVoxel ) {

  xVoxelRef this = (xVoxelRef) calloc (1,sizeof(xVoxel));
  *oppVoxel = this;
}

void xVoxl_Delete ( xVoxelRef* ioppVoxel ) {

  xVoxelRef this = NULL;

  if ( NULL == ioppVoxel ) 
    return;

  this = *ioppVoxel;

  free ( this );
  *ioppVoxel = NULL;
}

void xVoxl_Copy ( xVoxelRef ipVoxDest, xVoxelRef ipVoxSrc ) {

  xVoxl_SetFloat ( ipVoxDest, 
       xVoxl_GetFloatX ( ipVoxSrc ),
       xVoxl_GetFloatY ( ipVoxSrc ),
       xVoxl_GetFloatZ ( ipVoxSrc ) );
}

char xVoxl_IsEqualInt ( xVoxelRef ipVox1, xVoxelRef ipVox2 ) {

  if ( xVoxl_GetX(ipVox1) == xVoxl_GetX(ipVox2)
       && xVoxl_GetY(ipVox1) == xVoxl_GetY(ipVox2)
       && xVoxl_GetZ(ipVox1) == xVoxl_GetZ(ipVox2) ) {
    return TRUE;
  } else {
    return FALSE;
  }
}

char xVoxl_IsEqualFloat ( xVoxelRef ipVox1, xVoxelRef ipVox2 ) {

  if ( xVoxl_GetFloatX(ipVox1) == xVoxl_GetFloatX(ipVox2)
       && xVoxl_GetFloatY(ipVox1) == xVoxl_GetFloatY(ipVox2)
       && xVoxl_GetFloatZ(ipVox1) == xVoxl_GetFloatZ(ipVox2) )
    return TRUE;
  else
    return FALSE;
}

inline
void xVoxl_Set ( xVoxelRef this, int x, int y, int z ) {
  
  this->mfX = x;
  this->mfY = y;
  this->mfZ = z;
}

inline
int xVoxl_GetX ( xVoxelRef this ) {

  return this->mfX;
}

inline
int xVoxl_GetY ( xVoxelRef this ) {

  return this->mfY;
}

inline
int xVoxl_GetZ ( xVoxelRef this ) {

  return this->mfZ;
}

inline
int xVoxl_GetRoundX ( xVoxelRef this ) {

  return rint( this->mfX );
}

inline
int xVoxl_GetRoundY ( xVoxelRef this ) {

  return rint( this->mfY );
}

inline
int xVoxl_GetRoundZ ( xVoxelRef this ) {

  return rint( this->mfZ );
}

inline 
int xVoxl_GetI ( xVoxelRef this ) {
  return xVoxl_GetX ( this );
}

inline 
int xVoxl_GetJ ( xVoxelRef this ) {
  return xVoxl_GetY ( this );
}

inline 
int xVoxl_GetK ( xVoxelRef this ) {
  return xVoxl_GetZ ( this );
}

inline 
void xVoxl_SetFloat ( xVoxelRef this, float x, float y, float z ) {

  this->mfX = x;
  this->mfY = y;
  this->mfZ = z;
}

inline void xVoxl_SetX ( xVoxelRef this, int x ) {

  this->mfX = (float)x;
}

inline void xVoxl_SetY ( xVoxelRef this, int y ) {

  this->mfY = (float)y;
}

inline void xVoxl_SetZ ( xVoxelRef this, int z ) {

  this->mfZ = (float)z;
}

inline void xVoxl_SetFloatX ( xVoxelRef this, float x ) {

  this->mfX = x;
}

inline void xVoxl_SetFloatY ( xVoxelRef this, float y ) {

  this->mfY = y;
}

inline void xVoxl_SetFloatZ ( xVoxelRef this, float z ) {

  this->mfZ = z;
}

inline 
float xVoxl_GetFloatX ( xVoxelRef this ) {

  return this->mfX;
}

inline 
float xVoxl_GetFloatY ( xVoxelRef this ) {

  return this->mfY;
}

inline
float xVoxl_GetFloatZ ( xVoxelRef this ) {

  return this->mfZ;
}

tBoolean xVoxl_IncrementUntilLimit ( xVoxelRef this, int inLimit ) {

  if( this->mfX < inLimit ) {
    this->mfX += 1;
    return TRUE;
  } else if( this->mfY < inLimit ) {
    this->mfX = 0;
    this->mfY += 1;
    return TRUE;
  } else if( this->mfZ < inLimit ) {
    this->mfX = 0;
    this->mfY = 0;
    this->mfZ += 1;
    return TRUE;
  } else {
    return FALSE;
  }
}

tBoolean xVoxl_IncrementUntilLimits ( xVoxelRef this, int inXLimit, 
              int inYLimit, int inZLimit ) {

  if( this->mfX < inXLimit ) {
    this->mfX += 1;
    return TRUE;
  } else if( this->mfY < inYLimit ) {
    this->mfX = 0;
    this->mfY += 1;
    return TRUE;
  } else if( this->mfZ < inZLimit ) {
    this->mfX = 0;
    this->mfY = 0;
    this->mfZ += 1;
    return TRUE;
  } else {
    return FALSE;
  }
}

int xVoxl_ExpandToIndex ( xVoxelRef this, int inDimensionX, int inDimensionY ) {

  return ( this->mfZ * inDimensionX * inDimensionY ) +
    ( this->mfY * inDimensionX  ) + this->mfX;
}

void
xVoxl_PrintDebug ( xVoxelRef this ) {

  DebugPrint( ("Voxel: %d, %d, %d\n", 
    xVoxl_GetX(this), xVoxl_GetY(this), xVoxl_GetZ(this) ) );
}
