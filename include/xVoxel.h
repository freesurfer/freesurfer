#ifndef xxVoxl_H
#define xxVoxl_H

#include "xTypes.h"

#ifdef IRIX
#define inline
#endif

typedef struct {

  float mfX, mfY, mfZ;

} xVoxel, *xVoxelRef;


void xVoxl_New    ( xVoxelRef* oppVoxel );
void xVoxl_Delete ( xVoxelRef* ioppVoxel );

void xVoxl_Copy    ( xVoxelRef ipVoxDest, xVoxelRef ipVoxSrc );

char xVoxl_IsEqualInt   ( xVoxelRef ipVox1, xVoxelRef ipVox2 );
char xVoxl_IsEqualFloat ( xVoxelRef ipVox1, xVoxelRef ipVox2 );

inline void xVoxl_Set      ( xVoxelRef this,   int x,   int y,   int z );
inline void xVoxl_SetFloat ( xVoxelRef this, float x, float y, float z );

inline void xVoxl_SetX ( xVoxelRef this, int x );
inline void xVoxl_SetY ( xVoxelRef this, int y );
inline void xVoxl_SetZ ( xVoxelRef this, int z );

inline void xVoxl_SetFloatX ( xVoxelRef this, float x );
inline void xVoxl_SetFloatY ( xVoxelRef this, float y );
inline void xVoxl_SetFloatZ ( xVoxelRef this, float z );

inline int xVoxl_GetX ( xVoxelRef this );
inline int xVoxl_GetY ( xVoxelRef this );
inline int xVoxl_GetZ ( xVoxelRef this );

inline int xVoxl_GetRoundX ( xVoxelRef this );
inline int xVoxl_GetRoundY ( xVoxelRef this );
inline int xVoxl_GetRoundZ ( xVoxelRef this );

inline int xVoxl_GetI ( xVoxelRef this );
inline int xVoxl_GetJ ( xVoxelRef this );
inline int xVoxl_GetK ( xVoxelRef this );

inline float xVoxl_GetFloatX ( xVoxelRef this );
inline float xVoxl_GetFloatY ( xVoxelRef this );
inline float xVoxl_GetFloatZ ( xVoxelRef this );

tBoolean xVoxl_IncrementUntilLimit ( xVoxelRef this, float inLimit );
tBoolean xVoxl_IncrementWithMinUntilLimit ( xVoxelRef this, 
              float inMin, float inLimit );
tBoolean xVoxl_IncrementUntilLimits ( xVoxelRef this, float inXLimit, 
              float inYLimit, float inZLimit );

int xVoxl_ExpandToIndex ( xVoxelRef this, int inDimensionX, int inDimensionY );

void xVoxl_PrintDebug ( xVoxelRef  this );

/* easy way to expand the voxel struct into x/y/z coords */
#define xVoxl_ExpandFloat(v) xVoxl_GetFloatX(v), xVoxl_GetFloatY(v), xVoxl_GetFloatZ(v)
#define xVoxl_ExpandInt(v)   xVoxl_GetX(v),      xVoxl_GetY(v),      xVoxl_GetZ(v)
#define xVoxl_ExpandRint(v)   (int)rint(xVoxl_GetX(v)), (int)rint(xVoxl_GetY(v)), (int)rint(xVoxl_GetZ(v))



#endif







