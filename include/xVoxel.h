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

inline int xVoxl_GetX ( xVoxelRef this );
inline int xVoxl_GetY ( xVoxelRef this );
inline int xVoxl_GetZ ( xVoxelRef this );

inline int xVoxl_GetI ( xVoxelRef this );
inline int xVoxl_GetJ ( xVoxelRef this );
inline int xVoxl_GetK ( xVoxelRef this );

inline float xVoxl_GetFloatX ( xVoxelRef this );
inline float xVoxl_GetFloatY ( xVoxelRef this );
inline float xVoxl_GetFloatZ ( xVoxelRef this );

void xVoxl_PrintDebug ( xVoxelRef  this );

/* easy way to expand the voxel struct into x/y/z coords */
#define xVoxl_ExpandFloat(v) xVoxl_GetFloatX(v), xVoxl_GetFloatY(v), xVoxl_GetFloatZ(v)
#define xVoxl_ExpandInt(v)   xVoxl_GetX(v),      xVoxl_GetY(v),      xVoxl_GetZ(v)



#endif







