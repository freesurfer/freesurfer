#ifndef mriSurface_h
#define mriSurface_h

#include "mriTypes.h"
#include "mrisurf.h"
#include "mriTransform.h"
#include "xVoxel.h"

typedef enum {

  Surf_tErr_NoErr = 0,
  Surf_tErr_InvalidObject,
  Surf_tErr_InvalidParameter,
  Surf_tErr_InvalidSignature,
  Surf_tErr_AllocationFailed,
  Surf_tErr_ErrorLoadingSurface,
  Surf_tErr_ErrorLoadingVertexSet,
  Surf_tErr_ErrorAccesssingSurface,
  Surf_tErr_LastFace,
  Surf_tErr_LastVertex,
  Surf_tErr_InvalidErrorCode,
  Surf_knNumErrorCodes
} Surf_tErr;

typedef enum {

  Surf_tOrientation_X = 0,
  Surf_tOrientation_Y,
  Surf_tOrientation_Z,
  Surf_knNumOrientations
} Surf_tOrientation;

typedef enum { 

  Surf_tVertexSet_None = -1,
  Surf_tVertexSet_Main = 0,
  Surf_tVertexSet_Original,
  Surf_tVertexSet_Pial,
  Surf_knNumVertexSets
} Surf_tVertexSet;


#define Surf_kSignature 0x0934cba4

typedef struct {

  xVoxel mVoxel[Surf_knNumVertexSets][VERTICES_PER_FACE];

} Surf_tFace, *Surf_tFaceRef;

typedef struct {

  long mSignature;

  /* surface object */
  MRIS* mSurface;
  
  /* face list */
  Surf_tFace* maFaces; 

  /* for checking which faces to iterate over */
  float mfLongestEdge;

  /* transform object. a->b = client->surface space */
  mriTransformRef mTransform;

  /* load status */
  tBoolean mabVertexSetLoaded[ Surf_knNumVertexSets ];

  /* iterator state */
  Surf_tOrientation mIterOrientation;
  float             mfIterPlane;
  int               mnCurFace;
  int               mnCurVertex;

} mriSurface, *mriSurfaceRef;

/* transformer should be a->b client->surface coordinate system */
Surf_tErr Surf_New    ( mriSurfaceRef*  opSurface,
      char*           isFileName,
      mriTransformRef iTransform );
Surf_tErr Surf_Delete ( mriSurfaceRef* iopSurface );

/* loads a vertex set. if Main, will reload the entire MRIS. if original
   or pial, will just shift out the vertex sets */
Surf_tErr Surf_LoadVertexSet ( mriSurfaceRef   this,
             char*           isName,
             Surf_tVertexSet iSet );

/* return status of a vertex set */
Surf_tErr Surf_IsVertexSetLoaded ( mriSurfaceRef   this,
           Surf_tVertexSet iSet,
           tBoolean*       obIsLoaded );

/* for each face in the mris surface, create a Surf_tFace, and put the
   vertex values in it in client coord system */
Surf_tErr Surf_ConvertSurfaceToClientSpace_ ( mriSurfaceRef   this,
                Surf_tVertexSet iSet );
               

/* for iterating thru the vertices in the surface. first step is to set the 
   iteration start point by passing an orientation and plane. the plane should
   be a voxel in client space with the coord that is constant on the plane
   set and the other coords set to 0. then use GetNextAndNeightborVertex to 
   get two neighboring verticies. often we compare pairs of
   verticies in faces, so we access them as next and last. returns the vertex
   in the client's coord system. returns Surf_tErr_LastVertex when done, at
   which point the next vertices returned will be in a new face. when
   all faces are done, Surf_tErr_LastFace will be returned. */
Surf_tErr Surf_SetIteratorPosition       ( mriSurfaceRef   this,
             xVoxelRef       plane );
Surf_tErr Surf_GetNextAndNeighborVertex  ( mriSurfaceRef   this,
             Surf_tVertexSet iSet, 
             xVoxelRef       oNextVoxel,
             xVoxelRef       oNeighborVoxel );

/* iteration implemetation functions */
float Surf_GetVertexValue ( vertex_type*      iVertex,
          Surf_tVertexSet   iSet,
          Surf_tOrientation iOrientation );

void Surf_ConvertVertexToVoxel ( vertex_type*    iVertex,
         Surf_tVertexSet iSet,
         mriTransformRef iTransform,
         xVoxelRef       oVoxel );

Surf_tErr Surf_Verify ( mriSurfaceRef this );
char* Surf_GetErrorString ( Surf_tErr ieCode );

#endif
