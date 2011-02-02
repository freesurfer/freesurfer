/**
 * @file  mriSurface.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.18 $
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


#ifndef mriSurface_h
#define mriSurface_h

#include "mriTypes.h"
#include "mrisurf.h"
#include "mriTransform.h"
#include "xVoxel.h"
#include "colortab.h"

typedef enum
{

  Surf_tErr_NoErr = 0,
  Surf_tErr_InvalidObject,
  Surf_tErr_InvalidParameter,
  Surf_tErr_InvalidSignature,
  Surf_tErr_AllocationFailed,
  Surf_tErr_ErrorLoadingSurface,
  Surf_tErr_ErrorLoadingVertexSet,
  Surf_tErr_ErrorLoadingAnnotation,
  Surf_tErr_ErrorAccesssingSurface,
  Surf_tErr_LastFace,
  Surf_tErr_LastVertex,
  Surf_tErr_InvalidErrorCode,
  Surf_knNumErrorCodes
} Surf_tErr;

typedef enum
{

  Surf_tOrientation_X = 0,
  Surf_tOrientation_Y,
  Surf_tOrientation_Z,
  Surf_knNumOrientations
} Surf_tOrientation;

typedef enum
{

  Surf_tVertexSet_None = -1,
  Surf_tVertexSet_Main = 0,
  Surf_tVertexSet_Original,
  Surf_tVertexSet_Pial,
  Surf_knNumVertexSets
} Surf_tVertexSet;

typedef enum
{

  Surf_tValueSet_None = -1,
  Surf_tValueSet_Val = 0,

  Surf_knNumValueSets
} Surf_tValueSet;

#define Surf_kSignature 0x0934cba4

typedef struct
{

  xVoxel mVoxel[Surf_knNumVertexSets][VERTICES_PER_FACE];
  int    mnVertexIndex[Surf_knNumVertexSets][VERTICES_PER_FACE];

}
Surf_tFace, *Surf_tFaceRef;

typedef struct
{

  long mSignature;

  /* source file name */
  char* msFileName;

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
  tBoolean mabVertexSetConverted[ Surf_knNumVertexSets ];

  /* iterator state */
  Surf_tOrientation mIterOrientation;
  float             mfIterPlane;
  int               mnCurFace;
  int               mnCurVertex;

}
mriSurface, *mriSurfaceRef;

/* transformer should be a->b client->surface coordinate system */
Surf_tErr Surf_New    ( mriSurfaceRef*  opSurface,
                        char*           isFileName );
Surf_tErr Surf_Delete ( mriSurfaceRef* iopSurface );

/* Set the transform and precalc all the vertex positions in client
   space. */
Surf_tErr Surf_SetTransform ( mriSurfaceRef this,
                              mriTransformRef iTransform );

/* ==================================================================== IO */

/* loads a vertex set. if Main, will reload the entire MRIS. if original
   or pial, will just shift out the vertex sets */
Surf_tErr Surf_LoadVertexSet ( mriSurfaceRef   this,
                               char*           isName,
                               Surf_tVertexSet iSet );

/* writes the val fields to a separate file. */
Surf_tErr Surf_WriteValues ( mriSurfaceRef  this,
                             char*          isFileName );

/* return status of a vertex set */
Surf_tErr Surf_IsVertexSetLoaded ( mriSurfaceRef   this,
                                   Surf_tVertexSet iSet,
                                   tBoolean*       obIsLoaded );

/* read an annotation file */
Surf_tErr Surf_LoadAnnotation ( mriSurfaceRef this,
                                char*         isFileName );

/* Query and return a color table based on internal data, like if a
   annotation had an embedded color table. This allocates a new table
   and the caller is repsonsible for deleting it. */
Surf_tErr Surf_IsInternalColorTablePresent ( mriSurfaceRef this,
    tBoolean*     obPresent );
Surf_tErr Surf_NewColorTableFromInternal   ( mriSurfaceRef           this,
    COLOR_TABLE**           opTable );


/* ======================================================= Vertex iteration */

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
    int*            onNextIndex,
    xVoxelRef       oNeighborVoxel,
    int*            onNeighborIndex);

/* ==================================================== Geometry managament */

/* Copies the MRIS's VOL_GEOM field. */
Surf_tErr Surf_CopyGeometryInformation ( mriSurfaceRef this,
    VOL_GEOM*     ioVolumeGeometry );

Surf_tErr Surf_TransformToVolumeGeometry ( mriSurfaceRef this,
    VOL_GEOM*     iVolumeGeometry );

/* ========================================================== Vertex access */

/* get a vertex by index in voxel space */
Surf_tErr Surf_GetNthVertex ( mriSurfaceRef   this,
                              Surf_tVertexSet iSet,
                              int             inIndex,
                              xVoxelRef       oClientVoxel,
                              char*           osDescription );

/* find the closest vertex to the given location */
Surf_tErr Surf_GetClosestVertexVoxel ( mriSurfaceRef   this,
                                       Surf_tVertexSet iSet,
                                       xVoxelRef       iClientVoxel,
                                       xVoxelRef       oClientVoxel,
                                       char*           osDescription );

/* sets a vertex value within the MRIS structure. */
Surf_tErr Surf_SetVertexValue ( mriSurfaceRef   this,
                                Surf_tVertexSet iVertexSet,
                                Surf_tValueSet  iValueSet,
                                xVoxelRef       iClientVoxel,
                                float           iValue );
Surf_tErr Surf_GetVertexValue ( mriSurfaceRef   this,
                                Surf_tVertexSet iVertexSet,
                                Surf_tValueSet  iValueSet,
                                xVoxelRef       iClientVoxel,
                                float*          opValue );

/* gets an annotation value from a vertex, based on vno. */
Surf_tErr Surf_GetVertexAnnotationByIndex ( mriSurfaceRef   this,
    int             iVertexIndex,
    int*            oAnnotation );

/* get distance between two points in surface space. */
Surf_tErr Surf_GetDistance ( mriSurfaceRef this,
                             xVoxelRef     iClientVoxel1,
                             xVoxelRef     iClientVoxel2,
                             float*        ofDistance );

/* Sometimes it's easier to just give access to the surface. */
Surf_tErr Surf_GetMRIS ( mriSurfaceRef this,
                         MRIS**        opSurface );

Surf_tErr Surf_GetSurfaceSetName ( Surf_tVertexSet iSet,
                                   char*           osName );


Surf_tErr Surf_UsesRealRAS ( mriSurfaceRef this,
                             tBoolean*     obUseRealRAS );

/* ============================================================== Processing */

/* Most of these are just wrappers to functions in mrisurf.c */

Surf_tErr Surf_AverageVertexPositions ( mriSurfaceRef this,
                                        int           inNumAverages );

/* =============================================================== Internal  */

/* for each face in the mris surface, create a Surf_tFace, and put the
   vertex values in it in client coord system */
Surf_tErr Surf_ConvertSurfaceToClientSpace_ ( mriSurfaceRef   this,
    Surf_tVertexSet iSet );


/* helper functions */
float Surf_GetVertexCoord ( vertex_type*      iVertex,
                            Surf_tVertexSet   iSet,
                            Surf_tOrientation iOrientation );
void Surf_GetClosestVertex ( mriSurfaceRef   this,
                             Surf_tVertexSet iSet,
                             xVoxelRef       iSurfaceVoxel,
                             vertex_type**   opVertex,
                             int*            onIndex,
                             float*          ofDistance);
void Surf_ConvertVertexToVoxel ( vertex_type*    iVertex,
                                 Surf_tVertexSet iSet,
                                 mriTransformRef iTransform,
                                 xVoxelRef       oVoxel );
void Surf_ConvertVoxelToSurfaceSpace ( xVoxelRef       iVoxel,
                                       mriTransformRef iTransform,
                                       xVoxelRef       oSurfVox );

Surf_tErr Surf_Verify ( mriSurfaceRef this );
char* Surf_GetErrorString ( Surf_tErr ieCode );

#endif
