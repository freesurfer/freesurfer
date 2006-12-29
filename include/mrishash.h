/**
 * @file  mrishash.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.16 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef MRISHASH_H
#define MRISHASH_H

#include "mrisurf.h" // MRI_SURFACE

#if 0
/* BF - should be defined in mrisurf.h */
/* kt - wasn't defined? */
#ifndef CURRENT_VERTICES
#define CURRENT_VERTICES   1
#endif
#ifndef ORIGINAL_VERTICES
#define ORIGINAL_VERTICES  2
#endif
#ifndef CANONICAL_VERTICES
#define CANONICAL_VERTICES 3
#endif
#endif

typedef struct
{
  int     fno ;
}
MRIS_HASH_BIN, MHB ;

typedef struct
{
  MRIS_HASH_BIN  *bins ;
  int            max_bins ;
  int            nused ;
}
MRIS_HASH_BUCKET, MHBT ;

#define FIELD_OF_VIEW  400
#define VOXEL_RES      1.0
#define TABLE_SIZE     ((int)(FIELD_OF_VIEW / VOXEL_RES))

#define WORLD_TO_VOLUME(mht,x)   (((x)+FIELD_OF_VIEW/2)/((mht)->vres))
#define WORLD_TO_VOXEL(mht,x)    ((int)(WORLD_TO_VOLUME(mht,x)))
#define VOXEL_TO_WORLD(mht,x)    ((((x)*(mht)->vres)-FIELD_OF_VIEW/2))

typedef struct
{
  float            fov ;         /* maximum extent of surface */
  float            vres ;        /* resolution of discretization */
  int              nbuckets ;    /* total # of buckets */
  MRIS_HASH_BUCKET **buckets[TABLE_SIZE][TABLE_SIZE] ;
  int              which_vertices ;       /* ORIGINAL, CANONICAL, CURRENT */
}
MRIS_HASH_TABLE, MHT ;


MRIS_HASH_TABLE *MHTfillTable(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht) ;
MRIS_HASH_TABLE *MHTfillTableAtResolution(MRI_SURFACE *mris,
    MRIS_HASH_TABLE *mht,
    int which, float res) ;
MRIS_HASH_TABLE *MHTfillVertexTable(MRI_SURFACE *mris,
                                    MRIS_HASH_TABLE *mht,
                                    int which) ;
MRIS_HASH_TABLE *MHTfillVertexTableRes(MRI_SURFACE *mris,
                                       MRIS_HASH_TABLE *mht,
                                       int which,
                                       float res) ;
int             MHTfree(MRIS_HASH_TABLE **pmht) ;
int             MHTcheckFaces(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht) ;
int             MHTcheckSurface(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht) ;
int             MHTisFilled(MRIS_HASH_TABLE *mht,
                            MRI_SURFACE *mris, int fno,
                            float xw, float yw, float zw);
int             MHTisVoxelFilled(MRIS_HASH_TABLE *mht,
                                 MRI_SURFACE *mris,
                                 int vno,
                                 int xv, int yv, int zv) ;
int             MHTisVectorFilled(MRIS_HASH_TABLE *mht,
                                  MRI_SURFACE *mris,
                                  int vno,
                                  float dx, float dy, float dz) ;
int             MHTaddAllFaces(MRIS_HASH_TABLE *mht,
                               MRI_SURFACE *mris,
                               VERTEX *v) ;
int             MHTremoveAllFaces(MRIS_HASH_TABLE *mht,
                                  MRI_SURFACE *mris,
                                  VERTEX *v) ;
MHBT            *MHTgetBucket(MRIS_HASH_TABLE *mht,
                              float x, float y, float z) ;

VERTEX          *MHTfindClosestVertex(MRIS_HASH_TABLE *mht,
                                      MRI_SURFACE *mris, VERTEX *v) ;
VERTEX          *MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht,
    MRI_SURFACE *mris, VERTEX *v, int which) ;
int             *MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht,
    MRI_SURFACE *mris,
    int vno, float max_dist,
    int *pvnum);
int MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht,
                           MRI_SURFACE *mris,
                           VERTEX *v,
                           float *min_dist);
VERTEX *MHTfindClosestVertexInTable(MRIS_HASH_TABLE *mht,
                                    MRI_SURFACE *mris,
                                    float x, float y, float z) ;
int MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris,int fno);


#endif
