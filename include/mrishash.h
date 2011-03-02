/**
 * @file  mrishash.h
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which 
 * need to compare vertices with one another or to a point.  See: 
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.25 $
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

//--------------------------------------------------------
// Structures, constants, macros. Include these only once
// In addition, include mrisurf.h
//--------------------------------------------------------
#ifndef MRISHASH_ONCE_H
#define MRISHASH_ONCE_H

// define the following to get a single inclusion of non-renamed
// functions
#define MRISHASH_VANILLA_FUNCS

//-----------------------------------------------------------
// [GW] REMOVED  #include "mrisurf.h"
// This is because with respect to mrisurf.h, mrishash.h is included in a
// very specific location in mrisurf.h, after definitions that mrishash,h
// needs from mrisurf.h, but before parts of mrisurf.h that needs defs
// from mrishash.h.
// So under normal conditions, if some caller wants to include mrishash,
// they *must* include mrisurf.h, and that should be enough.
//-----------------------------------------------------------

//--------------------------
typedef struct
{
  int     fno ;
} MRIS_HASH_BIN, MHB ;

//--------------------------
typedef struct
{
  MRIS_HASH_BIN  *bins ;
  int            max_bins ;
  int            nused ;
} MRIS_HASH_BUCKET, MHBT ;

//--------------------------
//-----------------------------------------------------------
// In mrishash, "voxel" means box in rectangular grid used
// to select spatial points into hash buckets. It does NOT refer to voxels
// in the sense of MRIs.
//
// WORLD   : surface space  (x,y,z)
// VOLUME  : WORLD rescaled and recentered at 200
// VOXEL   : VOLUME discretized (trunc VOLUME)
//
// Since the relationship between mrihash points and voxels is
// local to this unit, there's no special interpretation
// of x,y,z to worry about: eg: not necessarily RAS.
//-----------------------------------------------------------

// FIELD_OF_VIEW: Way more than needed even at 1mm resolution.
#define FIELD_OF_VIEW  400

// VOXEL_RES: Default value for MHT->vres for when caller doesn't set it.
#define VOXEL_RES      1.0

// TABLE_SIZE dimensions for array of hash buckets. As defined here
// TABLE_SIZE = 400.
#define TABLE_SIZE     ((int)(FIELD_OF_VIEW / VOXEL_RES))

#define WORLD_TO_VOLUME(mht,x)   (((x)+FIELD_OF_VIEW/2)/((mht)->vres))
#define WORLD_TO_VOXEL(mht,x)    ((int)(WORLD_TO_VOLUME(mht,x)))
#define VOXEL_TO_WORLD(mht,x)    ((((x)*(mht)->vres)-FIELD_OF_VIEW/2))

typedef enum {
    MHTFNO_FACE   = 0,
    MHTFNO_VERTEX = 1
} MHTFNO_t;

//--------------------------
typedef struct _mht
{
  float              vres ;     /* resolution of discretization */
  MHTFNO_t           fno_usage; /* 2007-03-20 GW Added: To enforce consistent 
                                   use of fno:  face number or vertex number */
  int                nbuckets ; /* total # of buckets */
  MRIS_HASH_BUCKET **buckets[TABLE_SIZE][TABLE_SIZE] ;
  int                which_vertices ;  /* ORIGINAL, CANONICAL, CURRENT */
  struct _mht       *mhts[MAX_SURFACES] ; // for MRI_SURFACE_ARRAYs
  MRI_SURFACE       *mris[MAX_SURFACES] ;
  int                ntables ;
} MRIS_HASH_TABLE, MHT ;

//------------------------------------------------
// GW wew functions post V1.27
//------------------------------------------------

MHBT * MHTgetBucketAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv);

// Ad hoc test functions
int MHT_gw_version(void);  // version of that unit
void MHTfindReportCounts(int * BucketsChecked, 
                         int * BucketsPresent, 
                         int * VtxNumByMHT);
int MHTtestIsMRISselfIntersecting(MRI_SURFACE *mris, float res);

#endif // END #ifndef MRISHASH_ONCE_H


//------------------------------------------------
// Surface --> MHT, store Face Numbers
//------------------------------------------------
MRIS_HASH_TABLE *MHTfillTable(MRI_SURFACE *mris, MRIS_HASH_TABLE *mht) ;

MRIS_HASH_TABLE *MHTfillTableAtResolution(MRI_SURFACE *mris, 
                                          MRIS_HASH_TABLE *mht,
                                          int which, 
                                          float res) ;

// Add/remove the faces of which vertex V is a part
int  MHTaddAllFaces(   MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v) ;
int  MHTremoveAllFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v) ;

//------------------------------------------------
// Surface --> MHT, store Vertex Numbers
//------------------------------------------------
MRIS_HASH_TABLE *MHTfillVertexTable(MRI_SURFACE *mris,
                                    MRIS_HASH_TABLE *mht,
                                    int which) ;
MRIS_HASH_TABLE *MHTfillVertexTableRes(MRI_SURFACE *mris,
                                       MRIS_HASH_TABLE *mht,
                                       int which,
                                       float res) ;

//------------------------------------------------
// Surface self-intersection (Uses MHT initialized with FACES)
//------------------------------------------------
int MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris,int fno);


int MHTisVectorFilled(MRIS_HASH_TABLE *mht,    MRI_SURFACE *mris,
                         int vno,  float dx, float dy, float dz) ;

//------------------------------------------------
// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//------------------------------------------------
//------- new generic find function ------------
int MHTfindClosestVertexGeneric(MRIS_HASH_TABLE *mht, 
                                MRI_SURFACE *mris,
                                double probex, double probey, double probez,
                                double in_max_distance_mm, 
                                int in_max_halfmhts,
                                VERTEX **pvtx, 
                                int *vtxnum, 
                                double *vtx_distance);

//------- original mrishash find functions ------------
VERTEX *MHTfindClosestVertex(MRIS_HASH_TABLE *mht, 
                             MRI_SURFACE *mris, 
                             VERTEX *v) ;
VERTEX *MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht, 
                                MRI_SURFACE *mris, 
                                VERTEX *v, 
                                int which) ;
VERTEX * MHTfindClosestVertexSetInDirection(MRIS_HASH_TABLE *mht, 
                                            MRI_SURFACE *mris, 
                                            VERTEX *v, 
                                            int which,
                                            double nx, double ny, double nz);
int    *MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht, 
                                        MRI_SURFACE *mris,
                                        int vno, 
                                        float max_dist, 
                                        int *pvnum);
int     MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht, 
                               MRI_SURFACE *mris, 
                               VERTEX *v, 
                               float *min_dist);
VERTEX *MHTfindClosestVertexInTable(MRIS_HASH_TABLE *mht, 
                                    MRI_SURFACE *mris,
                                    float x, float y, float z, int do_global_search) ;

//------------------------------------------------
//  Utility
//------------------------------------------------
// See also:
// MHBT * MHTgetBucketAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv);

MHBT * MHTgetBucket(MRIS_HASH_TABLE *mht, float x, float y, float z) ;
int    MHTfree(MRIS_HASH_TABLE **pmht) ;

//------------------------------------------------
// Diagnostic
//------------------------------------------------
int MHTcheckFaces(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht) ;
int MHTcheckSurface(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht);


//------------------------------------------------
// utilities for finding closest face
//------------------------------------------------
int MHTfindClosestFaceGeneric(MRIS_HASH_TABLE *mht, 
                              MRI_SURFACE *mris,
                              //---------- inputs --------------
                              double probex, double probey, double probez,
                              // How far to search: set one or both
                              double in_max_distance_mm, /* Use large number 
                                                            to ignore */
                              int    in_max_mhts,  /* Use -1 to ignore */
                              // only faces that projection is interior to (Use -1 to ignore )
                              int    project_into_face, 
                              //---------- outputs -------------
                              FACE **pface, 
                              int *pfno, 
                              double *pface_distance);
int mhtBruteForceClosestFace(MRI_SURFACE *mris, 
                             float x, float y, float z, 
                             int which,                  // which surface within mris to search
                             float *dmin);
