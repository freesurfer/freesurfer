/**
 * @file  mrishash_internals.h
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which 
 * need to compare vertices with one another or to a point.  See: 
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * Moved here by Bevin Brett
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2015/03/18 17:04:00 $
 *    $Revision: 1.26 $
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


#ifndef MRISHASH_INTERNALS_ONCE_H
#define MRISHASH_INTERNALS_ONCE_H

#include "mrishash.h"


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
  int            xsize, ysize, zsize ;
} MRIS_HASH_BUCKET, MHBT ;


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
//#define TABLE_SIZE     ((int)(FIELD_OF_VIEW / VOXEL_RES))
#define TABLE_SIZE     2000

#define WORLD_TO_VOLUME(mht,x)   (((x)+FIELD_OF_VIEW/2)/((mht)->vres))
#define WORLD_TO_VOXEL(mht,x)    ((int)(WORLD_TO_VOLUME(mht,x)))
#define VOXEL_TO_WORLD(mht,x)    ((((x)*(mht)->vres)-FIELD_OF_VIEW/2))

typedef enum {
    MHTFNO_FACE   = 0,
    MHTFNO_VERTEX = 1
} MHTFNO_t;

struct _mht 
{
  float              vres ;     /* resolution of discretization */
  MHTFNO_t           fno_usage; /* 2007-03-20 GW Added: To enforce consistent 
                                   use of fno:  face number or vertex number */
  int                nbuckets ; /* total # of buckets */
  MRIS_HASH_BUCKET **buckets[TABLE_SIZE][TABLE_SIZE] ;
  int                which_vertices ;  /* ORIGINAL, CANONICAL, CURRENT */
  struct _mht       *mhts[MAX_SURFACES] ; // for MRI_SURFACE_ARRAYs
  MRI_SURFACE const *mris[MAX_SURFACES] ;
  int                ntables ;
} ;


MHBT * MHTgetBucketAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv);
MHBT * MHTgetBucket       (MRIS_HASH_TABLE *mht, float x, float y, float z) ;

#endif
