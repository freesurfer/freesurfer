/**
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which 
 * need to compare vertices with one another or to a point.  See: 
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * Moved here by Bevin Brett
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


#pragma once

#include "mrishash.h"
#include "romp_support.h"


//--------------------------
typedef struct
{
    int fno ;
} MRIS_HASH_BIN, MHB ;


//--------------------------
typedef struct MRIS_HASH_BUCKET
{
#ifdef HAVE_OPENMP
    omp_lock_t     mutable bucket_lock;
#endif
    MRIS_HASH_BIN  * const bins ;
    int              const max_bins ;
    int                    nused ;
    int                    size, ysize, zsize ;
} MHBT ;


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

// VOXEL_RES: Default value for MHT->vres for when caller doesn't set it.
#define VOXEL_RES      1.0

// TABLE_SIZE dimensions for array of hash buckets. As defined here
#define TABLE_SIZE     2000
#define TABLE_CENTER   (int)(TABLE_SIZE / 2)

typedef struct mht_face_t {
    // for per-vertex information that should not be stored in the MRIS FACE
    float cx,cy,cz; // centroid
} MHT_FACE;


MHBT * MHTacqBucketAtVoxIx(MRIS_HASH_TABLE *mht, int  xv, int   yv, int   zv);
MHBT * MHTacqBucket       (MRIS_HASH_TABLE *mht, float x, float y,  float z );

void MHTrelBucket(MHBT**);
void MHTrelBucketC(MHBT const **);


