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
    omp_lock_t     bucket_lock;
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


typedef struct mht_face_t {
    // for per-vertex information that should not be stored in the MRIS FACE
    float cx,cy,cz; // centroid
} MHT_FACE;


MHBT * MHTacqBucketAtVoxIx(MRIS_HASH_TABLE *mht, int  xv, int   yv, int   zv);
MHBT * MHTacqBucket       (MRIS_HASH_TABLE *mht, float x, float y,  float z );

void MHTrelBucket(MHBT**);
void MHTrelBucketC(MHBT const **);


struct _mht 
{
    MRIS const *       const mris ;                            //
    float              vres ;                                 // Resolution of discretization
    MHTFNO_t           fno_usage;                             // To enforce consistent use of fno:  face number or vertex number
    int                nbuckets ;                             // Total # of buckets

#ifdef HAVE_OPENMP
    omp_lock_t         buckets_lock;
#endif
    MRIS_HASH_BUCKET **buckets_mustUseAcqRel[TABLE_SIZE][TABLE_SIZE] ;
    int                which_vertices ;                       // ORIGINAL, CANONICAL, CURRENT, etc.

    int                nfaces;
    MHT_FACE*          f;

    _mht(MRIS const * mris);
    virtual ~_mht();
};


struct MRIS_HASH_TABLE : public _mht {                      // Later we may make this private inheritance...

    MRIS_HASH_TABLE(MRIS const * mris) : _mht(mris) {}

  // Implement the traditional functions as virtual or static member functions
  // so they will all be appropriately changed once this
  // becomes a template class
  //
#define MHT_VIRTUAL                 virtual
#define MHT_ABSTRACT                = 0
#define MHT_STATIC_MEMBER           static
#define MHT_FUNCTION(NAME)          NAME
#define MHT_FUNCTION(NAME)          NAME
#define MHT_CONST_THIS_PARAMETER
#define MHT_CONST_THIS              const
#define MHT_THIS_PARAMETER_NOCOMMA
#define MHT_THIS_PARAMETER
#define MHT_MRIS_PARAMETER_NOCOMMA  MRIS const   *mris
#define MHT_MRIS_PARAMETER          MHT_MRIS_PARAMETER_NOCOMMA ,
#include "mrishash_traditional_functions.h"

} ;
