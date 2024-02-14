/**
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which
 * need to compare vertices with one another or to a point.  See:
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
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

#include <math.h>
#include <stdlib.h>

//----------------------------------------------------
// Includes that differ for linux vs GW BC compile
//----------------------------------------------------
#ifdef __BORLANDC__
#include <mem.h>
#include "../mrishash_include/mrisurf.h"
#include "mrishash_include/error.h"
//... above includes vanilla mrihash.h
#include "mrishash_include/mrishash_gwmisc.h"
#else
#include <string.h>

#include "diag.h"
#include "error.h"

#include "mrisurf.h"
#include "mrishash_SurfaceFromMRIS.h"
#include "mrishash_SurfaceFromMRISPV.h"

#include "tritri.h"

#endif

#include "face_barycentric_coords.h"
#include "mrishash_internals.h"

//==================================================================
// Local macros
//==================================================================

// GW_VERSION reported by MHTmrishash_gw_version()
// (This is basically a reality check for compile and link)
//
#define GW_VERSION 126

int MHT_gw_version()
{
    return GW_VERSION;  // <-- change this as needed
}



// __func__ in gcc appears not to be a predefined macro that can be
// tested and also glued to string literals. Instead we have to use
// printf-style

#define __MYFUNCTION__ __func__

#ifdef __FUNCTION__
#undef __MYFUNCTION__
#define __MYFUNCTION__ __FUNCTION__
#endif

#ifdef __FUNC__
#undef __MYFUNCTION__
#define __MYFUNCTION__ __FUNC__
#endif

//==================================================================
// Local structure types
//==================================================================

typedef struct
{
  int xv;
  int yv;
  int zv;
} VOXEL_COORD;

//----------------------------------------------------------------------
// VOXEL_LISTgw adopted from mrishash.c 1.27
// 10000 seems way excessive: suggests that a face can span 10000 voxels,
// which, even with vres()=1mm, means a triangle approx 100x200 mm!
// But better safe than sorry, I guess.
//----------------------------------------------------------------------
#define MAX_VOXELS 10000
typedef struct
{
  int nused;
  int voxels[MAX_VOXELS][3];
} VOXEL_LISTgw;

// Ad hoc points
typedef struct
{
  double x;
  double y;
  double z;
} Ptdbl_t;

typedef struct mht_triangle_t {
    Ptdbl_t corners[3];
} MHT_TRIANGLE;


// Support for optional parallelism
// since don't want to slow down the serial version with unnecessary locks
//
#ifdef HAVE_OPENMP
static int parallelLevel;
#endif

void MHT_maybeParallel_begin()
{
#ifdef HAVE_OPENMP
#pragma omp atomic
    parallelLevel++;
#endif
}

void MHT_maybeParallel_end()
{
#ifdef HAVE_OPENMP
#pragma omp atomic
    --parallelLevel;
#endif
}

static void checkThread0() 
{
#ifdef HAVE_OPENMP
  int tid = omp_get_thread_num();
  if (tid != 0) {
    fprintf(stderr, "lock or unlock, not thread 0, but claiming no parallelism\n");
    fprintf(stdout, "lock or unlock, not thread 0, but claiming no parallelism\n");
    fprintf(stdout, "You may need to surround the hash code with   MHT_maybeParallel_{begin,end}()\n");
    *(int*)(-1) = 0;
    exit(1);
  }
#endif
}
 

//==================================================================
// Static forward declarations
//==================================================================
// Note: Exposed API functions start with uppercase MHT prefix.
// Static (private) functions start with lowercase mht prefix.


// Utilities that do not require a MRIS_HASH_TABLE
//
static void mhtVoxelList_Init(VOXEL_LISTgw *voxlist)
{
    voxlist->nused = 0;
}

static void mhtVoxelList_Add(VOXEL_LISTgw *voxlist, int xv, int yv, int zv)
{
    for (int i = 0; i < voxlist->nused; i++)
        if (voxlist->voxels[i][0] == xv && voxlist->voxels[i][1] == yv && voxlist->voxels[i][2] == zv) return;

    if (voxlist->nused >= MAX_VOXELS) {
        fprintf(stderr, "%s(%d, %d, %d): complete list too big!\n", __MYFUNCTION__, xv, yv, zv);

        ErrorPrintf(ERROR_NOMEMORY, "%s(%d, %d, %d): complete list too big!", __MYFUNCTION__, xv, yv, zv);
        return;
    }

    int i = voxlist->nused;
    voxlist->voxels[i][0] = xv;
    voxlist->voxels[i][1] = yv;
    voxlist->voxels[i][2] = zv;
    voxlist->nused++;
}


static void mhtVoxelList_AddCoord(VOXEL_LISTgw *voxlist, VOXEL_COORD vc)
{
    mhtVoxelList_Add(voxlist, vc.xv, vc.yv, vc.zv);
}


static void mhtVoxelList_AddPath(VOXEL_LISTgw *voxlist, VOXEL_COORD oldvc, VOXEL_COORD newvc)
{
    int incx = -1, incy = -1, incz = -1;

    if (newvc.xv >= oldvc.xv) incx = 1;
    if (newvc.yv >= oldvc.yv) incy = 1;
    if (newvc.zv >= oldvc.zv) incz = 1;

    int count = -1;
    for (int xv = oldvc.xv; incx == 1 ? xv <= newvc.xv : xv >= newvc.xv; xv += incx) {
        for (int yv = oldvc.yv; incy == 1 ? yv <= newvc.yv : yv >= newvc.yv; yv += incy) {
            for (int zv = oldvc.zv; incz == 1 ? zv <= newvc.zv : zv >= newvc.zv; zv += incz) {
                count++;
                // ie: ignore the start voxel, assume that's already added
                // (which can be done with mhtVoxelList_AddCoord)
                if (count >= 1) {
                    mhtVoxelList_Add(voxlist, xv, yv, zv);
                }
            }
        }
    }
}


/*-------------------------------------------------
  mhtVoxelList_SampleTriangle
  Scans edges and interior of triangle, finding all mht voxels that are
  impinged upon, listing those in voxlist.
  -------------------------------------------------*/
static void mhtVoxelList_SampleTriangle(
    float mhtres, 
    Ptdbl_t const *vptin0, 
    Ptdbl_t const *vptin1, 
    Ptdbl_t const *vptin2, 
    VOXEL_LISTgw *voxlist)
{
  //------------------------------------
  const float SamplesPerMHTRes = 2.0;
  double mhtres_recip;
  Ptdbl_t vpta, vptb, vptc;  // pts opposite the Long, Medium and Short sides
  Ptdbl_t dif_b2a, dif_c2a, dif_b2c;
  Ptdbl_t delta_b2a, delta_c2a, delta_b2c;
  Ptdbl_t posn_b2a, posn_c2a, posn_b2c;
  VOXEL_COORD voxco_b2a, voxco_c2a, voxco_b2c;
  VOXEL_COORD oldvoxco_b2a, oldvoxco_c2a, oldvoxco_b2c;
  oldvoxco_b2a.xv = 0;
  oldvoxco_b2a.yv = 0;
  oldvoxco_b2a.zv = 0;
  oldvoxco_c2a.xv = 0;
  oldvoxco_c2a.yv = 0;
  oldvoxco_c2a.zv = 0;
  double LenSq_b2a, LenSq_c2a;  // Distances along sides
  double TempLenSq, TempLen, StepsReqdMain_dbl, StepsReqdRung_dbl;
  int mainstep, rungstep, mainsteps_reqd, rungsteps_reqd;
  double MainStepFrac, RungStepFrac;
  bool changed;
#define PTLENSQ(a) (a.x * a.x + a.y * a.y + a.z * a.z)
#define PTDIF(ans, a, b) \
  ans.x = a.x - b.x;     \
  ans.y = a.y - b.y;     \
  ans.z = a.z - b.z
#define PTADD(ans, a, b) \
  ans.x = a.x + b.x;     \
  ans.y = a.y + b.y;     \
  ans.z = a.z + b.z
#define PTINC(ans, a) \
  ans.x += a.x;       \
  ans.y += a.y;       \
  ans.z += a.z
#define PTMULTK(ans, b, K) \
  ans.x = b.x * K;         \
  ans.y = b.y * K;         \
  ans.z = b.z * K
#define PTABS(ans, a) \
  ans.x = abs(a);     \
  ans.y = abs(a.y);   \
  ans.z = abs(a.z)
#define W2VOL(ares, x) ((x / ares) + TABLE_CENTER)
#define W2VOX(ares, x) ((int)(W2VOL(ares, x)))
#define PTWORLD2VOXEL(ans, ares, pt) \
  ans.xv = W2VOX(ares, pt.x);        \
  ans.yv = W2VOX(ares, pt.y);        \
  ans.zv = W2VOX(ares, pt.z)
#define SAME_VOXEL(a, b) ((a.xv == b.xv) && (a.yv == b.yv) && (a.zv == b.zv))

// Following is slightly less than one to avoid voxel-edge overshoot on rungs
// (makes example that are on exact coords look better.)
#define RUNGRECIPFACTOR 0.9995

  mhtres_recip = 1;
  if (mhtres != 0) mhtres_recip = 1.0 / mhtres;  // for a little speed gain

  vpta = *vptin0;
  vptb = *vptin1;
  vptc = *vptin2;

  PTDIF(dif_b2a, vpta, vptb);
  PTDIF(dif_c2a, vpta, vptc);
  //  PTDIF(dif_b2c, vptc, vptb);  not needed

  LenSq_b2a = PTLENSQ(dif_b2a);
  LenSq_c2a = PTLENSQ(dif_c2a);
  //  LenSq_b2c = PTLENSQ(dif_b2c);  not needed

  //--------------------------------------------------
  // Calcs for "main" walk along edges b-->a and c->a,
  // (Walk rungs b-->c direction)
  //--------------------------------------------------
  TempLenSq = MAX(LenSq_b2a, LenSq_c2a);
  TempLen = sqrt(TempLenSq);
  StepsReqdMain_dbl = TempLen * SamplesPerMHTRes * mhtres_recip;
  mainsteps_reqd = ceil(StepsReqdMain_dbl);
  if (0 >= mainsteps_reqd)  // can't be less than!
    mainsteps_reqd = 1;
  MainStepFrac = 1.0 / mainsteps_reqd;

  PTMULTK(delta_b2a, dif_b2a, MainStepFrac);
  PTMULTK(delta_c2a, dif_c2a, MainStepFrac);

  // Set starting positions
  posn_b2a = vptb;
  posn_c2a = vptc;

  for (mainstep = 0; mainstep <= mainsteps_reqd; mainstep++) {
    if (mainstep >= 1) {
      PTINC(posn_b2a, delta_b2a);
      PTINC(posn_c2a, delta_c2a);
    }
    PTWORLD2VOXEL(voxco_b2a, mhtres, posn_b2a);
    PTWORLD2VOXEL(voxco_c2a, mhtres, posn_c2a);

    //---------------------------------------------
    // On zero iteration, add the two starting corners
    // to voxlist
    //---------------------------------------------
    changed = false;
    if (0 == mainstep) {
      mhtVoxelList_AddCoord(voxlist, voxco_b2a);
      mhtVoxelList_AddCoord(voxlist, voxco_c2a);
      changed = true;
    }
    else {
      //---------------------------------------------
      // If we crossed a boundary, add that "path" to voxlist
      //---------------------------------------------
      if (!SAME_VOXEL(oldvoxco_b2a, voxco_b2a)) {
        mhtVoxelList_AddPath(voxlist, oldvoxco_b2a, voxco_b2a);
        changed = true;
      }

      if (!SAME_VOXEL(oldvoxco_c2a, voxco_c2a)) {
        mhtVoxelList_AddPath(voxlist, oldvoxco_c2a, voxco_c2a);
        changed = true;
      }
    }

    if (changed) {
      oldvoxco_b2a = voxco_b2a;
      oldvoxco_c2a = voxco_c2a;
    }

    //----------------------------------------------------
    // Rung:
    // 1. Endpoints are current positions of voxco_b2a, and voxco_c2a.
    // 2. Proceed in b-->c direction, rung-ing across from point on
    //    line ba to point on line ca.
    // 3. Variables suffixed _b2c
    // 4. First "rung" is actual side bc
    //----------------------------------------------------

    if (SAME_VOXEL(voxco_c2a, voxco_b2a))
      goto rung_done; /* Nothing to be done if both ends
                         are already in same voxel */

    PTDIF(dif_b2c, posn_c2a, posn_b2a);
    TempLenSq = PTLENSQ(dif_b2c);
    TempLen = sqrt(TempLenSq);
    StepsReqdRung_dbl = TempLen * SamplesPerMHTRes * mhtres_recip;
    rungsteps_reqd = ceil(StepsReqdRung_dbl);

    if (rungsteps_reqd < 1) rungsteps_reqd = 1;

    RungStepFrac = RUNGRECIPFACTOR / rungsteps_reqd;
    PTMULTK(delta_b2c, dif_b2c, RungStepFrac);

    // Starting conditions
    posn_b2c = posn_b2a;
    voxco_b2c = voxco_b2a;
    oldvoxco_b2c = voxco_b2c;

    //----------------------------------------------------
    // Step along a rung (xxx_b2c)
    // No action on step 0 (ie: starting position) , as this voxel is
    // already in the mht as a result of the main stepping.
    // Constrastingly, we *do* need action on *last* step even though that
    // voxel is in mht, because we may need to "draw" from some
    // other voxel *to* that one
    //------------------------------------------------------------
    for (rungstep = 1; rungstep <= rungsteps_reqd; rungstep++) {
      PTINC(posn_b2c, delta_b2c);
      PTWORLD2VOXEL(voxco_b2c, mhtres, posn_b2c);

      //---------------------------------------------
      // If we crossed a boundary, add voxels in that
      // "path" to voxlist
      //---------------------------------------------
      if (!SAME_VOXEL(oldvoxco_b2c, voxco_b2c)) {
        mhtVoxelList_AddPath(voxlist, oldvoxco_b2c, voxco_b2c);
        oldvoxco_b2c = voxco_b2c;
      }
    }          // for rungstep
  rung_done:;  // semicolon avoids error message
  }            // for mainstep
}


// Utilities that do require a MRIS_HASH_TABLE
//
#define MHT_MAX_TOUCHING_FACES 10000

static void lockBucket(const MHBT *bucketc) {
#ifdef HAVE_OPENMP
    MHBT *bucket = (MHBT *)bucketc;
    if (parallelLevel) omp_set_lock(&bucket->bucket_lock); else checkThread0();
#endif
}
static void unlockBucket(const MHBT *bucketc) {
#ifdef HAVE_OPENMP
    MHBT *bucket = (MHBT *)bucketc;
    if (parallelLevel) omp_unset_lock(&bucket->bucket_lock); else checkThread0();
#endif
}

// A Bucket's ->bins and ->binsDebug 
//      starts          NULL,       0
//      each realloc    set,        incremented
//      upon free,      dangling,   0xDEADDEAD
// then the Bucket is freed
//
// But I am seeing one Bucket having problems when it bins are freed
//
static void reallocBins(MHBT* bucket, int atLeast) {
    int max_bins = MAX(1,bucket->max_bins);
    while (max_bins < atLeast) max_bins *= 2;
    if (max_bins <= bucket->max_bins) return;
  
    MHB* bins = (MHB *)realloc(bucket->bins, max_bins*sizeof(MHB));
    if (!bins)
        ErrorExit(ERROR_NO_MEMORY, "%s: could not allocate %d bins.\n", __MYFUNCTION__, bucket->max_bins);

    // The const is there to stop any other code from modifying these fields
    //
    *(MHB**)&bucket->bins     = bins;
    *(int *)&bucket->max_bins = max_bins;
}


static void freeBins(MHBT* bucket) {
    MHB* bins = bucket->bins;
    free(bins);
    *(MHB**)&bucket->bins     = NULL;
    *(int *)&bucket->max_bins = 0;
}


int FindBucketsChecked_Count;
int FindBucketsPresent_Count;
int VertexNumFoundByMHT; /* 2007-07-30 GW: Added to allow diagnostics even
                            with fallback-to-brute-force */

void MHTfindReportCounts(int *BucketsChecked, int *BucketsPresent, int *VtxNumByMHT)
{
  if (BucketsChecked) *BucketsChecked = FindBucketsChecked_Count;
  if (BucketsPresent) *BucketsPresent = FindBucketsPresent_Count;
  if (VtxNumByMHT) *VtxNumByMHT = VertexNumFoundByMHT;  // 2007-07-30 GW
}



struct MRIS_HASH_TABLE_NoSurface : public MRIS_HASH_TABLE {

#ifdef HAVE_OPENMP
    omp_lock_t mutable buckets_lock;
#endif
    int                 nbuckets ;                      // Total # of buckets
    MRIS_HASH_BUCKET **buckets_mustUseAcqRel[TABLE_SIZE][TABLE_SIZE] ;

    int                nfaces;
    MHT_FACE*          f;


    virtual MRIS_HASH_TABLE_NoSurface       * toMRIS_HASH_TABLE_NoSurface_Wkr()       { return this; }
    virtual MRIS_HASH_TABLE_NoSurface const * toMRIS_HASH_TABLE_NoSurface_Wkr() const { return this; }

    MRIS_HASH_TABLE_NoSurface(MHTFNO_t fno_usage, float vres, int which, int nfaces) 
      : MRIS_HASH_TABLE(fno_usage, vres, which), nbuckets(0), nfaces(0), f(nullptr) 
    {  
        bzero(&buckets_mustUseAcqRel, sizeof(buckets_mustUseAcqRel));
#ifdef HAVE_OPENMP
        omp_init_lock(&buckets_lock);
#endif
        this->nfaces = nfaces;
        f = (MHT_FACE*)calloc(nfaces, sizeof(MHT_FACE));
    }
    
    ~MRIS_HASH_TABLE_NoSurface();
   
    void init() 
    {
        for (int fno = 0; fno < nfaces; fno++) {

            float xt,yt,zt;
            computeFaceCentroid(which(), fno, &xt, &yt, &zt);

            MHT_FACE* face = &f[fno];
            face->cx = xt;
            face->cy = yt;
            face->cz = zt;
        }
    }
      
    virtual void computeFaceCentroid(int which, int fno, float *x, float *y, float *z) = 0;

    double WORLD_TO_VOLUME(double x) const { return (x / vres()) + TABLE_CENTER; }
    int    WORLD_TO_VOXEL (double x) const { return int(WORLD_TO_VOLUME(x));  }
    float  WORLD_TO_VOLUME(float  x) const { return (x / vres()) + TABLE_CENTER; }
    int    WORLD_TO_VOXEL (float  x) const { return int(WORLD_TO_VOLUME(x));  }

    void checkConstructedWithFaces   () const;
    void checkConstructedWithVertices() const;

    void  lockBuckets() const;
    void  unlockBuckets() const;
    
    MHBT* acqBucket       (float x, float y, float z) const;
    MHBT* acqBucket       (int xv, int yv, int zv) const;
    MHBT* acqBucketAtVoxIx(int xv, int yv, int zv) const;
    MHBT* makeAndAcqBucket(int xv, int yv, int zv);
    bool  existsBuckets2  (int xv, int yv) const;

    int mhtAddFaceOrVertexAtCoords   (float x, float y, float z, int forvnum);
    int mhtAddFaceOrVertexAtVoxIx    (int xv, int yv, int zv, int forvnum);
    int mhtRemoveFaceOrVertexAtVoxIx (int xv, int yv, int zv, int forvnum);

    void mhtFaceCentroid2xyz_float   (int fno, float *x, float *y, float *z);
        // Centroids are computed once and stored in the MHT_FACE
        // They *should* be updated when the face is moved - but I suspect that they are not!
};


MRIS_HASH_TABLE_NoSurface::~MRIS_HASH_TABLE_NoSurface() 
{
    for (int xv = 0; xv < TABLE_SIZE; xv++) {
        for (int yv = 0; yv < TABLE_SIZE; yv++) {
            if (!buckets_mustUseAcqRel[xv][yv]) continue;
            for (int zv = 0; zv < TABLE_SIZE; zv++) {
                MHBT* bucket = buckets_mustUseAcqRel[xv][yv][zv];
                if (!bucket) continue;
#ifdef HAVE_OPENMP
                omp_destroy_lock(&bucket->bucket_lock);
#endif
                if (bucket->bins) freeBins(bucket);
                ::free(bucket);
            }
            ::free(buckets_mustUseAcqRel[xv][yv]);
        }
    }

#ifdef HAVE_OPENMP
    omp_destroy_lock(&buckets_lock);
#endif
}


//  Constructors for MRIS_HASH_TABLE that deal with faces are different to those for vertices
//  They should have been a different type...
//
void MRIS_HASH_TABLE_NoSurface::checkConstructedWithVertices() const
{
    if (fno_usage() != MHTFNO_VERTEX) {
        ErrorExit(ERROR_BADPARM, "%s: mht not initialized for vertices\n", __MYFUNCTION__);
    }
}

void MRIS_HASH_TABLE_NoSurface::checkConstructedWithFaces() const
{
    if (fno_usage() != MHTFNO_FACE) {
        ErrorExit(ERROR_BADPARM, "%s: MRIS_HASH_TABLE_NoSurface not initialized for faces\n", __MYFUNCTION__);
    }
}


void MHTfree(MRIS_HASH_TABLE **pmht)
{
    auto mht = *pmht;
    if (!mht) return;

    *pmht = NULL;  // sets pointer to null to indicated free'ed
    delete mht;
}


void MRIS_HASH_TABLE_NoSurface::lockBuckets() const {
#ifdef HAVE_OPENMP
    if (parallelLevel) omp_set_lock(&buckets_lock); else checkThread0();
#endif
}
void MRIS_HASH_TABLE_NoSurface::unlockBuckets() const {
#ifdef HAVE_OPENMP
    if (parallelLevel) omp_unset_lock(&buckets_lock); else checkThread0();
#endif
}

MHBT* MRIS_HASH_TABLE_NoSurface::makeAndAcqBucket(int xv, int yv, int zv) 
{
  //-----------------------------------------------
  // Allocate space if needed
  //-----------------------------------------------
  // 1. Allocate a 1-D array at buckets_mustUseAcqRel[xv][yv]
  
  lockBuckets();
  
  if (!buckets_mustUseAcqRel[xv][yv]) {
    buckets_mustUseAcqRel[xv][yv] = (MHBT **)calloc(TABLE_SIZE, sizeof(MHBT *));
    if (!buckets_mustUseAcqRel[xv][yv]) ErrorExit(ERROR_NO_MEMORY, "%s: could not allocate slice.", __MYFUNCTION__);
  }
  
  // 2. Allocate a bucket at buckets_mustUseAcqRel[xv][yv][zv]
  MHBT *bucket = buckets_mustUseAcqRel[xv][yv][zv];
  
  
  if (!bucket) {
    buckets_mustUseAcqRel[xv][yv][zv] = bucket = (MHBT *)calloc(1, sizeof(MHBT));
    if (!bucket) ErrorExit(ERROR_NOMEMORY, "%s couldn't allocate bucket.\n", __MYFUNCTION__);
#ifdef HAVE_OPENMP
    omp_init_lock(&bucket->bucket_lock);
#endif
  }
  
  unlockBuckets();
  
  lockBucket(bucket);
  
  // 3. Allocate bins
  if (!bucket->max_bins) /* nothing in this bucket yet - allocate bins */
    reallocBins(bucket, 4);

  // returns with the bucket locked
  
  return bucket;
}


MHBT * MRIS_HASH_TABLE_NoSurface::acqBucket(float x, float y, float z) const
{
  int xv = WORLD_TO_VOXEL(x);
  int yv = WORLD_TO_VOXEL(y);
  int zv = WORLD_TO_VOXEL(z);

  return acqBucket(xv, yv, zv);
}


MHBT * MRIS_HASH_TABLE_NoSurface::acqBucketAtVoxIx(int xv, int yv, int zv) const
{
  return acqBucket(xv, yv, zv);
}


MHBT* MRIS_HASH_TABLE_NoSurface::acqBucket(int xv, int yv, int zv) const
{
  if (xv >= TABLE_SIZE || yv >= TABLE_SIZE || zv >= TABLE_SIZE || xv < 0 || yv < 0 || zv < 0) return (NULL);

  lockBuckets();

  MHBT* bucket = NULL;
  if (buckets_mustUseAcqRel[xv][yv]) 
      bucket = buckets_mustUseAcqRel[xv][yv][zv];

  unlockBuckets();
  
  // returns with the bucket, if any, locked
  //  
  if (bucket) lockBucket(bucket);
  
  return bucket;
}



static void relBucketC(MHBT const ** bucket)
{
  if (*bucket) {
    unlockBucket(*bucket);
    *bucket = NULL;
  }
}

static void relBucket(MHBT ** bucket)
{
  MHBT const * bucketC = *bucket;
  *bucket = NULL;
  relBucketC(&bucketC);
}

MHBT * MHTacqBucketAtVoxIx(MRIS_HASH_TABLE *mht, int  xv, int   yv, int   zv)
{
    return mht->toMRIS_HASH_TABLE_NoSurface()->acqBucketAtVoxIx(xv, yv, zv);
}

MHBT * MHTacqBucket(MRIS_HASH_TABLE *mht, float x, float y,  float z)
{
    return mht->toMRIS_HASH_TABLE_NoSurface()->acqBucket(x, y, z);
}

void MHTrelBucket (MHBT       ** bucket) { relBucket (bucket); }
void MHTrelBucketC(MHBT const ** bucket) { relBucketC(bucket); }


bool MRIS_HASH_TABLE_NoSurface::existsBuckets2(int xv, int yv) const
{
    bool result = false;
    if (xv >= TABLE_SIZE || yv >= TABLE_SIZE || xv < 0 || yv < 0) goto Done;
    if (!buckets_mustUseAcqRel[xv][yv]) goto Done;
    result = true;
Done:
    return result;
}

#define buckets_mustUseAcqRel SHOULD_NOT_ACCESS_BUCKETS_DIRECTLY


void MRIS_HASH_TABLE_NoSurface::mhtFaceCentroid2xyz_float(
    int fno, 
    float *px, float *py, float *pz)
{
    MHT_FACE const* face = &f[fno];
    *px = face->cx;
    *py = face->cy;
    *pz = face->cz;
}

/*------------------------------------------------------------
  mhtAddFaceOrVertexAtVoxIx  Was: mhtAddFacePosition
  Adds forvnum (Face or Vertex number) to mht, in bucket(xv,yv,zv)
  Coerces (xv,yv,zv) to be sane.
  -------------------------------------------------------------*/
int MRIS_HASH_TABLE_NoSurface::mhtAddFaceOrVertexAtVoxIx(int xv, int yv, int zv, int forvnum)
{
  if (xv < 0) xv = 0;
  if (xv >= TABLE_SIZE) xv = TABLE_SIZE - 1;
  if (yv < 0) yv = 0;
  if (yv >= TABLE_SIZE) yv = TABLE_SIZE - 1;
  if (zv < 0) zv = 0;
  if (zv >= TABLE_SIZE) zv = TABLE_SIZE - 1;

  {
    MHBT *bucket = makeAndAcqBucket(xv, yv, zv);

    //-------------------------------------------------------------
    // If this forvnum already is listed in this bucket then
    // nothing to do.
    //-------------------------------------------------------------
    MHB *bin = bucket->bins;

    int i;
    for (i = 0; i < bucket->nused; i++, bin++) {
      if (bin->fno == forvnum) goto done;
    }

    //-------------------------------------------------------------
    // Add forvnum to this bucket
    //-------------------------------------------------------------
    if (i == bucket->nused) /* forvnum not already listed at this bucket */
    {
      reallocBins(bucket, bucket->nused + 1);
      //----- add this face-position to this bucket ------
      bucket->bins[bucket->nused++].fno = forvnum;
    }

  done:
    relBucket(&bucket);
  }  

  return (NO_ERROR);
}


int MRIS_HASH_TABLE_NoSurface::mhtAddFaceOrVertexAtCoords(float x, float y, float z, int forvnum)
{
  int xv, yv, zv;
  xv = WORLD_TO_VOXEL( x);
  yv = WORLD_TO_VOXEL( y);
  zv = WORLD_TO_VOXEL( z);

  return mhtAddFaceOrVertexAtVoxIx(xv, yv, zv, forvnum);
}


// Reverse of mhtAddFaceOrVertexAtIndexes
int MRIS_HASH_TABLE_NoSurface::mhtRemoveFaceOrVertexAtVoxIx(int xv, int yv, int zv, int forvnum)
{
  int i;

  if (xv < 0) xv = 0;
  if (yv < 0) yv = 0;
  if (zv < 0) zv = 0;
  if (xv >= TABLE_SIZE) xv = TABLE_SIZE - 1;
  if (yv >= TABLE_SIZE) yv = TABLE_SIZE - 1;
  if (zv >= TABLE_SIZE) zv = TABLE_SIZE - 1;

  if (!existsBuckets2(xv,yv)) return (NO_ERROR);  // no bucket at such coordinates
  
  MHBT *bucket = acqBucket(xv,yv,zv);
  if (!bucket) return (NO_ERROR);  // no bucket at such coordinates

  MHB *bin = bucket->bins;
  for (i = 0; i < bucket->nused; i++, bin++) {
    if (bin->fno == forvnum) /* found face */
      break;
  }
  if (i < bucket->nused) /* found face bucket - remove it */
  {
    bucket->nused--;
    if (i < bucket->nused) /* not the last one in the list - compact list */
    {
      int nbytes = (bucket->nused - i) * sizeof(MHB);
      MHB *src_bin, *dst_bin;
      src_bin = &bucket->bins[i + 1];
      dst_bin = &bucket->bins[i];
      memmove(dst_bin, src_bin, nbytes);
    }
  }
  relBucket(&bucket);
  return (NO_ERROR);
}


// Now the algorithms that depend on the surface representation
//
template <class Surface, class Face, class Vertex>
struct MRIS_HASH_TABLE_IMPL : public MRIS_HASH_TABLE_NoSurface {

    static MRIS_HASH_TABLE_IMPL* newMHT(MHTFNO_t fno_usage, float vres, int which, Surface surface) 
    {
        auto mht = new MRIS_HASH_TABLE_IMPL(fno_usage, vres, surface, which);
        if (!mht) ErrorExit(ERROR_NO_MEMORY, "%s: could not allocate hash table.\n", __MYFUNCTION__);
        return mht;
    }

    static void mhtVertex2xyz(Vertex  const vtx, int which, float *x, float *y, float *z)
    {
      vtx.which_coords(which, x, y, z);
    }

    static void mhtVertex2xyz(Vertex  const vtx, int which, Ptdbl_t *pt)
    {
      float x,y,z;
      vtx.which_coords(which, &x, &y, &z);
      pt->x = x;
      pt->y = y;
      pt->z = z;
    }

    static void mhtVertex2xyz(Vertex  const vtx, int which, double *array3)
    {
      float x,y,z;
      vtx.which_coords(which, &x, &y, &z);
      array3[0] = x;
      array3[1] = y;
      array3[2] = z;
    }

    static void mhtComputeFaceCentroid(
        Face const face, int which, int fno,
        float *px, float *py, float* pz) 
    {
        float xt = 0.0, yt = 0.0, zt = 0.0;

        for (int n = 0; n < VERTICES_PER_FACE; n++) {
            float x = 0, y = 0, z = 0;
            mhtVertex2xyz(face.v(n), which, &x, &y, &z);
            xt += x;
            yt += y;
            zt += z;
        }

        xt /= VERTICES_PER_FACE;
        yt /= VERTICES_PER_FACE;
        zt /= VERTICES_PER_FACE;

        *px = xt;
        *py = yt;
        *pz = zt;
    }

    static void mhtComputeFaceCentroid(
        Surface const surface, int which, int fno,
        float *px, float *py, float* pz) 
    {
        mhtComputeFaceCentroid(surface.faces(fno), which, fno, px, py, pz);
    }

    static int BruteForceClosestFace(
        Surface const surface,
        float   x,
        float   y,
        float   z,
        int     which,  // which surface within mris to search
        float * dmin);
    
    Surface surface;

    MRIS_HASH_TABLE_IMPL(MHTFNO_t fno_usage, float vres, Surface surface, int which) 
      : MRIS_HASH_TABLE_NoSurface(fno_usage, vres, which, surface.nfaces()), surface(surface) 
    {
        init();
    }

    ~MRIS_HASH_TABLE_IMPL() {}

    virtual void computeFaceCentroid(int which, int fno, float *x, float *y, float *z) {
        mhtComputeFaceCentroid(surface, which, fno, x, y, z);
    }

    void checkConstructedWithFacesAndSurface   (Surface surface) const;
    void checkConstructedWithVerticesAndSurface(Surface surface) const;

    void captureFaceData();
    void captureVertexData();
    
    int mhtFaceToMHT                 (Face f, bool on);
    int mhtDoesTriangleVoxelListIntersect(
                                      MHT_TRIANGLE const    * const triangle, 
                                      VOXEL_LISTgw const    * const voxlistForTriangle,
                                      int                     const nFaceToIgnore,
                                      int const             * const fnoToIgnore,
                                      int                     const trace) const;


    int MHTexpandToTouchingFaces     (int   const fno,
                                      int   const fnoListCapacity,
                                      int * const fnoList,
                                      int   const trace) const;

    int MHTdoesFaceIntersect_new     (int fno, int const trace) const;
    
    int MHTdoesTriangleIntersect     (MHT_TRIANGLE const    * const triangle,
                                      int                     const nFaceToIgnore,
                                      int const             * const fnoToIgnore,
                                      int                     const trace) const;

	
    int mhtBruteForceClosestVertex(float x, float y, float z, int which, float *dmin);

    void mhtfindClosestVertexGenericInBucket(
                                        //---------- inputs --------------
                                        int xv,
                                        int yv,
                                        int zv,
                                        double probex,
                                        double probey,
                                        double probez,
                                        //---------- in/outs -------------
                                        int *MinDistVtxNum,
                                        double *MinDistSq);

    int mhtfindClosestFaceCentroidGenericInBucket(
                                              //---------- inputs --------------
                                              int xv,
                                              int yv,
                                              int zv,
                                              double probex,
                                              double probey,
                                              double probez,
                                              int project_into_face,
                                              //---------- in/outs -------------
                                              int *MinDistFaceNum,
                                              double *MinDistSq);

#define MHT_ONLY_VIRTUAL
#define MHT_VIRTUAL                 virtual
#define MHT_ABSTRACT                
#define MHT_STATIC_MEMBER           static
#define MHT_FUNCTION(NAME)          NAME
#define MHT_FUNCTION(NAME)          NAME
#define MHT_CONST_THIS_PARAMETER
#define MHT_CONST_THIS              const
#define MHT_THIS_PARAMETER_NOCOMMA
#define MHT_THIS_PARAMETER
#define MHT_MRIS_PARAMETER_NOCOMMA
#define MHT_MRIS_PARAMETER
#include "mrishash_traditional_functions.h"
#undef MHT_ONLY_VIRTUAL

};

template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::BruteForceClosestFace(
    Surface const surface,
    float   x,
    float   y,
    float   z,
    int     which,  // which surface within mris to search
    float * dmin)
{

    int   min_fno = -1;
    float min_dsq = 1e8;

    for (int fno = 0; fno < surface.nfaces(); fno++) {

        float tryx, tryy, tryz;
        mhtComputeFaceCentroid(surface, which, fno, &tryx, &tryy, &tryz);

        float const dx = tryx - x, dy = tryy - y, dz = tryz - z;

        float dsq = dx * dx + dy * dy + dz * dz;  // squared distance is fine for detecting min
        if (min_dsq <= dsq) continue;
        min_dsq = dsq;
        min_fno = fno;
    }
  
    if (dmin) *dmin = sqrt(min_dsq);

    return min_fno;
}


template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::captureFaceData()
{
    static int ncalls = 0, ncalls_limit = 1;
    ncalls++;

    // Capture data from caller and surface
    //    
    for (int fno = 0; fno < surface.nfaces(); fno++) {
        auto f = surface.faces(fno);
        if (f.ripflag()) continue;
        mhtFaceToMHT(f, true);
    }

    // Diagnostics
    //
    if ((Gdiag & DIAG_SHOW) && (ncalls == ncalls_limit)) {
        ncalls_limit *= 2;
        
        double  mean = 0.0, var = 0.0;
        int     max_nused = -1;

        for (int pass = 0; pass < 2; pass++) {
            int n = 0;
            for (int xv = 0; xv < TABLE_SIZE; xv++) {
                for (int yv = 0; yv < TABLE_SIZE; yv++) {
                    if (!existsBuckets2(xv,yv)) continue;
                    for (int zv = 0; zv < TABLE_SIZE; zv++) {

                        MHBT* bucket = acqBucket(xv,yv,zv);
                        if (!bucket) continue;

                        if (pass == 0) {
                            if (bucket->nused) {
                                mean += bucket->nused;
                                n++;
                            }
                            if (bucket->nused > max_nused) max_nused = bucket->nused;
                        } else {
                            double v = mean - bucket->nused;
                            var += v*v;
                        }
                        
                        relBucket(&bucket);
                    }
                }
            }
            if (n == 0) n = 1;
            if (pass == 0) 
                mean /= n;
            else
                var /= (n - 1);
        }

        if (0) {
            fprintf(stderr, "%s buckets: ncalls:%d mean = %2.1f +- %2.2f, max_nused = %d\n", 
                __MYFUNCTION__, ncalls, mean, sqrt(var), max_nused);
        }
    }
}
    

template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::captureVertexData()
{
    static int ncalls = 0;
    ncalls++;
    
    for (int vno = 0; vno < surface.nvertices(); vno++) {
        auto v = surface.vertices(vno);
        if (v.ripflag()) continue;
        float x, y, z;
        mhtVertex2xyz(v, which(), &x, &y, &z);
        mhtAddFaceOrVertexAtCoords(x, y, z, vno);
    }
}


//  Constructors for MRIS_HASH_TABLE that deal with faces are different to those for vertices
//  They should have been a different type...
//
template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::checkConstructedWithVerticesAndSurface(Surface surface) const
{
    if (this->surface != surface) ErrorExit(ERROR_BADPARM, "%s: mris is bad\n", __MYFUNCTION__);
    checkConstructedWithVertices();
}

template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::checkConstructedWithFacesAndSurface(Surface surface) const
{
    if (this->surface != surface) ErrorExit(ERROR_BADPARM, "%s: mris is bad\n", __MYFUNCTION__);
    checkConstructedWithFaces();
}

//  Appends or removes all faces of the vertex.
//  Returns: NO_ERROR
//
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::addAllFaces(int vno)
{
    Vertex v = surface.vertices(vno);
    for (int fi = 0; fi < v.num(); fi++) mhtFaceToMHT(v.f(fi), true);
    return (NO_ERROR);
}


template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::removeAllFaces(int vno)
{
    Vertex v = surface.vertices(vno);
    for (int fi = 0; fi < v.num(); fi++) mhtFaceToMHT(v.f(fi), false);
    return (NO_ERROR);
}


//  Adds face fno to mht. 
//  Calls mhtVoxelList_SampleTriangle to get a list of MHT Voxels (buckets) in which to list fno.
//
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtFaceToMHT(Face const face, bool const on)
{
    if (face.ripflag()) return (NO_ERROR);
    auto const fno = face.fno();

    Vertex const v0 = face.v(0);
    Vertex const v1 = face.v(1);
    Vertex const v2 = face.v(2);

    Ptdbl_t vpt0, vpt1, vpt2;
    mhtVertex2xyz(v0, which(), &vpt0);
    mhtVertex2xyz(v1, which(), &vpt1);
    mhtVertex2xyz(v2, which(), &vpt2);

    if (Gx >= 0) {
        double dist0 = sqrt(SQR(vpt0.x - Gx) + SQR(vpt0.y - Gy) + SQR(vpt0.z - Gz));
        double dist1 = sqrt(SQR(vpt1.x - Gx) + SQR(vpt1.y - Gy) + SQR(vpt1.z - Gz));
        double dist2 = sqrt(SQR(vpt2.x - Gx) + SQR(vpt2.y - Gy) + SQR(vpt2.z - Gz));
        if (dist0 < vres() || dist1 < vres() || dist2 < vres()) DiagBreak();
    }

    VOXEL_LISTgw voxlist;
    mhtVoxelList_Init(&voxlist);
    mhtVoxelList_SampleTriangle(vres(), &vpt0, &vpt1, &vpt2, &voxlist);

    for (int vlix = 0; vlix < voxlist.nused; vlix++) {

        int i = voxlist.voxels[vlix][0];
        int j = voxlist.voxels[vlix][1];
        int k = voxlist.voxels[vlix][2];

        if (on) mhtAddFaceOrVertexAtVoxIx   (i, j, k, fno);
        else    mhtRemoveFaceOrVertexAtVoxIx(i, j, k, fno);
    }

    return (NO_ERROR);
}


//=============================================================================
// MRIS_HASH_TABLE_IMPL that stores Vertex Numbers

//  Surface self-intersection (Uses MHT initialized with FACES)
//
//  If vertex vtxno is moved by (dx,dy,dz), will the
//  adjoining faces then intersect with some other faces in the surface?".
//  .Returns 1=intersect, else 0.
//
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::isVectorFilled(
    int   const vtxno, 
    float const dx, 
    float const dy, 
    float const dz) const
{
    Vertex const vtx = surface.vertices(vtxno);
    
    // Check whether each face adjoining vtxno will now intersect other faces
    //
    // There is two complications
    //
    // 1 - These faces are in the MHT in their pre-moved positions, and since only one vertex has moved, of course they might intersect!
    //            When the old code had these spurious hits, it then compared the faces based on the updated vertex positions
    //            and hence got the right answer
    //
    // 2 - The moved faces might intersect with each other, and their post-moved positions are not in the MHT, and so
    //     this intersection should be calculated without regard to the MHT, which is right because it would just report
    //     that they intersect at the shared vertex
    //            The old code updated the vertex position BUT it is comparing vno's to recognize abutting faces
    //            and was not comparing faces that share even one vertex!

    checkConstructedWithFaces();

    static int count; count++;
  
    int const trace = 0;
  
    // The following code moves x,y,z regardless of 'which()'
    // and that only makes sense if 'which()' is CURRENT_VERTICES
    //
    if (which() != CURRENT_VERTICES) {
        ErrorExit(ERROR_BADPARM, "%s: mht not loaded using CURRENT_VERTICES\n", __MYFUNCTION__);
    }

    if (trace) {
        fprintf(stderr, "The faces surrounding vertex %d are ", vtx.vno());
        for (int fi = 0; fi < vtx.num(); fi++) {
            auto face = vtx.f(fi);
            fprintf(stderr, " %d", face.fno());
        }
        fprintf(stderr, "\n");
    }

    float const moved_x = vtx.x() + dx;
    float const moved_y = vtx.y() + dy;
    float const moved_z = vtx.z() + dz;

    if (trace) {
        fprintf(stderr, "vertex:%d (%g,%g,%g) moved to (%g,%g,%g)\n",
            vtx.vno(), vtx.x(), vtx.y(), vtx.z(), moved_x, moved_y, moved_z);
    }
  
    // Try each changed faced in turn
    //
    for (int fi = 0; fi < vtx.num(); fi++) {
        Face const face = vtx.f(fi);
        
        MHT_TRIANGLE triangle;
        for (int corneri = 0; corneri < 3; corneri++) {
            Vertex corner = face.v(corneri);

            Ptdbl_t* point = &triangle.corners[corneri];
            
            if (corner.vno() == vtxno) {
                point->x = moved_x; point->y = moved_y; point->z = moved_z;
            } else {
                mhtVertex2xyz(corner, which(), point);
                if (trace) {
                    fprintf(stderr, "corner:%d vertex:%d stays at (%g,%g,%g)\n",
                        corneri, corner.vno(), point->x, point->y, point->z);
                }
            }
        }

        // Intersect with the non-changed faces that are not adjacent to this face
        // The changed faces share the vertex so will be in the touchingFnos
        //
        int touchingFnos[MHT_MAX_TOUCHING_FACES];
        int touchingFnosSize = MHTexpandToTouchingFaces(face.fno(), MHT_MAX_TOUCHING_FACES, touchingFnos, trace);
    
        int result = MHTdoesTriangleIntersect(&triangle, touchingFnosSize, touchingFnos, trace);
        if (trace) {
            fprintf(stderr, " MHT::doesFaceIntersect_new face:%d returns %d\n", face.fno(), result);
        }
        if (result) return result;
    }
  
    return 0;
}


template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::MHTdoesFaceIntersect_new(int const fno, int const trace) const
{
    checkConstructedWithFaces();
  
    auto face = surface.faces(fno);

    if (face.ripflag()) return 0;

    MHT_TRIANGLE triangle;
    for (int corneri = 0; corneri < 3; corneri++) {
        mhtVertex2xyz(face.v(corneri), which(), &triangle.corners[corneri]);
    }

    int touchingFnos[MHT_MAX_TOUCHING_FACES];
    int touchingFnosSize = MHTexpandToTouchingFaces(fno, MHT_MAX_TOUCHING_FACES, touchingFnos, trace);
  
    if (trace) {
        fprintf(stderr, "MHTdoesFaceIntersect_new ignoring ");
        for (int tfi = 0; tfi < touchingFnosSize; tfi++) {
            fprintf(stderr, " %d ", touchingFnos[tfi]);
        }
        fprintf(stderr, "\n");
    }
  
    return MHTdoesTriangleIntersect(&triangle, touchingFnosSize, touchingFnos, trace); 
}


template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::doesFaceIntersect(int fno) {
    return MHTdoesFaceIntersect_new(fno, false);
}

template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::MHTexpandToTouchingFaces(
    int   const fno,
    int   const fnoListCapacity,
    int * const fnoList,
    int   const trace) const
    //
    // Puts fno in the list, and all the faces that (by intersection rules) are touching fno
    // Returns the length of the list
{
  int size = 0;
  if (size == fnoListCapacity) { fprintf(stderr, "%s:%d exceeded fnoListCapacity\n",__FILE__,__LINE__); exit (1); }
  fnoList[size++] = fno;

  Face const face = surface.faces(fno);
  for (int vi = 0; vi < VERTICES_PER_FACE; vi++) {
    Vertex const vertex = face.v(vi);

    for (int fi = 0; fi < vertex.num(); fi++) {
        Face const face2 = vertex.f(fi);
      
        int const fno2 = face2.fno();
        if (fno2 == fno) continue;

        for (int fli = 0; fli < size; fli++) {
            if (fnoList[fli] == fno2) goto already_in_list;
        }
      
        if (size == fnoListCapacity) { fprintf(stderr, "%s:%d exceeded fnoListCapacity\n",__FILE__,__LINE__); exit (1); }
        fnoList[size++] = fno2;
      
        if (trace) fprintf(stderr, "MHTexpandToTouchingFaces adding %d\n", fno2);
 already_in_list:;
        } // each other face of the vertex
    } // each vertex of the face
    
    return size;
}


template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::MHTdoesTriangleIntersect(
    MHT_TRIANGLE    const * const triangle,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace) const
{
    // Find all the voxels that the triangle intersects
    //
    VOXEL_LISTgw voxlist;
    mhtVoxelList_Init(&voxlist);
    mhtVoxelList_SampleTriangle(vres(), &triangle->corners[0], &triangle->corners[1], &triangle->corners[2], &voxlist);

    int retval = 0;
    if (mhtDoesTriangleVoxelListIntersect(triangle, &voxlist, nFaceToIgnore, fnoToIgnore, trace)) retval = 1;
    return (retval);
}


#define MHT_MAX_FACES 10000
/*-------------------------------------------------------------------
  mhtDoesTriangleVoxelListIntersect
  Subsidiary function for MHTdoesTriangleIntersect, 
  triangle already analyzed into voxlist.
  ------------------------------------------------------------------*/
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtDoesTriangleVoxelListIntersect(
    MHT_TRIANGLE    const * const triangle, 
    VOXEL_LISTgw    const * const voxlist,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace) const
//------------------------------------------------------------------
{
  //------------------------------------------------------------------------
  // Iterate through all the voxels pertaining to the triangle, 
  // looking at those voxel buckets in mht to find faces to add to list in flist 
  // for later intersection check.
  //------------------------------------------------------------------------

  int facelist[MHT_MAX_FACES];
  int facelistSize = 0;

  int voxnum;
  for (voxnum = 0; voxnum < voxlist->nused; voxnum++) {
    int const (* voxel)[3] = &voxlist->voxels[voxnum];
    int const xv = (*voxel)[0];
    int const yv = (*voxel)[1];
    int const zv = (*voxel)[2];

    //----------------------------------------------------------
    // Get corresponding bucket from mht, if any.
    // There might not be...
    //----------------------------------------------------------
    if (!existsBuckets2(xv,yv)) continue;
    
    MHBT const * bucket = acqBucket(xv,yv,zv);
    if (!bucket) continue;

    // Consider all the faces that intersect this voxel
    //
    int binix;
    for (binix = 0; binix < bucket->nused; binix++) {
      MHB const * const bin = &bucket->bins[binix];
      int const fno = bin->fno;

      // Is this one of the faces to be ignored?
      // This is used to ignore
      //        the face that made the triangle
      //        the faces that are moving because a vertex is moving
      //
      int ignoreI;
      for (ignoreI = 0; ignoreI < nFaceToIgnore; ignoreI++) {
        if (fno == fnoToIgnore[ignoreI])
          goto skip_this_face;
      }
      
      // The test that was here has been moved into MHTexpandToTouchingFaces
      // because it is faster to decide up front which faces should be ignored

      // append it if not already present
      //
      int fli;
      for (fli = 0; fli < facelistSize; fli++) {
        if (facelist[fli] == fno) goto skip_this_face;
      }
      
      if (facelistSize >= MHT_MAX_FACES) {
        ErrorPrintf(ERROR_NO_MEMORY, "%s: MHT_MAX_FACES exceeded!", __MYFUNCTION__);
        exit (1);
      }
      
      facelist[facelistSize++] = fno;

    skip_this_face:;
    }                     // for binix
    relBucketC(&bucket);
  }                       // for voxnum

  // get vertices of 1st triangle
  //
  double v0[3], v1[3], v2[3];
  v0[0] = triangle->corners[0].x;  v0[1] = triangle->corners[0].y;  v0[2] = triangle->corners[0].z;
  v1[0] = triangle->corners[1].x;  v1[1] = triangle->corners[1].y;  v1[2] = triangle->corners[1].z;
  v2[0] = triangle->corners[2].x;  v2[1] = triangle->corners[2].y;  v2[2] = triangle->corners[2].z;

  // see if any face intersects
  //
  int fli;
  for (fli = 0; fli < facelistSize; fli++) {
    Face const face = surface.faces(facelist[fli]);
    
    double u0[3], u1[3], u2[3];   
    mhtVertex2xyz(face.v(0), which(), u0);
    mhtVertex2xyz(face.v(1), which(), u1);
    mhtVertex2xyz(face.v(2), which(), u2);

    int const intersect = tri_tri_intersect(v0, v1, v2, u0, u1, u2);
    
    if (intersect) {
      if (trace) {
        fprintf(stderr, "mhtDoesTriangleVoxelListIntersect thinks face:%d intersects\n", face.fno());
      }
      return (1);
    }
  }

  return (0);
}


//=================================================================
// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//

template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtfindClosestVertexGenericInBucket(
                                        //---------- inputs --------------
                                        int xv,
                                        int yv,
                                        int zv,
                                        double probex,
                                        double probey,
                                        double probez,
                                        //---------- in/outs -------------
                                        int     *MinDistVtxNum,
                                        double  *MinDistSq)
{
    FindBucketsChecked_Count++;

    auto bucket = acqBucketAtVoxIx(xv, yv, zv);
    if (!bucket) return;

    MHB* bin = bucket->bins;
    for (int vtxix = 0; vtxix < bucket->nused; vtxix++, bin++) {
        int const AVtxNum = bin->fno;
        auto AVtx = surface.vertices(AVtxNum);
        if (AVtx.ripflag()) continue ;

        float tryx, tryy, tryz;
        mhtVertex2xyz(AVtx, which(), &tryx, &tryy, &tryz);

        double ADistSq = SQR(tryx - probex) + SQR(tryy - probey) + SQR(tryz - probez);

        if (ADistSq >= *MinDistSq) continue;
        
        *MinDistSq      = ADistSq;
        *MinDistVtxNum  = AVtxNum;
    }
  
    relBucket(&bucket);
}

/*----------------------------------------------------------------
  mhtfindClosestFaceCentroidGenericInBucket

  find the face whose centroid is closest to the specified coordinate
  -----------------------------------------------------------------*/
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtfindClosestFaceCentroidGenericInBucket(
                                              //---------- inputs --------------
                                              int xv,
                                              int yv,
                                              int zv,
                                              double probex,
                                              double probey,
                                              double probez,
                                              int project_into_face,
                                              //---------- in/outs -------------
                                              int *MinDistFaceNum,
                                              double *MinDistSq)
{
    FindBucketsChecked_Count++;

    MHBT* bucket = acqBucketAtVoxIx(xv, yv, zv);
    if (!bucket) return NO_ERROR;

    FindBucketsPresent_Count++;

    MHB *bin = bucket->bins;
  
    for (int faceix = 0; faceix < bucket->nused; faceix++, bin++) {

        int const fno = bin->fno;

        if (fno == Gdiag_no) DiagBreak();

        float tryx, tryy, tryz;
        mhtFaceCentroid2xyz_float(fno, &tryx, &tryy, &tryz);

        double lambda[3];
        if (project_into_face > 0 &&
            face_barycentric_coords_template(
                surface, fno, which(), probex, probey, probez, &lambda[0], &lambda[1], &lambda[2]) < 0)
            continue;
    
        double ADistSq = SQR(tryx - probex) + SQR(tryy - probey) + SQR(tryz - probez);

        if (ADistSq >= *MinDistSq) continue;
        *MinDistSq      = ADistSq;
        *MinDistFaceNum = fno;
    }    // for faceix

    MHTrelBucket(&bucket);  

    return NO_ERROR;
}

/*
 ----------------------------------------------------------------
 findClosestVertexGeneric
 Generic find routine satisfying all other finds
 Inputs:
 -------
 probex, probey, probez: The point relative to which to search

 n_max_distance_mm: Furthest vertex distance to accept. Set to a large
 number to ignore, use only in_max_mhts.
 in_max_mhts       : Furthest voxels to search
 0 : One half mht, so    2 x 2 x 2
 1 : One mht, so         3 x 3 x 3
 2 : Two mht             5 x 5 x 5
 -1 : ignore, use only n_max_distance_mm
 (If both the above are used, the tigher one prevails.)
 Outputs:
 --------
 int *vtxnum,              Closest vertex number
 double *vtx_distance      Distance to closest vertex

 Design Notes:
 -------------
 1. I've incorporated ability to search near voxels first, and bail early if
 remaining voxels are further away.

 2. Elaborateness of calcs for measuring whether other voxels are still useful
 to search could be increased, but the tradeoffs for that are quite different
 depending on whether voxel size (mht->vres()) is large or small.

 3. Will not return a vertex that has been found, but which is outside the
 distance that the inspected voxels cover exhaustively.
 Eg: If in_max_mhts = 0 that means inspect 2 x 2 x 2, and covers all vertices
 within 0.5 mhtres. If closest vertex found is at say 0.75 then function
 returns no result because there could be a closer vertex in the 3 x 3 x 3 box.

 4. Version to return *list* of voxels within range is not implemented,
 because I could find no code that used that.
 -----------------------------------------------------------------
*/

template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::findClosestVertexGeneric(
                                //---------- inputs --------------
                                double probex,
                                double probey,
                                double probez,
                                // How far to search: set one or both
                                double in_max_distance_mm, /* Use large number
                                                              to ignore */
                                int in_max_mhts,           /* Use -1 to ignore */
                                //---------- outputs -------------
                                int *vtxnum,
                                double *vtx_distance)
{  
  const int max_mhts_MAX = 5;
  double mhtres, max_distance_mm, tempdbl;
  int max_mhts;
  double probex_vol, probey_vol, probez_vol;     // probex etc translated to volume space
  int probex_vox, probey_vox, probez_vox;        // probex_vol etc to voxel index
  double probex_mod, probey_mod, probez_mod;     // probex_vol remainder (posn of probe within voxel)
  int near8offsetx, near8offsety, near8offsetz;  // probe?_mod to -1..0
  int xv, yv, zv, xvi, yvi, zvi;                 // voxel indices
  int MinDistVtxNum;
  double MinDistSq, MinDistTemp;
  double RemainingVoxelDistance;
  int RVox, WallJump;
  bool isWall;
  unsigned char central27[3][3][3];  // Indexes 0..2 stand for -1..+1
  //----------------------------------

  mhtres = vres();

  //--------------------------------------------------
  // Initialize instrumentation
  //--------------------------------------------------
  FindBucketsChecked_Count = 0;
  FindBucketsPresent_Count = 0;
  VertexNumFoundByMHT = -1;  // -1 = "none found"  2007-07-30 GW

  //--------------------------------------------------
  // Figure how far afield to search
  //--------------------------------------------------
  if (-1 == in_max_mhts) {  // use in_max_distance_mm
    max_distance_mm = in_max_distance_mm;
    if ((max_distance_mm * 2) <= mhtres) {
      max_mhts = 0;
    }
    else {
      max_mhts = ceil(max_distance_mm / mhtres);
    }
  }
  else {
    max_mhts = in_max_mhts;
    max_distance_mm = in_max_distance_mm;

    // How far does max_mhts cover in mm?
    if (max_mhts >= 1) {
      tempdbl = max_mhts * mhtres;
    }
    else {
      tempdbl = 0.5 * mhtres;
    }
    // Must not include points beyond the exhaustive coverage distance
    // of the chosen max_mhts....
    if (max_distance_mm > tempdbl) max_distance_mm = tempdbl;
  }
  // Safety limit
  if (max_mhts > max_mhts_MAX) max_mhts = max_mhts_MAX;

  // printf("\nmax_distance_mm=%f\n",max_distance_mm);

  //--------------------------------------------------
  // Initialize mins
  //--------------------------------------------------
  MinDistSq = 1e6;
  MinDistVtxNum = -1;

  //--------------------------------------------------
  // Translate probe point to voxel-space coord and indexes
  //--------------------------------------------------
  probex_vol = WORLD_TO_VOLUME( probex);
  probey_vol = WORLD_TO_VOLUME( probey);
  probez_vol = WORLD_TO_VOLUME( probez);

  // (Note: In following (int) truncs toward zero, but that's OK because
  // range of probex_vol is all positive, centered at TABLE_CENTER)
  probex_vox = (int)probex_vol;
  probey_vox = (int)probey_vol;
  probez_vox = (int)probez_vol;

  probex_mod = probex_vol - (double)probex_vox;
  probey_mod = probey_vol - (double)probey_vox;
  probez_mod = probez_vol - (double)probez_vox;

  near8offsetx = (probex_mod <= 0.5) ? -1 : 0;
  near8offsety = (probey_mod <= 0.5) ? -1 : 0;
  near8offsetz = (probez_mod <= 0.5) ? -1 : 0;

  //--------------------------------------------------
  // Initialize checklist for central 27 voxels
  //--------------------------------------------------
  memset(central27, 0, sizeof(central27));
  /*
    for (    xv = 0; xv <= 2; xv++)
    for (  yv = 0; yv <= 2; yv++)
    for (zv = 0; zv <= 2; zv++)
    central27[xv][yv][zv] = 0;
  */

  //--------------------------------------------------
  // Look in home vertex and closest 7 voxels. Without measuring
  // exact probe?_mod position, these "closest 7" could have vertices as close
  // as zero mm, and they completely cover the region out to 0.5 mht voxels.
  // Note: Could possibly recode to abort early based on exact measurements,
  // which might be worthwhile for large voxel.
  //--------------------------------------------------
  for (xvi = 0; xvi <= 1; xvi++) {
    xv = xvi + near8offsetx;
    for (yvi = 0; yvi <= 1; yvi++) {
      yv = yvi + near8offsety;
      for (zvi = 0; zvi <= 1; zvi++) {
        zv = zvi + near8offsetz;

        mhtfindClosestVertexGenericInBucket(probex_vox + xv,
                                            probey_vox + yv,
                                            probez_vox + zv,
                                            probex,
                                            probey,
                                            probez,
                                            &MinDistVtxNum,
                                            &MinDistSq);

        central27[xv + 1][yv + 1][zv + 1] = 1;
      }
    }
  }

  if (max_mhts == 0) goto done;  // stop if caller restricts us to "nearest 8", regardless of whether vertex found

  RemainingVoxelDistance = 0.5 * mhtres;                     // all other voxels contain space at least this far away
  if (max_distance_mm <= RemainingVoxelDistance) goto done;  // Stop if caller restricts us to less than this

  //---------------------------------------------------------------------------
  // We can stop now if found vertex's distance is < 0.5 mhtres, because all
  // other voxels are at at least that far away
  //---------------------------------------------------------------------------

  if (MinDistVtxNum >= 0)                       // if a vertex was found...
  {                                             // not NULL if one has been found)
    MinDistTemp = sqrt(MinDistSq);              // take sqrt
    if (MinDistTemp <= RemainingVoxelDistance)  // if less than all remaining space, we can stop
      goto done;
  }

  //--------------------------------------------------
  // Continue with rest of central 27
  //--------------------------------------------------
  for (xv = -1; xv <= 1; xv++) {
    for (yv = -1; yv <= 1; yv++) {
      for (zv = -1; zv <= 1; zv++) {
        if (!central27[xv + 1][yv + 1][zv + 1])  // skip ones already done
        {
          mhtfindClosestVertexGenericInBucket(probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex,
                                              probey,
                                              probez,
                                              &MinDistVtxNum,
                                              &MinDistSq);

          central27[xv + 1][yv + 1][zv + 1] = 1;
        }
      }
    }
  }

  //--------------------------------------------------
  // Continue with further-away voxels
  // "RVox" roughly "radius in voxel units"
  //  max_mhts     RVox      box dims
  //      0        prev      2 x 2 x 2
  //      1        prev      3 x 3 x 3
  //      2         2        5 x 5 x 5
  //      3         3        7 x 7 x 7
  //      4         4        9 x 9 x 9 etc -- obviously getting slow!
  //--------------------------------------------------

  for (RVox = 2; RVox <= max_mhts; RVox++) {
    //----------------------------------------------------------------------
    // We can stop now if found vertex's distance is<(1.0 x (RVox-1) x mhtres),
    // because all other voxels are at at least that far away
    //----------------------------------------------------------------------
    RemainingVoxelDistance = (mhtres * (RVox - 1));

    if (MinDistVtxNum >= 0) {  // not NULL if one has been found)
      MinDistTemp = sqrt(MinDistSq);
      if (MinDistTemp <= RemainingVoxelDistance) goto done;
    }
    if (max_distance_mm <= RemainingVoxelDistance) goto done;

    //-------------------------------------------------
    // Inspect "shell" of voxels
    //-------------------------------------------------
    WallJump = RVox + RVox;  // jump from one side to the other across the empty middle
    for (xv = -RVox; xv <= RVox; xv++) {
      for (yv = -RVox; yv <= RVox; yv++) {
        isWall = ((xv == -RVox) || (xv == RVox)) || ((yv == -RVox) || (yv == RVox));
        for (zv = -RVox; zv <= RVox; zv = isWall ? zv + 1 : zv + WallJump) {
          mhtfindClosestVertexGenericInBucket(probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex,
                                              probey,
                                              probez,
                                              &MinDistVtxNum,
                                              &MinDistSq);
        }  // zv
      }    // yv
    }      // xv
  }        // RVox

done:
  MinDistTemp = 1e3;

  //--------------------------------------------
  // Enforce not returning a vertex if it's outside
  // max_distance_mm.
  //--------------------------------------------
  if (MinDistVtxNum >= 0) {
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp > max_distance_mm) {  // Legit vertex not found, so set "not-found" values
      MinDistVtxNum = -1;
      MinDistTemp = 1e3;
    }
  }

  // Copy to output
  if (vtxnum) *vtxnum = MinDistVtxNum;
  if (vtx_distance) *vtx_distance = MinDistTemp;

  return NO_ERROR;
}


/*
 ----------------------------------------------------------------
 findClosestFaceGeneric
 Generic find routine satisfying all other finds
 Inputs:
 -------
 probex, probey, probez: The point relative to which to search

 n_max_distance_mm: Furthest vertex distance to accept. Set to a large
 number to ignore, use only in_max_mhts.
 in_max_mhts       : Furthest voxels to search
 0 : One half mht, so    2 x 2 x 2
 1 : One mht, so         3 x 3 x 3
 2 : Two mht             5 x 5 x 5
 -1 : ignore, use only n_max_distance_mm
 (If both the above are used, the tigher one prevails.)
 Outputs:
 --------
 int *facenum,              Closest face number -1 if none
 double *face_distance      Distance to closest face

 Design Notes:
 -------------
 1. I've incorporated ability to search near voxels first, and bail early if
 remaining voxels are further away.

 2. Elaborateness of calcs for measuring whether other voxels are still useful
 to search could be increased, but the tradeoffs for that are quite different
 depending on whether voxel size (mht->vres()) is large or small.

 3. Will not return a vertex that has been found, but which is outside the
 distance that the inspected voxels cover exhaustively.
 Eg: If in_max_mhts = 0 that means inspect 2 x 2 x 2, and covers all vertices
 within 0.5 mhtres. If closest vertex found is at say 0.75 then function
 returns no result because there could be a closer vertex in the 3 x 3 x 3 box.

 4. Version to return *list* of voxels within range is not implemented,
 because I could find no code that used that.
 -----------------------------------------------------------------
*/

template <class Surface, class Face, class Vertex>
void MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::findClosestFaceNoGeneric(
    //---------- inputs --------------
    double probex,
    double probey,
    double probez,
    // How far to search: set one or both
    double in_max_distance_mm, /* Use large number
                                to ignore */
    int in_max_mhts,           /* Use -1 to ignore */
    int project_into_face,
    //---------- outputs -------------
    int *pfno,
    double *pface_distance)
{
  //  const int max_mhts_MAX = 5;
  double mhtres, max_distance_mm, tempdbl;
  int max_mhts;
  double probex_vol, probey_vol, probez_vol;     // probex etc translated to volume space
  int probex_vox, probey_vox, probez_vox;        // probex_vol etc to voxel index
  double probex_mod, probey_mod, probez_mod;     // probex_vol remainder (posn of probe within voxel)
  int near8offsetx, near8offsety, near8offsetz;  // probe?_mod to -1..0
  int xv, yv, zv, xvi, yvi, zvi;                 // voxel indices
  int MinDistFaceNum;
  double MinDistSq, MinDistTemp;
  double RemainingVoxelDistance;
  int RVox, WallJump;
  bool isWall;
  unsigned char central27[3][3][3];  // Indexes 0..2 stand for -1..+1
  //----------------------------------

  mhtres = vres();

  //--------------------------------------------------
  // Initialize instrumentation
  //--------------------------------------------------
  FindBucketsChecked_Count = 0;
  FindBucketsPresent_Count = 0;
  VertexNumFoundByMHT = -1;  // -1 = "none found"  2007-07-30 GW

  //--------------------------------------------------
  // Figure how far afield to search
  //--------------------------------------------------
  if (-1 == in_max_mhts) {  // use in_max_distance_mm
    max_distance_mm = in_max_distance_mm;
    if ((max_distance_mm * 2) <= mhtres) {
      max_mhts = 0;
    }
    else {
      max_mhts = ceil(max_distance_mm / mhtres);
    }
  }
  else {
    max_mhts = in_max_mhts;
    max_distance_mm = in_max_distance_mm;

    // How far does max_mhts cover in mm?
    if (max_mhts >= 1) {
      tempdbl = max_mhts * mhtres;
    }
    else {
      tempdbl = 0.5 * mhtres;
    }
    // Must not include points beyond the exhaustive coverage distance
    // of the chosen max_mhts....
    if (max_distance_mm > tempdbl) max_distance_mm = tempdbl;
  }
#if 0  // disabled by BRF
  // Safety limit
  if (max_mhts > max_mhts_MAX)
    max_mhts = max_mhts_MAX;
#endif

  // printf("\nmax_distance_mm=%f\n",max_distance_mm);

  //--------------------------------------------------
  // Initialize mins
  //--------------------------------------------------
  MinDistSq = 1e6;
  MinDistFaceNum = -1;

  //--------------------------------------------------
  // Translate probe point to voxel-space coord and indexes
  //--------------------------------------------------
  probex_vol = WORLD_TO_VOLUME( probex);
  probey_vol = WORLD_TO_VOLUME( probey);
  probez_vol = WORLD_TO_VOLUME( probez);

  // (Note: In following (int) truncs toward zero, but that's OK because
  // range of probex_vol is all positive, centered at TABLE_CENTER)
  probex_vox = (int)probex_vol;
  probey_vox = (int)probey_vol;
  probez_vox = (int)probez_vol;

  probex_mod = probex_vol - (double)probex_vox;
  probey_mod = probey_vol - (double)probey_vox;
  probez_mod = probez_vol - (double)probez_vox;

  near8offsetx = (probex_mod <= 0.5) ? -1 : 0;
  near8offsety = (probey_mod <= 0.5) ? -1 : 0;
  near8offsetz = (probez_mod <= 0.5) ? -1 : 0;

  //--------------------------------------------------
  // Initialize checklist for central 27 voxels
  //--------------------------------------------------
  memset(central27, 0, sizeof(central27));

  //--------------------------------------------------
  // Look in home vertex and closest 7 voxels. Without measuring
  // exact probe?_mod position, these "closest 7" could have vertices as close
  // as zero mm, and they completely cover the region out to 0.5 mht voxels.
  // Note: Could possibly recode to abort early based on exact measurements,
  // which might be worthwhile for large voxel.
  //--------------------------------------------------
  for (xvi = 0; xvi <= 1; xvi++) {
    xv = xvi + near8offsetx;
    for (yvi = 0; yvi <= 1; yvi++) {
      yv = yvi + near8offsety;
      for (zvi = 0; zvi <= 1; zvi++) {
        zv = zvi + near8offsetz;

        mhtfindClosestFaceCentroidGenericInBucket(probex_vox + xv,
                                                  probey_vox + yv,
                                                  probez_vox + zv,
                                                  probex,
                                                  probey,
                                                  probez,
                                                  project_into_face,
                                                  &MinDistFaceNum,
                                                  &MinDistSq);

        central27[xv + 1][yv + 1][zv + 1] = 1;
      }
    }
  }

  if (max_mhts == 0) goto done;  // stop if caller restricts us to "nearest 8", regardless of whether face was found

  RemainingVoxelDistance = 0.5 * mhtres;                     // all other voxels contain space at least this far away
  if (max_distance_mm <= RemainingVoxelDistance) goto done;  // Stop if caller restricts us to less than this

  //---------------------------------------------------------------------------
  // We can stop now if found vertex's distance is < 0.5 mhtres, because all
  // other voxels are at at least that far away
  //---------------------------------------------------------------------------

  if (MinDistFaceNum >= 0)                      // if a face was found...
  {                                             // not NULL if one has been found)
    MinDistTemp = sqrt(MinDistSq);              // take sqrt
    if (MinDistTemp <= RemainingVoxelDistance)  // if less than all remaining space, we can stop
      goto done;
  }

  //--------------------------------------------------
  // Continue with rest of central 27
  //--------------------------------------------------
  for (xv = -1; xv <= 1; xv++) {
    for (yv = -1; yv <= 1; yv++) {
      for (zv = -1; zv <= 1; zv++) {
        if (!central27[xv + 1][yv + 1][zv + 1])  // skip ones already done
        {
          mhtfindClosestFaceCentroidGenericInBucket(probex_vox + xv,
                                                    probey_vox + yv,
                                                    probez_vox + zv,
                                                    probex,
                                                    probey,
                                                    probez,
                                                    project_into_face,
                                                    &MinDistFaceNum,
                                                    &MinDistSq);

          central27[xv + 1][yv + 1][zv + 1] = 1;
        }
      }
    }
  }

  //--------------------------------------------------
  // Continue with further-away voxels
  // "RVox" roughly "radius in voxel units"
  //  max_mhts     RVox      box dims
  //      0        prev      2 x 2 x 2
  //      1        prev      3 x 3 x 3
  //      2         2        5 x 5 x 5
  //      3         3        7 x 7 x 7
  //      4         4        9 x 9 x 9 etc -- obviously getting slow!
  //--------------------------------------------------

  for (RVox = 2; RVox <= max_mhts; RVox++) {
    //----------------------------------------------------------------------
    // We can stop now if found vertex's distance is<(1.0 x (RVox-1) x mhtres),
    // because all other voxels are at at least that far away
    //----------------------------------------------------------------------
    RemainingVoxelDistance = (mhtres * (RVox - 1));

    if (MinDistFaceNum >= 0) {
      MinDistTemp = sqrt(MinDistSq);
      if (MinDistTemp <= RemainingVoxelDistance) goto done;
    }
    if (max_distance_mm <= RemainingVoxelDistance) goto done;

    //-------------------------------------------------
    // Inspect "shell" of voxels
    //-------------------------------------------------
    WallJump = RVox + RVox;  // jump from one side to the other across the empty middle
    for (xv = -RVox; xv <= RVox; xv++) {
      for (yv = -RVox; yv <= RVox; yv++) {
        isWall = ((xv == -RVox) || (xv == RVox)) || ((yv == -RVox) || (yv == RVox));
        for (zv = -RVox; zv <= RVox; zv = isWall ? zv + 1 : zv + WallJump) {
          mhtfindClosestFaceCentroidGenericInBucket(probex_vox + xv,
                                                    probey_vox + yv,
                                                    probez_vox + zv,
                                                    probex,
                                                    probey,
                                                    probez,
                                                    project_into_face,
                                                    &MinDistFaceNum,
                                                    &MinDistSq);
        }  // zv
      }    // yv
    }      // xv
  }        // RVox

done:
  MinDistTemp = 1e3;

  //--------------------------------------------
  // Enforce not returning a vertex if it's outside
  // max_distance_mm.
  //--------------------------------------------
  if (MinDistFaceNum >= 0) {
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp > max_distance_mm) {  // Legit vertex not found, so set "not-found" values
      MinDistFaceNum = -1;
      MinDistTemp = 1e3;
    }
  }

  // Copy to output
  if (pfno) *pfno = MinDistFaceNum;
  if (pface_distance) *pface_distance = MinDistTemp;
}


/*---------------------------------------------------------------------
  findClosestSetVertexNo

  Returns vertex in mht and mris that is closest to which-selected ###xyz.

  -----------------------------------------------------------------------*/
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::findClosestSetVertexNo(float x, float y, float z)
{
  int vno;
  findClosestVertexGeneric(   x,
                              y,
                              z,
                              1000,
                              3,  // max_mhts: search out to 7 x 7 x 7
                              &vno,
                              NULL);

  if (vno < 0) {  // did not find a vertex, so use brute-force
    vno = mhtBruteForceClosestVertex(x, y, z, which(), NULL);
  }

  return vno;
}


/*--------------------------------------------------------------------
  MHTfindClosestVertexNo()
  Returns vertex number and distance from vertex v to closest vertex
  in mris & mht.
  --------------------------------------------------------------------*/
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::findClosestVertexNoXYZ(
                               float x, float y, float z, 
                               float *min_dist) {

  double min_dist_dbl;
  int vno, rslt;
  rslt = findClosestVertexGeneric(   x,
                                     y,
                                     z,
                                     1000,
                                     1,  // max_mhts: search out to 3 x 3 x 3
                                     &vno,
                                     &min_dist_dbl);

  *min_dist = min_dist_dbl;
  //----------------------------------------------------------
  // [GW] Fixup for "no vertex found". V1.27 function returns
  // zero for no vertex found. I think this is a bug because
  // there is a legit vertex zero I think.
  // However, for now just duplicate existing output
  //----------------------------------------------------------
  if (-1 == rslt) {
    // brf - should always find one that is closest
    vno = mhtBruteForceClosestVertex(x, y, z, which(), min_dist);
  }

  return vno;
}

/*---------------------------------------------------------------
  findVnoOfClosestVertexInTable
  Returns vertex from mris & mht that's closest to provided coordinates.
  ---------------------------------------------------------------*/
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::findVnoOfClosestVertexInTable(
    float x, float y, float z, int do_global_search)
{
  int vno;
  findClosestVertexGeneric(   x,
                              y,
                              z,
                              1000,
                              1,  // max_mhts: search out to 7x7x7 (was 1, BRF)
                              &vno,
                              NULL);
                              
  if ((vno < 0) && do_global_search) { // do more local search first
    for (int i = 2; i <= 4; i++) {
      findClosestVertexGeneric(   x,
                                  y,
                                  z,
                                  1000,
                                  i,  // max_mhts: search out to 7x7x7 (was 1, BRF)
                                  &vno,
                                  NULL);
      if (vno >= 0)  // found it
        break;
    }

    if (vno < 0) { // did not find a vertex, so use brute-force (BRF)
      vno = mhtBruteForceClosestVertex(x, y, z, which(), NULL);
    }
  }

  return vno;
}

#define MAX_VERTICES 50000
template <class Surface, class Face, class Vertex>
int MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtBruteForceClosestVertex(
    float   x,
    float   y,
    float   z,
    int     which,  // which surface within mris to search
    float * dmin)
{
    int   min_vno = -1;
    float min_dsq = 1e8;

    for (int vno = 0; vno < surface.nvertices(); vno++) {
        auto vtx = surface.vertices(vno);
        if (vtx.ripflag()) continue;
            
        float tryx, tryy, tryz;
        mhtVertex2xyz(vtx, which, &tryx, &tryy, &tryz);

        float dx  = tryx - x, dy = tryy - y, dz = tryz - z;
        float dsq = dx * dx + dy * dy + dz * dz;  // squared distance is fine for detecting min

        if (dsq >= min_dsq) continue;
        min_dsq = dsq;
        min_vno = vno;
    }
  
    if (dmin != NULL) *dmin = sqrt(min_dsq);

    return min_vno;
}


// 
//
using namespace Minimal_Surface_MRIS;
//using namespace SurfaceFromMRIS::XYZPositionConsequences;
//using namespace SurfaceFromMRISPV::XYZPositionConsequences;


#define CONSTRUCTORS(ARG,NS) \
MRIS_HASH_TABLE* MHTcreateVertexTable_Resolution(ARG mris, int which, float res)                                    \
{                                                                                                                   \
    NS::Surface surface(mris);                                                                                      \
    auto mht = MRIS_HASH_TABLE_IMPL<NS::Surface,NS::Face,NS::Vertex>::newMHT(MHTFNO_VERTEX, res, which, surface);   \
    mht->captureVertexData();                                                                                       \
    return (mht);                                                                                                   \
}                                                                                                                   \
                                                                                                                    \
MRIS_HASH_TABLE* MHTcreateVertexTable(ARG mris, int which)                                                          \
{                                                                                                                   \
    return MHTcreateVertexTable_Resolution(mris, which, VOXEL_RES);                                                 \
}                                                                                                                   \
                                                                                                                    \
MRIS_HASH_TABLE* MHTcreateFaceTable_Resolution(                                                                     \
    ARG         mris,                                                                                               \
    int         which,                                                                                              \
    float       res)                                                                                                \
{                                                                                                                   \
    NS::Surface surface(mris);                                                                                      \
    auto mht = MRIS_HASH_TABLE_IMPL<NS::Surface,NS::Face,NS::Vertex>::newMHT(MHTFNO_FACE, res, which, surface);     \
    mht->captureFaceData();                                                                                         \
    return mht;                                                                                                     \
}                                                                                                                   \
                                                                                                                    \
MRIS_HASH_TABLE* MHTcreateFaceTable(ARG mris)                                                                       \
{                                                                                                                   \
    return MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, VOXEL_RES);                                        \
}                                                                                                                   \
// end of macro

CONSTRUCTORS(MRIS *,                                               Minimal_Surface_MRIS)
CONSTRUCTORS(Minimal_Surface_MRIS::Surface,                   Minimal_Surface_MRIS)
CONSTRUCTORS(SurfaceFromMRIS::XYZPositionConsequences::Surface,    SurfaceFromMRIS::XYZPositionConsequences)
CONSTRUCTORS(SurfaceFromMRISPV::XYZPositionConsequences::Surface,  SurfaceFromMRISPV::XYZPositionConsequences)

#undef CONSTRUCTORS


int MRIS_HASH_TABLE::BruteForceClosestFace(
    MRIS *mris,
    float   x,
    float   y,
    float   z,
    int     which,  // which surface within mris to search
    float * dmin)
{
    return 
        MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::BruteForceClosestFace(
            Surface(mris),
            x, y, z, which, dmin);
}


//=================================================================
// Diagnostic
//=================================================================

/*--------------------------------------------------------
  MHTtestIsMRISselfIntersecting

  Basically a copy of mrisurf IsMRISselfIntersecting();
  -------------------------------------------------------*/
int MRIS_HASH_TABLE::testIsMRISselfIntersecting(MRIS *mris, float res)
{
    MRIS_HASH_TABLE * mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, res);
    int rslt = 0;
    for (int fno = 0; fno < mris->nfaces; fno++) {
        if (mht->doesFaceIntersect(fno)) {
            rslt = 1;
            goto done;
        }
    }
done:
    MHTfree(&mht);
    return rslt;
}


int MHTBruteForceClosestFace(MRIS* mris,  
                             float x, float y, float z, 
                             int which,                  // which surface within mris to search
                             float *dmin)
{ return MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::BruteForceClosestFace(mris,x,y,z,which,dmin); }


int MHTtestIsMRISselfIntersecting(MRIS* mris,  float res)
{ return MRIS_HASH_TABLE::testIsMRISselfIntersecting(mris, res); }

// Simple properties
//
// Given that the mris is stored in the MHT, it is very unclear why some of the following functions require it to be passed in again!
//
int  MHTwhich(MRIS_HASH_TABLE const * mht) { return mht->which(); }


// Add/remove the faces of which vertex vno is a part
//
int  MHTaddAllFaces   (MRIS_HASH_TABLE* mht, MRIS* mris, int vno) 
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->addAllFaces(vno); }

int  MHTremoveAllFaces(MRIS_HASH_TABLE* mht, MRIS* mris, int vno) 
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->removeAllFaces(vno); }


// Surface self-intersection (Uses MHT initialized with FACES)
//
int MHTdoesFaceIntersect(MRIS_HASH_TABLE* mht, MRIS* mris, int fno) 
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();
  return mht->doesFaceIntersect(fno); }


int MHTisVectorFilled(MRIS_HASH_TABLE const* mht,  int vno, 
                                                float dx, float dy, float dz)
{ return mht->isVectorFilled(vno,dx,dy,dz); }


int MHTfindClosestVertexGeneric(MRIS_HASH_TABLE* mht,
                                double probex, double probey, double probez,
                                double in_max_distance_mm, 
                                int in_max_halfmhts,
                                int *vtxnum, 
                                double *vtx_distance) 
{ return mht->findClosestVertexGeneric(
                                probex, probey, probez,
                                in_max_distance_mm, 
                                in_max_halfmhts,
                                vtxnum, 
                                vtx_distance); }

int MHTfindClosestVertexNoXYZ(MRIS_HASH_TABLE* mht,
                               MRIS* mris, 
                               float x, float y, float z, 
                               float *min_dist) 
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findClosestVertexNoXYZ(x,y,z,min_dist); }

                             
int MHTfindClosestSetVertexNo(MRIS_HASH_TABLE* mht,
                                MRIS* mris,
                                float x, float y, float z) 

{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findClosestSetVertexNo(x,y,z); }


                                            
int MHTfindVnoOfClosestVertexInTable(MRIS_HASH_TABLE* mht,
                                MRIS* mris,
                                float x, float y, float z, int do_global_search) 
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithVertices();
  return mht->findVnoOfClosestVertexInTable(x,y,z,do_global_search); }



// utilities for finding closest face
//
void MHTfindClosestFaceNoGeneric(MRIS_HASH_TABLE* mht,
                              MRIS* mris, 
                              //---------- inputs --------------
                              double probex, double probey, double probez,
                              // How far to search: set one or both
                              double in_max_distance_mm, /* Use large number 
                                                            to ignore */
                              int    in_max_mhts,  /* Use -1 to ignore */
                              // only faces that projection is interior to (Use -1 to ignore )
                              int    project_into_face, 
                              //---------- outputs -------------
                              int *pfno, 
                              double *pface_distance)
{ mht->toMRIS_HASH_TABLE_NoSurface()->checkConstructedWithFaces();    // seen to fail when Vertices
  mht->findClosestFaceNoGeneric(probex, probey, probez,
                              in_max_distance_mm,
                              in_max_mhts,
                              project_into_face,
                              pfno, 
                              pface_distance); }



// This had to be moved to mrisurf_sphere_interp.cpp because something weird was happening 
// with the includes.
//int MHTfindClosestFaceSph(MRIS *surf, MRIS_HASH_TABLE *mht, SphericalInterpolator *si, 
//			  double *cxyz, double *w=NULL, int debug=0);

// These need to be tidied up

int MHTfindClosestVertexNo2(
    MRIS_HASH_TABLE *mht,   // the result will be a vno in the mris used to create this
    MRIS *mris,             // this mris should be the one used to create the mht
    MRIS *mris_for_v,       // this may not be the one for mht nor the one for the result
    VERTEX const *v,        // but this is looked up in it 
    float *pmin_dist)
{
    int vno = v - mris_for_v->vertices;
    cheapAssert(0 <= vno && vno < mris_for_v->nvertices);
    Surface surface(mris_for_v);
    float x,y,z;
    MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtVertex2xyz(surface.vertices(vno), mht->which(), &x, &y, &z);
    return mht->findClosestVertexNoXYZ(x,y,z, pmin_dist);
}

VERTEX * MHTfindClosestVertex2(
    MRIS_HASH_TABLE *mht, 
    MRIS *mris, 
    MRIS *mris_for_v,       // this may not be the one for mht nor the one for the result
    VERTEX const *v,        // but this must be in it 
    float *pmin_dist )
{
    int closest_vno = MHTfindClosestVertexNo2(mht, mris, mris_for_v, v, pmin_dist);
    return (closest_vno < 0) ? nullptr : &mris->vertices[closest_vno];
}


VERTEX* MHTfindClosestVertexSet2(
    MRIS_HASH_TABLE *mht,
    MRIS *mris,
    MRIS *mris_for_v,
    VERTEX const *v)
{
    int vno = v - mris_for_v->vertices;
    cheapAssert(0 <= vno && vno < mris_for_v->nvertices);
    Surface surface(mris_for_v);
    float x,y,z;
    MRIS_HASH_TABLE_IMPL<Surface,Face,Vertex>::mhtVertex2xyz(surface.vertices(vno), mht->which(), &x, &y, &z);
    int closest_vno =  mht->findClosestSetVertexNo(x,y,z);
    return (closest_vno < 0) ? nullptr : &mris->vertices[closest_vno];
}


VERTEX* MHTfindClosestVertexInTable(
    MRIS_HASH_TABLE *mht,
    MRIS *mris,
    float x, float y, float z, int do_global_search)
{
    int closest_vno =  mht->findVnoOfClosestVertexInTable(x, y, z, do_global_search);
    return (closest_vno < 0) ? nullptr : &mris->vertices[closest_vno];
}


void MHTfindClosestFaceGeneric(
    MRIS_HASH_TABLE *mht,
    MRIS *mris, 
    //---------- inputs --------------
    double probex, double probey, double probez,
    // How far to search: set one or both
    double in_max_distance_mm,  /* Use large number to ignore */
    int    in_max_mhts,         /* Use -1 to ignore */
    // only faces that projection is interior to (Use -1 to ignore )
    int    project_into_face, 
    //---------- outputs -------------
    FACE **pface, 
    int *pfno, 
    double *pface_distance)
{
    //checkConstructedWithFacesAndSurface(mris);    // seen to fail when VerticesAndSurface
    mht->findClosestFaceNoGeneric(probex,probey,probez, in_max_distance_mm,in_max_mhts,project_into_face,pfno,pface_distance);
    *pface = (*pfno < 0) ? nullptr : &mris->faces[*pfno];
}
