/**
 * @file  mrishash.c
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which
 * need to compare vertices with one another or to a point.  See:
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2016/06/15 17:49:47 $
 *    $Revision: 1.53 $
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
#include "tritri.h"

#endif

#include "mrishash_internals.h"

//==================================================================
// Local macros
//==================================================================

// GW_VERSION reported by MHTmrishash_gw_version()
// (This is basically a reality check for compile and link)
#define GW_VERSION 126

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
// which, even with vres=1mm, means a triangle approx 100x200 mm!
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


//==================================================================
// Static forward declarations
//==================================================================
// Note: Exposed API functions start with uppercase MHT prefix.
// Static (private) functions start with lowercase mht prefix.

static void mhtComputeFaceCentroid(
                MRIS const *mris, int which, int fno, 
                float *x, float *y, float *z);

static void mhtStoreFaceCentroids(MHT* mht, 
                MRIS const *mris, int which);

static void mhtFaceCentroid2xyz_float(MHT const *mht, 
                int fno, 
                float *x, float *y, float *z);
    // Centroids are computed once and stored in the MHT_FACE
    // They *should* be updated when the face is moved - but I suspect that they are not!

static void mhtVertex2xyz_float    (VERTEX const *vtx, int which, float  *x, float  *y, float  *z);
static void mhtVertex2xyz_double   (VERTEX const *vtx, int which, double *x, double *y, double *z);
static void mhtVertex2Ptxyz_double (VERTEX const *vtx, int which, Ptdbl_t *pt);
static void mhtVertex2array3_double(VERTEX const *vtx, int which, double  *array3);

static int mhtAddFaceOrVertexAtCoords(MRIS_HASH_TABLE *mht, float x, float y, float z, int forvnum);
static int mhtAddFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv, int forvnum);
static int mhtRemoveFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv, int forvnum);

static int mhtFaceToMHT(MRIS_HASH_TABLE *mht, MRIS const *mris, int fno, int on);
static int mhtVoxelList_Init(VOXEL_LISTgw *voxlist);
static int mhtVoxelList_SampleTriangle(
    float mhtres, Ptdbl_t const *vpt0, Ptdbl_t const *vpt1, Ptdbl_t const *vpt2, VOXEL_LISTgw *voxlist);
static int mhtVoxelList_Add(VOXEL_LISTgw *voxlist, int xv, int yv, int zv);
static int mhtVoxelList_AddCoord(VOXEL_LISTgw *voxlist, VOXEL_COORD vc);
static int mhtVoxelList_AddPath (VOXEL_LISTgw *voxlist, VOXEL_COORD oldvc, VOXEL_COORD newvc);

#define MHT_MAX_TOUCHING_FACES 10000
static int MHTexpandToTouchingFaces(
    MRIS_HASH_TABLE const * const mht,
    int   const fno,
    int   const fnoListCapacity,
    int * const fnoList,
    int   const trace);

static int MHTdoesFaceIntersect_old(MRIS_HASH_TABLE const *mht, MRIS const *mris, int fno, int const trace);
static int MHTdoesFaceIntersect_new(MRIS_HASH_TABLE const *mht, MRIS const *mris, int fno, int const trace);
    
static int MHTdoesTriangleIntersect(
    MRIS_HASH_TABLE const * const mht, 
    MHT_TRIANGLE const    * const triangle,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace);

static int mhtDoesTriangleVoxelListIntersect(
    MRIS_HASH_TABLE const * const mht, 
    MHT_TRIANGLE const    * const triangle, 
    VOXEL_LISTgw const    * const voxlistForTriangle,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace);
	
// DELETE THIS
static int mhtDoesFaceVoxelListIntersect(
	MRIS_HASH_TABLE const *mht, 
        MRIS const *mris, VOXEL_LISTgw *voxlist, int fno, int const trace);

int mhtBruteForceClosestVertex(
        MRIS const *mris, float x, float y, float z, int which, float *dmin);

//--------- test -----------
static int checkFace(MRIS_HASH_TABLE *mht, MRIS const *mris, int fno1);


// Primitives that have to be correct...
//
static MRIS_HASH_TABLE* newMHT(MRIS const   *mris)
{
  MRIS_HASH_TABLE* mht = (MRIS_HASH_TABLE *)calloc(1, sizeof(MRIS_HASH_TABLE));
  if (!mht) {
    ErrorExit(ERROR_NO_MEMORY, "%s: could not allocate hash table.\n", __MYFUNCTION__);
  }
  mht->mris = mris;
#ifdef HAVE_OPENMP
  omp_init_lock(&mht->buckets_lock);
#endif
  return mht;
}

static void freeBins(MHBT* bucket);
void MHTfree(MRIS_HASH_TABLE **pmht)
{
  MRIS_HASH_TABLE* mht = *pmht;
  if (!mht) return;

  *pmht = NULL;  // sets pointer to null to signal free'ed

  int xv, yv, zv;
  for (xv = 0; xv < TABLE_SIZE; xv++) {
    for (yv = 0; yv < TABLE_SIZE; yv++) {
      if (!mht->buckets_mustUseAcqRel[xv][yv]) continue;
      for (zv = 0; zv < TABLE_SIZE; zv++) {
        MHBT *bucket = mht->buckets_mustUseAcqRel[xv][yv][zv];
        if (bucket) {
#ifdef HAVE_OPENMP
          omp_destroy_lock(&bucket->bucket_lock);
#endif
          if (bucket->bins) freeBins(bucket);
          free(bucket);
        }
      }
      free(mht->buckets_mustUseAcqRel[xv][yv]);
    }
  }

#ifdef HAVE_OPENMP
  omp_destroy_lock(&mht->buckets_lock);
#endif
  
  free(mht);
}

static void checkThread0() 
{
#ifdef HAVE_OPENMP
  int tid = omp_get_thread_num();
  if (tid != 0) {
    fprintf(stderr, "lock or unlock, not thread 0, but claiming no parallelism\n");
    *(int*)(-1) = 0;
    exit(1);
  }
#endif
}
 
static void lockBuckets(const MRIS_HASH_TABLE *mhtc) {
    MRIS_HASH_TABLE *mht = (MRIS_HASH_TABLE *)mhtc;
#ifdef HAVE_OPENMP
    if (parallelLevel) omp_set_lock(&mht->buckets_lock); else checkThread0();
#endif
}
static void unlockBuckets(const MRIS_HASH_TABLE *mhtc) {
    MRIS_HASH_TABLE *mht = (MRIS_HASH_TABLE *)mhtc;
#ifdef HAVE_OPENMP
    if (parallelLevel) omp_unset_lock(&mht->buckets_lock); else checkThread0();
#endif
}

static void lockBucket(const MHBT *bucketc) {
    MHBT *bucket = (MHBT *)bucketc;
#ifdef HAVE_OPENMP
    if (parallelLevel) omp_set_lock(&bucket->bucket_lock); else checkThread0();
#endif
}
static void unlockBucket(const MHBT *bucketc) {
    MHBT *bucket = (MHBT *)bucketc;
#ifdef HAVE_OPENMP
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

static MHBT* makeAndAcqBucket(MRIS_HASH_TABLE *mht, int xv, int yv, int zv) 
{
  //-----------------------------------------------
  // Allocate space if needed
  //-----------------------------------------------
  // 1. Allocate a 1-D array at buckets_mustUseAcqRel[xv][yv]
  
  lockBuckets(mht);
  
  if (!mht->buckets_mustUseAcqRel[xv][yv]) {
    mht->buckets_mustUseAcqRel[xv][yv] = (MHBT **)calloc(TABLE_SIZE, sizeof(MHBT *));
    if (!mht->buckets_mustUseAcqRel[xv][yv]) ErrorExit(ERROR_NO_MEMORY, "%s: could not allocate slice.", __MYFUNCTION__);
  }
  
  // 2. Allocate a bucket at buckets_mustUseAcqRel[xv][yv][zv]
  MHBT *bucket = mht->buckets_mustUseAcqRel[xv][yv][zv];
  
  
  if (!bucket) {
    mht->buckets_mustUseAcqRel[xv][yv][zv] = bucket = (MHBT *)calloc(1, sizeof(MHBT));
    if (!bucket) ErrorExit(ERROR_NOMEMORY, "%s couldn't allocate bucket.\n", __MYFUNCTION__);
#ifdef HAVE_OPENMP
    omp_init_lock(&bucket->bucket_lock);
#endif
  }
  
  unlockBuckets(mht);
  
  lockBucket(bucket);
  
  // 3. Allocate bins
  if (!bucket->max_bins) /* nothing in this bucket yet - allocate bins */
    reallocBins(bucket, 4);

  // returns with the bucket locked
  
  return bucket;
}

static MHBT *acqBucket(MRIS_HASH_TABLE const *mht, int xv, int yv, int zv)
{
  if (xv >= TABLE_SIZE || yv >= TABLE_SIZE || zv >= TABLE_SIZE || xv < 0 || yv < 0 || zv < 0) return (NULL);
  if (!mht) return (NULL);

  lockBuckets(mht);

  MHBT * bucket = NULL;
  if (mht->buckets_mustUseAcqRel[xv][yv]) 
      bucket = mht->buckets_mustUseAcqRel[xv][yv][zv];

  unlockBuckets(mht);
  
  // returns with the bucket, if any, locked
  //  
  if (bucket) lockBucket(bucket);
  
  return (bucket);
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

void MHTrelBucket (MHBT       ** bucket) { relBucket (bucket); }
void MHTrelBucketC(MHBT const ** bucket) { relBucketC(bucket); }

static bool existsBuckets2(MRIS_HASH_TABLE const *mht, int xv, int yv)
{
  bool result = false;
  if (xv >= TABLE_SIZE || yv >= TABLE_SIZE || xv < 0 || yv < 0) goto Done;
  if (!mht) goto Done;
  if (!mht->buckets_mustUseAcqRel[xv][yv]) goto Done;
  result = true;
Done:
  return result;
}


#define buckets_mustUseAcqRel SHOULD_NOT_ACCESS_BUCKETS_DIRECTLY

//=============================================================================
// Surface --> MHT, store Face Numbers
//=============================================================================

int MHTwhich(MRIS_HASH_TABLE* mht) {
  return mht->which_vertices;
}

MRIS_HASH_TABLE* MHTcreateFaceTable(
    MRIS const   *mris)
{
  return MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, VOXEL_RES);
}

MRIS_HASH_TABLE *MHTcreateFaceTable_Resolution(
    MRIS const *mris, 
    int   which, 
    float res) {

  static int ncalls = 0;
  ncalls++;

  MRIS_HASH_TABLE* mht = newMHT(mris);
  
  mhtStoreFaceCentroids(mht, mris, which);

  //--------------------------------------
  // Capture data from caller and surface
  //--------------------------------------
  mht->vres = res;
  mht->which_vertices = which;
  mht->fno_usage = MHTFNO_FACE;

  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE const* f = &mris->faces[fno];
    if (f->ripflag) continue;
    if (fno == Gdiag_no) DiagBreak();
    mhtFaceToMHT(mht, mris, fno, 1);
  }

  //-------------------------------------------
  // Diagnostics
  //-------------------------------------------
  if ((Gdiag & DIAG_SHOW) && !ncalls) {
    double mean, var, v, n;
    int mx;

    n = mean = 0.0;
    mx = -1;
    int xv,yv,zv;
    for (xv = 0; xv < TABLE_SIZE; xv++) {
      for (yv = 0; yv < TABLE_SIZE; yv++) {
        for (zv = 0; zv < TABLE_SIZE; zv++) {

          MHBT* bucket = acqBucket(mht,xv,yv,zv);
          if (!bucket) continue;
                    
          if (bucket->nused) {
            mean += bucket->nused;
            n++;
          }
          if (bucket->nused > mx) mx = bucket->nused;
          
          relBucket(&bucket);
        }
      }
    }
    mean /= n;
    var = 0.0;
    for (xv = 0; xv < TABLE_SIZE; xv++) {
      for (yv = 0; yv < TABLE_SIZE; yv++) {
        if (!existsBuckets2(mht,xv,yv)) continue;
        for (zv = 0; zv < TABLE_SIZE; zv++) {
          MHBT* bucket = acqBucket(mht,xv,yv,zv);
          if (bucket && bucket->nused) {
            v = mean - bucket->nused;
            var += v * v;
          }
          relBucket(&bucket);
        }
      }
    }
    var /= (n - 1);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%s buckets: mean = %2.1f +- %2.2f, max = %d\n", __MYFUNCTION__, mean, sqrt(var), mx);
  }

  return (mht);
}

/*------------------------------------------------------
  MHTaddAllFaces
  Appends to mht all faces in mris of which VERTEX const v is a part.
  Returns: NO_ERROR unless mht was initialized with vertices rather
  than faces.
  -------------------------------------------------------*/
//------------------------------------
int MHTaddAllFaces(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX_TOPOLOGY const *v)
//------------------------------------
{
  int fi;

  if ((!mht) || (mht->fno_usage != MHTFNO_FACE)) {
    ErrorExit(ERROR_BADPARM, "%s: mht not initialized for faces\n", __MYFUNCTION__);
  }

  for (fi = 0; fi < v->num; fi++) mhtFaceToMHT(mht, mris, v->f[fi], 1);
  return (NO_ERROR);
}

/*------------------------------------------------------
  MHTaddAllFaces
  Appends to mht all faces in mris of which VERTEX const v is a part.
  Returns: NO_ERROR unless mht was initialized with vertices rather
  than faces.
  -------------------------------------------------------*/
//------------------------------------
int MHTremoveAllFaces(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX_TOPOLOGY const *v)
//------------------------------------
{
  int fno;

  if ((!mht) || (mht->fno_usage != MHTFNO_FACE)) {
    ErrorExit(ERROR_BADPARM, "%s: mht not initialized for faces\n", __MYFUNCTION__);
  }

  for (fno = 0; fno < v->num; fno++) mhtFaceToMHT(mht, mris, v->f[fno], 0);
  return (NO_ERROR);
}

/*-------------------------------------------------
  mhtFaceToMHT  (was mhtHatchFace)
  Adds face fno to mht. Calls mhtVoxelList_SampleTriangle to get a list
  of MHT Voxels (buckets) in which to list fno.
  -------------------------------------------------*/
static int mhtFaceToMHT(MRIS_HASH_TABLE *mht, MRIS const *mris, int fno, int on)
{
  //------------------------------------
  FACE const *face;
  VERTEX const *v0, *v1, *v2;
  Ptdbl_t vpt0, vpt1, vpt2;
  VOXEL_LISTgw voxlist;
  int vlix, i, j, k;

  mhtVoxelList_Init(&voxlist);

  face = &mris->faces[fno];
  if (face->ripflag) return (NO_ERROR);

  v0 = &mris->vertices[face->v[0]];
  v1 = &mris->vertices[face->v[1]];
  v2 = &mris->vertices[face->v[2]];

  if (face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no) DiagBreak();
  mhtVertex2Ptxyz_double(v0, mht->which_vertices, &vpt0);
  mhtVertex2Ptxyz_double(v1, mht->which_vertices, &vpt1);
  mhtVertex2Ptxyz_double(v2, mht->which_vertices, &vpt2);

  if (Gx >= 0) {
    double dist0, dist1, dist2;
    dist0 = sqrt(SQR(vpt0.x - Gx) + SQR(vpt0.y - Gy) + SQR(vpt0.z - Gz));
    dist1 = sqrt(SQR(vpt1.x - Gx) + SQR(vpt1.y - Gy) + SQR(vpt1.z - Gz));
    dist2 = sqrt(SQR(vpt2.x - Gx) + SQR(vpt2.y - Gy) + SQR(vpt2.z - Gz));
    if (dist0 < mht->vres || dist1 < mht->vres || dist2 < mht->vres) DiagBreak();
  }
  mhtVoxelList_SampleTriangle(mht->vres, &vpt0, &vpt1, &vpt2, &voxlist);

  for (vlix = 0; vlix < voxlist.nused; vlix++) {
    i = voxlist.voxels[vlix][0];
    j = voxlist.voxels[vlix][1];
    k = voxlist.voxels[vlix][2];

    if (on)
      mhtAddFaceOrVertexAtVoxIx(mht, i, j, k, fno);
    else
      mhtRemoveFaceOrVertexAtVoxIx(mht, i, j, k, fno);
  }
  return (NO_ERROR);
}

/*-------------------------------------------------
  mhtVoxelList_SampleTriangle
  Scans edges and interior of triangle, finding all mht voxels that are
  impinged upon, listing those in voxlist.
  -------------------------------------------------*/
static int mhtVoxelList_SampleTriangle(
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
#define W2VOL(ares, x) (((x) + FIELD_OF_VIEW / 2) / (ares))
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
  return NO_ERROR;
}

//=============================================================================
// Surface --> MHT, store Vertex Numbers
//=============================================================================

//---------------------------------------------------------
MRIS_HASH_TABLE *MHTcreateVertexTable(MRIS const *mris, int which)
//---------------------------------------------------------
{
  return (MHTcreateVertexTable_Resolution(mris, which, VOXEL_RES));
}

//---------------------------------------------------------
MRIS_HASH_TABLE *MHTcreateVertexTable_Resolution(MRIS const *mris, int which, float res)
//---------------------------------------------------------
{
  int vno;
  int xv, yv, zv;
  float x = 0.0, y = 0.0, z = 0.0;
  VERTEX const *v;
  static int ncalls = 0;

  //-----------------------------
  // Allocation and initialization
  //-----------------------------

  MRIS_HASH_TABLE* mht = newMHT(mris);

  mhtStoreFaceCentroids(mht, mris, which);

  //--------------------------------------
  // Capture data from caller and surface
  //--------------------------------------
  mht->vres = res;
  mht->which_vertices = which;
  mht->fno_usage = MHTFNO_VERTEX;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    if (v->ripflag) continue;

    mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);
    mhtAddFaceOrVertexAtCoords(mht, x, y, z, vno);
  }

  //-------------------------------------------
  // Diagnostics
  // Not modified from 1.27
  //-------------------------------------------
  if ((Gdiag & DIAG_SHOW) && !ncalls) {
    double mean, var, v, n;
    int mx;

    n = mean = 0.0;
    mx = -1;
    for (xv = 0; xv < TABLE_SIZE; xv++) {
      for (yv = 0; yv < TABLE_SIZE; yv++) {
        if (!existsBuckets2(mht,xv,yv)) continue;
        for (zv = 0; zv < TABLE_SIZE; zv++) {
          MHBT const *bucket = acqBucket(mht,xv,yv,zv);
          if (!bucket) continue;
          if (bucket->nused) {
            mean += bucket->nused;
            n++;
          }
          if (bucket->nused > mx) mx = bucket->nused;
          relBucketC(&bucket);
        }
      }
    }
    mean /= n;
    var = 0.0;
    for (xv = 0; xv < TABLE_SIZE; xv++) {
      for (yv = 0; yv < TABLE_SIZE; yv++) {
        if (!existsBuckets2(mht,xv,yv)) continue;
        for (zv = 0; zv < TABLE_SIZE; zv++) {
          MHBT const *bucket = acqBucket(mht,xv,yv,zv);
          if (!bucket) continue;
          if (bucket->nused) {
            v = mean - bucket->nused;
            var += v*v;
          }
          relBucketC(&bucket);
        }
      }
    }
    var /= (n - 1);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%s buckets: mean = %2.1f +- %2.2f, max = %d\n", __MYFUNCTION__, mean, sqrt(var), mx);
  }
  ncalls++;
  return (mht);
}

//=================================================================
// AddFaceOrVertexXxxx functions
//=================================================================

/*------------------------------------------------------------
  mhtAddFaceOrVertexAtVoxIx  Was: mhtAddFacePosition
  Adds forvnum (Face or Vertex number) to mht, in bucket(xv,yv,zv)
  Coerces (xv,yv,zv) to be sane.
  -------------------------------------------------------------*/
static int mhtAddFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv, int forvnum)
{
  if (xv < 0) xv = 0;  if (xv >= TABLE_SIZE) xv = TABLE_SIZE - 1;
  if (yv < 0) yv = 0;  if (yv >= TABLE_SIZE) yv = TABLE_SIZE - 1;
  if (zv < 0) zv = 0;  if (zv >= TABLE_SIZE) zv = TABLE_SIZE - 1;

  // [GW] (removed xv,yv,zv check here, as it's impossible to trigger)

  {
    MHBT *bucket = makeAndAcqBucket(mht, xv, yv, zv);

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

/*------------------------------------------------------------
  mhtAddFaceOrVertexAtCoords
  Converts x, y, z to bucket indexes, calls mhtAddFaceOrVertexAtIndexes
  Adds forvnum (Face or Vertex number) to mht, in bucket(xv,yv,zv)
  -------------------------------------------------------------*/
static int mhtAddFaceOrVertexAtCoords(MRIS_HASH_TABLE *mht, float x, float y, float z, int forvnum)
{
  //---------------------------------------
  int xv, yv, zv;

  xv = WORLD_TO_VOXEL(mht, x);
  yv = WORLD_TO_VOXEL(mht, y);
  zv = WORLD_TO_VOXEL(mht, z);

  return mhtAddFaceOrVertexAtVoxIx(mht, xv, yv, zv, forvnum);
}

/*------------------------------------------------------------
  mhtRemoveFaceOrVertexAtVoxIx (was mhtRemoveFacePosition)
  Reverse of mhtAddFaceOrVertexAtIndexes
  -------------------------------------------------------------*/
static int mhtRemoveFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv, int forvnum)
{
  //---------------------------------------
  int i;

  if (xv < 0) xv = 0;
  if (xv >= TABLE_SIZE) xv = TABLE_SIZE - 1;
  if (yv < 0) yv = 0;
  if (yv >= TABLE_SIZE) yv = TABLE_SIZE - 1;
  if (zv < 0) zv = 0;
  if (zv >= TABLE_SIZE) zv = TABLE_SIZE - 1;

  if (!existsBuckets2(mht,xv,yv)) return (NO_ERROR);  // no bucket at such coordinates
  
  MHBT *bucket = acqBucket(mht,xv,yv,zv);
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

//=============================================================================
// Surface self-intersection (Uses MHT initialized with FACES)
//=============================================================================

/*-----------------------------------------------------
  MHTisVectorFilled
  Obscurely named, but means "If vertex vtxno is moved by (dx,dy,dz), will the
  adjoining faces then intersect with some other faces in the surface?".
  .Returns 1=intersect, else 0.
  Notes:
  1. MUST be used with mht prepared for CURRENT_VERTICES
  ------------------------------------------------------*/
int MHTisVectorFilled(
    MRIS_HASH_TABLE const * const mht, 
    MRIS const     * const mris, 
    int   const vtxno, 
    float const dx, 
    float const dy, 
    float const dz)
{
  static int count;
  count++;
  int const trace = 0;
  
  //----------------------------------------------------
  // Sanity check
  //----------------------------------------------------
  if (!mht)  ErrorExit(ERROR_BADPARM, "%s: mht  is NULL\n", __MYFUNCTION__);
  if (!mris) ErrorExit(ERROR_BADPARM, "%s: mris is NULL\n", __MYFUNCTION__);
  if (mht->mris != mris) ErrorExit(ERROR_BADPARM, "%s: mris is wrong\n", __MYFUNCTION__);
  if (mht->fno_usage != MHTFNO_FACE) {
             ErrorExit(ERROR_BADPARM, "%s: mht not initialized for vertices\n", __MYFUNCTION__);
  }
  if (mht->which_vertices != CURRENT_VERTICES) {
             ErrorExit(ERROR_BADPARM, "%s: mht not loaded using CURRENT_VERTICES\n", __MYFUNCTION__);
  }

  // This code moves x,y,z regardless of 'which'
  // and that only makes sense if 'which' is CURRENT_VERTICES
  //
  if (mht->which_vertices != CURRENT_VERTICES) {
    fprintf(stderr, "%s:%d MHTisVectorFilled requires a CURRENT_VERTICES MHT\n", __FILE__, __LINE__);
    exit(1);
  }

  //-------------------------------------------
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
  // 
  //-------------------------------------------

  if (trace) {
    VERTEX_TOPOLOGY const * vtxt = &mris->vertices_topology[vtxno];
    fprintf(stderr, "The faces surrounding vertex %d are ", vtxno);
    int fi;
    for (fi = 0; fi < vtxt->num; fi++) {
      int fno = vtxt->f[fi];
      fprintf(stderr, " %d", fno);
    }
    fprintf(stderr, "\n");
  }

  int result = 0;
  
  VERTEX_TOPOLOGY const * const vtxt = &mris->vertices_topology[vtxno];
  VERTEX          const * const vtx  = &mris->vertices         [vtxno];
  float const moved_x = vtx->x + dx;
  float const moved_y = vtx->y + dy;
  float const moved_z = vtx->z + dz;

  if (trace) {
    fprintf(stderr, "vertex:%d (%g,%g,%g) moved to (%g,%g,%g)\n",
        vtxno, vtx->x, vtx->y, vtx->z, moved_x, moved_y, moved_z);
  }
  
  // Try each changed faced in turn
  //
  int fi;
  for (fi = 0; fi < vtxt->num; fi++) {
    int const fno = vtxt->f[fi];
    FACE const * face = &mris->faces[fno];

    MHT_TRIANGLE triangle;
    int corner;
    for (corner = 0; corner < 3; corner++) {
      Ptdbl_t* point = &triangle.corners[corner];
      int cornerVno = face->v[corner];
      if (cornerVno == vtxno) {
        point->x = moved_x; point->y = moved_y; point->z = moved_z;
      } else {
        VERTEX const * v = &mris->vertices[cornerVno];
        mhtVertex2Ptxyz_double(v, mht->which_vertices, point);
        if (trace) {
          fprintf(stderr, "corner:%d vertex:%d (%g,%g,%g)\n",
            corner, cornerVno, point->x, point->y, point->z);
        }
      }
    }

    // Intersect with the non-changed faces that are not adjacent to this face
    // The changed faces share the vertex so will be in the touchingFnos
    //
    int touchingFnos[MHT_MAX_TOUCHING_FACES];
    int touchingFnosSize = MHTexpandToTouchingFaces(mht, fno, MHT_MAX_TOUCHING_FACES, touchingFnos, trace);
    
    result = MHTdoesTriangleIntersect(mht, &triangle, touchingFnosSize, touchingFnos, trace);
    if (trace) {
      fprintf(stderr, " MHTdoesFaceIntersect_new face:%d returns %d\n", fno, result);
    }
    if (result) break;
  }
  
  return result;
}

/*-------------------------------------------------------------
  MHTdoesFaceIntersect
  Does a particular face intersect with any other faces in surface mris,
  which was previously hashed into mht using MHTcreateFaceTablexxx
  Returns 1 for "yes -- intersection detected", else NO_ERROR.
  -------------------------------------------------------------*/
static int MHTdoesFaceIntersect_old(MRIS_HASH_TABLE const *mht, MRIS const *mris, int fno, int const trace)
{
  //------------------------------------------------------------
  VERTEX const *v0, *v1, *v2;
  Ptdbl_t vpt0, vpt1, vpt2;
  FACE const *face;
  VOXEL_LISTgw voxlist;

  int retval = 0;

  //----------------------------------------------------
  // Sanity check
  //----------------------------------------------------
  if (!mht) ErrorExit(ERROR_BADPARM, "%s: mht is NULL\n", __MYFUNCTION__);
  if (!mris) ErrorExit(ERROR_BADPARM, "%s: mris is NULL\n", __MYFUNCTION__);
  if (mht->fno_usage != MHTFNO_FACE) {
    ErrorExit(ERROR_BADPARM, "%s: mht not initialized for vertices\n", __MYFUNCTION__);
  }

  //--------------------------------
  // Create a voxlist for face fno
  //--------------------------------
  mhtVoxelList_Init(&voxlist);
  face = &mris->faces[fno];

  if (face->ripflag) return (NO_ERROR);

  if (fno == Gdiag_no) DiagBreak();

  v0 = &mris->vertices[face->v[0]];
  v1 = &mris->vertices[face->v[1]];
  v2 = &mris->vertices[face->v[2]];

  mhtVertex2Ptxyz_double(v0, mht->which_vertices, &vpt0);
  mhtVertex2Ptxyz_double(v1, mht->which_vertices, &vpt1);
  mhtVertex2Ptxyz_double(v2, mht->which_vertices, &vpt2);

  mhtVoxelList_SampleTriangle(mht->vres, &vpt0, &vpt1, &vpt2, &voxlist);

  //--------------------------------
  // Any intersections?
  //--------------------------------
  if (mhtDoesFaceVoxelListIntersect(mht, mris, &voxlist, fno, trace)) retval = 1;

  return (retval);
}

static int MHTdoesFaceIntersect_new(MRIS_HASH_TABLE const * const mht, MRIS const * const mris, int const fno, int const trace)
{
  if (!mht ) ErrorExit(ERROR_BADPARM, "%s: mht is NULL\n",  __MYFUNCTION__);
  if (!mris) ErrorExit(ERROR_BADPARM, "%s: mris is NULL\n", __MYFUNCTION__);
  if (mht->fno_usage != MHTFNO_FACE) {
    ErrorExit(ERROR_BADPARM, "%s: mht not initialized for vertices\n", __MYFUNCTION__);
  }
  if (fno == Gdiag_no) DiagBreak();

  FACE const *face = &mris->faces[fno];

  if (face->ripflag) return 0;

  MHT_TRIANGLE triangle;

  int corner;
  for (corner = 0; corner < 3; corner++) {
    mhtVertex2Ptxyz_double(&mris->vertices[face->v[corner]], mht->which_vertices, &triangle.corners[corner]);
  }

  int touchingFnos[MHT_MAX_TOUCHING_FACES];
  int touchingFnosSize = MHTexpandToTouchingFaces(mht, fno, MHT_MAX_TOUCHING_FACES, touchingFnos, trace);
  
  if (trace) {
    int tfi;
    fprintf(stderr, "MHTdoesFaceIntersect_new ignoring ");
    for (tfi = 0; tfi < touchingFnosSize; tfi++) {
      fprintf(stderr, " %d ", touchingFnos[tfi]);
    }
    fprintf(stderr, "\n");
  }
  
  return MHTdoesTriangleIntersect(mht, &triangle, touchingFnosSize, touchingFnos, trace); 
}

int MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRIS const *mris, int fno) {
    static const bool do_old = false;
    static const bool do_new = true;
    static int count;
    count++;
    int old_result = !do_old ? -1 : MHTdoesFaceIntersect_old(mht, mris, fno, count == 0);
    int new_result = !do_new ? -1 : MHTdoesFaceIntersect_new(mht, mris, fno, count == 0);
    if (do_old && do_new && old_result != new_result) {
        fprintf(stderr, "%s:%d MHTdoesFaceIntersect results differ old:%d new:%d count:%d\n", 
            __FILE__, __LINE__,
            old_result, new_result, count);
        exit (1);
    }
    return do_old ? old_result : new_result;
}

static int MHTexpandToTouchingFaces(
    MRIS_HASH_TABLE const * const mht, 
    int   const fno,
    int   const fnoListCapacity,
    int * const fnoList,
    int   const trace)
    //
    // Puts fno in the list, and all the faces that (by intersection rules) are touching fno
    // Returns the length of the list
{
  MRIS const * const mris = mht->mris;    
  FACE const * const face = &mris->faces[fno];

  int size = 0;
  if (size == fnoListCapacity) { fprintf(stderr, "%s:%d exceeded fnoListCapacity\n",__FILE__,__LINE__); exit (1); }
  fnoList[size++] = fno;

  // Consider each vertex of the face
  //
  int vi;
  for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
    int const vno = face->v[vi];
    VERTEX_TOPOLOGY const * const vertex = &mris->vertices_topology[vno];

    #if VERTICES_PER_FACE != 3
    #error assumes 3 vertices per face
    #endif
    int const vno_prev = (0  < vi                 ) ? face->v[vi-1] : face->v[2];
    int const vno_next = (vi < VERTICES_PER_FACE-1) ? face->v[vi+1] : face->v[0];

    // Consider each other face that this is a vertex of
    //
    int fi;
    for (fi = 0; fi < vertex->num; fi++) {
      int const fno2 = vertex->f[fi];
      if (fno2 == fno) continue;
      FACE const * const face2 = &mris->faces[fno2];

      // Is this face already known to touch?
      //
      int fli;
      for (fli = 0; fli < size; fli++) {
        if (fnoList[fli] == fno2) goto already_in_list;
      }
      
      // Note: the previous code only checked for one vertex in common, but such faces CAN intersect...
      // Without this, the code goes into a loop ...
      //
      if (1) { 
        goto touches;
      } else {
        // Does this nearby face share two vertex?
        //
        int vi2;
        for (vi2 = 0; vi2 < VERTICES_PER_FACE; vi2++) {
          int const vno2 = face2->v[vi2];
          if (vno2 == vno_prev || vno2 == vno_next) goto touches;
        } 
      }
      
      // This face does not share an edge therefore can intersect
      //
      continue;

      // This is a new face does share an edge
 touches:     
      if (size == fnoListCapacity) { fprintf(stderr, "%s:%d exceeded fnoListCapacity\n",__FILE__,__LINE__); exit (1); }
      fnoList[size++] = fno2;
      
      if (trace) fprintf(stderr, "MHTexpandToTouchingFaces adding %d\n", fno2);
 already_in_list:;
    } // each other face of the vertex
  } // each vertex of the face
  
  return size;
}

static int MHTdoesTriangleIntersect(
    MRIS_HASH_TABLE const * const mht, 
    MHT_TRIANGLE    const * const triangle,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace) 
{
  int retval = 0;

  // Find all the voxels that the triangle intersects
  //
  VOXEL_LISTgw voxlist;
  mhtVoxelList_Init(&voxlist);
  mhtVoxelList_SampleTriangle(mht->vres, &triangle->corners[0], &triangle->corners[1], &triangle->corners[2], &voxlist);

  //--------------------------------
  // Any intersections?
  //--------------------------------
  if (mhtDoesTriangleVoxelListIntersect(mht, triangle, &voxlist, nFaceToIgnore, fnoToIgnore, trace)) retval = 1;

  return (retval);
}

#define MHT_MAX_FACES 10000
/*-------------------------------------------------------------------
  mhtDoesFaceVoxelListIntersect
  Subsidiary function for MHTdoesFaceIntersect, given face fno already
  analyzed into voxlist.
  ------------------------------------------------------------------*/
static int mhtDoesFaceVoxelListIntersect(
    MRIS_HASH_TABLE const * const mht, 
    MRIS     const * const mris, 
    VOXEL_LISTgw          * const voxlist, 
    int                     const fno, 
    int                     const trace)
//------------------------------------------------------------------
{
  int xv, yv, zv, intersect, voxnum;
  int binix, faceix, facetestno;
  int facelist[MHT_MAX_FACES], nfaces;
  FACE const *facein, *facetest;
  VERTEX const *avtx;
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3];

  //----------------------------------------------------
  // Pertaining to cross-checking vertices of fno versus
  // test face
  //----------------------------------------------------
  int facein_vtx_ix, facetest_vtx_ix;    // indices of vertices in face's list
  int facein_vtx_vno, facetest_vtx_vno;  // Fno of vertices in mris
  VERTEX const *facein_vtx, *facetest_vtx;     // The vertices themselves

  //----------------------------------------------------
  facein = &mris->faces[fno];

  if (fno == Gdiag_no) DiagBreak();

  //------------------------------------------------------------------------
  // Iterate through all the voxels pertaining to face fno (ie:
  // the voxels listed in voxlist), looking at those voxel buckets
  // in mht to find triangles to add to list in flist for later intersection
  // check
  //------------------------------------------------------------------------
  nfaces = 0;
  for (voxnum = 0; voxnum < voxlist->nused; voxnum++) {
    xv = voxlist->voxels[voxnum][0];
    yv = voxlist->voxels[voxnum][1];
    zv = voxlist->voxels[voxnum][2];

    //----------------------------------------------------------
    // Get corresponding bucket from mht. This *should* always
    // succeed, given that the same info was just used to put xv,yv,zv
    // into voxlist as was used to put faces into mht buckets.
    //----------------------------------------------------------
    if (!existsBuckets2(mht,xv,yv)) return (0);

    MHBT *bucket = acqBucket(mht,xv,yv,zv);
    if (!bucket) continue;

    MHB *bin = bucket->bins;

    //-------------------------------------------------
    // Following is one iteration per face
    // (ie: binix is essentially facetestix)
    //-------------------------------------------------
    for (binix = 0; binix < bucket->nused; binix++, bin++) {
      facetestno = bin->fno;  // get the one face from this bin

      const char* reason = NULL;
      
      if (facetestno == fno) { // no point comparing to same face!
        reason = "facetestno == fno";
        goto skip_this_facetest;
      }

      if (trace) {
        fprintf(stderr, "mhtDoesFaceVoxelListIntersect considering face:%d\n", facetestno);
      }
      
      facetest = &mris->faces[facetestno];

      //-------------------------------------------------------------
      // Tests: Several tests to see whether we should skip the
      // actual intersection check for this face, eg because facetest
      // adjoins facein.
      //------------------------------------------------------------
      for (facein_vtx_ix = 0; facein_vtx_ix < VERTICES_PER_FACE; facein_vtx_ix++) {
        facein_vtx_vno = facein->v[facein_vtx_ix];
        facein_vtx = &(mris->vertices[facein_vtx_vno]);

        for (facetest_vtx_ix = 0; facetest_vtx_ix < VERTICES_PER_FACE; facetest_vtx_ix++) {
          facetest_vtx_vno = facetest->v[facetest_vtx_ix];
          facetest_vtx = &(mris->vertices[facetest_vtx_vno]);

          // Do the faces share a vertex?
          if (facein_vtx_vno == facetest_vtx_vno) {
            reason = "facein_vtx_vno == facetest_vtx_vno";
            goto skip_this_facetest;
          }
          
        }      // for facetest_vtx_ix
      }        // for facein_vtx_ix

      //-----------------------------------------------------
      // If we get to here, there was no reason to reject facetest,
      // so append it to faces list
      //-----------------------------------------------------
      for (faceix = 0; faceix < nfaces; faceix++) {
        if (facelist[faceix] == facetestno) goto skip_this_facetest;
      }
      
      
      if (nfaces >= MHT_MAX_FACES) {
        ErrorPrintf(ERROR_NO_MEMORY, "%s: MHT_MAX_FACES exceeded!", __MYFUNCTION__);
        exit (1);
      }
      
      if (trace) {
        fprintf(stderr, "mhtDoesFaceVoxelListIntersect further considering face:%d\n", facetestno);
      }

      facelist[nfaces++] = facetestno;
      reason = NULL;
      
    skip_this_facetest:
      if (trace && reason) fprintf(stderr, "mhtDoesFaceVoxelListIntersect discard face:%d because %s\n", facetestno, reason);
    }                     // for binix

    relBucket(&bucket);
  }                       // for voxnum

  //-----------------------------------------------------------------------
  // Now that we have a list of faces in facelist run the actual geometry
  // intersection test
  // Notes:
  // 1. Code up to v1.27 used CURRENT_VERTICES, which meant that using this
  // code with an mht that was created with MHTcreateFaceTable_Resolution and some
  // other value for "which" would lead to incorrect results.
  // Code below instead calls mhtVertex2array3_double, and selects vertex
  // set based on mht->which_vertices
  //-----------------------------------------------------------------------

  /* set vertices of 1st triangle */
  avtx = &(mris->vertices[facein->v[0]]);
  mhtVertex2array3_double(avtx, mht->which_vertices, v0);

  avtx = &(mris->vertices[facein->v[1]]);
  mhtVertex2array3_double(avtx, mht->which_vertices, v1);

  avtx = &(mris->vertices[facein->v[2]]);
  mhtVertex2array3_double(avtx, mht->which_vertices, v2);

  for (faceix = 0; faceix < nfaces; faceix++) {
    /* set vertices of 2nd triangle */
    facetest = &mris->faces[facelist[faceix]];

    avtx = &(mris->vertices[facetest->v[0]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u0);

    avtx = &(mris->vertices[facetest->v[1]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u1);

    avtx = &(mris->vertices[facetest->v[2]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u2);

    intersect = tri_tri_intersect(v0, v1, v2, u0, u1, u2);

    if (trace && (facelist[faceix] == 1)) {
        fprintf(stderr, "tri_tri_intersect face:%d face:%d returned %d\n", fno, facelist[faceix], intersect);
        fprintf(stderr, "  vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g)\n",
            facein->v[0],v0[0],v0[1],v0[2],
            facein->v[1],v1[0],v1[1],v1[2],
            facein->v[2],v2[0],v2[1],v2[2]);
        fprintf(stderr, "  vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g)\n",
            facetest->v[0],u0[0],u0[1],u0[2],
            facetest->v[1],u1[0],u1[1],u1[2],
            facetest->v[2],u2[0],u2[1],u2[2]);
    }
    
    if (intersect) return (1);
  }
  return (0);
}

/*-------------------------------------------------------------------
  mhtDoesTriangleVoxelListIntersect
  Subsidiary function for MHTdoesTriangleIntersect, 
  triangle already analyzed into voxlist.
  ------------------------------------------------------------------*/
static int mhtDoesTriangleVoxelListIntersect(
    MRIS_HASH_TABLE const * const mht, 
    MHT_TRIANGLE    const * const triangle, 
    VOXEL_LISTgw    const * const voxlist,
    int                     const nFaceToIgnore,
    int const             * const fnoToIgnore,
    int                     const trace)
//------------------------------------------------------------------
{
  MRIS const * const mris = mht->mris;
  
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
    if (!existsBuckets2(mht,xv,yv)) continue;
    
    MHBT const * bucket = acqBucket(mht,xv,yv,zv);
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

  // Find the first face that intersects the triangle
  //
  // Notes:
  // 1. Code up to v1.27 used CURRENT_VERTICES, which meant that using this
  // code with an mht that was created with MHTcreateFaceTable_Resolution and some
  // other value for "which" would lead to incorrect results.
  // Code below instead calls mhtVertex2array3_double, and selects vertex
  // set based on mht->which_vertices
  //-----------------------------------------------------------------------

  double v0[3], v1[3], v2[3];

  // set vertices of 1st triangle
  v0[0] = triangle->corners[0].x;  v0[1] = triangle->corners[0].y;  v0[2] = triangle->corners[0].z;
  v1[0] = triangle->corners[1].x;  v1[1] = triangle->corners[1].y;  v1[2] = triangle->corners[1].z;
  v2[0] = triangle->corners[2].x;  v2[1] = triangle->corners[2].y;  v2[2] = triangle->corners[2].z;
  
  // see if any face intersects
  int fli;
  for (fli = 0; fli < facelistSize; fli++) {
    int const fno = facelist[fli];
    FACE const * const face = &mris->faces[fno];
    
    int const vno0 = face->v[0];
    int const vno1 = face->v[1];
    int const vno2 = face->v[2];
    
    // set vertices of 2nd triangle
    double u0[3], u1[3], u2[3];
    
    mhtVertex2array3_double(&mris->vertices[vno0], mht->which_vertices, u0);
    mhtVertex2array3_double(&mris->vertices[vno1], mht->which_vertices, u1);
    mhtVertex2array3_double(&mris->vertices[vno2], mht->which_vertices, u2);

    int const intersect = tri_tri_intersect(v0, v1, v2, u0, u1, u2);
    
    if (trace && (fno == 1)) {
        fprintf(stderr, "tri_tri_intersect triangle face:%d returned %d\n", fno, intersect);
        fprintf(stderr, "  vertex:unk (%g,%g,%g), vertex:unk (%g,%g,%g), vertex:unk (%g,%g,%g)\n", 
            v0[0],v0[1],v0[2],
            v1[0],v1[1],v1[2],
            v2[0],v2[1],v2[2]);
        fprintf(stderr, "  vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g), vertex:%d (%g,%g,%g)\n",
            vno0,u0[0],u0[1],u0[2],
            vno1,u1[0],u1[1],u1[2],
            vno2,u2[0],u2[1],u2[2]);
    }

    if (intersect) {
      if (trace) {
        fprintf(stderr, "mhtDoesTriangleVoxelListIntersect thinks face:%d intersects\n", fno);
      }
      return (1);
    }
  }

  return (0);
}

//=================================================================
// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//=================================================================

//---------------------------------------------
void mhtFindCommonSanityCheck(MRIS_HASH_TABLE *mht, MRIS const *mris)
{
  //---------------------------------------------
  if (!mht) ErrorExit(ERROR_BADPARM, "%s: mht is NULL\n", __MYFUNCTION__);
  if (!mris) ErrorExit(ERROR_BADPARM, "%s: mris is NULL\n", __MYFUNCTION__);
  if ((!mht) || (mht->fno_usage != MHTFNO_VERTEX)) {
    ErrorExit(ERROR_BADPARM, "%s: mht not initialized for vertices\n", __MYFUNCTION__);
  }
}

//------------------------------------------
// Simple instrumentation
//------------------------------------------
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

/*----------------------------------------------------------------
  mhtfindClosestVertexGenericInBucket

  2007-07-27 GW: Prior to this, ADistSq was based on CURRENT_VERTICES,
  when it should have been based on the vertices specified by
  mht->which_vertices. Search for [2007-07-27 GW]

  -----------------------------------------------------------------*/
int mhtfindClosestVertexGenericInBucket(MRIS_HASH_TABLE *mht,
                                        MRIS const *mris,
                                        //---------- inputs --------------
                                        int xv,
                                        int yv,
                                        int zv,
                                        double probex,
                                        double probey,
                                        double probez,
                                        //---------- in/outs -------------
                                        VERTEX **MinDistVtx,
                                        int *MinDistVtxNum,
                                        double *MinDistSq)
{
  int vtxix, rslt;
  VERTEX *AVtx;
  int AVtxNum;
  double ADistSq;
  float tryx = 0.0, tryy = 0.0, tryz = 0.0;  // Added [2007-07-27 GW]

  //----------------------------------
  rslt = NO_ERROR;

  FindBucketsChecked_Count++;

  MHB *bin;
  MHBT* bucket = MHTacqBucketAtVoxIx(mht, xv, yv, zv);
  if (!bucket) goto done;

  FindBucketsPresent_Count++;

  //-----------------------------------------
  // Iterate through vertices in this bucket
  //-----------------------------------------
  bin = bucket->bins;
  for (vtxix = 0; vtxix < bucket->nused; vtxix++, bin++) {
    AVtxNum = bin->fno;

    AVtx = &mris->vertices[AVtxNum];
    if (AVtx->ripflag)
      continue ;

    if (AVtxNum == Gdiag_no) DiagBreak();

    //    Replaced below [2007-07-27 GW]
    //    ADistSq = SQR(AVtx->x - probex)
    //      + SQR(AVtx->y - probey)
    //      + SQR(AVtx->z - probez) ;

    //----- New [2007-07-27 GW] -----
    mhtVertex2xyz_float(AVtx, mht->which_vertices, &tryx, &tryy, &tryz);  // Added [2007-07-27 GW]

    ADistSq = SQR(tryx - probex) + SQR(tryy - probey) + SQR(tryz - probez);
    //----- end new -----

    if (ADistSq < *MinDistSq) {
      *MinDistSq = ADistSq;
      *MinDistVtxNum = AVtxNum;
      *MinDistVtx = AVtx;
    }  // if
  }    // for vtxix
  
  MHTrelBucket(&bucket);
done:
  return rslt;
}

/*----------------------------------------------------------------
  mhtfindClosestFaceCentroidGenericInBucket

  find the face whose centroid is closest to the specified coordinate
  -----------------------------------------------------------------*/
int mhtfindClosestFaceCentroidGenericInBucket(MRIS_HASH_TABLE *mht,
                                              MRIS const *mris,
                                              //---------- inputs --------------
                                              int xv,
                                              int yv,
                                              int zv,
                                              double probex,
                                              double probey,
                                              double probez,
                                              int project_into_face,
                                              //---------- in/outs -------------
                                              FACE **MinDistFace,
                                              int *MinDistFaceNum,
                                              double *MinDistSq)
{
  int rslt = NO_ERROR;

  FindBucketsChecked_Count++;

  MHB *bin;
  MHBT *bucket = MHTacqBucketAtVoxIx(mht, xv, yv, zv);
  if (!bucket) goto done;

  FindBucketsPresent_Count++;

  //-----------------------------------------
  // Iterate through vertices in this bucket
  //-----------------------------------------
  bin = bucket->bins;
  int faceix;
  for (faceix = 0; faceix < bucket->nused; faceix++, bin++) {

    int const fno = bin->fno;

    if (fno == Gdiag_no) DiagBreak();

    FACE* face = &mris->faces[fno];
    float tryx = 0.0, tryy = 0.0, tryz = 0.0;
    mhtFaceCentroid2xyz_float(mht, fno, &tryx, &tryy, &tryz);

    double lambda[3];
    if (project_into_face > 0 &&
        face_barycentric_coords(
            mris, fno, mht->which_vertices, probex, probey, probez, &lambda[0], &lambda[1], &lambda[2]) < 0)
      continue;
    
    double ADistSq = SQR(tryx - probex) + SQR(tryy - probey) + SQR(tryz - probez);

    if (ADistSq < *MinDistSq) {
      *MinDistSq      = ADistSq;
      *MinDistFaceNum = fno;
      *MinDistFace    = face;
    }  // if
  }    // for faceix

  MHTrelBucket(&bucket);  
done:
  return rslt;
}

/*
 ----------------------------------------------------------------
 MHTfindClosestVertexGeneric
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
 VERTEX **pvtx,            Closest vertex
 int *vtxnum,              Closest vertex number
 double *vtx_distance      Distance to closest vertex

 Design Notes:
 -------------
 1. I've incorporated ability to search near voxels first, and bail early if
 remaining voxels are further away.

 2. Elaborateness of calcs for measuring whether other voxels are still useful
 to search could be increased, but the tradeoffs for that are quite different
 depending on whether voxel size (mht->vres) is large or small.

 3. Will not return a vertex that has been found, but which is outside the
 distance that the inspected voxels cover exhaustively.
 Eg: If in_max_mhts = 0 that means inspect 2 x 2 x 2, and covers all vertices
 within 0.5 mhtres. If closest vertex found is at say 0.75 then function
 returns no result because there could be a closer vertex in the 3 x 3 x 3 box.

 4. Version to return *list* of voxels within range is not implemented,
 because I could find no code that used that.
 -----------------------------------------------------------------
*/

int MHTfindClosestVertexGeneric(MRIS_HASH_TABLE *mht,
                                MRIS const *mris,
                                //---------- inputs --------------
                                double probex,
                                double probey,
                                double probez,
                                // How far to search: set one or both
                                double in_max_distance_mm, /* Use large number
                                                              to ignore */
                                int in_max_mhts,           /* Use -1 to ignore */
                                //---------- outputs -------------
                                VERTEX **pvtx,
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
  VERTEX *MinDistVtx;
  int MinDistVtxNum;
  double MinDistSq, MinDistTemp;
  double RemainingVoxelDistance;
  int RVox, WallJump;
  bool isWall;
  unsigned char central27[3][3][3];  // Indexes 0..2 stand for -1..+1
  //----------------------------------

  mhtFindCommonSanityCheck(mht, mris);
  mhtres = mht->vres;

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
  MinDistVtx = NULL;
  MinDistVtxNum = -1;

  //--------------------------------------------------
  // Translate probe point to voxel-space coord and indexes
  //--------------------------------------------------
  probex_vol = WORLD_TO_VOLUME(mht, probex);
  probey_vol = WORLD_TO_VOLUME(mht, probey);
  probez_vol = WORLD_TO_VOLUME(mht, probez);

  // (Note: In following (int) truncs toward zero, but that's OK because
  // range of probex_vol is all positive, centered at FIELD_OF_VIEW/2)
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

        mhtfindClosestVertexGenericInBucket(mht,
                                            mris,
                                            probex_vox + xv,
                                            probey_vox + yv,
                                            probez_vox + zv,
                                            probex,
                                            probey,
                                            probez,
                                            &MinDistVtx,
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

  if (MinDistVtx)                               // if a vertex was found...
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
          mhtfindClosestVertexGenericInBucket(mht,
                                              mris,
                                              probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex,
                                              probey,
                                              probez,
                                              &MinDistVtx,
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

    if (MinDistVtx) {  // not NULL if one has been found)
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
          mhtfindClosestVertexGenericInBucket(mht,
                                              mris,
                                              probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex,
                                              probey,
                                              probez,
                                              &MinDistVtx,
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
  if (MinDistVtx) {
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp > max_distance_mm) {  // Legit vertex not found, so set "not-found" values
      MinDistVtx = NULL;
      MinDistVtxNum = -1;
      MinDistTemp = 1e3;
    }
  }

  // Copy to output
  if (pvtx) *pvtx = MinDistVtx;
  if (vtxnum) *vtxnum = MinDistVtxNum;
  if (vtx_distance) *vtx_distance = MinDistTemp;

  // 2007-07-30 GW added additional "instrumentation"

  return NO_ERROR;
}
/*
 ----------------------------------------------------------------
 MHTfindClosestFaceGeneric
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
 VERTEX **pface,            Closest face
 int *facenum,              Closest face number
 double *face_distance      Distance to closest face

 Design Notes:
 -------------
 1. I've incorporated ability to search near voxels first, and bail early if
 remaining voxels are further away.

 2. Elaborateness of calcs for measuring whether other voxels are still useful
 to search could be increased, but the tradeoffs for that are quite different
 depending on whether voxel size (mht->vres) is large or small.

 3. Will not return a vertex that has been found, but which is outside the
 distance that the inspected voxels cover exhaustively.
 Eg: If in_max_mhts = 0 that means inspect 2 x 2 x 2, and covers all vertices
 within 0.5 mhtres. If closest vertex found is at say 0.75 then function
 returns no result because there could be a closer vertex in the 3 x 3 x 3 box.

 4. Version to return *list* of voxels within range is not implemented,
 because I could find no code that used that.
 -----------------------------------------------------------------
*/

int MHTfindClosestFaceGeneric(MRIS_HASH_TABLE *mht,
                              MRIS const *mris,
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
                              FACE **pface,
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
  FACE *MinDistFace;
  int MinDistFaceNum;
  double MinDistSq, MinDistTemp;
  double RemainingVoxelDistance;
  int RVox, WallJump;
  bool isWall;
  unsigned char central27[3][3][3];  // Indexes 0..2 stand for -1..+1
  //----------------------------------

  //  mhtFindCommonSanityCheck(mht, mris);
  mhtres = mht->vres;

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
  MinDistFace = NULL;
  MinDistFaceNum = -1;

  //--------------------------------------------------
  // Translate probe point to voxel-space coord and indexes
  //--------------------------------------------------
  probex_vol = WORLD_TO_VOLUME(mht, probex);
  probey_vol = WORLD_TO_VOLUME(mht, probey);
  probez_vol = WORLD_TO_VOLUME(mht, probez);

  // (Note: In following (int) truncs toward zero, but that's OK because
  // range of probex_vol is all positive, centered at FIELD_OF_VIEW/2)
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

        mhtfindClosestFaceCentroidGenericInBucket(mht,
                                                  mris,
                                                  probex_vox + xv,
                                                  probey_vox + yv,
                                                  probez_vox + zv,
                                                  probex,
                                                  probey,
                                                  probez,
                                                  project_into_face,
                                                  &MinDistFace,
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

  if (MinDistFace)                              // if a face was found...
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
          mhtfindClosestFaceCentroidGenericInBucket(mht,
                                                    mris,
                                                    probex_vox + xv,
                                                    probey_vox + yv,
                                                    probez_vox + zv,
                                                    probex,
                                                    probey,
                                                    probez,
                                                    project_into_face,
                                                    &MinDistFace,
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

    if (MinDistFace) {  // not NULL if one has been found)
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
          mhtfindClosestFaceCentroidGenericInBucket(mht,
                                                    mris,
                                                    probex_vox + xv,
                                                    probey_vox + yv,
                                                    probez_vox + zv,
                                                    probex,
                                                    probey,
                                                    probez,
                                                    project_into_face,
                                                    &MinDistFace,
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
  if (MinDistFace) {
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp > max_distance_mm) {  // Legit vertex not found, so set "not-found" values
      MinDistFace = NULL;
      MinDistFaceNum = -1;
      MinDistTemp = 1e3;
    }
  }

  // Copy to output
  if (pface) *pface = MinDistFace;
  if (pfno) *pfno = MinDistFaceNum;
  if (pface_distance) *pface_distance = MinDistTemp;

  return NO_ERROR;
}

/*----------------------------------------------------------------
  MHTfindClosestVertex
  Returns VERTEX const *, closest vertex to v in mris & mht (or NULL).
  -----------------------------------------------------------------*/
VERTEX *MHTfindClosestVertex(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX const *v)
{
  //------------------------------------------------------
  float x = 0.0, y = 0.0, z = 0.0;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);

  VERTEX *vtx;
  MHTfindClosestVertexGeneric(mht,
                              mris,
                              x,
                              y,
                              z,
                              1000,
                              1,  // max_mhts: search out to 3 x 3 x 3
                              &vtx,
                              NULL,
                              NULL);

  return vtx;
}

/*----------------------------------------------------------------
  MHTfindClosestFaceToVertex
  Returns index of face whose centroid is closest to the specified vertex
  -----------------------------------------------------------------*/
FACE const *MHTfindClosestFaceToVertex(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX const *v)
{
  //------------------------------------------------------
  float x = 0.0, y = 0.0, z = 0.0;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);

  FACE *face;
  MHTfindClosestFaceGeneric(mht,
                            mris,
                            x,
                            y,
                            z,
                            1000,
                            1,   // max_mhts: search out to 3 x 3 x 3
                            -1,  // don't force it to project into the face
                            &face,
                            NULL,
                            NULL);

  return face;
}

/*---------------------------------------------------------------------
  MHTfindClosestVertexSet.

  updated by BRF 6/10/2016. I don't think this is called anywhere so I'm going
  to revert to the indended usage and have it use the vertex locations specified
  by the which parameter instead of  v->[xyz].

  Comment revised 2007-07-25:
  Returns vertex in mht and mris that is closest to v->x, v->y, v->z.
  The which parameter is actually ignored. However, it suggests that
  MHTfindClosestVertexSet allows some choice of either the v coordinate
  to use as reference point, or of the surface to be searched.
  However, neither is true: The choice of surface coords is established
  when the surface is hashed into the mht, and in MHTfindClosestVertexSet
  the which parameter is not used to select the coordinates within v.
  Because this all might be misleading, GW's version of this code
  calls ErrorExit when which doesn't match mht->which_vertices, as that
  implies a misunderstanding.

  Code revision 2007-07-25 GW: I had misunderstood the v1.27 code regarding
  which coords of v to use as the reference point. See updated lines below
  marked [GW 2007-07-25]

  -----------------------------------------------------------------------*/
VERTEX *MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX const *v, int which_ignored)
{
  //------------------------------------------------------
  float x = 0.0, y = 0.0, z = 0.0;

  //---------------------------------
  // Sanity checks
  //---------------------------------
  mhtFindCommonSanityCheck(mht, mris);

  if (mht->which_vertices != which_ignored)
    ErrorExit(ERROR_BADPARM, "%s called with mismatched 'which' parameter\n", __MYFUNCTION__);

  mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);

  VERTEX *vtx = NULL;
  MHTfindClosestVertexGeneric(mht,
                              mris,
                              x,
                              y,
                              z,
                              1000,
                              3,  // max_mhts: search out to 7 x 7 x 7
                              &vtx,
                              NULL,
                              NULL);

  if (!vtx)  // did not find a vertex, so use brute-force
  {
    int vnum = mhtBruteForceClosestVertex(mris, x, y, z, mht->which_vertices, NULL);
    vtx = &mris->vertices[vnum];
  }

  return vtx;
}

/*--------------------------------------------------------------------
  MHTfindClosestVertexNo()
  Returns vertex number and distance from vertex v to closest vertex
  in mris & mht.
  --------------------------------------------------------------------*/
int MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX const *v, float *min_dist)
{
  //------------------------------------------------------
  double x = 0.0, y = 0.0, z = 0.0;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_double(v, mht->which_vertices, &x, &y, &z);
  
  return MHTfindClosestVertexNoXYZ(mht, mris, (float)x, (float)y, (float)z, min_dist);
}

int MHTfindClosestVertexNoXYZ(MRIS_HASH_TABLE *mht, 
                               MRIS const *mris, 
                               float x, float y, float z, 
                               float *min_dist) {
  double min_dist_dbl;
  int vtxnum, rslt;
  rslt = MHTfindClosestVertexGeneric(mht,
                                     mris,
                                     x,
                                     y,
                                     z,
                                     1000,
                                     1,  // max_mhts: search out to 3 x 3 x 3
                                     NULL,
                                     &vtxnum,
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
    rslt = mhtBruteForceClosestVertex(mris, x, y, z, mht->which_vertices, min_dist);
  }
  return vtxnum;
}

/*---------------------------------------------------------------
  MHTfindClosestVertexInTable
  Returns vertex from mris & mht that's closest to provided coordinates.
  ---------------------------------------------------------------*/
VERTEX *MHTfindClosestVertexInTable(
    MRIS_HASH_TABLE *mht, MRIS const *mris, float x, float y, float z, int do_global_search)
{
  //------------------------------------------------------
  int i;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  VERTEX *vtx;
  MHTfindClosestVertexGeneric(mht,
                              mris,
                              x,
                              y,
                              z,
                              1000,
                              1,  // max_mhts: search out to 7x7x7 (was 1, BRF)
                              &vtx,
                              NULL,
                              NULL);
  if (!vtx && do_global_search)  // do more local search first
  {
    for (i = 2; i <= 4; i++) {
      MHTfindClosestVertexGeneric(mht,
                                  mris,
                                  x,
                                  y,
                                  z,
                                  1000,
                                  i,  // max_mhts: search out to 7x7x7 (was 1, BRF)
                                  &vtx,
                                  NULL,
                                  NULL);
      if (vtx)  // found it
        break;
    }
    if (!vtx)  // did not find a vertex, so use brute-force (BRF)
    {
      int vnum = mhtBruteForceClosestVertex(mris, x, y, z, mht->which_vertices, NULL);
      vtx = &mris->vertices[vnum];
    }
  }

  return vtx;
}

#define MAX_VERTICES 50000
/*------------------------------------------------------------
  MHTgetAllVerticesWithinDistance
  Returns a list of vertex numbers from mris & mht that are within
  max_dist of vertex vno.

  Allocates separate list pvnum which must be freed by caller.

  As of 2007-03-20, there don't appear to be any callers of this function
  in all the FS source code, so I'm not bothering to implement its
  functionality in the generic Find function, however it would be easy to do.
  ------------------------------------------------------------*/
int *MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht, MRIS const *mris, int vno, float max_dist, int *pvnum)
{
  //------------------------------------------------------
  ErrorExit(ERROR_UNSUPPORTED, "%s: not implemented.\n", __MYFUNCTION__);

  return NULL;
}

/*------------------------------------------------------------
  mhtBruteForceClosestVertex
  Finds closest vertex by exhaustive search of surface. This is
  useful as a fallback for other functions. (As of 2007-07-30,
  it's being introduced only to MHTfindClosestVertexSet so far.)
  ------------------------------------------------------------*/
int mhtBruteForceClosestVertex(MRIS const *mris,
                               float x,
                               float y,
                               float z,
                               int which,  // which surface within mris to search
                               float *dmin)
{
  int vno, min_v = -1;
  VERTEX const *vtx;
  float dsq, min_dsq;  //  Work with squares, avoid square root operation
  float tryx = 0.0, tryy = 0.0, tryz = 0.0, dx, dy, dz;

  min_dsq = 1e8;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vtx = &mris->vertices[vno];
    if (vtx->ripflag) continue;

    //-----------------------------------
    // For maintainability probably better to call mhtVertex2xyz_float,
    // and probably would get inlined anyway,
    // but due to frequency of call, implemented locally.
    //-----------------------------------

    switch (which) {
      case ORIGINAL_VERTICES:
        tryx = vtx->origx;
        tryy = vtx->origy;
        tryz = vtx->origz;
        break;
      case GOOD_VERTICES:
        break;
      case TMP_VERTICES:
        break;
      case CANONICAL_VERTICES:
        tryx = vtx->cx;
        tryy = vtx->cy;
        tryz = vtx->cz;
        break;
      case CURRENT_VERTICES:
        tryx = vtx->x;
        tryy = vtx->y;
        tryz = vtx->z;
        break;
      case INFLATED_VERTICES:
        break;
      case FLATTENED_VERTICES:
        tryx = vtx->fx;
        tryy = vtx->fy;
        tryz = vtx->fz;
        break;
      case PIAL_VERTICES:
        tryx = vtx->pialx;
        tryy = vtx->pialy;
        tryz = vtx->pialz;
        break;
      case TMP2_VERTICES:
        break;
      case WHITE_VERTICES:
        tryx = vtx->whitex;
        tryy = vtx->whitey;
        tryz = vtx->whitez;
        break;
    }

    dx = tryx - x;
    dy = tryy - y;
    dz = tryz - z;

    dsq = dx * dx + dy * dy + dz * dz;  // squared distance is fine for detecting min
    if (dsq < min_dsq) {
      min_dsq = dsq;
      min_v = vno;
    }
  }
  if (dmin != NULL) *dmin = sqrt(min_dsq);

  return (min_v);
}
/*------------------------------------------------------------
  mhtBruteForceClosestFace
  Finds closest vertex by exhaustive search of surface. This is
  useful as a fallback for other functions. (As of 2007-07-30,
  it's being introduced only to MHTfindClosestVertexSet so far.)
  ------------------------------------------------------------*/
int mhtBruteForceClosestFace(MRIS const *mris,
                             float x,
                             float y,
                             float z,
                             int which,  // which surface within mris to search
                             float *dmin)
{
  int   min_fno = -1;
  float min_dsq = 1e8;

  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {

    float tryx, tryy, tryz;
    mhtComputeFaceCentroid(mris, which, fno, &tryx, &tryy, &tryz);

    float const dx = tryx - x;
    float const dy = tryy - y;
    float const dz = tryz - z;

    float dsq = dx * dx + dy * dy + dz * dz;  // squared distance is fine for detecting min
    if (min_dsq > dsq) {
      min_dsq = dsq;
      min_fno = fno;
    }
  }
  
  if (dmin) *dmin = sqrt(min_dsq);

  return (min_fno);
}

static void mhtFaceCentroid2xyz_float(
    MHT const *mht, int fno, 
    float *px, float *py, float *pz)
{
  MHT_FACE const* face = &mht->f[fno];
  *px = face->cx;
  *py = face->cy;
  *pz = face->cz;
}

//---------------------------------------------
static void mhtVertex2xyz_float(VERTEX const *vtx, int which, float *x, float *y, float *z)
{
  //---------------------------------------------
  if (!vtx) return;

  switch (which) {
    case ORIGINAL_VERTICES:
      *x = vtx->origx;
      *y = vtx->origy;
      *z = vtx->origz;
      break;
    case GOOD_VERTICES:
      break;
    case TMP_VERTICES:
      break;
    case CANONICAL_VERTICES:
      *x = vtx->cx;
      *y = vtx->cy;
      *z = vtx->cz;
      break;
    case CURRENT_VERTICES:
      *x = vtx->x;
      *y = vtx->y;
      *z = vtx->z;
      break;
    case INFLATED_VERTICES:
      break;
    case FLATTENED_VERTICES:
      *x = vtx->fx;
      *y = vtx->fy;
      *z = 0;
      break;
    case PIAL_VERTICES:
      *x = vtx->pialx;
      *y = vtx->pialy;
      *z = vtx->pialz;
      break;
    case TMP2_VERTICES:
      break;
    case WHITE_VERTICES:
      *x = vtx->whitex;
      *y = vtx->whitey;
      *z = vtx->whitez;
      break;
  }
  return;
}
//---------------------------------------------
static void mhtVertex2xyz_double(VERTEX const *vtx, int which, double *x, double *y, double *z)
{
  //---------------------------------------------
  if (!vtx) return;

  switch (which) {
    case ORIGINAL_VERTICES:
      *x = vtx->origx;
      *y = vtx->origy;
      *z = vtx->origz;
      break;
    case GOOD_VERTICES:
      break;
    case TMP_VERTICES:
      break;
    case CANONICAL_VERTICES:
      *x = vtx->cx;
      *y = vtx->cy;
      *z = vtx->cz;
      break;
    case CURRENT_VERTICES:
      *x = vtx->x;
      *y = vtx->y;
      *z = vtx->z;
      break;
    case INFLATED_VERTICES:
      break;
    case FLATTENED_VERTICES:
      *x = vtx->fx;
      *y = vtx->fy;
      *z = 0;
      break;
    case PIAL_VERTICES:
      *x = vtx->pialx;
      *y = vtx->pialy;
      *z = vtx->pialz;
      break;
    case TMP2_VERTICES:
      break;
    case WHITE_VERTICES:
      *x = vtx->whitex;
      *y = vtx->whitey;
      *z = vtx->whitez;
      break;
  }
  return;
}
//---------------------------------------------
static void mhtVertex2Ptxyz_double(VERTEX const *vtx, int which, Ptdbl_t *pt)
{
  //---------------------------------------------
  if ((!vtx) || (!pt)) return;

  switch (which) {
    case ORIGINAL_VERTICES:
      pt->x = vtx->origx;
      pt->y = vtx->origy;
      pt->z = vtx->origz;
      break;
    case GOOD_VERTICES:
      break;
    case TMP_VERTICES:
      break;
    case CANONICAL_VERTICES:
      pt->x = vtx->cx;
      pt->y = vtx->cy;
      pt->z = vtx->cz;
      break;
    case CURRENT_VERTICES:
      pt->x = vtx->x;
      pt->y = vtx->y;
      pt->z = vtx->z;
      break;
    case INFLATED_VERTICES:
      break;
    case FLATTENED_VERTICES:
      pt->x = vtx->fx;
      pt->y = vtx->fy;
      pt->z = 0;
      break;
    case PIAL_VERTICES:
      pt->x = vtx->pialx;
      pt->y = vtx->pialy;
      pt->z = vtx->pialz;
      break;
    case TMP2_VERTICES:
      break;
    case WHITE_VERTICES:
      pt->x = vtx->whitex;
      pt->y = vtx->whitey;
      pt->z = vtx->whitez;
      break;
  }
  return;
}

//---------------------------------------------
static void mhtVertex2array3_double(VERTEX const *vtx, int which, double *array3)
{
  //---------------------------------------------
  if (!vtx) return;

  switch (which) {
    case ORIGINAL_VERTICES:
      array3[0] = vtx->origx;
      array3[1] = vtx->origy;
      array3[2] = vtx->origz;
      break;
    case GOOD_VERTICES:
      break;
    case TMP_VERTICES:
      break;
    case CANONICAL_VERTICES:
      array3[0] = vtx->cx;
      array3[1] = vtx->cy;
      array3[2] = vtx->cz;
      break;
    case CURRENT_VERTICES:
      array3[0] = vtx->x;
      array3[1] = vtx->y;
      array3[2] = vtx->z;
      break;
    case INFLATED_VERTICES:
      break;
    case FLATTENED_VERTICES:
      array3[0] = vtx->fx;
      array3[1] = vtx->fy;
      array3[2] = 0;
      break;
    case PIAL_VERTICES:
      array3[0] = vtx->pialx;
      array3[1] = vtx->pialy;
      array3[2] = vtx->pialz;
      break;
    case TMP2_VERTICES:
      break;
    case WHITE_VERTICES:
      array3[0] = vtx->whitex;
      array3[1] = vtx->whitey;
      array3[2] = vtx->whitez;
      break;
  }
  return;
}

/*-----------------------------------------------------------------*/
MHBT *MHTacqBucket(MRIS_HASH_TABLE *mht, float x, float y, float z)
{
  //-------------------------------------------------------------------
  int xv, yv, zv;

  xv = WORLD_TO_VOXEL(mht, x);
  yv = WORLD_TO_VOXEL(mht, y);
  zv = WORLD_TO_VOXEL(mht, z);

  return acqBucket(mht, xv, yv, zv);
}

/*-----------------------------------------------------------------*/
MHBT *MHTacqBucketAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv)
{
  return acqBucket(mht, xv, yv, zv);
}

/*------------------------------------------------
  MH_gw_version
  Confidence check that correct version of code is
  compiled in;
  ------------------------------------------------*/
int MHT_gw_version(void)
{
  //-------------------------------
  return GW_VERSION;  // <-- change this as needed
}

//=================================================================
// VOXEL_LIST
//=================================================================

//--------------------------------------------------
static int mhtVoxelList_Init(VOXEL_LISTgw *voxlist)
{
  //--------------------------------------------------
  voxlist->nused = 0;
  return (NO_ERROR);
}

//--------------------------------------------------
static int mhtVoxelList_Add(VOXEL_LISTgw *voxlist, int xv, int yv, int zv)
{
  //--------------------------------------------------
  int i;

  for (i = 0; i < voxlist->nused; i++)
    if (voxlist->voxels[i][0] == xv && voxlist->voxels[i][1] == yv && voxlist->voxels[i][2] == zv) return (NO_ERROR);

  if (voxlist->nused >= MAX_VOXELS) {
    fprintf(stderr, "%s(%d, %d, %d): complete list too big!\n", __MYFUNCTION__, xv, yv, zv);

    ErrorPrintf(ERROR_NOMEMORY, "%s(%d, %d, %d): complete list too big!", __MYFUNCTION__, xv, yv, zv);
    return (ERROR_NOMEMORY);
  }
  //-----------------------------
  // Actually append it
  //-----------------------------
  i = voxlist->nused;
  voxlist->voxels[i][0] = xv;
  voxlist->voxels[i][1] = yv;
  voxlist->voxels[i][2] = zv;
  voxlist->nused++;
  return (NO_ERROR);
}

//--------------------------------------------------
static int mhtVoxelList_AddCoord(VOXEL_LISTgw *voxlist, VOXEL_COORD vc)
{
  //--------------------------------------------------
  return mhtVoxelList_Add(voxlist, vc.xv, vc.yv, vc.zv);
}

//--------------------------------------------------
static int mhtVoxelList_AddPath(VOXEL_LISTgw *voxlist, VOXEL_COORD oldvc, VOXEL_COORD newvc)
{
  //--------------------------------------------------
  int xv, yv, zv, count;
  int incx = -1, incy = -1, incz = -1;

  if (newvc.xv >= oldvc.xv) incx = 1;
  if (newvc.yv >= oldvc.yv) incy = 1;
  if (newvc.zv >= oldvc.zv) incz = 1;

  count = -1;
  for (xv = oldvc.xv; incx == 1 ? xv <= newvc.xv : xv >= newvc.xv; xv += incx) {
    for (yv = oldvc.yv; incy == 1 ? yv <= newvc.yv : yv >= newvc.yv; yv += incy) {
      for (zv = oldvc.zv; incz == 1 ? zv <= newvc.zv : zv >= newvc.zv; zv += incz) {
        count++;
        // ie: ignore the start voxel, assume that's already added
        // (which can be done with mhtVoxelList_AddCoord)
        if (count >= 1) {
          mhtVoxelList_Add(voxlist, xv, yv, zv);
        }
      }
    }
  }
  return NO_ERROR;
}

//=================================================================
// Diagnostic
//=================================================================

/*--------------------------------------------------------
  MHTtestIsMRISselfIntersecting

  Basically a copy of mrisurf IsMRISselfIntersecting();
  -------------------------------------------------------*/
int MHTtestIsMRISselfIntersecting(MRIS const *mris, float res)
//--------------------------------------------------------
{
  int fno, rslt;

  rslt = 0;

  MRIS_HASH_TABLE * mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, res);
  for (fno = 0; fno < mris->nfaces; fno++) {
    if (MHTdoesFaceIntersect(mht, mris, fno)) {
      rslt = 1;
      goto done;
    }
  }
done:
  MHTfree(&mht);
  return rslt;
}

/*--------------------------------------------------------
  MHTcheckFaces
  Unchanged from mrishash.c 1.27
  -------------------------------------------------------*/
int MHTcheckFaces(MRIS const *mris, MRIS_HASH_TABLE *mht)
//--------------------------------------------------------
{
  static int ncalls = 0;

  if (ncalls++ >= 0) {
#if 0
    checkFace(mht, mris, 144) ;
    checkFace(mht, mris, 185) ;
    checkFace(mht, mris, 16960) ;
    checkFace(mht, mris, 18168) ;
    checkFace(mht, mris, 39705) ;
    checkFace(mht, mris, 32319) ;
    checkAllVertexFaces(mht, mris, 15300) ;
    checkAllVertexFaces(mht, mris, 4303) ;
    checkAllVertexFaces(mht, mris, 35701) ;
    checkAllVertexFaces(mht, mris, 4632) ;
    checkAllVertexFaces(mht, mris, 1573) ;
#endif
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  MHTcheckSurface
  Unchanged from mrishash.c 1.27
  However, subsidiary checkFace is disabled
  ------------------------------------------------------*/
int MHTcheckSurface(MRIS const *mris, MRIS_HASH_TABLE *mht)
//--------------------------------------------------------
{
  int fno, alloced = 0;

  if (!mht) {
    mht = MHTcreateFaceTable(mris);
    alloced = 1;
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    if (!(fno % (mris->nfaces / 10))) DiagHeartbeat((float)fno / (float)mris->nfaces);

    checkFace(mht, mris, fno);
  }
  if (alloced) MHTfree(&mht);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  checkFace
  Conducts a brute-force test of self-intersection relative to
  the one provided face.
  I think the intent here was to conduct brute-force intersections test,
  and if diag flags set, compare that to MHT results.
  However, the only output produced is printed or to file, this
  func doesn't return result of intersection test.
  I have disabled this whole routine because it doesn't look useful,
  and it calls
  local functions that I've refactored away.
  ------------------------------------------------------*/
static int checkFace(MRIS_HASH_TABLE *mht, MRIS const *mris, int fno1)
//--------------------------------------------------------
{
#if 0
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;
  int    fno2, filled, n1, n2, nbr ;
  FACE const   *f1, *f2 ;

  if (fno1 >= mris->nfaces)
    return(NO_ERROR) ;

  f1 = &mris->faces[fno1] ;

  /* fill vertices of 1st triangle */
  v0[0] = (double)mris->vertices[f1->v[0]].x ;
  v0[1] = (double)mris->vertices[f1->v[0]].y ;
  v0[2] = (double)mris->vertices[f1->v[0]].z ;
  v1[0] = (double)mris->vertices[f1->v[1]].x ;
  v1[1] = (double)mris->vertices[f1->v[1]].y ;
  v1[2] = (double)mris->vertices[f1->v[1]].z ;
  v2[0] = (double)mris->vertices[f1->v[2]].x ;
  v2[1] = (double)mris->vertices[f1->v[2]].y ;
  v2[2] = (double)mris->vertices[f1->v[2]].z ;
  for (fno2 = 0 ; fno2 < mris->nfaces ; fno2++)
  {
    f2 = &mris->faces[fno2] ;

    nbr = 0 ;
    for (n1 = 0 ; !nbr && n1 < VERTICES_PER_FACE ; n1++)
      for (n2 = 0 ; !nbr && n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f1->v[n1] == f2->v[n2])
          nbr = 1 ;  /* they share a vertex - don't count it as filled */
      }
    if (nbr)
      continue ;
    u0[0] = (double)mris->vertices[f2->v[0]].x ;
    u0[1] = (double)mris->vertices[f2->v[0]].y ;
    u0[2] = (double)mris->vertices[f2->v[0]].z ;
    u1[0] = (double)mris->vertices[f2->v[1]].x ;
    u1[1] = (double)mris->vertices[f2->v[1]].y ;
    u1[2] = (double)mris->vertices[f2->v[1]].z ;
    u2[0] = (double)mris->vertices[f2->v[2]].x ;
    u2[1] = (double)mris->vertices[f2->v[2]].y ;
    u2[2] = (double)mris->vertices[f2->v[2]].z ;
    filled = tri_tri_intersect(v0,v1,v2,u0,u1,u2) ;
    if (filled && (Gdiag & DIAG_SHOW))
    {
      int    intersect, n ;
      VERTEX const *v ;

      fprintf(stderr,
              "face %d (%d,%d,%d) intersects with face %d (%d,%d,%d)!!!\n",
              fno1, f1->v[0],f1->v[1],f1->v[2],
              fno2, f2->v[0],f2->v[1],f2->v[2]) ;
      if (Gdiag & DIAG_WRITE)
        MRISwrite(mris, "bad") ;
      DiagBreak() ;
      intersect = MHTdoesFaceIntersect(mht, mris, fno1) ;
      if (!intersect)
        mhtHatchFace(mht, mris, fno1, 1) ;
      intersect = MHTdoesFaceIntersect(mht, mris, fno2) ;
      if (!intersect)
        mhtHatchFace(mht, mris, fno2, 1) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        v = &mris->vertices[f1->v[n]] ;
        intersect = MHTisVectorFilled(mht, mris, f1->v[n],v->dx,v->dy,v->dz);
        v = &mris->vertices[f2->v[n]] ;
      }
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
        v = &mris->vertices[f2->v[n]] ;
        intersect = MHTisVectorFilled(mht, mris, f2->v[n],v->dx,v->dy,v->dz);
        v = &mris->vertices[f2->v[n]] ;
      }
    }
  }
#endif
  return (NO_ERROR);
}

VERTEX *MHTfindClosestVertexSetInDirection(
    MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX const *v, int which, double nx, double ny, double nz)
{
  VERTEX *v_closest, *vn;
  double dx, dy, dz, dot, dist, min_dist;
  int vno;

  v_closest = MHTfindClosestVertexSet(mht, mris, v, which);

  if (v_closest) {
    switch (which) {
      case FLATTENED_VERTICES:
        dx = v_closest->fx - v->x;
        dy = v_closest->fx - v->y;
        dz = 0 - v->z;
        break;
      case PIAL_VERTICES:
        dx = v_closest->pialx - v->x;
        dy = v_closest->pialy - v->y;
        dz = v_closest->pialz - v->z;
        break;
      case WHITE_VERTICES:
        dx = v_closest->whitex - v->x;
        dy = v_closest->whitey - v->y;
        dz = v_closest->whitez - v->z;
        break;
      default:
        dx = dy = dz = 0;
        ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MHTfindClosestVertexSet: unsupported which %d", which));
    }
    dot = dx * nx + dy * ny + dz * nz;
    if (dot > 0) return (v_closest);
  }

  min_dist = 1e10;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vn = &mris->vertices[vno];
    if (vn->ripflag) continue;
    switch (which) {
      case FLATTENED_VERTICES:
        dx = vn->fx - v->x;
        dy = vn->fy - v->y;
        dz = vn->fz - v->z;
        break;
      case PIAL_VERTICES:
        dx = vn->pialx - v->x;
        dy = vn->pialy - v->y;
        dz = vn->pialz - v->z;
        break;
      case WHITE_VERTICES:
        dx = vn->whitex - v->x;
        dy = vn->whitey - v->y;
        dz = vn->whitez - v->z;
        break;
      default:
        ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MHTfindClosestVertexSet: unsupported which %d", which));
    }
    dot = dx * nx + dy * ny + dz * nz;
    if (dot < 0) continue;
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    if (dist < min_dist) {
      min_dist = dist;
      v_closest = vn;
    }
  }

  return (v_closest);
}

static void mhtComputeFaceCentroid(
    MRIS const *mris, int which, int fno,
    float *px, float *py, float* pz) 
{
  float xt, yt, zt;
  xt = yt = zt = 0.0;
    
  FACE const * face = &mris->faces[fno];
   
  int n;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    float x = 0, y = 0, z = 0;
    mhtVertex2xyz_float(&mris->vertices[face->v[n]], which, &x, &y, &z);
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

static void mhtStoreFaceCentroids(MHT* mht, MRIS const *mris, int which)
{
  if (mris != mht->mris || mht->nfaces || mht->f) {
    fprintf(stderr, "%s:%d wrong initial state\n", __FILE__, __LINE__);
    exit(1);
  }

  mht->nfaces = mris->nfaces;
  mht->f = (MHT_FACE*)calloc(mht->nfaces, sizeof(MHT_FACE));
  
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    
    float xt,yt,zt;
    mhtComputeFaceCentroid(mris, which, fno, &xt, &yt, &zt);
    
    MHT_FACE *face = &mht->f[fno];
    face->cx = xt;
    face->cy = yt;
    face->cz = zt;
  }
}
