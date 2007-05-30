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
 *    $Date: 2007/05/30 15:35:25 $
 *    $Revision: 1.33 $
 *
 * Copyright (C) 2007,
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

/*-------------------------------------------------------
  GW Design Notes
  ---------------

  1. float, double, Real.
  float and double are official C types.
  Real is typedef'ed indirectly in MINC's basic.h as double.

  mrishash.c 1.27 is inconsistent as to using double or Real, and I'm guessing
  the choice of one or the other was not selected deliberately.
  Since I prefer to get rid of unneeded dependencies, I have replaced Real with
  double.

  ---------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>

//----------------------------------------------------
// Includes that differ for linux vs GW BC compile
//----------------------------------------------------
#ifdef __BORLANDC__
#include <mem.h>
#include "mrishash_include/error.h"
#include "../mrishash_include/mrisurf.h"
//... above includes vanilla mrihash.h
#include "mrishash_include/mrishash_gwmisc.h"
#else
#include <string.h>

#include "diag.h"
#include "error.h"

#include "mrisurf.h"
#include "tritri.h"

#endif

#include "mrishash.h"

//==================================================================
// Local macros
//==================================================================

// GW_VERSION reported by MHTmrishash_gw_version()
#define GW_VERSION 125

// __func__ in gcc appears not to be a predefined macro that can be
// tested and also glued to string literals. Instead we have to use
// printf-style

#define __MYFUNCTION__ __func__

#ifdef  __FUNCTION__
#undef  __MYFUNCTION__
#define __MYFUNCTION__  __FUNCTION__
#endif

#ifdef  __FUNC__
#undef  __MYFUNCTION__
#define __MYFUNCTION__  __FUNC__
#endif

//==================================================================
// Local structure types
//==================================================================
#undef bool
#undef true
#undef false
typedef enum {
  false = 0,
  true  = 1
} bool;

typedef struct
{
  int xv;
  int yv;
  int zv;
}
VOXEL_COORD;

//----------------------------------------------------------------------
// VOXEL_LISTgw adopted from mrishash.c 1.27
// 10000 seems way excessive: suggests that a face can span 10000 voxels,
// which, even with vres=1mm, means a triangle approx 100x200 mm!
// But better safe than sorry, I guess.
//----------------------------------------------------------------------
#define MAX_VOXELS  10000
typedef struct
{
  int nused ;
  int voxels[MAX_VOXELS][3] ;
}
VOXEL_LISTgw ;

// Ad hoc points
typedef struct
{
  double x;
  double y;
  double z;
}
Ptdbl_t;

//==================================================================
// Static forward declarations
//==================================================================
// Note: Exposed API functions start with uppercase MHT prefix.
// Static (private) functions start with lowercase mht prefix.

static void mhtVertex2xyz_float (VERTEX * vtx,
                                 int which,
                                 float  *x, float  *y, float  *z);
static void mhtVertex2xyz_double(VERTEX * vtx,
                                 int which,
                                 double *x, double *y, double *z);
static void mhtVertex2Ptxyz_double(VERTEX * vtx, int which, Ptdbl_t *pt);
static void mhtVertex2array3_double(VERTEX * vtx, int which, double *array3);

static int mhtAddFaceOrVertexAtCoords(MRIS_HASH_TABLE *mht,
                                      float x, float y, float z,
                                      int forvnum);
static int   mhtAddFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht,
                                       int xv, int yv, int zv,
                                       int forvnum);
static int   mhtRemoveFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht,
                                          int xv, int yv, int zv,
                                          int forvnum);

static int mhtFaceToMHT(MRIS_HASH_TABLE *mht,
                        MRI_SURFACE *mris,
                        int fno,
                        int on);
static int   mhtVoxelList_Init(VOXEL_LISTgw *voxlist);
static int   mhtVoxelList_SampleFace(float mhtres,
                                     Ptdbl_t *vpt0,
                                     Ptdbl_t *vpt1,
                                     Ptdbl_t *vpt2,
                                     int fno,
                                     VOXEL_LISTgw *voxlist);
static int     mhtVoxelList_Add(VOXEL_LISTgw *voxlist,
                                int xv, int yv, int zv,
                                int fno);
static int     mhtVoxelList_AddCoord(VOXEL_LISTgw *voxlist,
                                     VOXEL_COORD vc,
                                     int fno);
static int     mhtVoxelList_AddPath(VOXEL_LISTgw *voxlist,
                                    VOXEL_COORD oldvc,
                                    VOXEL_COORD newvc,
                                    int fno );

static int mhtDoesFaceVoxelListIntersect(MRIS_HASH_TABLE *mht,
                                         MRI_SURFACE *mris,
                                         VOXEL_LISTgw *voxlist,
                                         int fno);

//--------- test -----------
static int checkFace(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno1);

//=============================================================================
// Surface --> MHT, store Face Numbers
//=============================================================================

//------------------------------------
MRIS_HASH_TABLE *MHTfillTable(MRI_SURFACE *mris, MRIS_HASH_TABLE *mht)
{
//------------------------------------
  return(MHTfillTableAtResolution(mris, mht, CURRENT_VERTICES, VOXEL_RES)) ;
}

/*------------------------------------------------------
  MHTfillTableAtResolution
  Returns new MHT (freeing old mht if it exists), and loads it
  with FACE NUMBERS of faces in mris.
  Args:
  which     Which vertices to use: CURRENT_VERTICES, WHITE_VERTICES etc
  res       Resolution [mm]
  -------------------------------------------------------*/
//------------------------------------
MRIS_HASH_TABLE * MHTfillTableAtResolution(
  MRI_SURFACE *mris,MRIS_HASH_TABLE *mht,
  int which, float res)
//------------------------------------
{
  int     fno ;
  FACE    *f ;
  int     xv, yv, zv ;
  MHBT    *bucket ;
  static int ncalls = 0 ;

  //-----------------------------
  // Allocation and initialization
  //-----------------------------
  if (mht)    /* free old one */
    MHTfree(&mht) ;
  mht = (MRIS_HASH_TABLE *)calloc(1, sizeof(MRIS_HASH_TABLE)) ;
  if (!mht)
  {
    ErrorExit(ERROR_NO_MEMORY,
              "%s: could not allocate hash table.\n",
              __MYFUNCTION__) ;
  }

  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
      mht->buckets[xv][yv] = NULL ;
    }
  }

  //--------------------------------------
  // Capture data from caller and surface
  //--------------------------------------
  mht->vres            = res ;
  mht->which_vertices  = which ;
  mht->fno_usage       = MHTFNO_FACE ;

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    mhtFaceToMHT(mht, mris, fno, 1) ;
  }

  //-------------------------------------------
  // Diagnostics
  //-------------------------------------------
  if ((Gdiag & DIAG_SHOW) && !ncalls)
  {
    double mean, var, v, n ;
    int    mx ;

    n = mean = 0.0 ;
    mx = -1 ;
    for (xv = 0 ; xv < TABLE_SIZE ; xv++)
    {
      for (yv = 0 ; yv < TABLE_SIZE ; yv++)
      {
        for (zv = 0 ; zv < TABLE_SIZE ; zv++)
        {
          if (!mht->buckets[xv][yv] || !mht->buckets[xv][yv][zv])
            continue ;
          bucket = mht->buckets[xv][yv][zv] ;
          if (bucket->nused)
          {
            mean += bucket->nused ;
            n++ ;
          }
          if (bucket->nused > mx)
            mx = bucket->nused ;
        }
      }
    }
    mean /= n ;
    var = 0.0 ;
    for (xv = 0 ; xv < TABLE_SIZE ; xv++)
    {
      for (yv = 0 ; yv < TABLE_SIZE ; yv++)
      {
        if (!mht->buckets[xv][yv])
          continue ;
        for (zv = 0 ; zv < TABLE_SIZE ; zv++)
        {
          bucket = mht->buckets[xv][yv][zv] ;
          if (bucket && bucket->nused)
          {
            v = mean - bucket->nused ;
            var += v*v ;
          }
        }
      }
    }
    var /= (n-1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%s buckets: mean = %2.1f +- %2.2f, max = %d\n",
              __MYFUNCTION__, mean, sqrt(var), mx) ;
  }
  ncalls++ ;
  return(mht) ;
}

/*------------------------------------------------------
  MHTaddAllFaces
  Appends to mht all faces in mris of which VERTEX v is a part.
  Returns: NO_ERROR unless mht was initialized with vertices rather
  than faces.
  -------------------------------------------------------*/
//------------------------------------
int MHTaddAllFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v)
//------------------------------------
{
  int     fno ;

  if ( (!mht) || (mht->fno_usage != MHTFNO_FACE) )
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not initialized for faces\n",
              __MYFUNCTION__) ;
  }

  for (fno = 0 ; fno < v->num ; fno++)
    mhtFaceToMHT(mht, mris, v->f[fno], 1) ;
  return(NO_ERROR) ;
}

/*------------------------------------------------------
  MHTaddAllFaces
  Appends to mht all faces in mris of which VERTEX v is a part.
  Returns: NO_ERROR unless mht was initialized with vertices rather
  than faces.
  -------------------------------------------------------*/
//------------------------------------
int MHTremoveAllFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v)
//------------------------------------
{
  int     fno ;

  if ( (!mht) || (mht->fno_usage != MHTFNO_FACE) )
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not initialized for faces\n",
              __MYFUNCTION__) ;
  }

  for (fno = 0 ; fno < v->num ; fno++)
    mhtFaceToMHT(mht, mris, v->f[fno], 0) ;
  return(NO_ERROR) ;
}

/*-------------------------------------------------
  mhtFaceToMHT  (was mhtHatchFace)
  Adds face fno to mht. Calls mhtVoxelList_SampleFace to get a list
  of MHT Voxels (buckets) in which to list fno.
  -------------------------------------------------*/
static int mhtFaceToMHT(MRIS_HASH_TABLE *mht,
                        MRI_SURFACE *mris,
                        int fno, int on)
{
//------------------------------------
  FACE   *face ;
  VERTEX *v0, *v1, *v2 ;
  Ptdbl_t vpt0, vpt1, vpt2;
  VOXEL_LISTgw voxlist;
  int vlix, i, j, k;

  mhtVoxelList_Init(&voxlist);

  face = &mris->faces[fno] ;
  if (face->ripflag)
    return(NO_ERROR) ;

  v0 = &mris->vertices[face->v[0]] ;
  v1 = &mris->vertices[face->v[1]] ;
  v2 = &mris->vertices[face->v[2]] ;

  mhtVertex2Ptxyz_double(v0, mht->which_vertices, &vpt0);
  mhtVertex2Ptxyz_double(v1, mht->which_vertices, &vpt1);
  mhtVertex2Ptxyz_double(v2, mht->which_vertices, &vpt2);

  mhtVoxelList_SampleFace(mht->vres, &vpt0, &vpt1, &vpt2, fno, &voxlist);

  for (vlix = 0; vlix < voxlist.nused; vlix++)
  {
    i = voxlist.voxels[vlix][0];
    j = voxlist.voxels[vlix][1];
    k = voxlist.voxels[vlix][2];

    if (on)
      mhtAddFaceOrVertexAtVoxIx(   mht, i, j, k, fno);
    else
      mhtRemoveFaceOrVertexAtVoxIx(mht, i, j, k, fno);
  }
  return(NO_ERROR) ;
}

/*-------------------------------------------------
  mhtVoxelList_SampleFace
  Scans edges and interior of triangle, finding all mht voxels that are
  impinged upon, listing those in voxlist.
  -------------------------------------------------*/
static int mhtVoxelList_SampleFace(float mhtres,
                                   Ptdbl_t *vptin0,
                                   Ptdbl_t *vptin1,
                                   Ptdbl_t *vptin2,
                                   int fno,
                                   VOXEL_LISTgw *voxlist)
{
//------------------------------------
  const float SamplesPerMHTRes = 2.0;
  double mhtres_recip;
  Ptdbl_t vpta, vptb, vptc; // pts opposite the Long, Medium and Short sides
  Ptdbl_t       dif_b2a,    dif_c2a,   dif_b2c;
  Ptdbl_t     delta_b2a,  delta_c2a,  delta_b2c;
  Ptdbl_t      posn_b2a,   posn_c2a,   posn_b2c;
  VOXEL_COORD voxco_b2a,    voxco_c2a,    voxco_b2c;
  VOXEL_COORD oldvoxco_b2a, oldvoxco_c2a, oldvoxco_b2c;
  oldvoxco_b2a.xv=0; oldvoxco_b2a.yv=0; oldvoxco_b2a.zv=0;
  oldvoxco_c2a.xv=0; oldvoxco_c2a.yv=0; oldvoxco_c2a.zv=0;
  double LenSq_b2a, LenSq_c2a;    // Distances along sides
  double TempLenSq, TempLen, StepsReqdMain_dbl, StepsReqdRung_dbl;
  int    mainstep, rungstep, mainsteps_reqd, rungsteps_reqd;
  double MainStepFrac, RungStepFrac;
  bool changed;
#define PTLENSQ(a) ( a.x*a.x + a.y*a.y +a.z*a.z )
#define PTDIF(ans, a, b)  ans.x  = a.x - b.x;   ans.y  = a.y - b.y;   ans.z  = a.z - b.z
#define PTADD(ans, a, b)  ans.x  = a.x + b.x;   ans.y  = a.y + b.y;   ans.z  = a.z + b.z
#define PTINC(ans, a)     ans.x += a.x;         ans.y += a.y;         ans.z += a.z
#define PTMULTK(ans,b,K)  ans.x  = b.x * K;     ans.y  = b.y * K;     ans.z  = b.z * K
#define PTABS(ans, a)     ans.x  = abs(a);      ans.y  = abs(a.y);    ans.z  = abs(a.z)
#define W2VOL(ares, x)    ( ((x)+FIELD_OF_VIEW/2)/(ares) )
#define W2VOX(ares, x)    ( (int)(W2VOL(ares,x)) )
#define PTWORLD2VOXEL(ans, ares, pt) ans.xv = W2VOX(ares, pt.x); ans.yv = W2VOX(ares, pt.y); ans.zv = W2VOX(ares, pt.z)
#define SAME_VOXEL(a, b) ( (a.xv == b.xv) && (a.yv == b.yv) && (a.zv == b.zv) )

// Following is slightly less than one to avoid voxel-edge overshoot on rungs
// (makes example that are on exact coords look better.)
#define RUNGRECIPFACTOR 0.9995

  mhtres_recip = 1;
  if (mhtres != 0)
    mhtres_recip = 1.0/mhtres; // for a little speed gain

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
  TempLenSq         = MAX(LenSq_b2a, LenSq_c2a);
  TempLen           = sqrt(TempLenSq);
  StepsReqdMain_dbl = TempLen * SamplesPerMHTRes * mhtres_recip;
  mainsteps_reqd    = ceil(StepsReqdMain_dbl);
  if (0 >= mainsteps_reqd)    // can't be less than!
    mainsteps_reqd = 1;
  MainStepFrac      = 1.0/mainsteps_reqd;

  PTMULTK(delta_b2a, dif_b2a, MainStepFrac);
  PTMULTK(delta_c2a, dif_c2a, MainStepFrac);

  // Set starting positions
  posn_b2a = vptb;
  posn_c2a = vptc;

  for (mainstep = 0; mainstep <= mainsteps_reqd; mainstep++)
  {
    if (mainstep >= 1)
    {
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
    if (0 == mainstep)
    {
      mhtVoxelList_AddCoord(voxlist, voxco_b2a, fno);
      mhtVoxelList_AddCoord(voxlist, voxco_c2a, fno);
      changed = true;
    }
    else
    {
      //---------------------------------------------
      // If we crossed a boundary, add that "path" to voxlist
      //---------------------------------------------
      if (!SAME_VOXEL(oldvoxco_b2a, voxco_b2a) )
      {
        mhtVoxelList_AddPath(voxlist, oldvoxco_b2a, voxco_b2a, fno);
        changed = true;
      }

      if (!SAME_VOXEL(oldvoxco_c2a, voxco_c2a) )
      {
        mhtVoxelList_AddPath(voxlist, oldvoxco_c2a, voxco_c2a, fno);
        changed = true;
      }
    }

    if (changed)
    {
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

    if (SAME_VOXEL(voxco_c2a, voxco_b2a) )
      goto rung_done; /* Nothing to be done if both ends
                         are already in same voxel */

    PTDIF(dif_b2c, posn_c2a, posn_b2a);
    TempLenSq         = PTLENSQ(dif_b2c);
    TempLen           = sqrt(TempLenSq);
    StepsReqdRung_dbl =  TempLen * SamplesPerMHTRes * mhtres_recip;
    rungsteps_reqd     = ceil(StepsReqdRung_dbl);

    if (rungsteps_reqd < 1)
      rungsteps_reqd = 1;


    RungStepFrac      = RUNGRECIPFACTOR/rungsteps_reqd;
    PTMULTK(delta_b2c, dif_b2c, RungStepFrac);

    // Starting conditions
    posn_b2c      = posn_b2a;
    voxco_b2c     = voxco_b2a;
    oldvoxco_b2c  = voxco_b2c;

    //----------------------------------------------------
    // Step along a rung (xxx_b2c)
    // No action on step 0 (ie: starting position) , as this voxel is
    // already in the mht as a result of the main stepping.
    // Constrastingly, we *do* need action on *last* step even though that
    // voxel is in mht, because we may need to "draw" from some
    // other voxel *to* that one
    //------------------------------------------------------------
    for (rungstep = 1; rungstep <= rungsteps_reqd; rungstep++)
    {
      PTINC(posn_b2c, delta_b2c);
      PTWORLD2VOXEL(voxco_b2c, mhtres, posn_b2c);
      //---------------------------------------------
      // If we crossed a boundary, add voxels in that
      // "path" to voxlist
      //---------------------------------------------
      if (!SAME_VOXEL(oldvoxco_b2c, voxco_b2c) )
      {
        mhtVoxelList_AddPath(voxlist, oldvoxco_b2c, voxco_b2c, fno);
        oldvoxco_b2c = voxco_b2c;
      }
    } // for rungstep
  rung_done:
    ; // semicolon avoids error message
  } // for mainstep
  return NO_ERROR;
}

//=============================================================================
// Surface --> MHT, store Vertex Numbers
//=============================================================================

//---------------------------------------------------------
MRIS_HASH_TABLE *MHTfillVertexTable(
  MRI_SURFACE *mris,MRIS_HASH_TABLE *mht, int which)
//---------------------------------------------------------
{
  return(MHTfillVertexTableRes(mris, mht, which, VOXEL_RES)) ;
}

//---------------------------------------------------------
MRIS_HASH_TABLE *MHTfillVertexTableRes(
  MRI_SURFACE *mris,MRIS_HASH_TABLE *mht, int which, float res)
//---------------------------------------------------------
{
  int     vno ;
  int     xv, yv, zv ;
  float   x=0.0, y=0.0, z=0.0;
  MHBT    *bucket ;
  VERTEX  *v ;
  static int ncalls = 0 ;

  //-----------------------------
  // Allocation and initialization
  //-----------------------------
  if (mht)    /* free old one */
    MHTfree(&mht) ;
  mht = (MRIS_HASH_TABLE *)calloc(1, sizeof(MRIS_HASH_TABLE)) ;
  if (!mht)
    ErrorExit(ERROR_NO_MEMORY,
              "%s: could not allocate hash table.\n", __MYFUNCTION__ ) ;

  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
      mht->buckets[xv][yv] = NULL ;
    }
  }

  //--------------------------------------
  // Capture data from caller and surface
  //--------------------------------------
  mht->vres            = res ;
  mht->which_vertices  = which ;
  mht->fno_usage       = MHTFNO_VERTEX ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);
    mhtAddFaceOrVertexAtCoords(mht, x, y, z, vno);
  }

  //-------------------------------------------
  // Diagnostics
  // Not modified from 1.27
  //-------------------------------------------
  if ((Gdiag & DIAG_SHOW) && !ncalls)
  {
    double mean, var, v, n ;
    int    mx ;

    n = mean = 0.0 ;
    mx = -1 ;
    for (xv = 0 ; xv < TABLE_SIZE ; xv++)
    {
      for (yv = 0 ; yv < TABLE_SIZE ; yv++)
      {
        for (zv = 0 ; zv < TABLE_SIZE ; zv++)
        {
          if (!mht->buckets[xv][yv] || !mht->buckets[xv][yv][zv])
            continue ;
          bucket = mht->buckets[xv][yv][zv] ;
          if (bucket->nused)
          {
            mean += bucket->nused ;
            n++ ;
          }
          if (bucket->nused > mx)
            mx = bucket->nused ;
        }
      }
    }
    mean /= n ;
    var = 0.0 ;
    for (xv = 0 ; xv < TABLE_SIZE ; xv++)
    {
      for (yv = 0 ; yv < TABLE_SIZE ; yv++)
      {
        if (!mht->buckets[xv][yv])
          continue ;
        for (zv = 0 ; zv < TABLE_SIZE ; zv++)
        {
          bucket = mht->buckets[xv][yv][zv] ;
          if (bucket && bucket->nused)
          {
            v = mean - bucket->nused ;
            var += v*v ;
          }
        }
      }
    }
    var /= (n-1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%s buckets: mean = %2.1f +- %2.2f, max = %d\n",
              __MYFUNCTION__, mean, sqrt(var), mx) ;
  }
  ncalls++ ;
  return(mht) ;
}

//=================================================================
// AddFaceOrVertexXxxx functions
//=================================================================

/*------------------------------------------------------------
  mhtAddFaceOrVertexAtVoxIx  Was: mhtAddFacePosition
  Adds forvnum (Face or Vertex number) to mht, in bucket(xv,yv,zv)
  Coerces (xv,yv,zv) to be sane.
  -------------------------------------------------------------*/
static int mhtAddFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht,
                                     int xv, int yv, int zv,
                                     int forvnum)
{
//---------------------------------------
  int    i ;
  MHB   *bin ;
  MHBT  *bucket ;

  if (xv < 0)
    xv = 0 ;
  if (xv >= TABLE_SIZE)
    xv  = TABLE_SIZE-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= TABLE_SIZE)
    yv  = TABLE_SIZE-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= TABLE_SIZE)
    zv  = TABLE_SIZE-1 ;

  // [GW] (removed xv,yv,zv check here, as it's impossible to trigger)

  //-----------------------------------------------
  // Allocate space if needed
  //-----------------------------------------------
  // 1. Allocate a 1-D array at buckets[xv][yv]
  if (!mht->buckets[xv][yv])
  {
    mht->buckets[xv][yv] = (MHBT **)calloc(TABLE_SIZE, sizeof(MHBT *)) ;
    if (!mht->buckets[xv][yv])
      ErrorExit(ERROR_NO_MEMORY,
                "%s: could not allocate slice.",  __MYFUNCTION__) ;
  }
  // 2. Allocate a bucket at buckets[xv][yv][zv]
  bucket = mht->buckets[xv][yv][zv] ;
  if (!bucket)
  {
    mht->buckets[xv][yv][zv] = bucket = (MHBT *)calloc(1, sizeof(MHBT)) ;
    if (!bucket)
      ErrorExit(ERROR_NOMEMORY,
                "%s couldn't allocate bucket.\n",
                __MYFUNCTION__) ;
  }
  // 3. Allocate a bucket at buckets[xv][yv][zv]
  if (!bucket->max_bins)   /* nothing in this bucket yet - allocate bins */
  {
    bucket->max_bins = 4 ;
    bucket->bins = (MHB *)calloc(bucket->max_bins, sizeof(MHB));
    if (!bucket->bins)
      ErrorExit(ERROR_NO_MEMORY,
                "%s: could not allocate %d bins.\n",
                __MYFUNCTION__,  bucket->max_bins) ;
  }

  //-------------------------------------------------------------
  // If this forvnum already is listed in this bucket then
  // nothing to do.
  //-------------------------------------------------------------
  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    if (bin->fno == forvnum)
      goto done;
  }

  //-------------------------------------------------------------
  // Add forvnum to this bucket
  //-------------------------------------------------------------
  if (i >= bucket->nused)   /* forvnum not already listed at this bucket */
  {
    //------ Increase allocation if needed -------
    if (bucket->nused >= bucket->max_bins)
    {
      bin = bucket->bins ;
      bucket->max_bins *= 2 ;
      bucket->bins = (MHB *)calloc(bucket->max_bins, sizeof(MHB));
      if (!bucket->bins)
        ErrorExit(ERROR_NO_MEMORY,
                  "%s: could not allocate %d bins.\n",
                  __MYFUNCTION__,  bucket->max_bins) ;
      memmove(bucket->bins, bin, bucket->nused*sizeof(MHB)) ;
      free(bin) ;
      bin = &bucket->bins[i] ;
    }
    //----- add this face-position to this bucket ------
    bucket->nused++ ;
    bin->fno = forvnum ;
  }
 done:
  return(NO_ERROR) ;
}

/*------------------------------------------------------------
  mhtAddFaceOrVertexAtCoords
  Converts x, y, z to bucket indexes, calls mhtAddFaceOrVertexAtIndexes
  Adds forvnum (Face or Vertex number) to mht, in bucket(xv,yv,zv)
  -------------------------------------------------------------*/
static int mhtAddFaceOrVertexAtCoords(MRIS_HASH_TABLE *mht,
                                      float x, float y, float z,
                                      int forvnum)
{
//---------------------------------------
  int xv, yv, zv;

  xv = WORLD_TO_VOXEL(mht, x) ;
  yv = WORLD_TO_VOXEL(mht, y) ;
  zv = WORLD_TO_VOXEL(mht, z) ;

  return mhtAddFaceOrVertexAtVoxIx(mht, xv, yv, zv, forvnum);
}

/*------------------------------------------------------------
  mhtRemoveFaceOrVertexAtVoxIx (was mhtRemoveFacePosition)
  Reverse of mhtAddFaceOrVertexAtIndexes
  -------------------------------------------------------------*/
static int mhtRemoveFaceOrVertexAtVoxIx(MRIS_HASH_TABLE *mht,
                                        int xv, int yv, int zv,
                                        int forvnum)
{
//---------------------------------------
  int    i ;
  MHBT   *bucket ;
  MHB    *bin ;

  if (xv < 0)
    xv = 0 ;
  if (xv >= TABLE_SIZE)
    xv  = TABLE_SIZE-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= TABLE_SIZE)
    yv  = TABLE_SIZE-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= TABLE_SIZE)
    zv  = TABLE_SIZE-1 ;

  if (!mht->buckets[xv][yv])
    return(NO_ERROR) ;        // no bucket at such coordinates
  bucket = mht->buckets[xv][yv][zv] ;
  if (!bucket)
    return(NO_ERROR) ;       // no bucket at such coordinates

  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    if (bin->fno == forvnum)     /* found face */
      break ;
  }
  if (i < bucket->nused)    /* found face bucket - remove it */
  {
    bucket->nused-- ;
    if (i < bucket->nused)  /* not the last one in the list - compact list */
    {
      int nbytes = (bucket->nused-i) * sizeof(MHB) ;
      MHB *src_bin, *dst_bin ;
      src_bin = &bucket->bins[i+1] ;
      dst_bin = &bucket->bins[i] ;
      memmove(dst_bin, src_bin, nbytes) ;
    }
  }
  return(NO_ERROR) ;
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
int MHTisVectorFilled(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris,
                      int vtxno, float dx, float dy, float dz)
{
//--------------------------------------------------------------
  VERTEX   *vtx ;
  float    savex, savey, savez ;
  int      intersect, fno ;

  intersect = 0; // default result

  //----------------------------------------------------
  // Sanity check
  //----------------------------------------------------
  if (!mht ) ErrorExit(ERROR_BADPARM,  "%s: mht is NULL\n", __MYFUNCTION__) ;
  if (!mris) ErrorExit(ERROR_BADPARM,  "%s: mris is NULL\n", __MYFUNCTION__) ;
  if (mht->fno_usage != MHTFNO_FACE )
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not initialized for vertices\n",
              __MYFUNCTION__) ;
  }
  if (mht->which_vertices != CURRENT_VERTICES)
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not loaded using CURRENT_VERTICES\n",
              __MYFUNCTION__) ;
  }

  //----------------------------------------------------
  // Temporarily move the CURRENT_VERTICES
  //----------------------------------------------------
  vtx = &mris->vertices[vtxno] ;
  savex = vtx->x ;
  savey = vtx->y ;
  savez = vtx->z ;
  vtx->x += dx ;
  vtx->y += dy ;
  vtx->z += dz ;

  //-------------------------------------------
  // Check whether the faces adjoining current
  // vertex will now intersect
  //-------------------------------------------
  for (fno = 0 ; !intersect && (fno < vtx->num) ; fno++)
  {
    intersect = MHTdoesFaceIntersect(mht, mris, vtx->f[fno]) ;
  }

  //----------------------------------------------------
  // Restore CURRENT_VERTICES
  //----------------------------------------------------
  vtx->x = savex ;
  vtx->y = savey ;
  vtx->z = savez ;
  return(intersect) ;
}

/*-------------------------------------------------------------
  MHTdoesFaceIntersect
  Does a particular face intersect with any other faces in surface mris,
  which was previously hashed into mht using MHTfillTableXxx
  Returns 1 for "yes -- intersection detected", else NO_ERROR.
  -------------------------------------------------------------*/
int MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris,int fno)
{
//------------------------------------------------------------
  VERTEX     *v0, *v1, *v2 ;
  Ptdbl_t vpt0, vpt1, vpt2;
  FACE       *face ;
  VOXEL_LISTgw voxlist ;

  int retval = 0;

  //----------------------------------------------------
  // Sanity check
  //----------------------------------------------------
  if (!mht ) ErrorExit(ERROR_BADPARM,  "%s: mht is NULL\n", __MYFUNCTION__) ;
  if (!mris) ErrorExit(ERROR_BADPARM,  "%s: mris is NULL\n", __MYFUNCTION__) ;
  if (mht->fno_usage != MHTFNO_FACE)
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not initialized for vertices\n",
              __MYFUNCTION__) ;
  }

  //--------------------------------
  // Create a voxlist for face fno
  //--------------------------------
  mhtVoxelList_Init(&voxlist) ;
  face = &mris->faces[fno] ;

  if (face->ripflag)
    return(NO_ERROR) ;

  if (fno == Gdiag_no)
    DiagBreak() ;

  v0 = &mris->vertices[face->v[0]] ;
  v1 = &mris->vertices[face->v[1]] ;
  v2 = &mris->vertices[face->v[2]] ;

  mhtVertex2Ptxyz_double(v0, mht->which_vertices, &vpt0);
  mhtVertex2Ptxyz_double(v1, mht->which_vertices, &vpt1);
  mhtVertex2Ptxyz_double(v2, mht->which_vertices, &vpt2);

  mhtVoxelList_SampleFace(mht->vres, &vpt0, &vpt1, &vpt2, fno, &voxlist);

  //--------------------------------
  // Any intersections?
  //--------------------------------
  if (mhtDoesFaceVoxelListIntersect(mht, mris, &voxlist, fno))
    retval = 1;

  return(retval) ;
}

#define MHT_MAX_FACES 10000
/*-------------------------------------------------------------------
  mhtDoesFaceVoxelListIntersect
  Subsidiary function for MHTdoesFaceIntersect, given face fno already
  analyzed into voxlist.
  ------------------------------------------------------------------*/
static int mhtDoesFaceVoxelListIntersect(
  MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VOXEL_LISTgw *voxlist, int fno)
//------------------------------------------------------------------
{
  int    xv, yv, zv, n3, intersect, voxnum;
  int    binix, faceix, facetestno;
  int    facelist[MHT_MAX_FACES], nfaces;
  MHB    *bin ;
  MHBT   *bucket ;
  FACE   *facein, *facetest ;
  VERTEX * avtx;
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;

  //----------------------------------------------------
  // Pertaining to cross-checking vertices of fno versus
  // test face
  //----------------------------------------------------
  int facein_vtx_ix ,  facetest_vtx_ix; // indices of vertices in face's list
  int facein_vtx_vno,  facetest_vtx_vno;   // Fno of vertices in mris
  VERTEX  *facein_vtx    , *facetest_vtx;       // The vertices themselves
  VERTEX_INFO *vtxinfos, *vtxinfo;

  //----------------------------------------------------
  facein = &mris->faces[fno] ;

  if (fno == Gdiag_no)
    DiagBreak() ;

  //------------------------------------------------------------------------
  // Iterate through all the voxels pertaining to face fno (ie:
  // the voxels listed in voxlist), looking at those voxel buckets
  // in mht to find triangles to add to list in flist for later intersection
  // check
  //------------------------------------------------------------------------
  vtxinfos = (VERTEX_INFO *)mris->user_parms ;
  nfaces = 0;
  for (voxnum = 0 ; voxnum < voxlist->nused ; voxnum++)
  {
    xv = voxlist->voxels[voxnum][0] ;
    yv = voxlist->voxels[voxnum][1] ;
    zv = voxlist->voxels[voxnum][2] ;

    //----------------------------------------------------------
    // Get corresponding bucket from mht. This *should* always
    // succeed, given that the same info was just used to put xv,yv,zv
    // into voxlist as was used to put faces into mht buckets.
    //----------------------------------------------------------
    if (!mht->buckets[xv][yv])
      return(0) ;
    bucket = mht->buckets[xv][yv][zv] ;
    if (!bucket)
      continue ;

    bin = bucket->bins ;

    //-------------------------------------------------
    // Following is one iteration per face
    // (ie: binix is essentially facetestix)
    //-------------------------------------------------
    for (binix = 0 ; binix < bucket->nused ; binix++, bin++)
    {
      facetestno = bin->fno ;  // get the one face from this bin

      if (facetestno == fno)  // no point comparing to same face!
        goto skip_this_facetest;

      facetest   = &mris->faces[facetestno] ;

      //-------------------------------------------------------------
      // Tests: Several tests to see whether we should skip the
      // actual intersection check for this face, eg because facetest
      // adjoins facein or they are linked.
      //------------------------------------------------------------
      for (facein_vtx_ix = 0 ;
           facein_vtx_ix < VERTICES_PER_FACE ;
           facein_vtx_ix++)
      {
        facein_vtx_vno = facein->v[facein_vtx_ix];
        facein_vtx     = &(mris->vertices[facein_vtx_vno]);

        for (facetest_vtx_ix = 0 ;
             facetest_vtx_ix < VERTICES_PER_FACE ;
             facetest_vtx_ix++)
        {
          facetest_vtx_vno = facetest->v[facetest_vtx_ix];
          facetest_vtx     = &(mris->vertices[facetest_vtx_vno]);

          // Do the faces share a vertex?
          if (facein_vtx_vno == facetest_vtx_vno)
            goto skip_this_facetest;

          // Are they linked? Check facein's list of linked vertices
          if (facein_vtx->linked > 0)
          {
            vtxinfo  = &(vtxinfos[facein_vtx_vno]);
            for (n3 = 0 ; n3 < vtxinfo->nlinks ; n3++)
            {
              if (vtxinfo->linked_vno[n3] == facetest_vtx_vno)
                goto skip_this_facetest;
            }
          } // if

          // Are they linked? Check facetest's list of linked vertices
          if (facetest_vtx->linked > 0)
          {
            vtxinfo  = &(vtxinfos[facetest_vtx_vno]);
            for (n3 = 0 ; n3 < vtxinfo->nlinks ; n3++)
            {
              if (vtxinfo->linked_vno[n3] == facein_vtx_vno)
                goto skip_this_facetest;
            } // for
          } // if
        } // for facetest_vtx_ix
      } // for facein_vtx_ix

      //-----------------------------------------------------
      // If we get to here, there was no reason to reject facetest,
      // so append it to faces list
      //-----------------------------------------------------
      for (faceix = 0 ; faceix < nfaces ; faceix++)
      {
        if (facelist[faceix] == facetestno)
          goto skip_this_facetest;
      }
      if (nfaces >= MHT_MAX_FACES)
        ErrorPrintf(ERROR_NO_MEMORY,
                    "%s: MHT_MAX_FACES exceeded!",
                    __MYFUNCTION__) ;
      else
        facelist[nfaces++] = facetestno ;

    skip_this_facetest:
      ; // semicolon avoids error message
    } // for binix
  } // for voxnum

  //-----------------------------------------------------------------------
  // Now that we have a list of faces in facelist run the actual geometry
  // intersection test
  // Notes:
  // 1. Code up to v1.27 used CURRENT_VERTICES, which meant that using this
  // code with an mht that was created with MHTfillTableAtResolution and some
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

  for (faceix = 0 ; faceix < nfaces ; faceix++)
  {
    /* set vertices of 2nd triangle */
    facetest = &mris->faces[facelist[faceix]] ;

    avtx = &(mris->vertices[facetest->v[0]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u0);

    avtx = &(mris->vertices[facetest->v[1]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u1);

    avtx = &(mris->vertices[facetest->v[2]]);
    mhtVertex2array3_double(avtx, mht->which_vertices, u2);

    intersect = tri_tri_intersect(v0,v1,v2,  u0,u1,u2) ;
    if (intersect)
      return(1) ;
  }
  return(0) ;
}

//=================================================================
// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//=================================================================

//---------------------------------------------
void mhtFindCommonSanityCheck(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris)
{
//---------------------------------------------
  if (!mht ) ErrorExit(ERROR_BADPARM,  "%s: mht is NULL\n", __MYFUNCTION__) ;
  if (!mris) ErrorExit(ERROR_BADPARM,  "%s: mris is NULL\n", __MYFUNCTION__) ;
  if ( (!mht) || (mht->fno_usage != MHTFNO_VERTEX) )
  {
    ErrorExit(ERROR_BADPARM,
              "%s: mht not initialized for vertices\n",
              __MYFUNCTION__) ;
  }
}

//------------------------------------------
// Simple instrumentation
//------------------------------------------
int FindBucketsChecked_Count;
int FindBucketsPresent_Count;

void MHTfindReportCounts(int * BucketsChecked, int * BucketsPresent)
{
  *BucketsChecked = FindBucketsChecked_Count;
  *BucketsPresent = FindBucketsPresent_Count;
}

/*----------------------------------------------------------------
  mhtfindClosestVertexGenericInBucket

  -----------------------------------------------------------------*/
int mhtfindClosestVertexGenericInBucket(MRIS_HASH_TABLE *mht, 
                                        MRI_SURFACE *mris,
                                        //---------- inputs --------------
                                        int xv, int yv, int zv,
                                        double probex, 
                                        double probey, 
                                        double probez,
                                        //---------- in/outs -------------
                                        VERTEX **MinDistVtx, 
                                        int *MinDistVtxNum, 
                                        double *MinDistSq)
{
  int vtxix, rslt;
  VERTEX *AVtx ;
  int    AVtxNum;
  double ADistSq;
  MHB       *bin ;
  MHBT      *bucket ;
  //----------------------------------
  rslt = NO_ERROR;

  FindBucketsChecked_Count++;

  bucket = MHTgetBucketAtVoxIx(mht, xv, yv, zv);
  if (!bucket)
    goto done;

  FindBucketsPresent_Count++;

  //-----------------------------------------
  // Iterate through vertices in this bucket
  //-----------------------------------------
  bin = bucket->bins ;
  for (vtxix = 0 ; vtxix < bucket->nused ; vtxix++, bin++)
  {
    AVtxNum = bin->fno;
    AVtx = &mris->vertices[AVtxNum];

    if (AVtxNum == Gdiag_no)
      DiagBreak() ;

    ADistSq = SQR(AVtx->x - probex)
      + SQR(AVtx->y - probey)
      + SQR(AVtx->z - probez) ;

    if (ADistSq < *MinDistSq)
    {
      *MinDistSq     = ADistSq ;
      *MinDistVtxNum = AVtxNum;
      *MinDistVtx    = AVtx;
    } // if
  } // for vtxix
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
                                MRI_SURFACE *mris,
                                //---------- inputs --------------
                                double probex, double probey, double probez,
                                // How far to search: set one or both
                                double in_max_distance_mm, /* Use large number 
                                                              to ignore */
                                int    in_max_mhts,  /* Use -1 to ignore */
                                //---------- outputs -------------
                                VERTEX **pvtx, 
                                int *vtxnum, 
                                double *vtx_distance)
{
  const int max_mhts_MAX = 5;
  double mhtres, max_distance_mm, tempdbl;
  int    max_mhts;
  double probex_vol,
  probey_vol, probez_vol;// probex etc translated to volume space
  int    probex_vox,
  probey_vox, probez_vox;// probex_vol etc to voxel index
  double probex_mod,
  probey_mod, probez_mod;// probex_vol remainder (posn of probe within voxel)
  int    near8offsetx, near8offsety, near8offsetz; // probe?_mod to -1..0
  int   xv, yv, zv, xvi, yvi, zvi;                         // voxel indices
  VERTEX *MinDistVtx;
  int    MinDistVtxNum;
  double MinDistSq, MinDistTemp;
  double RemainingVoxelDistance;
  int    RVox, WallJump;
  bool   isWall;
  unsigned char  central27[3][3][3]; // Indexes 0..2 stand for -1..+1
  //----------------------------------

  mhtFindCommonSanityCheck(mht, mris);
  mhtres     = mht->vres;

  //--------------------------------------------------
  // Initialize instrumentation
  //--------------------------------------------------
  FindBucketsChecked_Count = 0;
  FindBucketsPresent_Count = 0;

  //--------------------------------------------------
  // Figure how far afield to search
  //--------------------------------------------------
  if (-1 == in_max_mhts)
  { // use in_max_distance_mm
    max_distance_mm = in_max_distance_mm;
    if ((max_distance_mm * 2) <= mhtres)
    {
      max_mhts = 0;
    }
    else
    {
      max_mhts        = ceil(max_distance_mm/mhtres);
    }
  }
  else
  {
    max_mhts        = in_max_mhts;
    max_distance_mm = in_max_distance_mm;

    // How far does max_mhts cover in mm?
    if (max_mhts >= 1)
    {
      tempdbl         = max_mhts * mhtres;
    }
    else
    {
      tempdbl         = 0.5 * mhtres;
    }
    // Must not include points beyond the exhaustive coverage distance
    // of the chosen max_mhts....
    if (max_distance_mm > tempdbl) max_distance_mm = tempdbl;
  }
  // Safety limit
  if (max_mhts > max_mhts_MAX)
  max_mhts = max_mhts_MAX;

  //--------------------------------------------------
  // Initialize mins
  //--------------------------------------------------
  MinDistSq     = 1e6;
  MinDistVtx    = NULL;
  MinDistVtxNum = -1;

  //--------------------------------------------------
  // Translate probe point to voxel-space coord and indexes
  //--------------------------------------------------
  probex_vol = WORLD_TO_VOLUME(mht, probex);
  probey_vol = WORLD_TO_VOLUME(mht, probey);
  probez_vol = WORLD_TO_VOLUME(mht, probez);

  // (Note: In following (int) truncs toward zero, but that's OK because
  // range of probex_vol is all positive, centered at FIELD_OF_VIEW/2)
  probex_vox = (int) probex_vol;
  probey_vox = (int) probey_vol;
  probez_vox = (int) probez_vol;

  probex_mod = probex_vol - (double) probex_vox;
  probey_mod = probey_vol - (double) probey_vox;
  probez_mod = probez_vol - (double) probez_vox;

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
  for (    xvi = 0; xvi <= 1; xvi++)
  {
    xv = xvi + near8offsetx;
    for (  yvi = 0; yvi <= 1; yvi++)
    {
      yv = yvi + near8offsety;
      for (zvi = 0; zvi <= 1; zvi++)
      {
        zv = zvi + near8offsetz;

        mhtfindClosestVertexGenericInBucket(mht, mris,
                                            probex_vox + xv,
                                            probey_vox + yv,
                                            probez_vox + zv,
                                            probex, probey, probez,
                                            &MinDistVtx,
                                            &MinDistVtxNum,
                                            &MinDistSq) ;

        central27[xv+1][yv+1][zv+1] = 1;
      }
    }
  }

  //---------------------------------------------------------------------------
  // We can stop now if found vertex's distance is < 0.5 mhtres, because all
  // other voxels are at at least that far away
  //---------------------------------------------------------------------------

  if (max_mhts == 0) goto done;

  RemainingVoxelDistance = 0.5 *  mhtres;
  if (max_distance_mm <= RemainingVoxelDistance)  goto done;

  if (MinDistVtx)
  { // not NULL if one has been found)
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp <= RemainingVoxelDistance)
      goto done;
  }

  //--------------------------------------------------
  // Continue with rest of central 27
  //--------------------------------------------------
  for (    xv = -1; xv <= 1; xv++)
  {
    for (  yv = -1; yv <= 1; yv++)
    {
      for (zv = -1; zv <= 1; zv++)
      {

        if (!central27[xv+1][yv+1][zv+1]) // skip ones already done
        {
          mhtfindClosestVertexGenericInBucket(mht, mris,
                                              probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex, probey, probez,
                                              &MinDistVtx,
                                              &MinDistVtxNum,
                                              &MinDistSq) ;

          central27[xv+1][yv+1][zv+1] = 1;
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

  for (RVox = 2; RVox <= max_mhts; RVox++)
  {
    //----------------------------------------------------------------------
    // We can stop now if found vertex's distance is<(1.0 x (RVox-1) x mhtres),
    // because all other voxels are at at least that far away
    //----------------------------------------------------------------------
    RemainingVoxelDistance = (mhtres * (RVox-1));

    if (MinDistVtx)
    { // not NULL if one has been found)
      MinDistTemp = sqrt(MinDistSq);
      if (MinDistTemp <= RemainingVoxelDistance)
        goto done;
    }
    if (max_distance_mm <= RemainingVoxelDistance)
    goto done;

    //-------------------------------------------------
    // Inspect "shell" of voxels
    //-------------------------------------------------
    WallJump =
    RVox + RVox; // jump from one side to the other across the empty middle
    for     (xv  = -RVox; xv <= RVox; xv++)
    {
      for   (yv  = -RVox; yv <= RVox; yv++)
      {
        isWall =  ( (xv == -RVox) || (xv == RVox) )
          || ( (yv == -RVox) || (yv == RVox) );
        for (zv  = -RVox; zv <= RVox; zv = isWall ? zv+1 : zv + WallJump)
        {

          mhtfindClosestVertexGenericInBucket(mht, mris,
                                              probex_vox + xv,
                                              probey_vox + yv,
                                              probez_vox + zv,
                                              probex, probey, probez,
                                              &MinDistVtx,
                                              &MinDistVtxNum,
                                              &MinDistSq) ;
        } // zv
      } // yv
    } // xv
  } // RVox

  done:
  MinDistTemp = 1e3;

  //--------------------------------------------
  // Enforce not returning a vertex if it's outside
  // max_distance_mm.
  //--------------------------------------------
  if (MinDistVtx)
  {
    MinDistTemp = sqrt(MinDistSq);
    if (MinDistTemp > max_distance_mm)
    {
      MinDistVtx     = NULL;
      MinDistVtxNum  = -1;
      MinDistTemp    = 1e3;
    }
  }

  // Copy to output
  if (pvtx)         *pvtx         = MinDistVtx;
  if (vtxnum)       *vtxnum       = MinDistVtxNum;
  if (vtx_distance) *vtx_distance = MinDistTemp;

  return NO_ERROR;
}

/*----------------------------------------------------------------
  MHTfindClosestVertex
  Returns VERTEX *, closest vertex to v in mris & mht (or NULL).
  -----------------------------------------------------------------*/
VERTEX * MHTfindClosestVertex(MRIS_HASH_TABLE *mht, 
                              MRI_SURFACE *mris, 
                              VERTEX *v)
{
//------------------------------------------------------
  VERTEX * vtx;
  int rslt;
  float x=0.0, y=0.0, z=0.0;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);

  rslt = MHTfindClosestVertexGeneric(mht, mris,
                                     x, y, z,
                                     1000,
                                     1,  // max_mhts: search out to 3 x 3 x 3
                                     &vtx, NULL, NULL);

  return vtx;
}

/*---------------------------------------------------------------------
  MHTfindClosestVertexSet.
  Returns vertex in mht and mris that is closest to v.
  Ostensibly allows choice of which surface vertex "set", but this is
  conceptually incorrect, as vertex set is determined at time MHT is created.
  (So GW revised this behavior to ErrorExit when which doesn't match
  mht->which_vertices.)
  -----------------------------------------------------------------------*/
VERTEX * MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht, 
                                 MRI_SURFACE *mris, 
                                 VERTEX *v, 
                                 int which)
{
//------------------------------------------------------
  VERTEX * vtx;
  int rslt;
  float x=0.0, y=0.0, z=0.0;
  //---------------------------------
  // Sanity checks
  //---------------------------------
  mhtFindCommonSanityCheck(mht, mris);

  if (mht->which_vertices != which )
    ErrorExit(ERROR_BADPARM,
              "%s called with mismatched 'which' parameter\n",
              __MYFUNCTION__) ;

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_float(v, mht->which_vertices, &x, &y, &z);

  rslt = MHTfindClosestVertexGeneric(mht, mris,
                                     x, y, z,
                                     1000,
                                     1, // max_mhts: search out to 3 x 3 x 3
                                     &vtx, NULL, NULL);
  return vtx;
}

/*--------------------------------------------------------------------
  MHTfindClosestVertexNo()
  Returns vertex number and distance from vertex v to closest vertex
  in mris & mht.
  --------------------------------------------------------------------*/
int MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht, 
                           MRI_SURFACE *mris, 
                           VERTEX *v, 
                           float *min_dist)
{
//------------------------------------------------------
  int vtxnum, rslt;
  double x=0.0, y=0.0, z=0.0, min_dist_dbl;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  mhtVertex2xyz_double(v, mht->which_vertices, &x, &y, &z);

  rslt = MHTfindClosestVertexGeneric(mht, mris,
                                     x, y, z,
                                     1000,
                                     1,  // max_mhts: search out to 3 x 3 x 3
                                     NULL, &vtxnum, &min_dist_dbl);

  *min_dist = min_dist_dbl;
  //----------------------------------------------------------
  // [GW] Fixup for "no vertex found". V1.27 function returns
  // zero for no vertex found. I think this is a bug because
  // there is a legit vertex zero I think.
  // However, for now just duplicate existing output
  //----------------------------------------------------------
  if (-1 == rslt)
    rslt = 0;
  return vtxnum;
}

/*---------------------------------------------------------------
  MHTfindClosestVertexInTable
  Returns vertex from mris & mht that's closest to provided coordinates.
  ---------------------------------------------------------------*/
VERTEX * MHTfindClosestVertexInTable(MRIS_HASH_TABLE *mht, 
                                     MRI_SURFACE *mris, 
                                     float x, float y, float z)
{
//------------------------------------------------------
  VERTEX * vtx;
  int rslt;

  mhtFindCommonSanityCheck(mht, mris);

  //---------------------------------
  // Generic find
  //---------------------------------
  rslt = MHTfindClosestVertexGeneric(mht, mris,
                                     x, y, z,
                                     1000,
                                     1, // max_mhts: search out to 3 x 3 x 3
                                     &vtx, NULL, NULL);
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
int * MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht,
                                      MRI_SURFACE *mris,
                                      int vno, float max_dist, int *pvnum)
{
//------------------------------------------------------
  ErrorExit(ERROR_UNSUPPORTED,  "%s: not implemented.\n", __MYFUNCTION__);

  return NULL;
}

//=================================================================
//  Utility
//=================================================================

/*---------------------------------------------------
  MHTfree
  Frees and nulls an MHT pointer.
  ----------------------------------------------------*/
//----------------------------------
int MHTfree(MRIS_HASH_TABLE **pmht)
//----------------------------------
{
  MRIS_HASH_TABLE  *mht ;
  int              xv, yv, zv ;

  if (!(*pmht)) // avoid crash if not initialized, or nulled previously
    return(NO_ERROR) ;

  mht = *pmht ;
  *pmht = NULL ;  // sets pointer to null to signal free'ed

  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
      if (!mht->buckets[xv][yv])
        continue ;
      for (zv = 0 ; zv < TABLE_SIZE ; zv++)
      {
        if (mht->buckets[xv][yv][zv])
        {
          if (mht->buckets[xv][yv][zv]->bins)
            free(mht->buckets[xv][yv][zv]->bins) ;
          free(mht->buckets[xv][yv][zv]) ;
        }
      }
      free(mht->buckets[xv][yv]) ;
    }
  }
  free(mht) ;
  return(NO_ERROR) ;
}

//---------------------------------------------
static void mhtVertex2xyz_float(VERTEX * vtx,
                                int which,
                                float *x, float *y, float *z)
{
//---------------------------------------------
  if (!vtx) return;

  switch (which)
  {
  case ORIGINAL_VERTICES  :
    *x = vtx->origx;
    *y = vtx->origy;
    *z = vtx->origz;
    break;
  case GOOD_VERTICES      :
    break;
  case TMP_VERTICES       :
    break;
  case CANONICAL_VERTICES :
    *x = vtx->cx;
    *y = vtx->cy;
    *z = vtx->cz;
    break;
  case CURRENT_VERTICES   :
    *x = vtx->x;
    *y = vtx->y;
    *z = vtx->z;
    break;
  case INFLATED_VERTICES  :
    break;
  case FLATTENED_VERTICES :
    break;
  case PIAL_VERTICES      :
    *x = vtx->pialx;
    *y = vtx->pialy;
    *z = vtx->pialz;
    break;
  case TMP2_VERTICES      :
    break;
  case WHITE_VERTICES     :
    *x = vtx->whitex;
    *y = vtx->whitey;
    *z = vtx->whitez;
    break;
  }
  return;
}
//---------------------------------------------
static void mhtVertex2xyz_double(VERTEX * vtx,
                                 int which,
                                 double *x, double *y, double *z)
{
//---------------------------------------------
  if (!vtx) return;

  switch (which)
  {
  case ORIGINAL_VERTICES  :
    *x = vtx->origx;
    *y = vtx->origy;
    *z = vtx->origz;
    break;
  case GOOD_VERTICES      :
    break;
  case TMP_VERTICES       :
    break;
  case CANONICAL_VERTICES :
    *x = vtx->cx;
    *y = vtx->cy;
    *z = vtx->cz;
    break;
  case CURRENT_VERTICES   :
    *x = vtx->x;
    *y = vtx->y;
    *z = vtx->z;
    break;
  case INFLATED_VERTICES  :
    break;
  case FLATTENED_VERTICES :
    break;
  case PIAL_VERTICES      :
    *x = vtx->pialx;
    *y = vtx->pialy;
    *z = vtx->pialz;
    break;
  case TMP2_VERTICES      :
    break;
  case WHITE_VERTICES     :
    *x = vtx->whitex;
    *y = vtx->whitey;
    *z = vtx->whitez;
    break;
  }
  return;
}
//---------------------------------------------
static void mhtVertex2Ptxyz_double(VERTEX * vtx, int which, Ptdbl_t *pt)
{
//---------------------------------------------
  if ((!vtx) || (!pt)) return;

  switch (which)
  {
  case ORIGINAL_VERTICES  :
    pt->x = vtx->origx;
    pt->y = vtx->origy;
    pt->z = vtx->origz;
    break;
  case GOOD_VERTICES      :
    break;
  case TMP_VERTICES       :
    break;
  case CANONICAL_VERTICES :
    pt->x = vtx->cx;
    pt->y = vtx->cy;
    pt->z = vtx->cz;
    break;
  case CURRENT_VERTICES   :
    pt->x = vtx->x;
    pt->y = vtx->y;
    pt->z = vtx->z;
    break;
  case INFLATED_VERTICES  :
    break;
  case FLATTENED_VERTICES :
    break;
  case PIAL_VERTICES      :
    pt->x = vtx->pialx;
    pt->y = vtx->pialy;
    pt->z = vtx->pialz;
    break;
  case TMP2_VERTICES      :
    break;
  case WHITE_VERTICES     :
    pt->x = vtx->whitex;
    pt->y = vtx->whitey;
    pt->z = vtx->whitez;
    break;
  }
  return;
}

//---------------------------------------------
static void mhtVertex2array3_double(VERTEX * vtx, int which, double *array3)
{
//---------------------------------------------
  if (!vtx) return;

  switch (which)
  {
  case ORIGINAL_VERTICES  :
    array3[0] = vtx->origx;
    array3[1] = vtx->origy;
    array3[2] = vtx->origz;
    break;
  case GOOD_VERTICES      :
    break;
  case TMP_VERTICES       :
    break;
  case CANONICAL_VERTICES :
    array3[0] = vtx->cx;
    array3[1] = vtx->cy;
    array3[2] = vtx->cz;
    break;
  case CURRENT_VERTICES   :
    array3[0] = vtx->x;
    array3[1] = vtx->y;
    array3[2] = vtx->z;
    break;
  case INFLATED_VERTICES  :
    break;
  case FLATTENED_VERTICES :
    break;
  case PIAL_VERTICES      :
    array3[0] = vtx->pialx;
    array3[1] = vtx->pialy;
    array3[2] = vtx->pialz;
    break;
  case TMP2_VERTICES      :
    break;
  case WHITE_VERTICES     :
    array3[0] = vtx->whitex;
    array3[1] = vtx->whitey;
    array3[2] = vtx->whitez;
    break;
  }
  return;
}

/*-----------------------------------------------------------------*/
MHBT * MHTgetBucket(MRIS_HASH_TABLE *mht, float x, float y, float z)
{
//-------------------------------------------------------------------
  int     xv, yv, zv ;

  xv = WORLD_TO_VOXEL(mht, x) ;
  yv = WORLD_TO_VOXEL(mht, y) ;
  zv = WORLD_TO_VOXEL(mht, z) ;

  return MHTgetBucketAtVoxIx(mht, xv, yv, zv);
}

/*-----------------------------------------------------------------*/
MHBT * MHTgetBucketAtVoxIx(MRIS_HASH_TABLE *mht, int xv, int yv, int zv)
{
//-------------------------------------------------------------------
  MHBT    *bucket ;

  if (  xv >= FIELD_OF_VIEW
        || yv >= FIELD_OF_VIEW
        || zv >= FIELD_OF_VIEW
        || xv <  0
        || yv <  0
        || zv <  0)
    return(NULL);

  if (!mht)
    return(NULL);

  if (!mht->buckets[xv][yv])
    return(NULL) ;

  bucket = mht->buckets[xv][yv][zv] ;
  return(bucket) ;
}

/*------------------------------------------------
  MHT_gw_version
  Confidence check that correct version of code is
  compiled in;
  ------------------------------------------------*/
int MHT_gw_version(void)
{
//-------------------------------
  return GW_VERSION;   // <-- change this as needed
}

//=================================================================
// VOXEL_LIST
//=================================================================

//--------------------------------------------------
static int mhtVoxelList_Init(VOXEL_LISTgw *voxlist)
{
//--------------------------------------------------
  voxlist->nused = 0 ;
  return(NO_ERROR) ;
}

//--------------------------------------------------
static int mhtVoxelList_Add(VOXEL_LISTgw *voxlist,
                            int xv, int yv, int zv, int fno)
{
//--------------------------------------------------
  int   i ;

  for (i = 0 ; i < voxlist->nused ; i++)
    if (voxlist->voxels[i][0] == xv &&
        voxlist->voxels[i][1] == yv &&
        voxlist->voxels[i][2] == zv)
      return(NO_ERROR) ;

  if (voxlist->nused >= MAX_VOXELS)
  {
    fprintf(stderr,
            "%s(%d, %d, %d, %d): complete list too big!\n",
            __MYFUNCTION__ ,
            xv,yv,zv, fno) ;

    ErrorPrintf(ERROR_NOMEMORY,
                "%s(%d, %d, %d, %d): complete list too big!",
                __MYFUNCTION__,
                xv,yv,zv, fno);
    return(ERROR_NOMEMORY);
  }
  //-----------------------------
  // Actually append it
  //-----------------------------
  i = voxlist->nused;
  voxlist->voxels[i][0] = xv ;
  voxlist->voxels[i][1] = yv ;
  voxlist->voxels[i][2] = zv ;
  voxlist->nused++ ;
  return(NO_ERROR) ;
}

//--------------------------------------------------
static int mhtVoxelList_AddCoord(VOXEL_LISTgw *voxlist,
                                 VOXEL_COORD vc, int fno)
{
//--------------------------------------------------
  return mhtVoxelList_Add(voxlist, vc.xv, vc.yv, vc.zv, fno);
}

//--------------------------------------------------
static int mhtVoxelList_AddPath(VOXEL_LISTgw *voxlist,
                                VOXEL_COORD oldvc, VOXEL_COORD newvc, int fno )
{
//--------------------------------------------------
  int xv, yv, zv, count;
  int incx = -1, incy = -1, incz = -1;

  if (newvc.xv >= oldvc.xv) incx = 1;
  if (newvc.yv >= oldvc.yv) incy = 1;
  if (newvc.zv >= oldvc.zv) incz = 1;

  count = -1;
  for (xv     = oldvc.xv;
       incx==1 ? xv <= newvc.xv: xv >= newvc.xv;
       xv += incx)
  {
    for (yv   = oldvc.yv;
         incy==1 ? yv <= newvc.yv: yv >= newvc.yv;
         yv += incy)
    {
      for (zv = oldvc.zv;
           incz==1 ? zv <= newvc.zv: zv >= newvc.zv;
           zv += incz)
      {
        count++;
        // ie: ignore the start voxel, assume that's already added
        // (which can be done with mhtVoxelList_AddCoord)
        if (count >= 1)
        {
          mhtVoxelList_Add(voxlist, xv, yv, zv, fno);
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
int MHTtestIsMRISselfIntersecting(MRI_SURFACE *mris, float res)
//--------------------------------------------------------
{
  MRIS_HASH_TABLE  *mht ;
  int fno, rslt ;

  rslt = 0;

//mht = MHTfillTable(mris, NULL) ;

  mht = MHTfillTableAtResolution(mris, NULL, CURRENT_VERTICES, res);
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (MHTdoesFaceIntersect(mht, mris, fno))
    {
      rslt = 1;
      goto done;
    }
  }
 done:
  MHTfree(&mht) ;
  return rslt ;
}


/*--------------------------------------------------------
  MHTcheckFaces
  Unchanged from mrishash.c 1.27
  -------------------------------------------------------*/
int MHTcheckFaces(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht)
//--------------------------------------------------------
{
  static int ncalls = 0 ;

  if (ncalls++ >= 0)
  {
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
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
  MHTcheckSurface
  Unchanged from mrishash.c 1.27
  However, subsidiary checkFace is disabled
  ------------------------------------------------------*/
int
MHTcheckSurface(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht)
//--------------------------------------------------------
{
  int    fno, alloced = 0 ;

  if (!mht)
  {
    mht = MHTfillTable(mris, mht) ;
    alloced = 1 ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (!(fno % (mris->nfaces/10)))
      DiagHeartbeat((float)fno / (float)mris->nfaces) ;

    checkFace(mht, mris, fno) ;
  }
  if (alloced)
    MHTfree(&mht) ;
  return(NO_ERROR) ;
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
static int
checkFace(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno1)
//--------------------------------------------------------
{
#if 0
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;
  int    fno2, filled, n1, n2, nbr ;
  FACE   *f1, *f2 ;

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
      VERTEX *v ;

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
  return(NO_ERROR) ;
}
