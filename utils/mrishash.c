#include <stdio.h>
#include <math.h>

#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "fio.h"
#include "mrisurf.h"
#include "proto.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "tritri.h"

/*---------------------------- STRUCTURES -----------------------------*/

#define MAX_VOXELS  10000  
typedef struct
{
  int nused ;
  int v[MAX_VOXELS][3] ;
} VOXEL_LIST ;

/*---------------------------- CONSTANTS ------------------------------*/

/*------------------------ STATIC PROTOTYPES --------------------------*/

static int mhtInitVoxelList(VOXEL_LIST *vl) ;
static int mhtAddVoxels(MRIS_HASH_TABLE *mht, VOXEL_LIST *vl, 
                        float xw, float yw, float zw) ;
static int mhtAddVoxel(VOXEL_LIST *vl, int xv, int yv, int zv) ;
static int mhtDoesFaceVoxelListIntersect(MRIS_HASH_TABLE *mht, 
                                         MRI_SURFACE *mris,
                                         VOXEL_LIST *vl, int fno) ;

static int mhtHatchFace(MRIS_HASH_TABLE *mht,MRI_SURFACE *mris,int fno,int on);

/*static int hash(MRIS_HASH_TABLE *mht, int z) ;*/
static int mhtAddFacePositions(MRIS_HASH_TABLE *mht, float x, float y, 
                              float z, int fno) ;
static int mhtRemoveFacePositions(MRIS_HASH_TABLE *mht, float x, float y, 
                                 float z, int fno) ;
static int mhtAddFacePosition(MRIS_HASH_TABLE *mht, int xv, int yv, 
                              int zv, int fno) ;
static int mhtRemoveFacePosition(MRIS_HASH_TABLE *mht, int xv, int yv, 
                                 int zv, int fno) ;
#if 0
static int mhtDoesFaceIntersect(MRIS_HASH_TABLE *mht,
                                MRI_SURFACE *mris,int fno);
#endif
#if 0
static int checkAllVertexFaces(MRIS_HASH_TABLE *mht,MRI_SURFACE *mris,int vno);
#endif
static int checkFace(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno1) ;


/*---------------------------- MACROS ---------------------------------*/


/*-------------------------- STATIC DATA ------------------------------*/

/*-------------------------- FUNCTIONS --------------------------------*/


MRIS_HASH_TABLE *
MHTfillVertexTable(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht, int which)
{
  return(MHTfillVertexTableRes(mris, mht, which, VOXEL_RES)) ;
}
MRIS_HASH_TABLE *
MHTfillVertexTableRes(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht, int which,
                      float res)
{
  int     vno ;
  int     xv, yv, zv ;
  MHBT    *bucket ;
  VERTEX  *v ;
  static int ncalls = 0 ;
  
  if (mht)    /* free old one */
    MHTfree(&mht) ;
  mht = (MRIS_HASH_TABLE *)calloc(1, sizeof(MRIS_HASH_TABLE)) ;
  if (!mht)
    ErrorExit(ERROR_NO_MEMORY, 
              "MHTfillVertexTable: could not allocate hash table.\n") ;

  mht->which_vertices = which ;
  mht->vres = res ;
  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
      mht->buckets[xv][yv] = NULL ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    switch (which)
    {
    case ORIGINAL_VERTICES:
      mhtAddFacePositions(mht, v->origx, v->origy, v->origz, vno) ;
      break ;
    case CANONICAL_VERTICES:
      mhtAddFacePositions(mht, v->cx, v->cy, v->cz, vno) ;
      break ;
    case CURRENT_VERTICES:
      mhtAddFacePositions(mht, v->x, v->y, v->z, vno) ;
      break ;
    case WHITE_VERTICES:
      mhtAddFacePositions(mht, v->whitex, v->whitey, v->whitez, vno) ;
      break ;
    case PIAL_VERTICES:
      mhtAddFacePositions(mht, v->pialx, v->pialy, v->pialz, vno) ;
      break ;
    }
  }

  if ((Gdiag & DIAG_SHOW) && !ncalls)
  {
    double mean, var, v, n ;
    int    mx ;

    n = mean = 0.0 ; mx = -1 ;
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
      fprintf(stderr, "buckets: mean = %2.1f +- %2.2f, max = %d\n",
              mean, sqrt(var), mx) ;
  }
  ncalls++ ;
  return(mht) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS_HASH_TABLE *
MHTfillTable(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht)
{
  return(MHTfillTableAtResolution(mris, mht, CURRENT_VERTICES, VOXEL_RES)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIS_HASH_TABLE *
MHTfillTableAtResolution(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht, int which,
                         float res)
{
  int     fno ;
  FACE    *f ;
  int     xv, yv, zv ;
  MHBT    *bucket ;
  static int ncalls = 0 ;

#if 0
  int     n, vno ;
  VERTEX  *v, *vn ;
  double  total_dist, nv, dist ;
#endif

  if (mht)    /* free old one */
    MHTfree(&mht) ;
  mht = (MRIS_HASH_TABLE *)calloc(1, sizeof(MRIS_HASH_TABLE)) ;
  if (!mht)
    ErrorExit(ERROR_NO_MEMORY, 
              "MHTfillTable: could not allocate hash table.\n") ;

  mht->which_vertices = which ;

  mht->vres = res ;
  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
#if 0
      for (zv = 0 ; zv < TABLE_SIZE ; zv++)
      {
        bucket = &mht->buckets[xv][yv][zv] ;
        bucket->nused = bucket->max_bins = 0 ;
        bucket->bins = NULL ;
      }
#else
      mht->buckets[xv][yv] = NULL ;
#endif
    }
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    mhtHatchFace(mht, mris, fno, 1) ;
  }

#if 0
  if ((Gdiag & DIAG_SHOW) && !(ncalls%10))
  {
    for (nv = total_dist = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        nv++ ;
        dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
        total_dist += dist ;
      }
    }
    total_dist /= nv ;
    fprintf(stderr, "average voxel spacing = %2.2f\n", total_dist) ;
  }
#endif
  if ((Gdiag & DIAG_SHOW) && !ncalls)
  {
    double mean, var, v, n ;
    int    mx ;

    n = mean = 0.0 ; mx = -1 ;
    for (xv = 0 ; xv < TABLE_SIZE ; xv++)
    {
      for (yv = 0 ; yv < TABLE_SIZE ; yv++)
      {
        for (zv = 0 ; zv < TABLE_SIZE ; zv++)
        {
#if 0
          bucket = &mht->buckets[xv][yv][zv] ;
#else
          if (!mht->buckets[xv][yv] || !mht->buckets[xv][yv][zv])
            continue ;
          bucket = mht->buckets[xv][yv][zv] ;
#endif
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
#if 0
          bucket = &mht->buckets[xv][yv][zv] ;
          if (bucket->nused)
#else
          bucket = mht->buckets[xv][yv][zv] ;
          if (bucket && bucket->nused)
#endif
          {
            v = mean - bucket->nused ;
            var += v*v ;
          }
        }
      }
    }
    var /= (n-1) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "buckets: mean = %2.1f +- %2.2f, max = %d\n",
              mean, sqrt(var), mx) ;
  }
  ncalls++ ;
  return(mht) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           check whether the position (xw,yw,zw) is filled by
           a face which shares no vertices with fno.
------------------------------------------------------*/
#if 0
int
MHTisFilled(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno,
            float xw, float yw, float zw)
{
  int   xv, yv, zv, fno2, i, n1, n2, intersect, nbr ;
  MHB   *bin ;
  MHBT  *bucket ;
  FACE  *f1, *f2 ;

  f1 = &mris->faces[fno] ;
  xv = WORLD_TO_VOXEL(mht,xw) ; 
  yv = WORLD_TO_VOXEL(mht,yw) ; 
  zv = WORLD_TO_VOXEL(mht,zw) ;

  if (xv < 0)
    xv = 0 ;
  if (xv >= TABLE_SIZE)
    xv = TABLE_SIZE-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= TABLE_SIZE)
    yv = TABLE_SIZE-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= TABLE_SIZE)
    zv = TABLE_SIZE-1 ;
  if (fno == Gdiag_no)
    DiagBreak() ;
  if (xv < 0 || yv < 0 || zv < 0 ||
      xv >= TABLE_SIZE || yv >= TABLE_SIZE || zv >= TABLE_SIZE ||
      fabs(xw) > FIELD_OF_VIEW ||
      fabs(yw) > FIELD_OF_VIEW ||
      fabs(zw) > FIELD_OF_VIEW)
  {
    fprintf(stderr, 
         "MHTisFilled: index out of bounds at face %d (%2.1f, %2.1f, %2.1f)\n",
          fno, xw, yw, zw) ;
    DiagBreak() ;
  }
  /*
    go through every entry in this bucket and see if any of them
    occupy the same voxel coordinate as the vertex in question. If so,
    check the face in the bucket to see if the vertex is part of if. If it
    is not, then the voxel is already filled and return 1.
  */
#if 0
  bucket = &mht->buckets[xv][yv][zv] ;
#else
  bucket = mht->buckets[xv][yv][zv] ;
#endif
  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    fno2 = bin->fno ; f2 = &mris->faces[fno2] ;
    nbr = 0 ;
    for (n1 = 0 ; !nbr && n1 < VERTICES_PER_FACE ; n1++)
      for (n2 = 0 ; !nbr && n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f1->v[n1] == f2->v[n2])
          nbr = 1 ;  /* they share a vertex - don't count it as filled */
      }
    if (!nbr)   /* use explicit intersection check */
    {
      double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;
      
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
      
      /* fill vertices of 2nd triangle */
      u0[0] = (double)mris->vertices[f2->v[0]].x ;
      u0[1] = (double)mris->vertices[f2->v[0]].y ;
      u0[2] = (double)mris->vertices[f2->v[0]].z ;
      u1[0] = (double)mris->vertices[f2->v[1]].x ;
      u1[1] = (double)mris->vertices[f2->v[1]].y ;
      u1[2] = (double)mris->vertices[f2->v[1]].z ;
      u2[0] = (double)mris->vertices[f2->v[2]].x ;
      u2[1] = (double)mris->vertices[f2->v[2]].y ;
      u2[2] = (double)mris->vertices[f2->v[2]].z ;
      intersect = tri_tri_intersect(v0,v1,v2,u0,u1,u2) ;
      if (intersect)
        return(1) ;
    }
  }
  return(0) ;
}
#endif
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           check whether 
------------------------------------------------------*/
int
MHTisVoxelFilled(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int vno,
            int xv, int yv, int zv)
{
  int   fno, i, n ;
  MHB   *bin ;
  MHBT  *bucket ;
  FACE  *f ;


  /*
    go through every entry in this bucket and see if any of them
    occupy the same voxel coordinate as the vertex in question. If so,
    check the face in the bucket to see if the vertex is part of if. If it
    is not, then the voxel is already filled and return 1.
  */
  bucket = &mht->buckets[xv][yv][zv] ;
  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    /* check to see if this face has vno in it */
    fno = bin->fno ; f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      if (f->v[n] == vno)   /* vertex is part of this face - ignore it */
        break ;
    }
    if (n >= VERTICES_PER_FACE)  /* vertex not in face - return filled */
      return(1) ;
  }

  return(0) ;   /* position is not filled with a face that v is not part of */
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTaddAllFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v)
{
  int     fno ;

  for (fno = 0 ; fno < v->num ; fno++)
    mhtHatchFace(mht, mris, v->f[fno], 1) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTremoveAllFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v)
{
  int     fno ;

  for (fno = 0 ; fno < v->num ; fno++)
    mhtHatchFace(mht, mris, v->f[fno], 0) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
        Description
        each face has 2 triangles defined by it:

       V0    b     V2
        o----------o
        |        /       
        |      /         
      a |    /            
        |  /             
        |/      
        o
       V1      b        V2        
------------------------------------------------------*/
#define SQRT_2        (1.4142136)
#define SAMPLE_DIST   (mht->vres/(2.0*SQRT_2))

static int
mhtHatchFace(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno, int on)
{
  Real   x, y, z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, dz, 
         cdx, cdy, cdz, alen, clen, delta_t0, delta_t1, len ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *face ;

  face = &mris->faces[fno] ;
  if (face->ripflag)
    return(NO_ERROR) ;
  
  v0 = &mris->vertices[face->v[0]] ;
  v1 = &mris->vertices[face->v[1]] ;
  v2 = &mris->vertices[face->v[2]] ;
  switch (mht->which_vertices)
  {
  case ORIGINAL_VERTICES:
    adx = v1->origx - v0->origx ; ady = v1->origy - v0->origy ; 
    adz = v1->origz - v0->origz ;
    cdx = v2->origx - v0->origx ; cdy = v2->origy - v0->origy ; 
    cdz = v2->origz - v0->origz ;
    break ;
  case CANONICAL_VERTICES:
    adx = v1->cx - v0->cx ; ady = v1->cy - v0->cy ; adz = v1->cz - v0->cz ;
    cdx = v2->cx - v0->cx ; cdy = v2->cy - v0->cy ; cdz = v2->cz - v0->cz ;
    break ;
  default:
  case CURRENT_VERTICES:
    adx = v1->x - v0->x ; ady = v1->y - v0->y ; adz = v1->z - v0->z ;
    cdx = v2->x - v0->x ; cdy = v2->y - v0->y ; cdz = v2->z - v0->z ;
    break ;
  case WHITE_VERTICES:
    adx = v1->whitex - v0->whitex ; ady = v1->whitey - v0->whitey ; 
		adz = v1->whitez - v0->whitez ;
    cdx = v2->whitex - v0->whitex ; cdy = v2->whitey - v0->whitey ; 
		cdz = v2->whitez - v0->whitez ;
    break ;
  case PIAL_VERTICES:
    adx = v1->pialx - v0->pialx ; ady = v1->pialy - v0->pialy ; adz = v1->pialz - v0->pialz ;
    cdx = v2->pialx - v0->pialx ; cdy = v2->pialy - v0->pialy ; cdz = v2->pialz - v0->pialz ;
    break ;
  }
  alen = sqrt(SQR(adx)+SQR(ady)+SQR(adz)) ;
  clen = sqrt(SQR(cdx)+SQR(cdy)+SQR(cdz)) ;
  
  /*
    sample along legs of the triangle making sure the maximum spacing
    between samples (along the longer leg) is SAMPLE_DIST.
  */
  
  /*
    move along v0->v1 and v3->v2 lines and draw in crossing line to fill face
    t0 parameterizes lines from v0->v1 and v0->v2 
  */
  if (FZERO(alen) && FZERO(clen))
    delta_t0 = 0.99 ;
  else
    delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen ) ;
  if (FZERO(delta_t0))
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,
                 "mrisFillFace: face %d has infinite leg (%d, %d)\n",
                 fno, alen, clen)) ;
  
  if (delta_t0 >= 1.0)
    delta_t0 = 0.99 ;
  
  /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
  for (t0 = 0 ; t0 <= 1.0f ; t0 += delta_t0)
    {
      /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
      switch (mht->which_vertices)
  {
  default:
  case CURRENT_VERTICES:
    xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
    break ;
  case ORIGINAL_VERTICES:
    xa = v0->origx + t0*adx ; ya = v0->origy + t0*ady ; 
    za = v0->origz + t0*adz ;
    xc = v0->origx + t0*cdx ; yc = v0->origy + t0*cdy ; 
    zc = v0->origz + t0*cdz ;
    break ;
  case CANONICAL_VERTICES:
    xa = v0->cx + t0*adx ; ya = v0->cy + t0*ady ; za = v0->cz + t0*adz ;
    xc = v0->cx + t0*cdx ; yc = v0->cy + t0*cdy ; zc = v0->cz + t0*cdz ;
    break ;
  }      
    dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    if (FZERO(len))
      delta_t1 = 0.99 ;
    else
    {
      delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f)
        delta_t1 = 0.99 ;
    }
    
    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
    {
      /* compute a point on the line connecting a and c */
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      if (on)
        mhtAddFacePositions(mht, x, y, z, fno) ;
      else
        mhtRemoveFacePositions(mht, x, y, z, fno) ;

    }
    /* compute last point on line */
    t1 = 1.0f ;
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    if (on)
      mhtAddFacePositions(mht, x, y, z, fno) ;
    else
      mhtRemoveFacePositions(mht, x, y, z, fno) ;
  }
  
  /* compute last line on the a and c lines resp. */
  t0 = 1.0f ;
  switch (mht->which_vertices)
  {
  default:
  case CURRENT_VERTICES:
    xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
    break ;
  case ORIGINAL_VERTICES:
    xa = v0->origx + t0*adx ; ya = v0->origy + t0*ady; za = v0->origz + t0*adz;
    xc = v0->origx + t0*cdx ; yc = v0->origy + t0*cdy; zc = v0->origz + t0*cdz;
    break ;
  case CANONICAL_VERTICES:
    xa = v0->cx + t0*adx ; ya = v0->cy + t0*ady ; za = v0->cz + t0*adz ;
    xc = v0->cx + t0*cdx ; yc = v0->cy + t0*cdy ; zc = v0->cz + t0*cdz ;
    break ;
  case WHITE_VERTICES:
    xa = v0->whitex + t0*adx ; ya = v0->whitey + t0*ady ; za = v0->whitez + t0*adz ;
    xc = v0->whitex + t0*cdx ; yc = v0->whitey + t0*cdy ; zc = v0->whitez + t0*cdz ;
    break ;
  case PIAL_VERTICES:
    xa = v0->pialx + t0*adx ; ya = v0->pialy + t0*ady ; za = v0->pialz + t0*adz ;
    xc = v0->pialx + t0*cdx ; yc = v0->pialy + t0*cdy ; zc = v0->pialz + t0*cdz ;
    break ;
  }
  dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
  len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
  if (FZERO(len))
    delta_t1 = 0.99 ;
  else
  {
    delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
    if (delta_t1 >= 1.0f)
      delta_t1 = 0.99 ;
  }
  
  /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
  for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
  {
    /* compute a point on the line connecting a and c */
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    if (on)
      mhtAddFacePositions(mht, x, y, z, fno) ;
    else
      mhtRemoveFacePositions(mht, x, y, z, fno) ;
  }
  /* compute last point on line */
  t1 = 1.0f ;
  x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
  if (on)
    mhtAddFacePositions(mht, x, y, z, fno) ;
  else
    mhtRemoveFacePositions(mht, x, y, z, fno) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtAddFacePositions(MRIS_HASH_TABLE *mht,float xw, float yw, float zw, int fno)
{
  int    xv, yv, zv ;
  float  x, y, z ;
  
  x = WORLD_TO_VOLUME(mht, xw) ; 
  y = WORLD_TO_VOLUME(mht, yw) ; 
  z = WORLD_TO_VOLUME(mht, zw) ;
  xv = (int)x ; yv = (int)y ; zv = (int)z ;
  mhtAddFacePosition(mht, xv, yv, zv, fno) ;

  if ((ceil(x) - x) < 0.5)
    mhtAddFacePosition(mht, xv+1, yv, zv, fno) ;
  else
    mhtAddFacePosition(mht, xv-1, yv, zv, fno) ;

  if ((ceil(y) - y) < 0.5)
    mhtAddFacePosition(mht, xv, yv+1, zv, fno) ;
  else
    mhtAddFacePosition(mht, xv, yv-1, zv, fno) ;
    
  if ((ceil(z) - z) < 0.5)
    mhtAddFacePosition(mht, xv, yv, zv+1, fno) ;
  else
    mhtAddFacePosition(mht, xv, yv, zv-1, fno) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtRemoveFacePositions(MRIS_HASH_TABLE *mht,float xw,float yw,float zw,int fno)
{
  int    xv, yv, zv ;
  float  x, y, z ;
  
  x = WORLD_TO_VOLUME(mht, xw) ; 
  y = WORLD_TO_VOLUME(mht, yw) ; 
  z = WORLD_TO_VOLUME(mht, zw) ;
  xv = (int)x ; yv = (int)y ; zv = (int)z ;
  mhtRemoveFacePosition(mht, xv, yv, zv, fno) ;

  if ((ceil(x) - x) < 0.5)
    mhtRemoveFacePosition(mht, xv+1, yv, zv, fno) ;
  if ((x - floor(x)) < 0.5)
    mhtRemoveFacePosition(mht, xv-1, yv, zv, fno) ;

  if ((ceil(y) - y) < 0.5)
    mhtRemoveFacePosition(mht, xv, yv+1, zv, fno) ;
  if ((y - floor(y)) < 0.5)
    mhtRemoveFacePosition(mht, xv, yv-1, zv, fno) ;
    
  if ((ceil(z) - z) < 0.5)
    mhtRemoveFacePosition(mht, xv, yv, zv+1, fno) ;
  if ((z - floor(z)) < 0.5)
    mhtRemoveFacePosition(mht, xv, yv, zv-1, fno) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtAddFacePosition(MRIS_HASH_TABLE *mht, int xv, int yv, int zv, int fno)
{
  int    i ;
  MHB   *bin ;
  MHBT  *bucket ;

  if (xv < 0)
    xv = 0 ;
  if (xv >= TABLE_SIZE)
    xv = TABLE_SIZE-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= TABLE_SIZE)
    yv = TABLE_SIZE-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= TABLE_SIZE)
    zv = TABLE_SIZE-1 ;
  if (xv < 0 || yv < 0 || zv < 0 ||
      xv >= TABLE_SIZE || yv >= TABLE_SIZE || zv >= TABLE_SIZE)
  {
    fprintf(stderr, 
            "mhtAddFacePosition: index out of bounds at face %d"
            " (%d, %d, %d)\n",
            fno, xv, yv, zv) ;
    DiagBreak() ;
  }
  if (!mht->buckets[xv][yv])
  {
    mht->buckets[xv][yv] = (MHBT **)calloc(TABLE_SIZE, sizeof(MHBT *)) ;
    if (!mht->buckets[xv][yv])
      ErrorExit(ERROR_NO_MEMORY, 
                "mhtAddFacePosition: could not allocate slice.") ;
  }
  bucket = mht->buckets[xv][yv][zv] ;
  if (!bucket)
  {
    mht->buckets[xv][yv][zv] = bucket = (MHBT *)calloc(1, sizeof(MHBT)) ;
    if (!bucket)
      ErrorExit(ERROR_NOMEMORY, "couldn't allocate bucket.\n") ;
  }
  if (!bucket->max_bins)   /* nothing in this bucket yet - allocate bins */
  {
    bucket->max_bins = 4 ;
    bucket->bins = (MHB *)calloc(bucket->max_bins, sizeof(MHB));
    if (!bucket->bins)
      ErrorExit(ERROR_NO_MEMORY, 
                "mhtAddFacePosition: could not allocate %d bins.\n",
                bucket->max_bins) ;
  }

  /* see if this face already is listed as occupying this position in this
     bin. If so, don't add it again, otherwise add it to the list.
  */
  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    if (bin->fno == fno)
      break ;
  }
  if (i >= bucket->nused)   /* face not already listed at this position */
  {
    if (bucket->nused >= bucket->max_bins)   /* allocate some more */
    {
      bin = bucket->bins ;
      bucket->max_bins *= 2 ;
      bucket->bins = (MHB *)calloc(bucket->max_bins, sizeof(MHB));
      if (!bucket->bins)
        ErrorExit(ERROR_NO_MEMORY, 
                  "mhtAddFacePosition: could not allocate %d bins.\n",
                  bucket->max_bins) ;
      memmove(bucket->bins, bin, bucket->nused*sizeof(MHB)) ;
      free(bin) ;
      bin = &bucket->bins[i] ;
    }
    /* add this face-position to this bucket */
    bucket->nused++ ;
    bin->fno = fno ; 
  }
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtRemoveFacePosition(MRIS_HASH_TABLE *mht, int xv, int yv, int zv,int fno)
{
  int    i ;
  MHBT   *bucket ;
  MHB    *bin ;

  if (xv < 0)
    xv = 0 ;
  if (xv >= TABLE_SIZE)
    xv = TABLE_SIZE-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= TABLE_SIZE)
    yv = TABLE_SIZE-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= TABLE_SIZE)
    zv = TABLE_SIZE-1 ;
  if (xv < 0 || yv < 0 || zv < 0 ||
      xv >= TABLE_SIZE || yv >= TABLE_SIZE || zv >= TABLE_SIZE)
  {
    fprintf(stderr, 
            "mhtRemoveFacePosition: index out of bounds at face %d"
            " (%d, %d, %df)\n",
            fno, xv, yv, zv) ;
    DiagBreak() ;
  }

#if 0
  bucket = &mht->buckets[xv][yv][zv] ;
#else
  if (!mht->buckets[xv][yv])
    return(NO_ERROR) ;
  bucket = mht->buckets[xv][yv][zv] ;
  if (!bucket)
    return(NO_ERROR) ;
#endif
  bin = bucket->bins ;
  for (i = 0 ; i < bucket->nused ; i++, bin++)
  {
    if (bin->fno == fno)     /* found face */
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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTisVectorFilled(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, 
                 int vno, float dx, float dy, float dz)
{
  VERTEX   *v ;
  float    ox, oy, oz ;
  int      intersect, fno ;

  v = &mris->vertices[vno] ;
  ox = v->x ; oy = v->y ; oz = v->z ;
  v->x += dx ; v->y += dy ; v->z += dz ;
  for (intersect = fno = 0 ; !intersect && fno < v->num ; fno++)
  {
    intersect = MHTdoesFaceIntersect(mht, mris, v->f[fno]) ;
  }
  v->x = ox ; v->y = oy ; v->z = oz ;
  return(intersect) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris,int fno)
{
  Real       x, y, z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, 
             dz, cdx, cdy, cdz, alen, clen, delta_t0, delta_t1, len ;
  VERTEX     *v0, *v1, *v2 ;
  FACE       *face ;
  VOXEL_LIST vl ;

  mhtInitVoxelList(&vl) ;
  face = &mris->faces[fno] ;
  if (face->ripflag)
    return(NO_ERROR) ;

  if (fno == Gdiag_no)
    DiagBreak() ;
  v0 = &mris->vertices[face->v[0]] ;
  v1 = &mris->vertices[face->v[1]] ;
  v2 = &mris->vertices[face->v[2]] ;
  adx = v1->x - v0->x ; ady = v1->y - v0->y ; adz = v1->z - v0->z ;
  alen = sqrt(SQR(adx)+SQR(ady)+SQR(adz)) ;
  cdx = v2->x - v0->x ; cdy = v2->y - v0->y ; cdz = v2->z - v0->z ;
  clen = sqrt(SQR(cdx)+SQR(cdy)+SQR(cdz)) ;
  
  /*
    sample along legs of the triangle making sure the maximum spacing
    between samples (along the longer leg) is SAMPLE_DIST.
  */
  
  /*
    move along v0->v1 and v3->v2 lines and draw in crossing line to fill face
    t0 parameterizes lines from v0->v1 and v0->v2 
  */
  if (FZERO(alen) && FZERO(clen))
    delta_t0 = 0.99 ;
  else
    delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen ) ;
  if (FZERO(delta_t0))
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,
                 "mrisFillFace: face %d has infinite leg (%d, %d)\n",
                 fno, alen, clen)) ;
  
  if (delta_t0 >= 1.0)
    delta_t0 = 0.99 ;
  
  /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
  for (t0 = 0 ; t0 <= 1.0f ; t0 += delta_t0)
  {
    /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
    xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
    xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
    dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
    len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
    if (FZERO(len))
      delta_t1 = 0.99 ;
    else
    {
      delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f)
        delta_t1 = 0.99 ;
    }
    
    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
    {
      /* compute a point on the line connecting a and c */
      x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
      mhtAddVoxels(mht, &vl, x, y, z) ;
    }
    /* compute last point on line */
    t1 = 1.0f ;
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    mhtAddVoxels(mht, &vl, x, y, z) ;
  }
  
  /* compute last line on the a and c lines resp. */
  t0 = 1.0f ;
  xa = v0->x + t0*adx ; ya = v0->y + t0*ady ; za = v0->z + t0*adz ;
  xc = v0->x + t0*cdx ; yc = v0->y + t0*cdy ; zc = v0->z + t0*cdz ;
  dx = xc-xa ; dy = yc-ya ; dz = zc-za ;
  len = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
  if (FZERO(len))
    delta_t1 = 0.99 ;
  else
  {
    delta_t1 = SAMPLE_DIST / len ;  /* sample at SAMPLE_DIST intervals */
    if (delta_t1 >= 1.0f)
      delta_t1 = 0.99 ;
  }
  
  /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
  for (t1 = 0 ; t1 <= 1.0f ; t1 += delta_t1)
  {
    /* compute a point on the line connecting a and c */
    x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
    mhtAddVoxels(mht, &vl, x, y, z) ;
  }
  /* compute last point on line */
  t1 = 1.0f ;
  x = xa + t1*dx ; y = ya + t1*dy ; z = za + t1*dz ;
  mhtAddVoxels(mht, &vl, x, y, z) ;


  if (mhtDoesFaceVoxelListIntersect(mht, mris, &vl, fno))
    return(1) ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
static int
checkAllVertexFaces(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int vno)
{
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;
  int    fno1, fno2, filled, n1, n2, nbr ;
  FACE   *f1, *f2 ;
  VERTEX *v ;
  
  v = &mris->vertices[vno] ;
  for (fno1 = 0 ; fno1 < v->num ; fno1++)
  {
    f1 = &mris->faces[v->f[fno1]] ;

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
      /* fill vertices of 2nd triangle */
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
      if (filled)
      {
        fprintf(stderr, "face %d intersects with face %d!!!!\n", 
                v->f[fno1], fno2) ;
        DiagBreak() ;
      }
    }
  }

  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTcheckFaces(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht)
{
  static int ncalls = 0 ;

  if (ncalls++ >= 0)
  {
    checkFace(mht, mris, 144) ;
    checkFace(mht, mris, 185) ;
#if 0
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
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTcheckSurface(MRI_SURFACE *mris,MRIS_HASH_TABLE *mht)
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
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MHTfree(MRIS_HASH_TABLE **pmht)
{
  MRIS_HASH_TABLE  *mht ;
  int              xv, yv, zv ;

  mht = *pmht ;
  *pmht = NULL ;

  for (xv = 0 ; xv < TABLE_SIZE ; xv++)
  {
    for (yv = 0 ; yv < TABLE_SIZE ; yv++)
    {
      if (!mht->buckets[xv][yv])
        continue ;
      for (zv = 0 ; zv < TABLE_SIZE ; zv++)
      {
#if 0
        if (mht->buckets[xv][yv][zv].bins)
#else
        if (mht->buckets[xv][yv][zv])
        {
          if (mht->buckets[xv][yv][zv]->bins)
            free(mht->buckets[xv][yv][zv]->bins) ;
          free(mht->buckets[xv][yv][zv]) ;
        }
#endif
      }
      free(mht->buckets[xv][yv]) ;
    }
  }
  free(mht) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
checkFace(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, int fno1)
{
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
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtInitVoxelList(VOXEL_LIST *vl)
{
  vl->nused = 0 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtAddVoxels(MRIS_HASH_TABLE *mht,VOXEL_LIST *vl, float xw, float yw, float zw)
{
  int   xv, yv, zv ;
  float x, y, z ;

  x = WORLD_TO_VOLUME(mht, xw) ; 
  y = WORLD_TO_VOLUME(mht, yw) ; 
  z = WORLD_TO_VOLUME(mht, zw) ;
  xv = (int)x ; yv = (int)y ; zv = (int)z ;
  mhtAddVoxel(vl, xv, yv, zv) ;
  if ((ceil(x) - x) < 0.5)
    mhtAddVoxel(vl, xv+1, yv, zv) ;
  if ((x - floor(x)) < 0.5)
    mhtAddVoxel(vl, xv-1, yv, zv) ;

  if ((ceil(y) - y) < 0.5)
    mhtAddVoxel(vl, xv, yv+1, zv) ;
  if ((y - floor(y)) < 0.5)
    mhtAddVoxel(vl, xv, yv-1, zv) ;
    
  if ((ceil(z) - z) < 0.5)
    mhtAddVoxel(vl, xv, yv, zv+1) ;
  if ((z - floor(z)) < 0.5)
    mhtAddVoxel(vl, xv, yv, zv-1) ;
    
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mhtAddVoxel(VOXEL_LIST *vl, int xv, int yv, int zv)
{
  int   i ;

  for (i = 0 ; i < vl->nused ; i++)
    if (vl->v[i][0] == xv && vl->v[i][1] == yv && vl->v[i][2] == zv)
      return(NO_ERROR) ;
  if (vl->nused >= MAX_VOXELS)
  {
    fprintf(stderr, "mhtAddVoxels(%d, %d, %d): complete list too big!\n",xv,yv,zv) ;
    ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, 
                                 "mhtAddVoxels(%d, %d, %d): complete list too big!",
                                 xv,yv,zv)) ;
  }
  vl->v[vl->nused][0] = xv ;
  vl->v[vl->nused][1] = yv ;
  vl->v[vl->nused][2] = zv ; vl->nused++ ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MHT_MAX_FACES 10000

static int
mhtDoesFaceVoxelListIntersect(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, 
                              VOXEL_LIST *vl, int fno)
{
  int    xv, yv, zv, fno2, i,j,n1,n2, intersect, nbr, vno, 
         flist[MHT_MAX_FACES], nfaces;
  MHB    *bin ;
  MHBT   *bucket ;
  FACE   *f1, *f2 ;
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3] ;
  

  f1 = &mris->faces[fno] ;

  if (fno == Gdiag_no)
    DiagBreak() ;

  for (nfaces = 0, vno = 0 ; vno < vl->nused ; vno++)
  {
    /*
      go through every entry in this bucket and see if any of them
      occupy the same voxel coordinate as the vertex in question. If so,
      check the face in the bucket to see if the vertex is part of if. If it
      is not, then the voxel is already filled and return 1.
    */
    xv = vl->v[vno][0] ; yv = vl->v[vno][1] ; zv = vl->v[vno][2] ;
#if 0
    bucket = &mht->buckets[xv][yv][zv] ;
#else
    if (!mht->buckets[xv][yv])
      return(0) ;
    bucket = mht->buckets[xv][yv][zv] ;
    if (!bucket)
			continue ;
#endif
    bin = bucket->bins ;
    for (i = 0 ; i < bucket->nused ; i++, bin++)
    {
      fno2 = bin->fno ; f2 = &mris->faces[fno2] ;
      nbr = 0 ;
      for (n1 = 0 ; !nbr && n1 < VERTICES_PER_FACE ; n1++)
        for (n2 = 0 ; !nbr && n2 < VERTICES_PER_FACE ; n2++)
        {
          if (f1->v[n1] == f2->v[n2])
            nbr = 1 ;  /* they share a vertex - don't count it as filled */
        }
      if (!nbr)   /* add it to list, if not already there */
      {
        for (j = 0 ; j < nfaces ; j++)
          if (flist[j] == fno2)
          {
            nbr = 1 ;
            break ;
          }
        if (!nbr)
        {
          if (nfaces >= MHT_MAX_FACES)
            ErrorPrintf(ERROR_NO_MEMORY, "mhtDoesFaceVoxelListIntersect: MHT_MAX_FACES exceeded!") ;
          else
            flist[nfaces++] = fno2 ;
        }
      }
    }
  }

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
    
  for (i = 0 ; i < nfaces ; i++)
  {
    /* fill vertices of 2nd triangle */
    f2 = &mris->faces[flist[i]] ;
    u0[0] = (double)mris->vertices[f2->v[0]].x ;
    u0[1] = (double)mris->vertices[f2->v[0]].y ;
    u0[2] = (double)mris->vertices[f2->v[0]].z ;
    u1[0] = (double)mris->vertices[f2->v[1]].x ;
    u1[1] = (double)mris->vertices[f2->v[1]].y ;
    u1[2] = (double)mris->vertices[f2->v[1]].z ;
    u2[0] = (double)mris->vertices[f2->v[2]].x ;
    u2[1] = (double)mris->vertices[f2->v[2]].y ;
    u2[2] = (double)mris->vertices[f2->v[2]].z ;
    intersect = tri_tri_intersect(v0,v1,v2,u0,u1,u2) ;
    if (intersect)
      return(1) ;
  }
  return(0) ;
}
/*-----------------------------------------------------------------*/
MHBT *
MHTgetBucket(MRIS_HASH_TABLE *mht, float x, float y, float z)
{
  int     xv, yv, zv ;
  MHBT    *bucket ;

  xv = WORLD_TO_VOXEL(mht,x) ; 
  yv = WORLD_TO_VOXEL(mht,y) ; 
  zv = WORLD_TO_VOXEL(mht,z) ;
  if (xv >= FIELD_OF_VIEW || yv >= FIELD_OF_VIEW || zv >= FIELD_OF_VIEW ||
      xv < 0 || yv < 0 || zv < 0)
    return(NULL) ;
  else if (!mht->buckets[xv][yv])
    return(NULL) ;

  bucket = mht->buckets[xv][yv][zv] ;

  return(bucket) ;
}
/*--------------------------------------------------------------------*/
VERTEX *
MHTfindClosestVertex(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v)
{
  VERTEX    *vmin, *vdst ;
  int       i, xk, yk, zk ;
  double    dist, min_dist ;
  MHB       *bin ;
  MHBT      *bucket ;
  float     x, y, z ;

  min_dist = 10000000 ; vmin = NULL ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
      {      
        x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->x)+xk) ;
        y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->y)+yk) ;
        z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->z)+zk) ;
        bucket = MHTgetBucket(mht, x, y, z) ;
        if (!bucket)
          continue ;
        bin = bucket->bins ; 
        for (i = 0 ; i < bucket->nused ; i++, bin++)
        {
          vdst = &mris->vertices[bin->fno] ;

          if (bin->fno == Gdiag_no)
            DiagBreak() ;
          dist = sqrt(SQR(vdst->x-v->x)+SQR(vdst->y-v->y)+SQR(vdst->z-v->z)) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            vmin = vdst ;
          }
        }
      }
    }
  }

  return(vmin) ;
}
VERTEX *
MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, VERTEX *v, int which)
{
  VERTEX    *vmin, *vdst ;
  int       i, xk, yk, zk ;
  double    dist, min_dist ;
  MHB       *bin ;
  MHBT      *bucket ;
  float     x, y, z, tx, ty, tz ;
  tx=0;
  ty=0;
  tz=0;

  min_dist = 10000000 ; vmin = NULL ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
      {      
        x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->x)+xk) ;
        y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->y)+yk) ;
        z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->z)+zk) ;
        bucket = MHTgetBucket(mht, x, y, z) ;
        if (!bucket)
          continue ;
        bin = bucket->bins ; 
        for (i = 0 ; i < bucket->nused ; i++, bin++)
        {
          vdst = &mris->vertices[bin->fno] ;

          if (bin->fno == Gdiag_no)
            DiagBreak() ;
					switch (which)
					{
					case CURRENT_VERTICES:
						tx = vdst->x ; ty = vdst->y  ; tz = vdst->z ; 
						break ;
					default:
					case ORIGINAL_VERTICES:
						tx = vdst->origx ; ty = vdst->origy  ; tz = vdst->origz ; 
						break ;
					case WHITE_VERTICES:
						tx = vdst->whitex ; ty = vdst->whitey  ; tz = vdst->whitez ; 
						break ;
					case PIAL_VERTICES:
						tx = vdst->pialx ; ty = vdst->pialy  ; tz = vdst->pialz ; 
						break ;
					}
          dist = sqrt(SQR(tx-v->x)+SQR(ty-v->y)+SQR(tz-v->z)) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            vmin = vdst ;
          }
        }
      }
    }
  }

  return(vmin) ;
}
/*--------------------------------------------------------------------
  MHTfindClosestVertexNo() - basically the same as findClosestVertex
  except it returns the vertex number instead of a pointer to the
  vertex. Also, passes the minimum distance found
  --------------------------------------------------------------------*/
int MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, 
         VERTEX *v, float *min_dist)
{
  VERTEX    *vmin, *vdst ;
  int       i, vtxno, vtxno_min ;
  double    dist;
  MHB       *bin ;
  MHBT      *bucket ;

  /* get bucket closest to the vertex */
  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  if (!bucket) 
    return(-1) ;

  *min_dist = 10000000 ; 
  vmin = NULL ;
  vtxno_min = 0;

  /* go through each bin in the bucket */
  for (i = 0 ; i < bucket->nused ; i++){
    bin = &(bucket->bins[i]);
    vtxno = bin->fno;
    vdst = &mris->vertices[vtxno] ;
    if(vdst->ripflag) continue;
    dist = sqrt(SQR(vdst->x-v->x)+SQR(vdst->y-v->y)+SQR(vdst->z-v->z)) ;
    if (dist < *min_dist) {
      vtxno_min = vtxno;
      *min_dist = dist ;
      vmin = vdst ;
    }
  }

  return(vtxno_min) ;
}
#define MAX_VERTICES 50000
int *
MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, 
                                int vno, float max_dist, int *pvnum)
{
  int       vertices[MAX_VERTICES], *returned_vertices ;
  VERTEX    *vdst, *v ;
  int       i, xk, yk, zk, vnum ;
  double    dist ;
  MHB       *bin ;
  MHBT      *bucket ;
  float     x, y, z ;

  v = &mris->vertices[vno] ;
  for (vnum = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
      {      
        switch (mht->which_vertices)
        {
        default:
        case CURRENT_VERTICES:
          x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->x)+xk) ;
          y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->y)+yk) ;
          z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->z)+zk) ;
          break ;
        case CANONICAL_VERTICES:
          x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->cx)+xk) ;
          y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->cy)+yk) ;
          z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->cz)+zk) ;
          break ;
        case WHITE_VERTICES:
          x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->whitex)+xk) ;
          y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->whitey)+yk) ;
          z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->whitez)+zk) ;
          break ;
        case PIAL_VERTICES:
          x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->pialx)+xk) ;
          y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->pialy)+yk) ;
          z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->pialz)+zk) ;
          break ;
        case ORIGINAL_VERTICES:
          x = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->origx)+xk) ;
          y = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->origy)+yk) ;
          z = VOXEL_TO_WORLD(mht, WORLD_TO_VOLUME(mht, v->origz)+zk) ;
          break ;
        }
        bucket = MHTgetBucket(mht, x, y, z) ;
        if (!bucket)
          continue ;
        bin = bucket->bins ; 
        for (i = 0 ; i < bucket->nused ; i++, bin++)
        {
          vdst = &mris->vertices[bin->fno] ;
          if (vdst->ripflag)
            continue ;  /* already processed */

          if (bin->fno == Gdiag_no)
            DiagBreak() ;
          switch (mht->which_vertices)
          {
          default:
          case CURRENT_VERTICES:
            dist = sqrt(SQR(vdst->x-v->x)+SQR(vdst->y-v->y)+SQR(vdst->z-v->z)) ;
            break ;
          case WHITE_VERTICES:
            dist = sqrt(SQR(vdst->whitex-v->whitex)+SQR(vdst->whitey-v->whitey)+SQR(vdst->whitez-v->whitez)) ;
            break ;
          case PIAL_VERTICES:
            dist = sqrt(SQR(vdst->pialx-v->pialx)+SQR(vdst->pialy-v->pialy)+SQR(vdst->pialz-v->pialz)) ;
            break ;
          case CANONICAL_VERTICES:
            dist = 
              sqrt(SQR(vdst->cx-v->cx)+SQR(vdst->cy-v->cy)+SQR(vdst->cz-v->cz));
            break ;
          case ORIGINAL_VERTICES:
            dist = 
              sqrt(SQR(vdst->origx-v->origx)+
                   SQR(vdst->origy-v->origy)+SQR(vdst->origz-v->origz));
            break ;
          }
          if (dist < max_dist)
          {
            if (vnum >= MAX_VERTICES)
              ErrorExit(ERROR_NOMEMORY, 
                        "MHTgetAllVerticesWithinDistance: couldn't"
                        "fit vertices!") ;
            vdst->ripflag = 1 ;
            vertices[vnum++] = bin->fno ;
          }
        }
      }
    }
  }

  returned_vertices = (int *)calloc(vnum, sizeof(int)) ;
  if (!returned_vertices)
    ErrorExit(ERROR_NOMEMORY, 
              "MHTgetAllVerticesWithinDistance: couldn't allocate %d int array",
              vnum) ;

  memmove(returned_vertices, vertices, vnum*sizeof(int)) ;
  *pvnum = vnum;
  for (i = 0 ; i < vnum ; i++)
    mris->vertices[vertices[i]].ripflag = 0 ;
  return(returned_vertices) ;
}


/*
	find vertex in hash table that is closest to input coordinates */
VERTEX *
MHTfindClosestVertexInTable(MRIS_HASH_TABLE *mht, MRI_SURFACE *mris, float x, float y, float z)
{
  VERTEX    *vmin, *vdst ;
  int       i, xk, yk, zk ;
  double    dist, min_dist ;
  MHB       *bin ;
  MHBT      *bucket ;

  min_dist = 10000000 ; vmin = NULL ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
      {      
        bucket = MHTgetBucket(mht, x+xk, y+yk, z+zk) ;
        if (!bucket)
          continue ;
        bin = bucket->bins ; 
        for (i = 0 ; i < bucket->nused ; i++, bin++)
        {
          vdst = &mris->vertices[bin->fno] ;

          if (bin->fno == Gdiag_no)
            DiagBreak() ;
          dist = sqrt(SQR(vdst->x-x)+SQR(vdst->y-y)+SQR(vdst->z-z)) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            vmin = vdst ;
          }
        }
      }
    }
  }

  return(vmin) ;
}
