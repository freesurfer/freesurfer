/*
 *       FILE NAME:   mripolv.c
 *
 *       DESCRIPTION: utilities for calculating and filtering
 *                    MRI data based on planes of least variance
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        1/8/97
 *
*/
/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>

#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "filter.h"
#include "box.h"
#include "region.h"
#include "nr.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x) == 7)&&((y)==15)&&((z)==15))
#define MAXLEN 256


/*-----------------------------------------------------
                    GLOBAL DATA
-------------------------------------------------------*/

/* vertices of an icosohedron (sp?) */
float ic_x_vertices[NVERTICES] =
{
   0,       0.8944, 0.2764, -0.7236, -0.7236, 0.2764, -0.4253,
  -0.8507, -0.4253, 0.1625, -0.2629,  0.5257, 0.6882,  0.1625,
   0.6882, -0.2629,-0.5878,  0,      -0.9511,-0.9511, -0.5878, 1.00
} ;

float ic_y_vertices[NVERTICES] =
{
  0,     0,      0.8507, 0.5257, -0.5257, -0.8507, -0.3090,   
  0,     0.3090,-0.5000,-0.8090,  0,      -0.5000,  0.5000,
  0.5000,0.8090, 0.8090, 1.0000,  0.3090, -0.3090, -0.8090, 0.00
} ;

float ic_z_vertices[NVERTICES] =
{
  1.0000, 0.4472, 0.4472, 0.4472, 0.4472, 0.4472, 0.8507,
  0.5257, 0.8507, 0.8507, 0.5257, 0.8507, 0.5257, 0.8507,
  0.5257, 0.5257, 0,      0,      0,      0,      0, 0.00
} ;

static int vertices_initialized = 0 ;
static float e1_x_v[NVERTICES] ;
static float e1_y_v[NVERTICES] ;
static float e1_z_v[NVERTICES] ;
static float e2_x_v[NVERTICES] ;
static float e2_y_v[NVERTICES] ;
static float e2_z_v[NVERTICES] ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int compare_sort_array(const void *pc1, const void *pc2) ;
static void init_basis_vectors(void) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIpolvMean(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi  ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase, total ;
  BUFTYPE  *pdst, *pptr ;

  init_basis_vectors() ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
  n = wsize*wsize ;
  for (z = whalf ; z < depth ; z++)
  {
    DiagHeartbeat((float)(z-whalf) / (float)(depth-whalf-1)) ;
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, whalf, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, whalf, y, z) ; /* ptr to normal vectors */
      for (x = whalf ; x < width ; x++)
      {
        vertex = *pptr++ ;
        e1_x = e1_x_v[vertex] ;  /* basis vectors for plane */
        e1_y = e1_y_v[vertex] ;
        e1_z = e1_z_v[vertex] ;
        e2_x = e2_x_v[vertex] ;
        e2_y = e2_y_v[vertex] ;
        e2_z = e2_z_v[vertex] ;

        /* 
           calculate the mean in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        /* now find the values in this plane */
        total = 0 ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          xbase = (float)x + (float)yk * e2_x ;
          ybase = (float)y + (float)yk * e2_y ;
          zbase = (float)z + (float)yk * e2_z ;
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk*e1_x) ;
            xi = pxi[xi] ;
            yi = nint(ybase + xk*e1_y) ;
            yi = pyi[yi] ;
            zi = nint(zbase + xk*e1_z) ;
            zi = pzi[zi] ;
            total += (float)MRIvox(mri_src, xi, yi, zi) ;
          }
        }
        *pdst++ = total / n ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIpolvMedian(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  BUFTYPE  *pdst, *pptr, plane_vals[MAXLEN], *pvals ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  init_basis_vectors() ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
  n = wsize*wsize ;
  for (z = whalf ; z < depth ; z++)
  {
    DiagHeartbeat((float)(z-whalf) / (float)(depth-whalf-1)) ;
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, whalf, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, whalf, y, z) ; /* ptr to normal vectors */
      for (x = whalf ; x < width ; x++)
      {
        vertex = *pptr++ ;
        e1_x = e1_x_v[vertex] ;  /* basis vectors for plane */
        e1_y = e1_y_v[vertex] ;
        e1_z = e1_z_v[vertex] ;
        e2_x = e2_x_v[vertex] ;
        e2_y = e2_y_v[vertex] ;
        e2_z = e2_z_v[vertex] ;

        /* 
           calculate the median in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        pvals = plane_vals ;
        /* now find the values in this plane */
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          xbase = (float)x + (float)yk * e2_x ;
          ybase = (float)y + (float)yk * e2_y ;
          zbase = (float)z + (float)yk * e2_z ;
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk*e1_x) ;
            xi = pxi[xi] ;
            yi = nint(ybase + xk*e1_y) ;
            yi = pyi[yi] ;
            zi = nint(zbase + xk*e1_z) ;
            zi = pzi[zi] ;
            *pvals++ = (float)MRIvox(mri_src, xi, yi, zi) ;
          }
        }
        qsort(plane_vals, n, sizeof(BUFTYPE), compare_sort_array) ;
        *pdst++ = plane_vals[n/2] ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIpolvMeanRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize,
                  MRI_REGION *region)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n,
           vertex, x0, y0, z0, xi, yi, zi, *pxi, *pyi, *pzi ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase, total ;
  BUFTYPE  *pdst, *pptr ;

  init_basis_vectors() ;
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, 
                 "MRIpolvMeanRegion: unsupported src type")) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;
  whalf = (wsize-1)/2 ;

  if (mri_dst && mri_dst->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, 
                 "MRIpolvMeanRegion: unsupported dst type")) ;

  if (!mri_dst)
  {
    int  w, h, d ;
    
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_dst = MRIalloc(w, h, d, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  n = wsize * wsize ;
  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y-y0, z-z0) ;    /* ptr to destination */
      pptr = &MRIvox(mri_polv, x0, y, z) ;   /* ptr to normal vectors */
      for (x = x0 ; x < width ; x++)
      {
        vertex = *pptr++ ;
        e1_x = e1_x_v[vertex] ;  /* basis vectors for plane */
        e1_y = e1_y_v[vertex] ;
        e1_z = e1_z_v[vertex] ;
        e2_x = e2_x_v[vertex] ;
        e2_y = e2_y_v[vertex] ;
        e2_z = e2_z_v[vertex] ;

        /* 
           calculate the mean in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        /* now find the values in this plane */
        total = 0 ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          xbase = (float)x + (float)yk * e2_x ;
          ybase = (float)y + (float)yk * e2_y ;
          zbase = (float)z + (float)yk * e2_z ;
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk*e1_x) ;
            xi = pxi[xi] ;
            yi = nint(ybase + xk*e1_y) ;
            yi = pyi[yi] ;
            zi = nint(zbase + xk*e1_z) ;
            zi = pzi[zi] ;
            total += (float)MRIvox(mri_src, xi, yi, zi) ;
          }
        }
        *pdst++ = (BUFTYPE)((float)total / (float)n) ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIpolvMedianRegion(MRI *mri_src, MRI *mri_dst,MRI *mri_polv,int wsize, 
                    MRI_REGION *region)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n, vertex,
           x0, y0, z0, xi, yi, zi, *pxi, *pyi, *pzi, median_index ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  BUFTYPE  *pdst, *pptr, plane_vals[MAXLEN], *pvals ;

  init_basis_vectors() ;
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, 
                 "MRIpolvMedianRegion: unsupported src type")) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;
  whalf = (wsize-1)/2 ;

  if (mri_dst && (mri_dst->type != MRI_UCHAR))
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, 
                 "MRIpolvMedianRegion: unsupported dst type")) ;
  if (!mri_dst)
  {
    int  w, h, d ;
    
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_dst = MRIalloc(w, h, d, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  n = wsize*wsize ;
  median_index = n/2 ;
  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y-y0, z-z0) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, x0, y, z) ;      /* ptr to normal vectors */
      for (x = x0 ; x < width ; x++)
      {
if (DEBUG_POINT(x,y,z))
  DiagBreak() ;

        vertex = *pptr++ ;
        e1_x = e1_x_v[vertex] ;  /* basis vectors for plane */
        e1_y = e1_y_v[vertex] ;
        e1_z = e1_z_v[vertex] ;
        e2_x = e2_x_v[vertex] ;
        e2_y = e2_y_v[vertex] ;
        e2_z = e2_z_v[vertex] ;

        /* 
           calculate the median in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        pvals = plane_vals ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          xbase = (float)x + (float)yk * e2_x ;
          ybase = (float)y + (float)yk * e2_y ;
          zbase = (float)z + (float)yk * e2_z ;
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk*e1_x) ;
            xi = pxi[xi] ;
            yi = nint(ybase + xk*e1_y) ;
            yi = pyi[yi] ;
            zi = nint(zbase + xk*e1_z) ;
            zi = pzi[zi] ;
            *pvals++ = (int)MRIvox(mri_src, xi, yi, zi) ;
          }
        }
        qsort(plane_vals, n, sizeof(BUFTYPE), compare_sort_array) ;
        *pdst++ = plane_vals[median_index] ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIextractPlanes(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int x, int y, 
                 int z, int wsize)
{
  float    a, b, c ;
  int      vertex, whalf, xk, yk, zk, pno, xp, yp, xoff ;
  BUFTYPE  *psrc, val ;

  init_basis_vectors() ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
    mri_dst = MRIalloc(wsize*wsize, wsize, 1, MRI_UCHAR) ;

  vertex = (int)MRIvox(mri_polv, x, y, z) ;
  a = ic_x_vertices[vertex] ;   /* vector on unit sphere */
  b = ic_y_vertices[vertex] ;
  c = ic_z_vertices[vertex] ;
  
  /* calculate the variance of all planes orthogonal to (a,b,c) */
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      psrc = &MRIvox(mri_src, x-whalf, y+yk, z+zk) ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        /* pno is the index (#) of the plane we are in */
        
        /* find the index of this plane (projection on normal) */
        pno = whalf + (int)(zk*c + yk*b + xk*a + 0.5) ;
        pno = (pno < 0) ? 0 : (pno > wsize-1) ? wsize-1 : pno ;
        val = (float)*psrc++ ;
        xp = whalf + (int)(xk*a + 0.5) ;
        xp = (xp < 0) ? 0 : (xp > wsize-1) ? wsize-1 : xp ;
        yp = whalf + (int)(yk*b + 0.5) ;
        yp = (yp < 0) ? 0 : (yp > wsize-1) ? wsize-1 : yp ;

        /* lay out the planes one after the other horizontally */
        xoff = pno + wsize ;
        MRIvox(mri_dst, xoff+xp, yp, 0) = val ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/


MRI *
MRIplaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int      width, height, depth, x, y, z, whalf, vertex, xk, yk, zk, pno,
           mini, maxi ;
  float    min_var, max_var, a, b, c, total[MAXLEN], total_sq[MAXLEN],
           nv[MAXLEN], varv[MAXLEN], avgv[MAXLEN], val, total_var ;
  BUFTYPE  *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
  for (z = whalf ; z < depth ; z++)
  {
    DiagHeartbeat((float)(z-whalf) / (float)(depth-whalf-1)) ;
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, whalf, y, z) ;
      for (x = whalf ; x < width ; x++)
      {
        /*
          for this point (x,y,z), go through a set of directions on the unit
          sphere. For each direction, find all the planes orthogonal to that
          dir. within our window, and pick the direction in which the variance
          of all the planes is smallest. This will hopefully be the normal to 
          the cortical surface.
          */
        maxi = mini = -1 ;
        min_var = 100000.0f ;    /* minimum variance of set of planes */
        max_var = -100000.0f ;   /* maximum variance of set of planes */
        for (vertex = 0 ; vertex < NVERTICES ; vertex++)
        {
          a = ic_x_vertices[vertex] ;   /* vector on unit sphere */
          b = ic_y_vertices[vertex] ;
          c = ic_z_vertices[vertex] ;

          memset(total, 0, wsize*sizeof(float)) ;
          memset(total_sq, 0, wsize*sizeof(float)) ;
          memset(nv, 0, wsize*sizeof(float)) ;

          /* calculate the variance of all planes orthogonal to (a,b,c) */
          for (zk = -whalf ; zk <= whalf ; zk++)
          {
            for (yk = -whalf ; yk <= whalf ; yk++)
            {
              psrc = &MRIvox(mri_src, x-whalf, y+yk, z+zk) ;
              for (xk = -whalf ; xk <= whalf ; xk++)
              {
                /* pno is the index (#) of the plane we are in */

                /* find the index of this plane (projection on normal) */
                pno = whalf + nint(zk*c + yk*b + xk*a) ;
                pno = (pno < 0) ? 0 : (pno > wsize-1) ? wsize-1 : pno ;
                val = (float)*psrc++ ;
                total[pno] += val ;       /* sum of all values in this plane */
                total_sq[pno] += val*val ; 
                nv[pno]++ ;                /* # of values in this plane */
              }
            }
          }
          total_var = 0.0f ;
          for (pno = 0 ; pno < wsize ; pno++)
          {
            if (!nv[pno])
              varv[pno] = avgv[pno] = 0.0f ;
            else
            {
              avgv[pno] = total[pno]/nv[pno];  /* mean of this plane */
              varv[pno] = total_sq[pno]/nv[pno]-avgv[pno]*avgv[pno];
            }
            total_var += varv[pno];   /* total variance of set of planes */
          }
          total_var /= (float)wsize ;
          if (total_var>max_var) {max_var=total_var;maxi=vertex;}
          if (total_var<min_var) {min_var=total_var;mini=vertex;}
        }
        /* done - put vector components into output */
        *pdst++ = (BUFTYPE)mini;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIcentralPlaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int      width, height, depth, x, y, z, whalf, vertex, xk, yk,
           mini, maxi, xi, yi, zi, *pxi, *pyi, *pzi, x1, y1, z1, x0, y0,z0;
  float    min_var, max_var, total, total_sq, nv, varv, avgv, val,
           background_val, fmax ;
  BUFTYPE  *pdst, max_val ;
  float    xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z,
           *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z ;

  if (getenv("USE_CACHED_CPOLV"))
  {
    char fname[100] ;
    MRI  *mri_tmp ;

    /* try and read previously computed CPOLV file from disk */
    sprintf(fname, "%s/cpolv.mnc", mri_src->fname) ;
    mri_tmp = MRIread(fname) ;
    if (mri_tmp)
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stderr, "reading previously calculated cpolv %s\n", fname) ;
      mri_dst = MRIcopy(mri_tmp, NULL) ;
      MRIfree(&mri_tmp) ;
      return(mri_dst) ;
    }
  }

  init_basis_vectors() ;
    
  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  MRIvalRange(mri_src, &background_val, &fmax) ;
  background_val *= 0.2f ;  /* anything smaller than 20% of peak is bg */

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
  if (mri_src->roi.dx > 0)
  {
    x0 = MAX(whalf, mri_src->roi.x) ;
    y0 = MAX(whalf, mri_src->roi.y) ;
    z0 = MAX(whalf, mri_src->roi.z) ;
    x1 = MIN(x0 + mri_src->roi.dx - 1, width) ;
    y1 = MIN(y0 + mri_src->roi.dy - 1, height) ;
    z1 = MIN(z0 + mri_src->roi.dz - 1, depth) ;
  }
  else
  {
    x0 = y0 = z0 = whalf ;
    x1 = width-1 ;
    y1 = height-1 ;
    z1 = depth-1 ;
  }

  for (z = z0 ; z <= z1 ; z++)
  {
    DiagHeartbeat((float)(z-z0) / (float)z1) ;
    for (y = y0 ; y <= y1 ; y++)
    {
      pdst = &MRIvox(mri_dst, x0, y, z) ;
      for (x = x0 ; x <= x1 ; x++)
      {
        /*
          for this point (x,y,z), go through a set of directions on the unit
          sphere. For each direction, find all the planes orthogonal to that
          dir. within our window, and pick the direction in which the variance
          of all the planes is smallest. This will hopefully be the normal to 
          the cortical surface.
          */
        maxi = mini = -1 ;
        min_var = 100000.0f ;    /* minimum variance of central planes */
        max_var = -100000.0f ;   /* maximum variance of central planes */
        pe1_x = e1_x_v ; pe1_y = e1_y_v ; pe1_z = e1_z_v ;
        pe2_x = e2_x_v ; pe2_y = e2_y_v ; pe2_z = e2_z_v ;
        for (vertex = 0 ; vertex < NVERTICES ; vertex++)
        {
          e1_x = *pe1_x++ ;   /* first in-plane basis vector */
          e1_y = *pe1_y++ ;
          e1_z = *pe1_z++ ;
          e2_x = *pe2_x++ ;   /* second in-plane basis vector */
          e2_y = *pe2_y++ ;
          e2_z = *pe2_z++ ;

          total = total_sq = nv = 0.0f ;
          max_val = 0 ;
          /* now find the values in this plane */
          for (yk = -whalf ; yk <= whalf ; yk++)
          {
            xbase = (float)x + (float)yk * e2_x ;
            ybase = (float)y + (float)yk * e2_y ;
            zbase = (float)z + (float)yk * e2_z ;
            for (xk = -whalf ; xk <= whalf ; xk++)
            {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk*e1_x) ;
              xi = pxi[xi] ;
              yi = nint(ybase + xk*e1_y) ;
              yi = pyi[yi] ;
              zi = nint(zbase + xk*e1_z) ;
              zi = pzi[zi] ;
              val = (float)MRIvox(mri_src, xi, yi, zi) ;
              total += val ;       /* sum of all values in this plane */
              total_sq += val*val ;/* sum of squared values in this plane */
              nv++ ;               /* # of values in this plane */
            }
          }

          if (!nv)
              varv = avgv = 0.0f ;
          else
          {
            avgv = total/nv;  /* mean of this plane */
            varv = total_sq / nv - avgv * avgv ;
          }

          if (varv>max_var) {max_var=varv;maxi=vertex;}
          if (varv<min_var) {min_var=varv;mini=vertex;}
          if (FZERO(varv))  /* zero variance - won't find anything less */
            break ;
        }
        /* done - put vector components into output */
        *pdst++ = (BUFTYPE)mini;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
compare_sort_array(const void *pc1, const void *pc2)
{
  register BUFTYPE c1, c2 ;

  c1 = *(BUFTYPE *)pc1 ;
  c2 = *(BUFTYPE *)pc2 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return(1) ;
  else if (c1 < c2)
    return(-1) ;

  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
             1st take the cross-product of ez with a random vector
             to obtain one vector in the plane, then take the cross
             product of that vector with the normal (ez) to obtain
             the 2nd in-plane basis vector.
------------------------------------------------------*/
static void
init_basis_vectors(void)
{
  float vx, vy, vz, *px, *py, *pz, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, 
        *pe2_z, e3_x, e3_y, e3_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, len ;
  int   vertex ;

  if (vertices_initialized)
    return ;

  px = ic_x_vertices ; py = ic_y_vertices ; pz = ic_z_vertices ;
  pe1_x = e1_x_v ; pe1_y = e1_y_v ; pe1_z = e1_z_v ;
  pe2_x = e2_x_v ; pe2_y = e2_y_v ; pe2_z = e2_z_v ;
  for (vertex = 0 ; vertex < NVERTICES ; vertex++)
  {
    e3_x = *px++ ;   /* vector on unit sphere */
    e3_y = *py++ ;
    e3_z = *pz++ ;

    /* pick some other unit non-linearly dependent vector */
    vx = e3_y ;
    vy = e3_z ;
    vz = e3_x ;

    /* don't care about sign (right-hand rule) */
    e1_x = vy*e3_z - vz*e3_y ;
    e1_y = vz*e3_x - vx*e3_z ;
    e1_z = vx*e3_y - vy*e3_x ;
    len = sqrt(e1_x*e1_x+e1_y*e1_y+e1_z*e1_z) ;
    *pe1_x++ = e1_x / len ;   
    *pe1_y++ = e1_y / len ;
    *pe1_z++ = e1_z / len ;
    e1_x /= len ; e1_y /= len ; e1_z /= len ;

    e2_x = e1_y*e3_z - e1_z*e3_y ;  
    e2_y = e1_x*e3_z - e1_z*e3_x ;
    e2_z = e1_y*e3_x - e1_x*e3_y ;
    len = sqrt(e2_x*e2_x+e2_y*e2_y+e2_z*e2_z) ;
    *pe2_x++ = e2_x / len ;   
    *pe2_y++ = e2_y / len ;
    *pe2_z++ = e2_z / len ;
    e2_x /= len ; e2_y /= len ; e2_z /= len ;
#if 0
    DiagFprintf(0L, 
              "vertex %d: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f) "
              "x (%2.2f, %2.2f, %2.2f)\n",
              vertex, ic_x_vertices[vertex], ic_y_vertices[vertex], 
              ic_z_vertices[vertex], e1_x_v[vertex], e1_y_v[vertex], 
              e1_z_v[vertex],e2_x_v[vertex], e2_y_v[vertex], e2_z_v[vertex]) ;
    DiagFprintf(0L, "lengths: %2.3f, %2.3f, %2.3f\n",
                sqrt(e3_x*e3_x+e3_y*e3_y+e3_z*e3_z),
                sqrt(e1_x*e1_x+e1_y*e1_y+e1_z*e1_z),
                sqrt(e2_x*e2_x+e2_y*e2_y+e2_z*e2_z)) ;
#endif
  }
  vertices_initialized = 1 ;
}

