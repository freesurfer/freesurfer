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

#define DEBUG_POINT(x,y,z)  (((x) == 32)&&((y)==25)&&((z)==32))
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
  n = wsize*wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
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
MRIpolvNormalCurvature(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int      width, height, depth, x, y, z, whalf, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi  ;
  float    nx, ny, nz, mean, var, val, std, *pdst ;
  BUFTYPE  *pptr ;

  init_basis_vectors() ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width,height,depth,MRI_FLOAT,mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  n = wsize-1 ;  /* excludes central point */
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIFvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
      {
        vertex = *pptr++ ;
        nx = ic_x_vertices[vertex] ;  /* normal vector */
        ny = ic_y_vertices[vertex] ;
        nz = ic_z_vertices[vertex] ;

        /* 
           calculate the mean in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        /* now find the values in the normal direction */
        mean = var = 0.0f ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          if (!yk)
            continue ;  /* skip central point */
          xi = pxi[nint((float)x + (float)yk * nx)] ;
          yi = pyi[nint((float)y + (float)yk * ny)] ;
          zi = pzi[nint((float)z + (float)yk * nz)] ;
          val = (float)MRIvox(mri_src, xi, yi, zi) ;
          var += val * val ;
          mean += val ;
        }
        mean /= (float)n ;
        val = (float)MRIvox(mri_src, x, y, z) ;
        std = sqrt(var / (float)n - mean*mean) ;
        if (!FZERO(std))
          *pdst++ = (val - mean) / std ;
        else
          *pdst++ = 0.0f ;
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
MRIpolvZscore(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int      width, height, depth, x, y, z, whalf, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi  ;
  float    nx, ny, nz, mean, var, val, std, *pdst, xf, yf, zf ;
  BUFTYPE  *pptr ;

  init_basis_vectors() ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width,height,depth,MRI_FLOAT,mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  n = wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIFvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
      {
        vertex = *pptr++ ;
        nx = ic_x_vertices[vertex] ;  /* normal vector */
        ny = ic_y_vertices[vertex] ;
        nz = ic_z_vertices[vertex] ;

        /* 
           calculate the mean in the plane orthogonal to (a,b,c), 
           through the current origin (x,y,z).
           */
        /* now find the values in the normal direction */
        mean = var = 0.0f ;
        xf = (float)x - whalf*nx ;
        yf = (float)y - whalf*ny ;
        zf = (float)z - whalf*nz ;
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          xi = pxi[nint(xf)] ;
          yi = pyi[nint(yf)] ;
          zi = pzi[nint(zf)] ;
          val = (float)MRIvox(mri_src, xi, yi, zi) ;
          var += val * val ;
          mean += val ;
          xf += nx ; yf += ny ; zf += nz ;
        }
        mean /= (float)n ;
        val = (float)MRIvox(mri_src, x, y, z) ;
        std = sqrt(var / (float)n - mean*mean) ;
        if (!FZERO(std))
          *pdst++ = (val - mean) / std ;
        else
          *pdst++ = 0.0f ;
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
  n = wsize*wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
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
MRIpolvOrder(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, int thresh)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi, order ;
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
  n = wsize*wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
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

        /* find the 1st supra-threshold value in the array */
        pvals = plane_vals ;
        for (order = 0 ; order < n ; order++)
          if (*pvals++ > thresh)
            break ;
        *pdst++ = (BUFTYPE)order ;
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
MRIextractCpolvCoords(MRI *mri_src, int *px, int *py, int *pz, 
                          MRI *mri_polv, int x, int y,int z, int wsize)
{
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      vertex, whalf, xk, yk, xi, yi, zi ;

  init_basis_vectors() ;
  whalf = (wsize-1)/2 ;

  vertex = (int)MRIvox(mri_polv, x, y, z) ;
  e1_x = e1_x_v[vertex] ;   /* get basis vectors for plane */
  e1_y = e1_y_v[vertex] ;
  e1_z = e1_z_v[vertex] ;
  e2_x = e2_x_v[vertex] ;
  e2_y = e2_y_v[vertex] ;
  e2_z = e2_z_v[vertex] ;
  
  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      *px++ = xi ; *py++ = yi ; *pz++ = zi ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIextractCpolv(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int x, int y, 
                 int z, int wsize)
{
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      vertex, whalf, xk, yk, xi, yi, zi ;

  init_basis_vectors() ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }

  vertex = (int)MRIvox(mri_polv, x, y, z) ;
  e1_x = e1_x_v[vertex] ;   /* get basis vectors for plane */
  e1_y = e1_y_v[vertex] ;
  e1_z = e1_z_v[vertex] ;
  e2_x = e2_x_v[vertex] ;
  e2_y = e2_y_v[vertex] ;
  e2_z = e2_z_v[vertex] ;
  
  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
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
MRIextractPlane(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int x, int y, 
                 int z, int wsize)
{
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      vertex, whalf, xk, yk, xi, yi, zi ;

  init_basis_vectors() ;
  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }

  vertex = (int)MRIvox(mri_polv, x, y, z) ;
  e1_x = e1_x_v[vertex] ;   /* get basis vectors for plane */
  e1_y = e1_y_v[vertex] ;
  e1_z = e1_z_v[vertex] ;
  e2_x = e2_x_v[vertex] ;
  e2_y = e2_y_v[vertex] ;
  e2_z = e2_z_v[vertex] ;
  
  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
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
        if (MRIvox(mri_src, x, y, z) < 50)
          continue ;

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
  float    min_mean, min_var, max_var, total, total_sq, nv, varv, avgv, val,
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
      if (mri_tmp->width == mri_src->width && 
          mri_tmp->height == mri_src->height &&
          mri_tmp->depth == mri_src->depth)
      {
        mri_dst = MRIcopy(mri_tmp, NULL) ;
        MRIfree(&mri_tmp) ;
        return(mri_dst) ;
      }
      MRIfree(&mri_tmp) ;
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

#if 0
  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
#endif

  if (mri_src->roi.dx > 0)
  {
    x0 = MAX(0, mri_src->roi.x) ;
    y0 = MAX(0, mri_src->roi.y) ;
    z0 = MAX(0, mri_src->roi.z) ;
    x1 = MIN(x0 + mri_src->roi.dx - 1, width-1) ;
    y1 = MIN(y0 + mri_src->roi.dy - 1, height-1) ;
    z1 = MIN(z0 + mri_src->roi.dz - 1, depth-1) ;
  }
  else
  {
    x0 = y0 = z0 = 0 ;
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
        min_mean = 1000.0f ;     /* mean of minimum variance plane */
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

          if (varv>max_var) { max_var=varv ; maxi=vertex ;}
          if (varv<min_var) { min_var=varv ; mini=vertex; min_mean = avgv ;}
          if (FZERO(varv))  /* zero variance - won't find anything less */
            break ;
        }
#if 0
        /* done - put vector components into output */
        MRIseq_vox(mri_dst, x, y, z, 1) = (BUFTYPE)nint(min_mean) ;
#endif
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
    e3_x = *px ;   /* vector on unit sphere */
    e3_y = *py ;
    e3_z = *pz ;
/* 
   now we must scale the length of the vector so that it reaches the
   border of the next voxel. Thus, 'diagonal' vectors must be extended
   by sqrt(2) relative to those which lie along the cardinal axes.
*/
    vx = fabs(e3_x) ; vy = fabs(e3_y) ; vz = fabs(e3_z) ; /* use symmetry */
    if ((vx > vy) && (vx > vz))  /* scale using x component */
      len = 1.0f / vx ;
    else if (vy > vz)            /* scale using y component */
      len = 1.0f / vy ;
    else                         /* scale using z component */
      len = 1.0f / vz ;
    *px++ = e3_x * len ;
    *py++ = e3_y * len ;
    *pz++ = e3_z * len ;


    /* pick some other unit non-linearly dependent vector */
    vx = e3_y ;
    vy = e3_z ;
    vz = e3_x ;

    /* don't care about sign (right-hand rule) */
    e1_x = vy*e3_z - vz*e3_y ;
    e1_y = vz*e3_x - vx*e3_z ;
    e1_z = vx*e3_y - vy*e3_x ;

    vx = fabs(e1_x) ; vy = fabs(e1_y) ; vz = fabs(e1_z) ; /* use symmetry */
    if ((vx > vy) && (vx > vz))  /* scale using x component */
      len = 1.0f / vx ;
    else if (vy > vz)            /* scale using y component */
      len = 1.0f / vy ;
    else                         /* scale using z component */
      len = 1.0f / vz ;

    e1_x = *pe1_x++ = e1_x * len ;   
    e1_y = *pe1_y++ = e1_y * len ;
    e1_z = *pe1_z++ = e1_z * len ;
    len = sqrt(e1_x*e1_x+e1_y*e1_y+e1_z*e1_z) ;

    e2_x = e1_y*e3_z - e1_z*e3_y ;  
    e2_y = e1_x*e3_z - e1_z*e3_x ;
    e2_z = e1_y*e3_x - e1_x*e3_y ;
    vx = fabs(e2_x) ; vy = fabs(e2_y) ; vz = fabs(e2_z) ; /* use symmetry */
    if ((vx > vy) && (vx > vz))  /* scale using x component */
      len = 1.0f / vx ;
    else if (vy > vz)            /* scale using y component */
      len = 1.0f / vy ;
    else                         /* scale using z component */
      len = 1.0f / vz ;
    e2_x = *pe2_x++ = e2_x * len ;   
    e2_y = *pe2_y++ = e2_y * len ;
    e2_z = *pe2_z++ = e2_z * len ;
#if 1
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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIfindThinWMStrands(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int      width, height, depth, x, y, z, whalf, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi, thin, was_white, black_white, val ;
  float    nx, ny, nz, xf, yf, zf ;
  BUFTYPE  *pdst, *psrc ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  init_basis_vectors() ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  n = wsize*wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;  /* ptr to destination */
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      for (x = 0 ; x < width ; x++)
      {
        thin = 0 ;
        if (*psrc++) for (vertex = 0 ; !thin && vertex < NVERTICES ; vertex++)
        {
          was_white = -1 ;
          black_white = -wsize-1 ;
          nx = ic_x_vertices[vertex] ;  /* normal vector */
          ny = ic_y_vertices[vertex] ;
          nz = ic_z_vertices[vertex] ;
          
          xf = (float)x - wsize*nx ;
          yf = (float)y - wsize*ny ;
          zf = (float)z - wsize*nz ;
          for (yk = -wsize ; yk <= wsize ; yk++)
          {
            xi = pxi[nint(xf)] ;
            yi = pyi[nint(yf)] ;
            zi = pzi[nint(zf)] ;
            val = (float)MRIvox(mri_src, xi, yi, zi) ;
            if ((was_white > 0) && !val)   /* white to black transition */
            {
               /* check to see if we previously had a b-to-w transition */
              if (black_white >= -wsize)
              {
                thin = ((yk - black_white) <= wsize) ;
                break ;
              }
            }
            else if (!was_white && val > 0)  /* black to white transition */
              black_white = yk ;
            was_white = val > 0 ;
            
            xf += nx ; yf += ny ; zf += nz ;
          }
        }
          
        *pdst++ = thin ;
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
MRIorderThreshold(MRI *mri_src, MRI *mri_dst, MRI *mri_order, int num)
{
  int     width, height, depth, x, y, z, frame ;
  BUFTYPE *psrc, *pdst, *porder, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_src->type != MRI_UCHAR || mri_dst->type != MRI_UCHAR ||
      mri_order->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRIorderThresh: all input MRs must be byte")) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
        pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
        porder = &MRIseq_vox(mri_order, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++, psrc++)
        {
          if (*porder++ >= num)
            val = *psrc ;
          else
            val = 0 ;
          *pdst++ = val ;
        }
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
MRIpolvCount(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, 
                    int lo_lim, int hi_lim)
{
  int      width, height, depth, x, y, z, whalf, xk, yk, n, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi, order ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  BUFTYPE  *pdst, *pptr, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  whalf = (wsize-1)/2 ;

  init_basis_vectors() ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  n = wsize*wsize ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
      for (x = 0 ; x < width ; x++)
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

        /* now find the values in this plane */
        for (order = 0, yk = -whalf ; yk <= whalf ; yk++)
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
            val = MRIvox(mri_src, xi, yi, zi) ;
            if (val >= lo_lim && val <= hi_lim)
              order++ ;
          }
        }

        *pdst++ = (BUFTYPE)order ;
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

#define WHITE_LOW    90
#define GRAY_HI      95
#define WHITE_HI    130
#define WSIZE       9

#define PSLOPE       1.0
#define NSLOPE       1.0

MRI *
MRIwmfilter(MRI *mri_src, MRI *mri_polv, MRI *mri_dst)
{
  int      width, height, depth, x, y, z, whalf, vertex,xi,yi,zi,
           *pxi, *pyi, *pzi, i ;
  float    nx, ny, nz, dx, dy, dz, curv, dsq, val ;
  BUFTYPE  *pdst, *pptr, val0, /* *psrc, */gray_hi, white_low/*,mean, *pmean*/;
  MRI      *mri_curv /*, *mri_tmp*/ ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (DIAG_VERBOSE_ON)
    mri_curv = MRIalloc(width, height, depth, MRI_FLOAT) ;
  /*  mri_tmp = MRIalloc(width, height, depth, MRI_UCHAR) ;*/

  whalf = (WSIZE-1)/2 ;

  init_basis_vectors() ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

#if 0
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_tmp, 0, y, z) ;  /* ptr to destination */
      psrc = &MRIvox(mri_src, 0, y, z) ;  /* ptr to source */
      for (x = 0 ; x < width ; x++)
      {
        val0 = *psrc++ ;
        if (val0 > WHITE_HI || val0 < WHITE_LOW)
          val0 = 0 ;
        *pdst++ = val0 ;
      }
    }
  }
#endif

  pxi = mri_src->xi ; pyi = mri_src->yi ; pzi = mri_src->zi ;
  for (z = 0 ; z < depth ; z++)
  {
    DiagHeartbeat((float)z / (float)(depth-1)) ;
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
#if 0
      pmean = &MRIseq_vox(mri_polv, 0, y, z, 1) ; /* ptr to means */
#endif
      for (x = 0 ; x < width ; x++)
      {
        if (x == 37 && y == 88 && z == 63)
          DiagBreak() ;
        vertex = *pptr++ ;
#if 0
        mean = *pmean++ ;
#endif
        nx = ic_x_vertices[vertex] ;
        ny = ic_y_vertices[vertex] ;
        nz = ic_z_vertices[vertex] ;
        val0 = MRIvox(mri_src, x, y, z) ;

        /* now compute the curvature in the normal direction */
        for (curv = 0.0f, i = -whalf ; i <= whalf ; i++)
        {
          if (!i)
            continue ;
          dx = (float)i * nx ; xi = pxi[x+nint(dx)] ;
          dy = (float)i * ny ; yi = pyi[y+nint(dy)] ;
          dz = (float)i * nz ; zi = pzi[z+nint(dz)] ;
          dsq = (dx*dx + dy*dy + dz*dz) ;
          val = (float)MRIvox(mri_src, xi, yi, zi) - (float)val0 ;
          curv += val / dsq ;
        }
        curv /= (WSIZE-1) ;
        if (DIAG_VERBOSE_ON)
          MRIFvox(mri_curv, x, y, z) = curv ;

        white_low = WHITE_LOW ; gray_hi = GRAY_HI ;
#if 0
        if (curv < -6.0f)  /* gyrus */
          white_low = WHITE_LOW-10 ;
        else if (curv > 6)
          gray_hi = GRAY_HI+10 ;
#else
        if (curv < 0.0f)  /* gyrus */
          white_low += nint(NSLOPE*curv) ;
        else if (curv > 0.0f)
          gray_hi += nint(PSLOPE*curv) ;
#endif
        /*        val0 = mean ;*/
        if (val0 > WHITE_HI)   /* too high to be white matter */
          val0 = 0 ;
        else if (val0 < white_low)  /* too low to be white matter */
          val0 = 0 ;
        else if (val0 < gray_hi)    /* ambiguous */
        {
          if ((gray_hi - val0) > (val0 - white_low))
            val0 = 0 ;
        }
        *pdst++ = val0 ;
      }
    }
  }

  if (DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_curv, "curv.mnc") ;
    MRIfree(&mri_curv) ;
  }
  /*  MRIfree(&mri_tmp) ;*/
  return(mri_dst) ;
}

