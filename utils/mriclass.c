/*
 *       FILE NAME:   mriclass.c
 *
 *       DESCRIPTION: utilities for MRI classification
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
#include "mriclass.h"
#include "matrix.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/


#define NCLASSES      4
#define BACKGROUND    0
#define GREY_MATTER   1
#define WHITE_MATTER  2
#define BRIGHT_MATTER 3

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static char *class_names[NCLASSES] =
{
  "BACKGROUND",
  "GREY MATTER",
  "WHITE MATTER",
  "BRIGHT MATTER"
} ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/


/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate an MRI structure, spacing the classifiers every
          scale pixels, with scale/2 at the left margin.

------------------------------------------------------*/
MRIC *
MRIclassAlloc(int width, int height, int depth, int scale, int nvars)
{
  MRIC *mric ;
  int  x, y, z ;

  mric = (MRIC *)calloc(1, sizeof(MRIC)) ;
  if (!mric)
    ErrorReturn(NULL, 
                (ERROR_NO_MEMORY, "MRIalloc: could not alloc struct")) ;

  mric->swidth = width ;
  mric->sheight = height ;
  mric->sdepth = depth ;

  width = nint((float)(width - scale/2) / (float)scale + 0.99f) ;
  height = nint((float)(height - scale/2) / (float)scale + 0.99f) ;
  depth = nint((float)(depth - scale/2) / (float)scale + 0.99f) ;
        
  mric->scale = scale ;
  mric->width = width ;
  mric->height = height ;
  mric->depth = depth ;
  mric->nvars = nvars ;

  mric->gcs = (GCLASSIFY ****)calloc(height, sizeof(GCLASSIFY ***)) ;
  if (!mric->gcs)
    ErrorExit(ERROR_NO_MEMORY, "MRIclassAlloc: could not allocate gcs") ;

  for (z = 0 ; z < depth ; z++)
  {
    mric->gcs[z] = (GCLASSIFY ***)calloc(height, sizeof(GCLASSIFY **)) ;
    if (!mric->gcs[z])
      ErrorExit(ERROR_NO_MEMORY, 
                "MRIclassAlloc: could not allocate gcs[%d]", z) ;

    for (y = 0 ; y < height ; y++)
    {
      mric->gcs[z][y] = (GCLASSIFY **)calloc(width, sizeof(GCLASSIFY *)) ;
      if (!mric->gcs[z][y])
        ErrorExit(ERROR_NO_MEMORY,
                  "GCalloc(%d,%d,%d,%d): could not allocate gc[%d][%d]",
                  width, height, depth, scale, y, z) ;

      for (x = 0 ; x < width ; x++)
      {
        mric->gcs[z][y][x] = GCalloc(NCLASSES, nvars, class_names) ;
        if (!mric->gcs[z][y][x])
          ErrorExit(ERROR_NO_MEMORY,
                    "GCalloc(%d,%d,%d,%d): could not allocate gc[%d][%d][%d]",
                    width, height, depth, scale, x, y, z) ;

      }
    }
  }

  return(mric) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRIclassFree(MRIC **pmric)
{
  MRIC  *mric ;
  int  x, y, z ;

  mric = *pmric ;
  *pmric = NULL ;


  for (z = 0 ; z < mric->depth ; z++)
  {
    if (mric->gcs[z])
    {
      for (y = 0 ; y < mric->height ; y++)
      {
        if (mric->gcs[z][y])
        {
          for (x = 0 ; x < mric->width ; x++)
          {
            if (mric->gcs[z][y][x])
              GCfree(&mric->gcs[z][y][x]) ;
          }
          free(mric->gcs[z][y]) ;
        }
      }
      free(mric->gcs[z]) ;
    }
  }

  free(mric->gcs) ;
  free(mric) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          allow each classifier to be trained on a part of 
          the space that it's neighbor is responsible (i.e. 
          use overlapping training regions). The overlap is
          defined to be MAX(1, scale/4) on each side of the
          region.
------------------------------------------------------*/
#define LO_LIM  50
#define HI_LIM  150

int
MRIclassTrain(MRIC *mric, MRI *mri_src, MRI *mri_norm, MRI *mri_target)
{
  MATRIX     *m_inputs[NCLASSES] ;
  GCLASSIFY  *gc, **pgc ;
  int        x, y, z, x0, y0, z0, x1, y1, z1, xm, ym, zm, 
             width, depth, height, scale, classno, nclasses, nobs[NCLASSES],
             swidth, sheight, sdepth, overlap, ninputs ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      *pnorm ;

  scale = mric->scale ;
  overlap = MAX(1, scale/4) ;
  ninputs = scale +2*overlap+1 ;
  ninputs = ninputs * ninputs * ninputs ;
  for (classno = 0 ; classno < NCLASSES ; classno++)
    m_inputs[classno] = MatrixAlloc(ninputs,mric->nvars,MATRIX_REAL);

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  scale = mric->scale ;
  nclasses = NCLASSES ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

  /* train each classifier, x,y,z are in classifier space */
  for (z = 0 ; z < depth ; z++)
  {
    z0 = MAX(0,z*scale - overlap) ;
    z1 = MIN(sdepth-1,z0+scale+2*overlap) ;
    for (y = 0 ; y < height ; y++)
    {
      y0 = MAX(0,y*scale-overlap) ;
      y1 = MIN(sheight-1,y0+scale+2*overlap) ;
      pgc = mric->gcs[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        gc = *pgc++ ;
        x0 = MAX(0,x*scale-overlap);
        x1 = MIN(swidth-1,x0+scale+2*overlap) ;
        
        memset(nobs, 0, NCLASSES*sizeof(nobs[0])) ;
        for (zm = z0 ; zm <= z1 ; zm++)
        {
          for (ym = y0 ; ym <= y1 ; ym++)
          {
            psrc = &MRIvox(mri_src, x0, ym, zm) ;
            ptarget = &MRIvox(mri_target, x0, ym, zm) ;
            pnorm = &MRIFvox(mri_norm, x0, ym, zm) ;
            for (xm = x0 ; xm <= x1 ; xm++)
            {
              src = *psrc++ ;
              target = *ptarget++ ;
              
              /* decide what class it is */
              if (target)
                classno = WHITE_MATTER ;
              else
              {
                if (src < LO_LIM)
                  classno = BACKGROUND ;
                else if (src > HI_LIM)
                  classno = BRIGHT_MATTER ;
                else
                  classno = GREY_MATTER ;
              }
              m_inputs[classno]->rptr[nobs[classno]+1][1] = src ;
              m_inputs[classno]->rptr[nobs[classno]+1][2] = *pnorm++ ;
              if (mric->nvars > 2)
              {
                m_inputs[classno]->rptr[nobs[classno]+1][3] = xm ;
                m_inputs[classno]->rptr[nobs[classno]+1][4] = ym ;
                m_inputs[classno]->rptr[nobs[classno]+1][5] = zm ;
              }
              nobs[classno]++ ;
            }
          }
        }

        /* now apply training vectors */
        for (classno = 0 ; classno < nclasses ; classno++)
        {
          m_inputs[classno]->rows = nobs[classno] ;
          GCtrain(gc, classno, m_inputs[classno]) ;
          m_inputs[classno]->rows = scale*scale*scale ;
        }
      }
    }
  }

  for (classno = 0 ; classno < NCLASSES ; classno++)
    MatrixFree(&m_inputs[classno]) ;

  return(ERROR_NONE) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define PRETTY_SURE   .90f

MRI *
MRIclassify(MRIC *mric, MRI *mri_src, MRI *mri_norm, MRI *mri_dst, float conf)
{
  MATRIX     *m_inputs ;
  GCLASSIFY  *gc, **pgc ;
  int        x, y, z, x0, y0, z0, x1, y1, z1, xm, ym, zm, 
             width, depth, height, scale, classno, nclasses,
             swidth, sheight, sdepth ;
  BUFTYPE    *psrc, src, *pdst ;
  float      *pnorm, prob ;

  if (mric->swidth != mri_src->width || mric->sheight != mri_src->height ||
      mric->sdepth != mri_src->depth)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIclassify: MRI does not match classifier"));

  if (conf < 0.0f || conf >= 1.0f)
    conf = PRETTY_SURE ;

  scale = mric->scale ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  m_inputs = MatrixAlloc(mric->nvars, 1, MATRIX_REAL) ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;
  scale = mric->scale ;
  nclasses = NCLASSES ;

  /* train each classifier, x,y,z are in classifier space */
  for (z = 0 ; z < depth ; z++)
  {
    z0 = MAX(0,z*scale) ;
    z1 = MIN(sdepth-1,z0+scale) ;
    for (y = 0 ; y < height ; y++)
    {
      y0 = MAX(0,y*scale) ;
      y1 = MIN(sheight-1,y0+scale) ;
      pgc = mric->gcs[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        x0 = MAX(0,x*scale);
        x1 = MIN(swidth-1,x0+scale) ;
        gc = *pgc++ ;

        for (zm = z0 ; zm <= z1 ; zm++)
        {
          for (ym = y0 ; ym <= y1 ; ym++)
          {
            psrc = &MRIvox(mri_src, x0, ym, zm) ;
            pdst = &MRIvox(mri_dst, x0, ym, zm) ;
            pnorm = &MRIFvox(mri_norm, x0, ym, zm) ;
            for (xm = x0 ; xm <= x1 ; xm++)
            {
              src = *psrc++ ;
              m_inputs->rptr[1][1] = src ;
              m_inputs->rptr[2][1] = *pnorm++ ;
              if (mric->nvars > 2)
              {
                m_inputs->rptr[3][1] = xm ;
                m_inputs->rptr[4][1] = ym ;
                m_inputs->rptr[5][1] = zm ;
              }

              /* now classify this observation */
              classno = GCclassify(gc, m_inputs, &prob) ;
              if (classno == WHITE_MATTER && prob > conf)
                *pdst++ = 255 ;
              else
                *pdst++ = src ;
            }
          }
        }

      }
    }
  }

  MatrixFree(&m_inputs) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRIclassToVoxel(MRIC *mric, int xc, int yc, int zc,
                int *pxv, int *pyv, int *pzv)
{
  int scale ;

  scale = mric->scale ;
  *pxv = xc*scale + scale/2 ;
  *pyv = yc*scale + scale/2 ;
  *pzv = zc*scale + scale/2 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRIvoxelToClass(MRIC *mric, int xv, int yv, int zv,
                int *pxc, int *pyc, int *pzc)
{
  int scale ;

  scale = mric->scale ;
  *pxc = nint((float)(xv - scale/2) / (float)scale + .99f) ;
  *pyc = nint((float)(yv - scale/2) / (float)scale + .99f) ;
  *pzc = nint((float)(zv - scale/2) / (float)scale + .99f) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRIclassSetTransform(MRIC *mric, Transform *transform, 
                     Transform *inverse_transform)
{
  mric->transform = transform ;
  mric->inverse_transform = inverse_transform ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRIC *
MRIclassRead(char *fname)
{
  MRIC  *mric ;
  FILE  *fp ;
  int   width, height, depth, swidth, sheight, sdepth, nvars, scale ;
  int   x, y, z ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NO_FILE,"MRIclassRead(%s): could not open file",fname));

  if (fscanf(fp, "%d %d %d %d %d %d %d %d\n", &scale, &width, &height, &depth,
             &swidth, &sheight, &sdepth, &nvars) != 8)
  {
    fclose(fp) ;
    ErrorReturn(NULL, 
                (ERROR_BADFILE, 
                 "MRIclassRead(%s): could scanf parms from file",fname));
  }

  mric = MRIclassAlloc(swidth, sheight, sdepth, scale, nvars) ;
  if (!mric)
  {
    fclose(fp) ;
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIclassRead(%s): mric allocation failed",
                 fname)) ;
  }

  for (z = 0 ; z < mric->depth ; z++)
  {
    for (y = 0 ; y < mric->height ; y++)
    {
      for (x = 0 ; x < mric->width ; x++)
      {
        GCasciiReadFrom(fp, mric->gcs[z][y][x]) ;
      }
    }
  }
  fclose(fp) ;
  return(mric) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRIclassWrite(MRIC *mric, char *fname)
{
  FILE  *fp ;
  int   x, y, z ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
              (ERROR_NO_FILE,"MRIclassWrite(%s): could not open file",fname));

  fprintf(fp, "%d %d %d %d %d %d %d %d\n",
          mric->scale,
          mric->width,
          mric->height,
          mric->depth,
          mric->swidth,
          mric->sheight,
          mric->sdepth,
          mric->nvars) ;

  for (z = 0 ; z < mric->depth ; z++)
  {
    for (y = 0 ; y < mric->height ; y++)
    {
      for (x = 0 ; x < mric->width ; x++)
      {
        GCasciiWriteInto(fp, mric->gcs[z][y][x]) ;
      }
    }
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

