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
int     MRIclassUpdate(MRIC *mric, MRI *mri_src,MRI *mri_norm,MRI *mri_target);
int     MRIclassFinish(MRIC *mric) ;
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
          Do iterative training of a classifier. First go through each
          training pair and compute the means for the inputs, then
          go through again and compute the covariance matrices.

------------------------------------------------------*/
int     
MRIclassTrainAll(MRIC *mric, char *training_file_name)
{
  char  source_fname[100], target_fname[100], line[300], *cp ;
  FILE  *fp ;
  int   fno, nfiles ;
  MRI   *mri_src, *mri_target, *mri_std, *mri_zscore, *mri_mean, *mris[10] ;

  /* first figure out the total # of files */
  fp = fopen(training_file_name, "r") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE,
                 "MRIclassTrainAll(%s): could not open file",
                 training_file_name)) ;

  nfiles = 0 ;
  while ((cp = fgetl(line, 299, fp)) != NULL)
    nfiles++ ;
  fprintf(stderr, "processing %d files\n", nfiles) ;
  rewind(fp) ;

  /* now calculate means */
  fprintf(stderr, "computing means...\n") ;
  fno = 0 ;
  while ((cp = fgetl(line, 299, fp)) != NULL)
  {
    sscanf(cp, "%s %s", source_fname, target_fname) ;
    fprintf(stderr, "file[%d]: %s --> %s\n", fno, source_fname, target_fname);
    mri_src = MRIread(source_fname) ;
    if (!mri_src)
    {
      fprintf(stderr, "could not read MR image %s\n", source_fname) ;
      continue ;
    }
    if (!fno)
    {
      mric->xstart = nint(mri_src->xstart) ;
      mric->ystart = nint(mri_src->ystart) ;
      mric->zstart = nint(mri_src->zstart) ;
    }
    mri_target = MRIread(target_fname) ;
    if (!mri_target)
    {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_src) ;
      continue ;
    }

    mri_mean = MRImean(mri_src, NULL, 3) ;
    mri_std = MRIstd(mri_src, NULL, mri_mean, 3) ;
    mri_zscore = MRInorm(mri_src, NULL, mri_mean, mri_std) ;
    mris[0] = mri_src ;
    mris[1] = mri_zscore ;
    MRIclassUpdateMeans(mric, mris, mri_target, 2) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    MRIfree(&mri_zscore) ;
    MRIfree(&mri_std) ;
    MRIfree(&mri_mean) ;
    fno++ ;
  }

  MRIclassComputeMeans(mric) ;  /* divide by # of observations */
  rewind(fp) ;

  /* now calculate covariances */
  fprintf(stderr, "computing covariance matrices...\n") ;
  fno = 0 ;
  while ((cp = fgetl(line, 299, fp)) != NULL)
  {
    sscanf(cp, "%s %s", source_fname, target_fname) ;
    fprintf(stderr, "file[%d]: %s --> %s\n", fno, source_fname, target_fname);
    mri_src = MRIread(source_fname) ;
    if (!mri_src)
    {
      fprintf(stderr, "could not read MR image %s\n", source_fname) ;
      continue ;
    }
    mri_target = MRIread(target_fname) ;
    if (!mri_target)
    {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_src) ;
      continue ;
    }

    mri_mean = MRImean(mri_src, NULL, 3) ;
    mri_std = MRIstd(mri_src, NULL, mri_mean, 3) ;
    mri_zscore = MRInorm(mri_src, NULL, mri_mean, mri_std) ;
    mris[0] = mri_src ;
    mris[1] = mri_zscore ;
    MRIclassUpdateCovariances(mric, mris, mri_target, 2) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    MRIfree(&mri_zscore) ;
    MRIfree(&mri_std) ;
    MRIfree(&mri_mean) ;
    fno++ ;
  }

  MRIclassComputeCovariances(mric) ;

  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Do iterative training of a classifier. First go through

          allow each classifier to be trained on a part of 
          the space that it's neighbor is responsible (i.e. 
          use overlapping training regions). The overlap is
          defined to be MAX(1, scale/4) on each side of the
          region.
------------------------------------------------------*/
int
MRIclassUpdate(MRIC *mric, MRI *mri_src,MRI *mri_norm,MRI *mri_target)
{
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
int     
MRIclassFinish(MRIC *mric)
{
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
#define LO_LIM  60
#define HI_LIM  150

int
MRIclassTrain(MRIC *mric, MRI *mri_src, MRI *mri_zscore, MRI *mri_target)
{
  MATRIX     *m_inputs[NCLASSES] ;
  GCLASSIFY  *gc, **pgc ;
  int        x, y, z, x0, y0, z0, x1, y1, z1, xm, ym, zm, xv, yv, zv,
             width, depth, height, scale, classno, nclasses, nobs[NCLASSES],
             swidth, sheight, sdepth, overlap, ninputs ;
  BUFTYPE    src, target ;
  Real       xrv, yrv, zrv ;

  scale = mric->scale ;
  overlap = MAX(1, scale/4) ;
  ninputs = scale +2*overlap+1 ;
  ninputs = ninputs * ninputs * ninputs ;
  for (classno = 0 ; classno < NCLASSES ; classno++)
    m_inputs[classno] = MatrixAlloc(ninputs,mric->nvars,MATRIX_REAL);

  MRItalairachVoxelToVoxel(mri_src, 0, 0, 0, &xrv, &yrv, &zrv) ;
  mric->xstart = (float)xrv ;
  mric->ystart = (float)yrv ;
  mric->zstart = (float)zrv ;
  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  nclasses = NCLASSES ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/
  /* train each classifier, x,y,z are in classifier space */
  for (z = 0 ; z < depth ; z++)
  {
    z0 = MAX(0,z*scale - overlap) ;
    z1 = MIN(sdepth-1,(z+1)*scale+overlap) ;
    for (y = 0 ; y < height ; y++)
    {
      y0 = MAX(0,y*scale-overlap) ;
      y1 = MIN(sheight-1,(y+1)*scale+overlap) ;
      pgc = mric->gcs[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        gc = *pgc++ ;
        x0 = MAX(0,x*scale-overlap);
        x1 = MIN(swidth-1,(x+1)*scale+overlap) ;
        
        memset(nobs, 0, NCLASSES*sizeof(nobs[0])) ;
        for (zm = z0 ; zm <= z1 ; zm++)
        {
          for (ym = y0 ; ym <= y1 ; ym++)
          {
            for (xm = x0 ; xm <= x1 ; xm++)
            {
              MRItalairachVoxelToVoxel(mri_src, (Real)xm, (Real)ym, (Real)zm,
                                       &xrv, &yrv, &zrv) ;
              xv = nint(xrv) ;
              yv = nint(yrv) ;
              zv = nint(zrv) ;
              if (xv < 0 || xv >= swidth ||
                  yv < 0 || yv >= sheight ||
                  zv < 0 || zv >= sdepth)
                continue ;
                                       
              src = MRIvox(mri_src, xv, yv, zv) ;
              target = MRIvox(mri_target, xv, yv, zv) ;
              
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
#if 0
if ((xv == 22 && yv == 27 && zv == 31) ||
    (xv == 22 && yv == 26 && zv == 31))
  fprintf(stderr, 
          "TRAIN: (%d, %d, %d) --> (%d, %d, %d) --> (%d, %d, %d) = (%d, %2.3f) = class %d\n",
      xm, ym, zm, x, y, z, xv, yv, zv, src, MRIFvox(mri_zscore, xv, yv, zv),
          classno) ;
#endif

              m_inputs[classno]->rptr[nobs[classno]+1][1] = src ;
              m_inputs[classno]->rptr[nobs[classno]+1][2] = 
                MRIFvox(mri_zscore, xv, yv, zv) ;
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
MRIclassify(MRIC *mric, MRI *mri_src, MRI *mri_zscore, MRI *mri_dst, 
            float conf, MRI *mri_probs, MRI *mri_classes)
{
  MATRIX     *m_inputs ;
  GCLASSIFY  *gc ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno, 
             nclasses, swidth, sheight, sdepth ;
  BUFTYPE    *psrc, src, *pdst, *pclasses ;
  float      *pzscore, prob, *pprobs = NULL ;
  Real       xt, yt, zt, xoff, yoff, zoff ;

#if 0
  if (mric->swidth != mri_src->width || mric->sheight != mri_src->height ||
      mric->sdepth != mri_src->depth)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIclassify: MRI does not match classifier"));
#endif

  if (conf < 0.0f || conf >= 1.0f)
    conf = PRETTY_SURE ;

  scale = mric->scale ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  m_inputs = MatrixAlloc(mric->nvars, 1, MATRIX_REAL) ;

  xoff = (Real)mri_src->xstart - (Real)mric->xstart ;
  yoff = (Real)mri_src->ystart - (Real)mric->ystart ;
  zoff = (Real)mri_src->zstart - (Real)mric->zstart ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;
  scale = mric->scale ;
  nclasses = NCLASSES ;

/* 
   x, y, and z are in the MR image space. To get the appropriate classifier
   for each spatial location we must convert to Talairach voxel coords.,
   then find the classifier assigned to that location in Talairach space.
   
   xc, yc, and zc are in classifier space, while xvt, yvt and zvt are
   in Talairach voxel coords.
*/
  for (z = 0 ; z < sdepth ; z++)
  {
    for (y = 0 ; y < sheight ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pzscore = &MRIFvox(mri_zscore, 0, y, z) ;
      if (mri_probs)
        pprobs = &MRIFvox(mri_probs, 0, y, z) ;
      else
        pprobs = NULL ;
      if (mri_classes)
        pclasses = &MRIvox(mri_classes, 0, y, z) ;
      else
        pclasses = NULL ;
      for (x = 0 ; x < swidth ; x++)
      {
        /* find the appropriate classifier for this location */
        MRIvoxelToTalairachVoxel(mri_src, (Real)x, (Real)y, (Real)z,
                                 &xt, &yt, &zt) ;
        xc = nint(((xt+xoff) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yt+yoff) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zt+zoff) - scale/2) / scale) ;
        if (zc < 0)
          zc = 0 ;
        else if (zc >= depth)
          zc = depth-1 ;
        gc = mric->gcs[zc][yc][xc] ;
        src = *psrc++ ;
        m_inputs->rptr[1][1] = src ;
        m_inputs->rptr[2][1] = *pzscore++ ;
        if (mric->nvars > 2)
        {
          m_inputs->rptr[3][1] = x ;
          m_inputs->rptr[4][1] = y ;
          m_inputs->rptr[5][1] = z ;
        }
        
        /* now classify this observation */
        classno = GCclassify(gc, m_inputs, &prob) ;
        if (pclasses)
          *pclasses++ = (BUFTYPE)classno ;
        if (pprobs)
          *pprobs++ = prob ;
        if (classno == WHITE_MATTER && prob > conf)
          *pdst++ = 255 ;
        else
          *pdst++ = src ;
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
  int   width, height, depth, swidth, sheight, sdepth, nvars, scale, x, y, z,
        xstart, ystart, zstart ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NO_FILE,"MRIclassRead(%s): could not open file",fname));

  if (fscanf(fp, "%d %d %d %d %d %d %d %d %d %d %d\n", 
             &scale, &width, &height, &depth, &swidth, &sheight, &sdepth, 
             &nvars, &xstart, &ystart, &zstart) != 11)
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

  mric->xstart = xstart ;
  mric->ystart = ystart ;
  mric->zstart = zstart ;
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

  fprintf(fp, "%d %d %d %d %d %d %d %d %d %d %d\n",
          mric->scale,
          mric->width,
          mric->height,
          mric->depth,
          mric->swidth,
          mric->sheight,
          mric->sdepth,
          mric->nvars,
          mric->xstart,
          mric->ystart,
          mric->zstart) ;

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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIclassThreshold(MRIC *mric, MRI *mri_probs, MRI *mri_classes,
                  MRI *mri_dst, float threshold)
{
  int      x, y, z, width, height, depth, class ;
  float    *pprobs, prob ;
  BUFTYPE  *pclasses, *pdst ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_classes, NULL) ;

  width = mri_classes->width ;
  height = mri_classes->height ;
  depth = mri_classes->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pprobs = &MRIFvox(mri_probs, 0, y, z) ;
      pclasses = &MRIvox(mri_classes, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        prob = *pprobs++ ;
        class = (int)*pclasses++ ;
        if (prob >= threshold && class == WHITE_MATTER)
          *pdst++ = 255 ;
        else if (class == WHITE_MATTER)
          *pdst++ = 128 ;
        else
          *pdst++ = 0 ;
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
GCLASSIFY *
MRIgetClassifier(MRIC *mric, MRI *mri, int xv, int yv, int zv)
{
  GCLASSIFY  *gc ;
  Real       xt, yt, zt ;
  int        width, depth, height, scale, xc, yc, zc ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  scale = mric->scale ;
  MRIvoxelToTalairachVoxel(mri, (Real)xv, (Real)yv, (Real)zv, &xt, &yt, &zt) ;
  xc = (int)((xt - scale/2) / scale) ;
  if (xc < 0)
    xc = 0 ;
  else if (xc >= width)
    xc = width - 1 ;
  yc = (int)((yt - scale/2) / scale) ;
  if (yc < 0)
    yc = 0 ;
  else if (yc >= height)
    yc = height-1 ;
  zc = (int)((zt - scale/2) / scale) ;
  if (zc < 0)
    zc = 0 ;
  else if (zc >= depth)
    zc = depth-1 ;
  gc = mric->gcs[zc][yc][xc] ;
#if 1
fprintf(stderr, 
        "classifier at (%d, %d, %d) --> (%d, %d, %d) is (%d, %d, %d)\n",
        xv, yv, zv, (int)xt, (int)yt, (int)zt, xc, yc, zc) ;
#endif

  return(gc) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          update the means using this set of images
------------------------------------------------------*/
int
MRIclassUpdateMeans(MRIC *mric, MRI *mris[], MRI *mri_target, int nimages)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno, 
             nclasses, swidth, sheight, sdepth, overlap ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      *pzscore ;
  Real       xt, yt, zt, xoff, yoff, zoff ;
  MRI        *mri_src ;

  mri_src = mris[0] ;

  xoff = (Real)mri_src->xstart - (Real)mric->xstart ;
  yoff = (Real)mri_src->ystart - (Real)mric->ystart ;
  zoff = (Real)mri_src->zstart - (Real)mric->zstart ;

  scale = mric->scale ;
  overlap = MAX(1, scale/4) ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  nclasses = NCLASSES ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/

/* 
   for each point in the image, find the its Talairach coordinates and
   therefore the classifier responsible for it, and update it's mean.
*/
  for (z = 0 ; z < sdepth ; z++)
  {
    for (y = 0 ; y < sheight ; y++)
    {
      psrc = &MRIvox(mris[0], 0, y, z) ;
      ptarget = &MRIvox(mri_target, 0, y, z) ;
      pzscore = &MRIFvox(mris[1], 0, y, z) ;
      for (x = 0 ; x < swidth ; x++)
      {
        /* find the appropriate classifier for this location */
        MRIvoxelToTalairachVoxel(mri_src, (Real)x, (Real)y, (Real)z, 
                                 &xt, &yt, &zt) ;
        xc = nint(((xt+xoff) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yt+yoff) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zt+zoff) - scale/2) / scale) ;
        if (zc < 0)
          zc = 0 ;
        else if (zc >= depth)
          zc = depth-1 ;

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
        
        gc = mric->gcs[zc][yc][xc] ;
        gcl = &gc->classes[classno] ;
        gcl->nobs++ ;
        gcl->m_u->rptr[1][1] += (float)src ;
        gcl->m_u->rptr[2][1] += *pzscore++ ;
        if (mric->nvars > 2)
        {
          gcl->m_u->rptr[3][1] += (float)xt ;
          gcl->m_u->rptr[4][1] += (float)yt ;
          gcl->m_u->rptr[5][1] += (float)zt ;
        }
      }
    }
  }

  return(ERROR_NONE) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           update the covariance estimates based on new observations
------------------------------------------------------*/
int
MRIclassUpdateCovariances(MRIC *mric, MRI *mris[],MRI *mri_target,int nimages)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno, 
             nclasses, swidth, sheight, sdepth, overlap, col, row ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      *pzscore, obs[6], covariance ;
  Real       xt, yt, zt, xoff, yoff, zoff ;
  MRI        *mri_src ;

  mri_src = mris[0] ;

  xoff = (Real)mri_src->xstart - (Real)mric->xstart ;
  yoff = (Real)mri_src->ystart - (Real)mric->ystart ;
  zoff = (Real)mri_src->zstart - (Real)mric->zstart ;

  scale = mric->scale ;
  overlap = MAX(1, scale/4) ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;
  nclasses = NCLASSES ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/

/* 
   for each point in the image, find the its Talairach coordinates and
   therefore the classifier responsible for it, and update it's mean.
*/
  for (z = 0 ; z < sdepth ; z++)
  {
    for (y = 0 ; y < sheight ; y++)
    {
      psrc = &MRIvox(mris[0], 0, y, z) ;
      ptarget = &MRIvox(mri_target, 0, y, z) ;
      pzscore = &MRIFvox(mris[1], 0, y, z) ;
      for (x = 0 ; x < swidth ; x++)
      {
        /* find the appropriate classifier for this location */
        MRIvoxelToTalairachVoxel(mri_src, (Real)x, (Real)y, (Real)z, 
                                 &xt, &yt, &zt) ;
        xc = nint(((xt+xoff) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yt+yoff) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zt+zoff) - scale/2) / scale) ;
        if (zc < 0)
          zc = 0 ;
        else if (zc >= depth)
          zc = depth-1 ;

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
        
        gc = mric->gcs[zc][yc][xc] ;
        gcl = &gc->classes[classno] ;
        obs[1] = (float)src - gcl->m_u->rptr[1][1] ;
        obs[2] = *pzscore++ - gcl->m_u->rptr[2][1] ;
        if (mric->nvars > 2)
        {
          obs[3] = (float)xt - gcl->m_u->rptr[3][1] ;
          obs[4] = (float)yt - gcl->m_u->rptr[4][1] ;
          obs[5] = (float)zt - gcl->m_u->rptr[5][1] ;
        }
        for (row = 1 ; row <= gcl->m_covariance->rows ; row++)
        {
          for (col = 1 ; col <= row ; col++)
          {
            covariance = obs[row] * obs[col] ;
            gcl->m_covariance->rptr[row][col] += covariance;
            gcl->m_covariance->rptr[col][row] += covariance;
          }
        }
      }
    }
  }

  return(ERROR_NONE) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           compute the means for each class
------------------------------------------------------*/
int
MRIclassComputeMeans(MRIC *mric)
{
  GCLASSIFY  *gc, **pgc ;
  GCLASS     *gcl ;
  int        x, y, z, width, depth, height, classno, nclasses, nobs, row ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/
  /* train each classifier, x,y,z are in classifier space */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pgc = mric->gcs[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        gc = *pgc++ ;
        nclasses = gc->nclasses ;
        for (classno = 0 ; classno < nclasses ; classno++)
        {
          gcl = &gc->classes[classno] ;
          nobs = gcl->nobs ;
          if (nobs)
          {
            for (row = 1 ; row <= gcl->m_u->rows ; row++)
              gcl->m_u->rptr[row][1] /= (float)nobs ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           compute the means for each class
------------------------------------------------------*/
int
MRIclassComputeCovariances(MRIC *mric)
{
  GCLASSIFY  *gc, **pgc ;
  GCLASS     *gcl ;
  int        x, y, z, width, depth, height, classno, nclasses, row,col ;
  float      nobs, covariance ;

  width = mric->width ;
  height = mric->height ;
  depth = mric->depth ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/
  /* train each classifier, x,y,z are in classifier space */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pgc = mric->gcs[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        gc = *pgc++ ;
        nclasses = gc->nclasses ;
        for (classno = 0 ; classno < nclasses ; classno++)
        {
          gcl = &gc->classes[classno] ;
          nobs = (float)(gcl->nobs-1) ;  /* ML estimate of covariance */
          if (nobs > 0.0f)
          {
            for (row = 1 ; row <= gcl->m_covariance->rows ; row++)
            {
              for (col = 1 ; col <= row ; col++)
              {
                covariance = gcl->m_covariance->rptr[row][col] / nobs ;
                gcl->m_covariance->rptr[row][col] = covariance ;
                gcl->m_covariance->rptr[col][row] = covariance ;
              }
            }
          }
          GCinit(gc, classno) ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

