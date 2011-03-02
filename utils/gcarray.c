/**
 * @file  gcarray.c
 * @brief utils for MRI classification using an array of Gaussian classifiers
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:44 $
 *    $Revision: 1.8 $
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
#include "gcarray.h"
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
int     GCarrayUpdate(GCARRAY *gcarray, MRI *mri_src,MRI *mri_norm,
                      MRI *mri_target);
int     GCarrayFinish(GCARRAY *gcarray) ;
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate an MRI structure, spacing the classifiers every
          scale pixels, with scale/2 at the left margin.

------------------------------------------------------*/
GCARRAY *
GCarrayAlloc(MRI *mri_template, int scale, int nvars)
{
  GCARRAY *gcarray ;
  int  x, y, z, width, height, depth ;
  double xw, yw, zw ;

  width = mri_template->width ;
  height = mri_template->height ;
  depth = mri_template->depth ;
  gcarray = (GCARRAY *)calloc(1, sizeof(GCARRAY)) ;
  if (!gcarray)
    ErrorReturn(NULL,
                (ERROR_NO_MEMORY, "GCarrayAlloc: could not alloc struct")) ;

  MRIvoxelToWorld(mri_template, 0, 0, 0, &xw, &yw, &zw) ;
  gcarray->xstart = xw ;
  gcarray->ystart = yw ;
  gcarray->zstart = zw ;
  gcarray->swidth = width ;
  gcarray->sheight = height ;
  gcarray->sdepth = depth ;

  width = nint((float)(width - scale/2) / (float)scale + 0.99f) ;
  height = nint((float)(height - scale/2) / (float)scale + 0.99f) ;
  depth = nint((float)(depth - scale/2) / (float)scale + 0.99f) ;

  gcarray->scale = scale ;
  gcarray->width = width ;
  gcarray->height = height ;
  gcarray->depth = depth ;
  gcarray->nvars = nvars ;

  gcarray->gcs = (GCLASSIFY ****)calloc(height, sizeof(GCLASSIFY ***)) ;
  if (!gcarray->gcs)
    ErrorExit(ERROR_NO_MEMORY, "GCarrayAlloc: could not allocate gcs") ;

  for (z = 0 ; z < depth ; z++)
  {
    gcarray->gcs[z] = (GCLASSIFY ***)calloc(height, sizeof(GCLASSIFY **)) ;
    if (!gcarray->gcs[z])
      ErrorExit(ERROR_NO_MEMORY,
                "GCarrayAlloc: could not allocate gcs[%d]", z) ;

    for (y = 0 ; y < height ; y++)
    {
      gcarray->gcs[z][y] = (GCLASSIFY **)calloc(width, sizeof(GCLASSIFY *)) ;
      if (!gcarray->gcs[z][y])
        ErrorExit(ERROR_NO_MEMORY,
                  "GCalloc(%d,%d,%d,%d): could not allocate gc[%d][%d]",
                  width, height, depth, scale, y, z) ;

      for (x = 0 ; x < width ; x++)
      {
        gcarray->gcs[z][y][x] = GCalloc(NCLASSES, nvars, class_names) ;
        if (!gcarray->gcs[z][y][x])
          ErrorExit(ERROR_NO_MEMORY,
                    "GCalloc(%d,%d,%d,%d): could not allocate gc[%d][%d][%d]",
                    width, height, depth, scale, x, y, z) ;

      }
    }
  }

  return(gcarray) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
GCarrayFree(GCARRAY **pgcarray)
{
  GCARRAY  *gcarray ;
  int  x, y, z ;

  gcarray = *pgcarray ;
  *pgcarray = NULL ;


  for (z = 0 ; z < gcarray->depth ; z++)
  {
    if (gcarray->gcs[z])
    {
      for (y = 0 ; y < gcarray->height ; y++)
      {
        if (gcarray->gcs[z][y])
        {
          for (x = 0 ; x < gcarray->width ; x++)
          {
            if (gcarray->gcs[z][y][x])
              GCfree(&gcarray->gcs[z][y][x]) ;
          }
          free(gcarray->gcs[z][y]) ;
        }
      }
      free(gcarray->gcs[z]) ;
    }
  }

  free(gcarray->gcs) ;
  free(gcarray) ;
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
GCARRAY *
GCarrayTrainAll(GCARRAY *gcarray,
                char *training_file_name,
                int scale,
                int ninputs)
{
  char  source_fname[100], target_fname[100], line[300], *cp ;
  FILE  *fp ;
  int   fno, nfiles ;
  MRI   *mri_src, *mri_target, *mri_std, *mri_zscore, *mri_mean, *mris[10] ;

  /* first figure out the total # of files */
  fp = fopen(training_file_name, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NO_FILE,
                 "GCarrayTrainAll(%s): could not open file",
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

    if (!fno && !gcarray)   /* allocate the GCARRAY */
      gcarray = GCarrayAlloc(mri_src, scale, ninputs) ;

    mri_mean = MRImean(mri_src, NULL, 3) ;
    mri_std = MRIstd(mri_src, NULL, mri_mean, 3) ;
    mri_zscore = MRIzScore(mri_src, NULL, mri_mean, mri_std) ;

    mris[0] = mri_src ;
    mris[1] = mri_zscore ;

    MRIfree(&mri_mean) ;
    MRIfree(&mri_std) ;

    mri_target = MRIread(target_fname) ;
    if (!mri_target)
    {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_src) ;
      MRIfree(&mri_zscore) ;
      continue ;
    }

    GCarrayUpdateMeans(gcarray, mris, mri_target, 2) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    MRIfree(&mri_zscore) ;
    fno++ ;
  }

  GCarrayComputeMeans(gcarray) ;  /* divide by # of observations */
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

    mri_mean = MRImean(mri_src, NULL, 3) ;
    mri_std = MRIstd(mri_src, NULL, mri_mean, 3) ;
    mri_zscore = MRIzScore(mri_src, NULL, mri_mean, mri_std) ;
    mris[0] = mri_src ;
    mris[1] = mri_zscore ;
    MRIfree(&mri_mean) ;
    MRIfree(&mri_std) ;

    mri_target = MRIread(target_fname) ;
    if (!mri_target)
    {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_zscore) ;
      MRIfree(&mri_src) ;
      continue ;
    }

    GCarrayUpdateCovariances(gcarray, mris, mri_target, 2) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    MRIfree(&mri_zscore) ;
    fno++ ;
  }

  GCarrayComputeCovariances(gcarray) ;

  fclose(fp) ;
  return(gcarray) ;
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
GCarrayUpdate(GCARRAY *gcarray, MRI *mri_src,MRI *mri_norm,MRI *mri_target)
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
GCarrayFinish(GCARRAY *gcarray)
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
GCarrayTrain(GCARRAY *gcarray, MRI *mri_src, MRI *mri_zscore, MRI *mri_target)
{
  MATRIX     *m_inputs[NCLASSES] ;
  GCLASSIFY  *gc, **pgc ;
  int        x, y, z, x0, y0, z0, x1, y1, z1, xm, ym, zm, xv, yv, zv,
  width, depth, height, scale, classno, nclasses, nobs[NCLASSES],
  swidth, sheight, sdepth, overlap, ninputs ;
  BUFTYPE    src, target ;
  double       xrv, yrv, zrv, xt, yt, zt ;

  scale = gcarray->scale ;
  overlap = MAX(1, scale/4) ;
  ninputs = scale +2*overlap+1 ;
  ninputs = ninputs * ninputs * ninputs ;
  for (classno = 0 ; classno < NCLASSES ; classno++)
    m_inputs[classno] = MatrixAlloc(ninputs,gcarray->nvars,MATRIX_REAL);

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;
  nclasses = NCLASSES ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

  /*
     the classifiers are distributed in Talairach space, whereas the
     input MR images are not (necessarily). Therefore we have to
     transform the voxel coordinates into Tal. space before selecting
     the training values.

     train each classifier, x,y,z are in classifier space, xm, ym, and
     zm are pixel coordinates in Talairach space
  */
  for (z = 0 ; z < depth ; z++)
  {
    z0 = MAX(0,z*scale - overlap) ;
    z1 = MIN(sdepth-1,(z+1)*scale+overlap) ;
    for (y = 0 ; y < height ; y++)
    {
      y0 = MAX(0,y*scale-overlap) ;
      y1 = MIN(sheight-1,(y+1)*scale+overlap) ;
      pgc = gcarray->gcs[z][y] ;
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
              /* convert to Talairach coords, then to voxel coords of src */
              /* convert from natural axes of coronal slices  Talairach axes */
              xt = (double)xm + gcarray->xstart ;   /* voxel to world */
              yt = (double)zm + gcarray->ystart ;
              zt = (double)-ym + gcarray->zstart ;
              MRItalairachToVoxel(mri_src, xt, yt, zt, &xrv, &yrv, &zrv) ;
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

              m_inputs[classno]->rptr[nobs[classno]+1][1] = src ;
              m_inputs[classno]->rptr[nobs[classno]+1][2] =
                MRIFvox(mri_zscore, xv, yv, zv) ;
              if (gcarray->nvars > 2)
              {
                m_inputs[classno]->rptr[nobs[classno]+1][3] = xt ;
                m_inputs[classno]->rptr[nobs[classno]+1][4] = yt ;
                m_inputs[classno]->rptr[nobs[classno]+1][5] = zt ;
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
GCarrayClassify(GCARRAY *gcarray, MRI *mri_src, MRI *mri_dst,
                float conf, MRI *mri_probs, MRI *mri_classes)
{
  MATRIX     *m_inputs ;
  GCLASSIFY  *gc ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno,
  nclasses, swidth, sheight, sdepth ;
  BUFTYPE    *psrc, src, *pdst, *pclasses ;
  float      *pzscore, prob, *pprobs = NULL ;
  double       xt, yt, zt, xv, yv, zv ;
  MRI        *mri_std, *mri_zscore, *mri_mean ;

  mri_mean = MRImean(mri_src, NULL, 3) ;
  mri_std = MRIstd(mri_src, mri_std, NULL, 3) ;
  mri_zscore = MRIzScore(mri_src, NULL, mri_mean, mri_std) ;
  MRIfree(&mri_std) ;
  MRIfree(&mri_mean) ;
  if (conf < 0.0f || conf >= 1.0f)
    conf = PRETTY_SURE ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  m_inputs = MatrixAlloc(gcarray->nvars, 1, MATRIX_REAL) ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;
  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;
  scale = gcarray->scale ;
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
        MRIvoxelToTalairach(mri_src,
                            (double)x, (double)y, (double)z,
                            &xt, &yt, &zt);

        /* convert from Talairach axes to natural axes of coronal slice data */
        xv = (xt - (double)gcarray->xstart) ;
        zv = (yt - (double)gcarray->ystart) ;
        yv = (-zt + (double)gcarray->zstart);
        xc = nint(((xv) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yv) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zv) - scale/2) / scale) ;
        if (zc < 0)
          zc = 0 ;
        else if (zc >= depth)
          zc = depth-1 ;
        gc = gcarray->gcs[zc][yc][xc] ;
        src = *psrc++ ;
        m_inputs->rptr[1][1] = src ;
        m_inputs->rptr[2][1] = *pzscore++ ;
        if (gcarray->nvars > 2)
        {
          m_inputs->rptr[3][1] = xt ;
          m_inputs->rptr[4][1] = yt ;
          m_inputs->rptr[5][1] = zt ;
        }

        /* now classify this observation */
        classno = GCclassify(gc, m_inputs, NULL, &prob) ;

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
  MRIfree(&mri_zscore) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
GCarrayToVoxel(GCARRAY *gcarray, int xc, int yc, int zc,
               int *pxv, int *pyv, int *pzv)
{
  int scale ;

  scale = gcarray->scale ;
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
GCarrayVoxelToClass(GCARRAY *gcarray, int xv, int yv, int zv,
                    int *pxc, int *pyc, int *pzc)
{
  int scale ;

  scale = gcarray->scale ;
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
GCarraySetTransform(GCARRAY *gcarray, Transform *transform,
                    Transform *inverse_transform)
{
  gcarray->transform = transform ;
  gcarray->inverse_transform = inverse_transform ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
GCARRAY *
GCarrayRead(char *fname)
{
  GCARRAY  *gcarray ;
  FILE  *fp ;
  int   width, height, depth, nvars, scale, x, y, z ;
  MRI   mri ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NO_FILE,"GCarrayRead(%s): could not open file",fname));

  if (fscanf(fp, "%d %d %d %d %d %d %d %d %f %f %f\n",
             &scale, &width, &height, &depth, &mri.width, &mri.height,
             &mri.depth, &nvars, &mri.xstart, &mri.ystart, &mri.zstart) != 11)
  {
    fclose(fp) ;
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "GCarrayRead(%s): could scanf parms from file",fname));
  }
  setDirectionCosine(&mri, MRI_CORONAL);
  mri.xsize = mri.ysize = mri.zsize = 1 ;
  gcarray = GCarrayAlloc(&mri, scale, nvars) ;
  if (!gcarray)
  {
    fclose(fp) ;
    ErrorReturn(NULL,
                (ERROR_BADPARM, "GCarrayRead(%s): gcarray allocation failed",
                 fname)) ;
  }
  gcarray->xstart = mri.xstart ;  /* a hack, I know, but what the hell... */
  gcarray->ystart = mri.ystart ;
  gcarray->zstart = mri.zstart ;

  for (z = 0 ; z < gcarray->depth ; z++)
  {
    for (y = 0 ; y < gcarray->height ; y++)
    {
      for (x = 0 ; x < gcarray->width ; x++)
      {
        GCasciiReadFrom(fp, gcarray->gcs[z][y][x]) ;
      }
    }
  }
  fclose(fp) ;
  return(gcarray) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
GCarrayWrite(GCARRAY *gcarray, char *fname)
{
  FILE  *fp ;
  int   x, y, z ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE,
                (ERROR_NO_FILE,"GCarrayWrite(%s): could not open file",fname));

  fprintf(fp, "%d %d %d %d %d %d %d %d %2.3f %2.3f %2.3f\n",
          gcarray->scale,
          gcarray->width,
          gcarray->height,
          gcarray->depth,
          gcarray->swidth,
          gcarray->sheight,
          gcarray->sdepth,
          gcarray->nvars,
          gcarray->xstart,
          gcarray->ystart,
          gcarray->zstart) ;

  for (z = 0 ; z < gcarray->depth ; z++)
  {
    for (y = 0 ; y < gcarray->height ; y++)
    {
      for (x = 0 ; x < gcarray->width ; x++)
      {
        GCasciiWriteInto(fp, gcarray->gcs[z][y][x]) ;
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
GCarrayThreshold(GCARRAY *gcarray, MRI *mri_probs, MRI *mri_classes,
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
#if 0
        else if (class == WHITE_MATTER)
          *pdst++ = 128 ;
#endif
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
MRIgetClassifier(GCARRAY *gcarray, MRI *mri, int x, int y, int z)
{
  GCLASSIFY  *gc ;
  double       xt, yt, zt ;
  int        width, depth, height, scale, xc, yc, zc, xv, yv, zv ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;
  scale = gcarray->scale ;

  /* find the appropriate classifier for this location */
  MRIvoxelToTalairach(mri, (double)x, (double)y, (double)z, &xt, &yt, &zt);

  /* convert from Talairach axes to natural axes of coronal slice data */
  xv = (xt - (double)gcarray->xstart) ;
  zv = (yt - (double)gcarray->ystart) ;
  yv = (-zt + (double)gcarray->zstart);
  xc = nint(((xv) - scale/2) / scale) ;
  if (xc < 0)
    xc = 0 ;
  else if (xc >= width)
    xc = width - 1 ;
  yc = nint(((yv) - scale/2) / scale) ;
  if (yc < 0)
    yc = 0 ;
  else if (yc >= height)
    yc = height-1 ;
  zc = nint(((zv) - scale/2) / scale) ;
  if (zc < 0)
    zc = 0 ;
  else if (zc >= depth)
    zc = depth-1 ;
  gc = gcarray->gcs[zc][yc][xc] ;
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
GCarrayUpdateMeans(GCARRAY *gcarray, MRI *mris[], MRI *mri_target, int nimages)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno,
  nclasses, swidth, sheight, sdepth, overlap ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      *pzscore ;
  double       xt, yt, zt, xv, yv, zv ;
  MRI        *mri_src ;

  mri_src = mris[0] ;

  scale = gcarray->scale ;
  overlap = MAX(1, scale/4) ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;
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
        MRIvoxelToTalairach(mri_src,
                            (double)x, (double)y, (double)z,
                            &xt, &yt, &zt);
        xv = (xt - (double)gcarray->xstart) ;
        zv = (yt - (double)gcarray->ystart) ;
        yv = (-zt + (double)gcarray->zstart);
        xc = nint(((xv) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yv) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zv) - scale/2) / scale) ;
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

        gc = gcarray->gcs[zc][yc][xc] ;
        gcl = &gc->classes[classno] ;
        gcl->nobs++ ;
        gcl->m_u->rptr[1][1] += (float)src ;
        gcl->m_u->rptr[2][1] += *pzscore++ ;
        if (gcarray->nvars > 2)
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
GCarrayUpdateCovariances(GCARRAY *gcarray,
                         MRI *mris[],
                         MRI *mri_target,
                         int nimages)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, xc, yc, zc, width, depth, height, scale, classno,
  nclasses, swidth, sheight, sdepth, overlap, col, row ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      *pzscore, obs[6], covariance ;
  double       xt, yt, zt, xv, yv, zv ;
  MRI        *mri_src ;

  mri_src = mris[0] ;

  scale = gcarray->scale ;
  overlap = MAX(1, scale/4) ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;
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
        MRIvoxelToTalairach(mri_src,
                            (double)x, (double)y, (double)z,
                            &xt, &yt, &zt);
        xv = (xt - (double)gcarray->xstart) ;
        zv = (yt - (double)gcarray->ystart) ;
        yv = (-zt + (double)gcarray->zstart);
        xc = nint(((xv) - scale/2) / scale) ;
        if (xc < 0)
          xc = 0 ;
        else if (xc >= width)
          xc = width - 1 ;
        yc = nint(((yv) - scale/2) / scale) ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= height)
          yc = height-1 ;
        zc = nint(((zv) - scale/2) / scale) ;
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

        gc = gcarray->gcs[zc][yc][xc] ;
        gcl = &gc->classes[classno] ;
        obs[1] = (float)src - gcl->m_u->rptr[1][1] ;
        obs[2] = *pzscore++ - gcl->m_u->rptr[2][1] ;
        if (gcarray->nvars > 2)
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
GCarrayComputeMeans(GCARRAY *gcarray)
{
  GCLASSIFY  *gc, **pgc ;
  GCLASS     *gcl ;
  int        x, y, z, width, depth, height, classno, nclasses, nobs, row ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;

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
      pgc = gcarray->gcs[z][y] ;
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
GCarrayComputeCovariances(GCARRAY *gcarray)
{
  GCLASSIFY  *gc, **pgc ;
  GCLASS     *gcl ;
  int        x, y, z, width, depth, height, classno, nclasses, row,col ;
  float      nobs, covariance ;

  width = gcarray->width ;
  height = gcarray->height ;
  depth = gcarray->depth ;

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
      pgc = gcarray->gcs[z][y] ;
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

