/*
 *       FILE NAME:   mriclass.c
 *
 *       DESCRIPTION: utilities for MRI classification using
 *                    a variety of classifiers
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
#include "backprop.h"
#include "artmap.h"
#include "gclass.h"
#include "mriclass.h"
#include "utils.h"
#include "region.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/


/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/
static char *gaussian_class_names[GAUSSIAN_NCLASSES] =
{
  "BACKGROUND",
  "GREY MATTER",
  "WHITE MATTER",
  "BRIGHT MATTER"
} ;

static int  total = 0, buffered = 0, total_computed = 0 ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/


/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
int
MRICfree(MRIC **pmric)
{
  MRIC  *mric ;
  int   round ;

  mric = *pmric ;
  *pmric = NULL ;
  if (mric)
  {
    for (round = 0 ; round < mric->nrounds ; round++)
    {
      switch (mric->type[round])
      {
      case CLASSIFIER_GAUSSIAN:
        GCfree(&mric->classifier[round].gc) ;
        break ;
      default:
        break ;
      }
    }
    free(mric) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *
MRICalloc(int nrounds, int *types, int *features, void *parms)
{
  MRIC  *mric ;
  int   f, ninputs, round ;

  mric = (MRIC *)calloc(1, sizeof(MRIC)) ;
  if (!mric)
    ErrorExit(ERROR_NO_MEMORY, "MRICalloc(%d): could not allocate struct",
              nrounds) ;

  mric->nrounds = nrounds ;
  for (round = 0 ; round < nrounds ; round++)
  {
    for (ninputs = 0, f = 0x001 ; f <= MAX_FEATURE ; f<<= 1)
      if (f & features[round])
        ninputs++ ;
    
    if (ninputs < 1 || ninputs > MAX_INPUTS)
      ErrorReturn(NULL, 
                  (ERROR_BADPARM, "MRICalloc(%d): bad # of inputs %d", round, 
                         ninputs)) ;
    
    mric->type[round] = types[round] ;
    mric->ninputs[round] = ninputs ;
    mric->features[round] = features[round] ;
    switch (types[round])
    {
    case CLASSIFIER_GAUSSIAN:
      mric->classifier[round].gc = 
        GCalloc(GAUSSIAN_NCLASSES, ninputs, gaussian_class_names);
      if (!mric->classifier[round].gc)
      {
        free(mric) ;
        return(NULL) ;
      }
      break ;
    default:
      MRICfree(&mric) ;
      ErrorReturn(NULL, 
                  (ERROR_UNSUPPORTED, "MRICalloc: classifier %d not supported",
                   types[round])) ;
      break ;
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
MRICtrain(MRIC *mric, char *file_name, char *prior_fname)
{
  char       source_fname[100], target_fname[100], line[300], *cp ;
  FILE       *fp ;
  int        fno, nfiles, round ;
  MRI        *mri_src, *mri_target ;
  MRI_REGION bounding_box ;

  if (prior_fname)
  {
    FileNameAbsolute(prior_fname, mric->prior_fname) ;
    mric->mri_priors = MRIread(prior_fname) ;
    if (!mric->mri_priors)
      ErrorReturn(ERROR_NO_FILE, 
                  (ERROR_NO_FILE, "MRICtrain: could not load prior file '%s'",
                  prior_fname)) ;
  }
  else
    mric->mri_priors = NULL ;

  for (round = 0 ; round < mric->nrounds ; round++)
  {
    /* first figure out the total # of files */
    fp = fopen(file_name, "r") ;
    if (!fp)
      ErrorReturn(ERROR_NO_FILE, 
                  (ERROR_NO_FILE, "MRICtrain(%s): could not open file",
                   file_name)) ;
    
    nfiles = 0 ;
    while ((cp = fgetl(line, 299, fp)) != NULL)
      nfiles++ ;
    fprintf(stderr, "processing %d files\n", nfiles) ;
    rewind(fp) ;
    
    /* now calculate statistics */
    fprintf(stderr, "computing classifier statistics...\n") ;
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
      
      MRIboundingBox(mri_src, DEFINITELY_BACKGROUND, &bounding_box) ;
      REGIONexpand(&bounding_box, &bounding_box, 5) ;
      MRIclipRegion(mri_src, &bounding_box, &bounding_box) ;
      REGIONcopy(&bounding_box, &mri_src->roi) ;
      MRICupdateStatistics(mric, round, mri_src, mri_target, &bounding_box) ;
      
      MRIfree(&mri_src) ;
      MRIfree(&mri_target) ;
      fno++ ;
    }
    
    MRICcomputeStatistics(mric, round) ;  /* divide by # of observations */
    fclose(fp) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *
MRICread(char *fname)
{
  MRIC  *mric ;
  int   ninputs, type[MAX_ROUNDS], features[MAX_ROUNDS], round, nrounds ;
  FILE  *fp ;
  char  prior_fname[100], line[100], *cp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NO_FILE,"MRICread(%s): could not open file", fname));

  prior_fname[0] = 0 ;
  cp = fgetl(line, 99, fp) ;
  if (sscanf(cp, "%d %s", &nrounds, prior_fname) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms",
                fname)) ;
  }
  for (round = 0 ; round < nrounds ; round++)
  {
    cp = fgetl(line, 99, fp) ;
    if (sscanf(cp, "%d %d 0x%x",&type[round], &ninputs,&features[round]) < 3)
    {
      fclose(fp) ;
      ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms",
                         fname)) ;
    }
  }

  mric = MRICalloc(nrounds, type, features, NULL) ;
  if (*prior_fname)
    mric->mri_priors = MRIread(prior_fname) ;
  for (round = 0 ; round < nrounds ; round++)
  {
    switch (type[round])
    {
    case CLASSIFIER_GAUSSIAN:
      GCasciiReadFrom(fp, mric->classifier[round].gc) ;
      break ;
    default:
      break ;
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
MRICwrite(MRIC *mric, char *fname)
{
  FILE  *fp ;
  int   round ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
              (ERROR_NO_FILE,"MRICwrite(%s): could not open file",fname));

  fprintf(fp, "%d ", mric->nrounds) ;
  if (mric->mri_priors)
    fprintf(fp, " %s", mric->prior_fname) ;
  fprintf(fp, "\n") ;

  for (round = 0 ; round < mric->nrounds ; round++)
    fprintf(fp, "%d %d 0x%x\n", mric->type[round], mric->ninputs[round], 
            mric->features[round]) ;

  for (round = 0 ; round < mric->nrounds ; round++)
  {
    switch (mric->type[round])
    {
    case CLASSIFIER_GAUSSIAN:
      GCasciiWriteInto(fp, mric->classifier[round].gc) ;
      break ;
    default:
      break ;
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
#define PRETTY_SURE              .90f
#define CLASS_SCALE              ((BUFTYPE)1)

MRI *
MRICclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, 
            float conf, MRI *mri_probs, MRI *mri_classes)
{
  MATRIX     *m_inputs, *m_priors ;
  GCLASSIFY  *gc ;
  int        x, y, z, width, depth, height, classno, nclasses, row, xt, yt,zt,
             round, x1, y1, z1, x0 ;
  BUFTYPE    *psrc, src, *pdst, *pclasses ;
  float      prob, *pprobs = NULL, inputs[MAX_INPUTS+1] ;
  Real       xrt, yrt, zrt ;
  MRI        *mri_priors, *mri_in ;
  MRI_REGION bounding_box ;

  MRIboundingBox(mri_src, DEFINITELY_BACKGROUND, &bounding_box) ;
  REGIONexpand(&bounding_box, &bounding_box, 5) ;
  MRIclipRegion(mri_src, &bounding_box, &bounding_box) ;
  REGIONcopy(&bounding_box, &mri_src->roi) ;
  x0 = bounding_box.x ;
  x1 = bounding_box.x + bounding_box.dx - 1 ;
  y1 = bounding_box.y + bounding_box.dy - 1 ;
  z1 = bounding_box.z + bounding_box.dz - 1 ;
  if (conf < 0.0f || conf >= 1.0f)
    conf = PRETTY_SURE ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  mri_priors = mric->mri_priors ;
  if (mri_priors)
    m_priors = MatrixAlloc(mric->classifier[0].gc->nclasses, 1, MATRIX_REAL) ;
  else
    m_priors = NULL ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  mri_in = mri_src ;
  for (round = 0 ; round < mric->nrounds ; round++)
  {
    gc = mric->classifier[round].gc ;
    nclasses = gc->nclasses ;
    m_inputs = MatrixAlloc(mric->ninputs[round], 1, MATRIX_REAL) ;
    
    for (z = bounding_box.z ; z <= z1 ; z++)
    {
      DiagHeartbeat((float)((z-bounding_box.z)+round*bounding_box.dz) / 
                    (float)(bounding_box.dz*mric->nrounds)) ;
      for (y = bounding_box.y ; y <= y1 ; y++)
      {
        psrc = &MRIvox(mri_src, x0, y, z) ;
        pdst = &MRIvox(mri_dst, x0, y, z) ;
        if (mri_probs)
          pprobs = &MRIFvox(mri_probs, x0, y, z) ;
        else
          pprobs = NULL ;
        if (mri_classes)
          pclasses = &MRIvox(mri_classes, x0, y, z) ;
        else
          pclasses = NULL ;
        for (x = x0 ; x <= x1 ; x++)
        {
          src = *psrc++ ;
          if (mri_priors)
          {
            MRIvoxelToVoxel(mri_src, mri_priors,
                            (Real)x, (Real)y, (Real)z,&xrt, &yrt,&zrt);
            xt = mri_priors->xi[nint(xrt)] ;
            yt = mri_priors->yi[nint(yrt)] ;
            zt = mri_priors->zi[nint(zrt)] ;
            for (classno = 0 ; classno < nclasses ; classno++)
              m_priors->rptr[classno+1][1] = 
                MRIFseq_vox(mri_priors, xt, yt, zt, classno) ;
          }
          MRICcomputeInputs(mri_in, x, y, z, inputs, mric->features[round]) ;
          for (row = 1 ; row <= mric->ninputs[round] ; row++)
            m_inputs->rptr[row][1] = inputs[row] ;
          
          /* now classify this observation */
          classno = GCclassify(gc, m_inputs, m_priors, &prob) ;
          
          if (pclasses)
            *pclasses++ = (BUFTYPE)classno*CLASS_SCALE ;
          if (pprobs)
            *pprobs++ = prob ;
          if (classno == WHITE_MATTER && prob > conf)
            *pdst++ = src ;
          else
            *pdst++ = 0 /*src*/ ;
        }
      }
    }
    MatrixFree(&m_inputs) ;
    mri_in = mri_classes ;  /* all subsequent rounds */
  }

  fprintf(stderr, "total = %d, buffered = %d, computed = %d\n",
          total, buffered, total_computed);
  fprintf(stderr, "efficiency = %2.3f%%\n", 
          100.0f*(float)(mri_src->width*mri_src->height*mri_src->depth) / 
          (float)total_computed) ;
  if (m_priors)
    MatrixFree(&m_priors) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define REGION_SIZE 8
int
MRICcomputeInputs(MRI *mri, int x,int y,int z,float *inputs,int features)
{
  static int         old_features = 0 ;
  static MRI_REGION  region = {0,0, 0, 0,0,0} ;
  static MRI  *mri_prev = NULL, *mri_zscore3 = NULL, *mri_zscore5 = NULL ,
              *mri_direction = NULL, *mri_mean3 = NULL, *mri_mean5 = NULL,
              *mri_cpolv = NULL, *mri_cpolv_mean3 = NULL, 
              *mri_cpolv_mean5 = NULL, *mri_cpolv_median3 = NULL,
              *mri_cpolv_median5 = NULL ;
  int         x0, y0, z0 ;
  MRI_REGION  rbig ;
  float       *in ;
  char        *cp ;

  if (mri != mri_prev)   /* reset counters */
    total = buffered = total_computed = 0 ;

  if ((features & FEATURE_CPOLV) && (!mri_cpolv || !MRImatch(mri,mri_cpolv)
                                     || (mri != mri_prev)))
  {
    if (mri_cpolv)
      MRIfree(&mri_cpolv) ;
    if (Gdiag & DIAG_HEARTBEAT)
      fprintf(stderr, "computing cpolv...\n") ;
    cp = getenv("cpolv") ;
    if (cp && !mri_cpolv)
    {
      mri_cpolv = MRIread(cp) ;
      if (!MRImatch(mri, mri_cpolv))
        MRIfree(&mri_cpolv) ;
    }

    if (!mri_cpolv)
      mri_cpolv = MRIplaneOfLeastVarianceNormal(mri, NULL, 5) ;
    if (Gdiag & DIAG_WRITE)
      MRIwrite(mri_cpolv, "cpolv.mnc") ;
  }
  total++ ;
/* 
   if the specified point is outside of the precomputed window,
   Update the window and compute a new set of input images.
   */
  if (!mri_prev || mri_prev != mri || (old_features != features) ||
      (REGIONinside(&region,x,y,z) == REGION_OUTSIDE))
  {
    MRI *mri_std, *mri_mean, *mri_region, *mri_grad ;

    old_features = features ;
    region.x = x ;
    region.y = y ;
    region.z = z ;
    region.dx = mri->width ;
    region.dy = mri->height ;
    region.dz = REGION_SIZE ;
    MRIclipRegion(mri, &region, &region) ;
    total_computed += (region.dx*region.dy*region.dz) ;
    mri_prev = mri ;
    if (mri_zscore3)
      MRIfree(&mri_zscore3) ;
    if (mri_zscore5)
      MRIfree(&mri_zscore5) ;
    if (mri_mean3)
      MRIfree(&mri_mean3) ;
    if (mri_mean5)
      MRIfree(&mri_mean5) ;
    if (mri_direction)
      MRIfree(&mri_direction) ;
    if (mri_cpolv_mean3)
      MRIfree(&mri_cpolv_mean3) ;
    if (mri_cpolv_mean5)
      MRIfree(&mri_cpolv_mean5) ;
    if (mri_cpolv_median3)
      MRIfree(&mri_cpolv_median3) ;
    if (mri_cpolv_median5)
      MRIfree(&mri_cpolv_median5) ;

    if (features & FEATURE_ZSCORE3)
    {
      static int first = 1 ;
      mri_mean = MRImeanRegion(mri, NULL, 3, &region) ;
      mri_std = MRIstdRegion(mri, NULL, mri_mean, 3, &region) ;
      mri_zscore3 = MRIzScoreRegion(mri, NULL, mri_mean, mri_std, &region);
      if (Gdiag & DIAG_WRITE && first)
      {
        first = 0 ;
        MRIwrite(mri_mean, "mean.mnc") ;
        MRIwrite(mri_std, "std.mnc") ;
        MRIwrite(mri_zscore3, "zscore.mnc") ;
      }
      MRIfree(&mri_mean) ;
      MRIfree(&mri_std) ;
    }
    if (features & FEATURE_ZSCORE5)
    {
      static int first = 1 ;
      mri_mean = MRImeanRegion(mri, NULL, 5, &region) ;
      mri_std = MRIstdRegion(mri, NULL, mri_mean, 5, &region) ;
      mri_zscore5 = MRIzScoreRegion(mri, NULL, mri_mean, mri_std, &region);
      if (Gdiag & DIAG_WRITE && first)
      {
        first = 0 ;
        MRIwrite(mri_mean, "mean.mnc") ;
        MRIwrite(mri_std, "std.mnc") ;
        MRIwrite(mri_zscore5, "zscore.mnc") ;
      }
      MRIfree(&mri_mean) ;
      MRIfree(&mri_std) ;
    }
    if (features & FEATURE_DIRECTION)
    {
      static int first = 1 ;
      MRI        *mri_tmp ;
      int        x0, y0, z0 ;

      REGIONexpand(&region, &rbig, 1) ; /* expand region by 1 voxel */
      MRIclipRegion(mri, &rbig, &rbig) ;
      mri_region = MRIextractRegion(mri, NULL, &rbig) ;
      mri_grad = MRIsobel(mri_region, NULL, NULL) ;
      mri_tmp = MRIdirectionMap(mri_grad, NULL, 3) ;
      x0 = region.x - rbig.x ;
      y0 = region.y - rbig.y ;
      z0 = region.z - rbig.z ;

      mri_direction = 
        MRIextract(mri_tmp, NULL, x0,y0,z0,region.dx,region.dy, region.dz) ;
      if (Gdiag & DIAG_WRITE && first)
      {
        first = 0 ;
        MRIwrite(mri_direction, "dir.mnc") ;
        MRIwrite(mri_region, "region.mnc") ;
        MRIwrite(mri_grad, "grad.mnc") ;
        MRIwrite(mri_tmp, "tmp.mnc") ;
      }
      MRIfree(&mri_region) ;
      MRIfree(&mri_grad) ;
      MRIfree(&mri_tmp) ;
    }
    if (features & FEATURE_MEAN3)
      mri_mean3 = MRImeanRegion(mri, NULL, 3, &region) ;
    if (features & FEATURE_MEAN5)
      mri_mean5 = MRImeanRegion(mri, NULL, 5, &region) ;
    if (features & FEATURE_CPOLV_MEAN3)
      mri_cpolv_mean3 = MRIpolvMeanRegion(mri, NULL, mri_cpolv, 3, &region) ;
    if (features & FEATURE_CPOLV_MEAN5)
      mri_cpolv_mean5 = MRIpolvMeanRegion(mri, NULL, mri_cpolv, 5, &region) ;
    if (features & FEATURE_CPOLV_MEDIAN3)
      mri_cpolv_median3 = MRIpolvMedianRegion(mri, NULL, mri_cpolv,3,&region);
    if (features & FEATURE_CPOLV_MEDIAN5)
    {
      static int first = 1 ;

      mri_cpolv_median5 = MRIpolvMedianRegion(mri, NULL, mri_cpolv,5,&region);
      if ((Gdiag & DIAG_WRITE) && first)
      {
        first = 0 ;
        MRIwrite(mri_cpolv_median5, "median5.mnc") ;
      }
    }
  }
  else
    buffered++ ;

  /* x0,y0,z0 are coordinates in region based images (not input mri or polv) */
  x0 = x-region.x ; y0 = y - region.y ; z0 = z - region.z ;
  in = &inputs[1] ;  /* inputs are 1-based because of NRC */
  if (features & FEATURE_INTENSITY)
    *in++ = (float)MRIvox(mri, x, y, z) ;
  if (features & FEATURE_ZSCORE3)
    *in++ = MRIFvox(mri_zscore3, x0, y0, z0) ;
  if (features & FEATURE_ZSCORE5)
    *in++ = MRIFvox(mri_zscore5, x0, y0, z0) ;
  if (features & FEATURE_DIRECTION)
    *in++ = MRIFvox(mri_direction, x0, y0, z0) ;
  if (features & FEATURE_MEAN3)
    *in++ = MRIFvox(mri_mean3, x0, y0, z0) ;
  if (features & FEATURE_MEAN5)
    *in++ = MRIFvox(mri_mean5, x0, y0, z0) ;
  if (features & FEATURE_CPOLV_MEAN3)
    *in++ = (float)MRIvox(mri_cpolv_mean3, x0, y0, z0) ;
  if (features & FEATURE_CPOLV_MEAN5)
    *in++ = (float)MRIvox(mri_cpolv_mean5, x0, y0, z0) ;
  if (features & FEATURE_CPOLV_MEDIAN3)
    *in++ = (float)MRIvox(mri_cpolv_median3, x0, y0, z0) ;
  if (features & FEATURE_CPOLV_MEDIAN5)
    *in++ = (float)MRIvox(mri_cpolv_median5, x0, y0, z0) ;
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Build an image of target values. Each point in the image
          contains the class # based on intensity thresholds either
          provided by the user, or defaults if parms are 0.
------------------------------------------------------*/
MRI *
MRICbuildTargetImage(MRI *mri_src, MRI *mri_target, MRI *mri_wm,
                     int lo_lim, int hi_lim)
{
  int     x, y, z, width, height, depth ;
  BUFTYPE *psrc, *pwm, *ptarget, src, wm, target ;

  if (lo_lim <= 0)
    lo_lim = LO_LIM ;
  if (hi_lim <= 0)
    hi_lim = HI_LIM ;

  if (!mri_target)
    mri_target = MRIclone(mri_wm, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      ptarget = &MRIvox(mri_target, 0, y, z) ;
      pwm = &MRIvox(mri_wm, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        src = *psrc++ ;
        wm = *pwm++ ;
        if (wm)
          target = WHITE_MATTER ;
        else if (src > hi_lim)
          target = BRIGHT_MATTER ;
        else if (src < lo_lim)
          target = BACKGROUND ;
        else
          target = GREY_MATTER ;
        *ptarget++ = target ;
      }
    }
  }
  return(mri_target) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Update the prior estimates using a frequency count.
------------------------------------------------------*/
#define COUNT_IMAGE   GAUSSIAN_NCLASSES

MRI *
MRICupdatePriors(MRI *mri_target, MRI *mri_priors, int scale)
{
  int      width, height, depth, x, y, z, class,w,h,d, xt, yt, zt ;
  BUFTYPE  *ptarget ;
  Real     xrt, yrt, zrt ;

  width = mri_target->width ;
  height = mri_target->height ;
  depth = mri_target->depth ;
  w = (int)(((float)width / (float)scale)+0.99f) ;
  h = (int)(((float)height / (float)scale)+0.99f) ;
  d = (int)(((float)depth / (float)scale)+0.99f) ;

/* 
   allocate one image for each class, plus one to keep track of the
   # of pixels mapped to that location.
   */
  if (!mri_priors)
  {
    mri_priors = 
      MRIallocSequence(w, h, d, MRI_FLOAT, GAUSSIAN_NCLASSES+1) ;
    MRIcopyHeader(mri_target, mri_priors) ;
    mri_priors->xsize *= scale ;
    mri_priors->ysize *= scale ;
    mri_priors->zsize *= scale ;
    mri_priors->linear_transform = mri_priors->inverse_linear_transform = NULL;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      ptarget = &MRIvox(mri_target, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        MRIvoxelToTalairachVoxel(mri_target, (Real)x, (Real)y, (Real)z,
                                 &xrt, &yrt, &zrt) ;
        xt = mri_priors->xi[nint(xrt/scale)] ;
        yt = mri_priors->yi[nint(yrt/scale)] ;
        zt = mri_priors->zi[nint(zrt/scale)] ;
        class = *ptarget++ ;
        MRIFseq_vox(mri_priors, xt, yt, zt, class) =
          MRIFseq_vox(mri_priors, xt, yt, zt, class) + 1.0f ;
        MRIFseq_vox(mri_priors, xt, yt, zt, COUNT_IMAGE) =
          MRIFseq_vox(mri_priors, xt, yt, zt, COUNT_IMAGE) + 1.0f ;
      }
    }
  }
  return(mri_priors) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Normalize the priors by dividing by the total # of
           times each pixel has been mapped.
------------------------------------------------------*/
int
MRInormalizePriors(MRI *mri_priors)
{
  int      width, height, depth, x, y, z, class ;
  float    *pnorm, norm ;

  width = mri_priors->width ;
  height = mri_priors->height ;
  depth = mri_priors->depth ;

/* 
   allocate one image for each class, plus one to keep track of the
   # of pixels mapped to that location.
   */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pnorm = &MRIFseq_vox(mri_priors, 0, y, z, COUNT_IMAGE) ;
      for (x = 0 ; x < width ; x++)
      {
        norm = *pnorm++ ;
        if (!FZERO(norm)) for (class = 0 ; class < GAUSSIAN_NCLASSES ; class++)
        {
          MRIFseq_vox(mri_priors, x, y, z, class) /= norm ;
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

------------------------------------------------------*/
int
MRICupdateStatistics(MRIC *mric, int round, MRI *mri_src, MRI *mri_target, 
                     MRI_REGION *box)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, classno, nclasses, width, height, depth, row, col,
             x1,y1,z1 ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      inputs[MAX_INPUTS+1], covariance ;

  nclasses = mric->classifier[round].gc->nclasses ;
  gc = mric->classifier[round].gc ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  x1 = box->x + box->dx - 1 ;
  y1 = box->y + box->dy - 1 ;
  z1 = box->z + box->dz - 1 ;
  for (z = box->z ; z <= z1 ; z++)
  {
    for (y = box->y ; y <= y1 ; y++)
    {
      psrc = &MRIvox(mri_src, box->x, y, z) ;
      ptarget = &MRIvox(mri_target, box->x, y, z) ;
      for (x = box->x ; x <= x1 ; x++)
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

        gcl = &gc->classes[classno] ;
        gcl->nobs++ ;

        MRICcomputeInputs(mri_src, x, y, z, inputs, mric->features[round]) ;
        for (row = 1 ; row <= mric->ninputs[round] ; row++)
          gcl->m_u->rptr[row][1] += inputs[row] ;
        for (row = 1 ; row <= gcl->m_covariance->rows ; row++)
        {
          for (col = 1 ; col <= row ; col++)
          {
            covariance = gcl->m_covariance->rptr[row][col] +
              inputs[row]*inputs[col] ;
            gcl->m_covariance->rptr[row][col] = covariance;
            gcl->m_covariance->rptr[col][row] = covariance;
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

------------------------------------------------------*/
int
MRICcomputeStatistics(MRIC *mric, int round)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        classno, nclasses, nobs, row, col ;
  float      mean_a, mean_b, covariance ;

  gc = mric->classifier[round].gc ;
  nclasses = gc->nclasses ;
  for (classno = 0 ; classno < nclasses ; classno++)
  {
    gcl = &gc->classes[classno] ;
    nobs = gcl->nobs ;
    if (nobs)
    {
      for (row = 1 ; row <= gcl->m_u->rows ; row++)
        gcl->m_u->rptr[row][1] /= (float)nobs ;
      for (row = 1 ; row <= gcl->m_covariance->rows ; row++)
      {
        mean_a = gcl->m_u->rptr[row][1] ;
        for (col = 1 ; col <= row ; col++)
        {
          mean_b = gcl->m_u->rptr[col][1] ;
          covariance = gcl->m_covariance->rptr[row][col] / nobs - mean_a*mean_b;
          gcl->m_covariance->rptr[row][col] = covariance ;
          gcl->m_covariance->rptr[col][row] = covariance ;
        }
      }
    }
    GCinit(gc, classno) ;
  }
  return(NO_ERROR) ;
}
char *
MRICclassName(MRIC *mric, int round, int classno)
{
  char *class_name = "unknown" ;

  switch (mric->type[round])
  {
  default:
  case CLASSIFIER_GAUSSIAN:
    class_name = gaussian_class_names[classno] ;
    break ;
  }

  return(class_name) ;
}

