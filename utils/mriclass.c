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

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define GAUSSIAN_NCLASSES    4
#define BACKGROUND           0
#define GREY_MATTER          1
#define WHITE_MATTER         2
#define BRIGHT_MATTER        3
#define LO_LIM               70
#define HI_LIM               150

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

  mric = *pmric ;
  *pmric = NULL ;
  if (mric)
  {
    switch (mric->type)
    {
    case CLASSIFIER_GAUSSIAN:
      GCfree(&mric->gc) ;
      break ;
    default:
      break ;
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
MRICalloc(int type, int ninputs, void *parms)
{
  MRIC  *mric ;

  if (ninputs < 1 || ninputs > MAX_INPUTS)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRICalloc(%d, %d): bad # of inputs", type, ninputs)) ;

  mric = (MRIC *)calloc(1, sizeof(MRIC)) ;
  if (!mric)
    ErrorExit(ERROR_NO_MEMORY, "MRICalloc(%d, %d): could not allocate struct",
              type, ninputs) ;

  mric->type = type ;
  mric->ninputs = ninputs ;
  switch (type)
  {
  case CLASSIFIER_GAUSSIAN:
    mric->gc = GCalloc(GAUSSIAN_NCLASSES, ninputs, gaussian_class_names);
    if (!mric->gc)
    {
      free(mric) ;
      return(NULL) ;
    }
    break ;
  default:
    MRICfree(&mric) ;
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "MRICalloc: classifier %d not supported",
                 type)) ;
    break ;
  }
  return(mric) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRICtrain(MRIC *mric, char *file_name)
{
  char  source_fname[100], target_fname[100], line[300], *cp ;
  FILE  *fp ;
  int   fno, nfiles ;
  MRI   *mri_src, *mri_target ;
  BOX   bounding_box ;

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

    mri_target = MRIread(target_fname) ;
    if (!mri_target)
    {
      fprintf(stderr, "could not read MR image %s\n", target_fname) ;
      MRIfree(&mri_src) ;
      continue ;
    }

    MRIboundingBox(mri_target, 1, &bounding_box) ;
    BoxExpand(&bounding_box, &bounding_box, 10, 10, 10) ;
    MRICupdateMeans(mric, mri_src, mri_target, &bounding_box) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    fno++ ;
  }

  MRICcomputeMeans(mric) ;  /* divide by # of observations */
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

    MRIboundingBox(mri_target, 1, &bounding_box) ;
    BoxExpand(&bounding_box, &bounding_box, 10, 10, 10) ;
    MRICupdateCovariances(mric, mri_src, mri_target, &bounding_box) ;

    MRIfree(&mri_src) ;
    MRIfree(&mri_target) ;
    fno++ ;
  }

  MRICcomputeCovariances(mric) ;

  fclose(fp) ;
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRICclassify(MRIC *mric, MRI *mri_src, float *pprob)
{
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *
MRICread(char *fname)
{
  MRIC  *mric ;
  int   ninputs, type ;
  FILE  *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NO_FILE,"MRICread(%s): could not open file", fname));

  if (fscanf(fp, "%d %d\n", &type, &ninputs) != 2)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms",
                fname)) ;
  }

  mric = MRICalloc(type, ninputs, NULL) ;
  switch (type)
  {
  case CLASSIFIER_GAUSSIAN:
    GCasciiReadFrom(fp, mric->gc) ;
    break ;
  default:
    break ;
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

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
              (ERROR_NO_FILE,"MRICwrite(%s): could not open file",fname));

  fprintf(fp, "%d %d\n", mric->type, mric->ninputs) ;
  switch (mric->type)
  {
  case CLASSIFIER_GAUSSIAN:
    GCasciiWriteInto(fp, mric->gc) ;
    break ;
  default:
    break ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define PRETTY_SURE   .90f
MRI *
MRICclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, 
            float conf, MRI *mri_probs, MRI *mri_classes)
{
  MATRIX     *m_inputs ;
  GCLASSIFY  *gc ;
  int        x, y, z, width, depth, height, classno, nclasses, row ;
  BUFTYPE    *psrc, src, *pdst, *pclasses ;
  float      prob, *pprobs = NULL, inputs[MAX_INPUTS+1] ;

  if (conf < 0.0f || conf >= 1.0f)
    conf = PRETTY_SURE ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  m_inputs = MatrixAlloc(mric->ninputs, 1, MATRIX_REAL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  gc = mric->gc ;
  nclasses = gc->nclasses ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      if (mri_probs)
        pprobs = &MRIFvox(mri_probs, 0, y, z) ;
      else
        pprobs = NULL ;
      if (mri_classes)
        pclasses = &MRIvox(mri_classes, 0, y, z) ;
      else
        pclasses = NULL ;
      for (x = 0 ; x < width ; x++)
      {
        src = *psrc++ ;
        MRICcomputeInputs(mri_src, x, y, z, inputs, mric->ninputs) ;
        for (row = 1 ; row <= mric->ninputs ; row++)
          m_inputs->rptr[row][1] = inputs[row] ;
        
        /* now classify this observation */
        classno = GCclassify(gc, m_inputs, &prob) ;

        if (pclasses)
          *pclasses++ = (BUFTYPE)classno ;
        if (pprobs)
          *pprobs++ = prob ;
        if (classno == WHITE_MATTER && prob > conf)
          *pdst++ = 255 ;
        else
          *pdst++ = 0 /*src*/ ;
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
MRICupdateMeans(MRIC *mric, MRI *mri_src, MRI *mri_target, BOX *box)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, classno, nclasses, width, height, depth, row ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      inputs[MAX_INPUTS+1] ;

  nclasses = mric->gc->nclasses ;
  gc = mric->gc ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

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
  for (z = box->z0 ; z <= box->z1 ; z++)
  {
    for (y = box->y0 ; y <= box->y1 ; y++)
    {
      psrc = &MRIvox(mri_src, box->x0, y, z) ;
      ptarget = &MRIvox(mri_target, box->x0, y, z) ;
      for (x = box->x0 ; x <= box->x1 ; x++)
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

        MRICcomputeInputs(mri_src, x, y, z, inputs, mric->ninputs) ;
        for (row = 1 ; row <= mric->ninputs ; row++)
          gcl->m_u->rptr[row][1] += inputs[row] ;
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
MRICcomputeMeans(MRIC *mric)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        classno, nclasses, nobs, row ;

/*
   the classifiers are distributed in Talairach space, whereas the
   input MR images are not (necessarily). Therefore we have to
   transform the voxel coordinates into Tal. space before selecting
   the training values.
*/
  gc = mric->gc ;
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
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRICupdateCovariances(MRIC *mric, MRI *mri_src, MRI *mri_target, BOX *box)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        x, y, z, width, depth, height, classno, 
             nclasses, col, row ;
  BUFTYPE    *psrc, *ptarget, src, target ;
  float      inputs[MAX_INPUTS+1], covariance ;

  gc = mric->gc ;
  nclasses = gc->nclasses ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = box->z0 ; z <= box->z1 ; z++)
  {
    for (y = box->y0 ; y <= box->y1 ; y++)
    {
      psrc = &MRIvox(mri_src, box->x0, y, z) ;
      ptarget = &MRIvox(mri_target, box->x0, y, z) ;
      for (x = box->x0 ; x <= box->x1 ; x++)
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

        MRICcomputeInputs(mri_src, x, y, z, inputs, mric->ninputs) ;
        for (row = 1 ; row <= mric->ninputs ; row++)  /* subtract means */
          inputs[row] -= gcl->m_u->rptr[row][1] ;
        for (row = 1 ; row <= gcl->m_covariance->rows ; row++)
        {
          for (col = 1 ; col <= row ; col++)
          {
            covariance = inputs[row] * inputs[col] ;
            gcl->m_covariance->rptr[row][col] += covariance;
            gcl->m_covariance->rptr[col][row] += covariance;
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
MRICcomputeCovariances(MRIC *mric)
{
  GCLASSIFY  *gc ;
  GCLASS     *gcl ;
  int        classno, nclasses, row,col ;
  float      nobs, covariance ;

  gc = mric->gc ;
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
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int
MRICcomputeInputs(MRI *mri, int x,int y,int z,float *inputs,int ninputs)
{
  int   i ;

  i = 1 ;
  inputs[i++] = (float)MRIvox(mri, x, y, z) ;
  if (ninputs > 1)
    inputs[i++] = MRIvoxelZscore(mri, x, y, z, 3) ;
  if (ninputs > 2)
    inputs[i++] = MRIvoxelZscore(mri, x, y, z, 5) ;
  if (ninputs > 3)
    inputs[i++] = MRIvoxelDirection(mri, x, y, z, 3) ;

  return(NO_ERROR) ;
}
