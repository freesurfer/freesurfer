/*
 *       FILE NAME:   mrisegment.c
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
#include "mrisegment.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/


#define VOX_INCREASE  1.1 // 1.25

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int mriSegmentReallocateVoxels(MRI_SEGMENTATION *mriseg, int sno,
                                      int max_voxels) ;
static int mriSegmentReallocateSegments(MRI_SEGMENTATION *mriseg, 
                                        int max_segments) ;
static int mriSegmentNew(MRI_SEGMENTATION *mriseg) ;
static int mriSegmentMerge(MRI_SEGMENTATION *mriseg, int s0, int s1,
                           MRI *mri_labeled) ;
static int mriComputeSegmentStatistics(MRI_SEGMENTATION *mriseg) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_SEGMENTS     8000
#define MAX_VOXELS         5
#define NBR_VOX          (3*3*3)

MRI_SEGMENTATION *
MRIsegment(MRI *mri, float low_val, float hi_val)
{
  MRI_SEGMENTATION  *mriseg ;
  MRI_SEGMENT       *mseg, *mseg2 ;
  int               x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk,
                    val, border_labels[NBR_VOX], nlabels, label, nvox ;
  MRI               *mri_labeled ;
  float             voxel_size ;

  voxel_size = mri->xsize * mri->ysize * mri->zsize ;
  width = mri->width ; height = mri->height ; depth = mri->depth ;
  mriseg = MRIsegmentAlloc(MAX_SEGMENTS, MAX_VOXELS) ;
  mriseg->mri = mri ;

  /* mri_labeled will contain the label number for each voxel (if it
     has been assigned to a segment).
  */
  mri_labeled = MRIalloc(width, height, depth, MRI_SHORT) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
        MRISvox(mri_labeled, x, y, z) = -1 ;
    }
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri, x, y, z) ;
        if (val >= low_val && val <= hi_val)
        {
          memset(border_labels, -1, NBR_VOX*sizeof(int)) ;
          for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
          {
            zi = z+zk ; 
            if ((zi < 0) || (zi >= depth))
              continue ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = y+yk ; 
              if ((yi < 0) || (yi >= height))
                continue ;
              for (xk = -1 ; xk <= 1 ; xk++, nvox++)
              {
#if 1
                if ((abs(xk) + abs(yk) + abs(zk)) > 1)
                  continue ;  /* only allow 4 (6 in 3-d) connectivity */
#endif

                xi = x+xk ; 
                if ((xi < 0) || (xi >= width))
                  continue ;
                label = MRISvox(mri_labeled, xi, yi, zi) ;
                if (label >= 0)
                  border_labels[nvox] = label ;
              }
            }
          }
          for (nlabels = nvox = 0 ; nvox < NBR_VOX ; nvox++)
          {
            label = border_labels[nvox] ;
            if ((label >= 0) && (!mriseg->segments[label].found))
            {
              mriseg->segments[label].found = 1 ;
              nlabels++ ;
            }
          }
          for (nvox = 0 ; nvox < NBR_VOX ; nvox++)
          {
            label = border_labels[nvox] ;
            if (label >= 0)
              mriseg->segments[label].found = 0 ; /* for next time */
          }
          label = 0 ;
          switch (nlabels)
          {
          case 0:          /* allocate a new segment */
            label = mriSegmentNew(mriseg) ;
            mseg = &mriseg->segments[label] ;
            if (DIAG_VERBOSE_ON)
              fprintf(stderr, "allocating new label %d (%d total)\n",
                      label, mriseg->nsegments) ;
            break ;
          case 1:          /* assign this voxel to the one that it borders */
            for (nvox = 0 ; nvox < NBR_VOX ; nvox++)
              if (border_labels[nvox] >= 0)
              {
                label = border_labels[nvox] ;
                break ;
              }
            mseg = &mriseg->segments[label] ;
            break ;
          default:         /* merge segments and assign to lowest number */
            mseg = NULL ;
            for (nvox = 0 ; nvox < NBR_VOX ; nvox++)
            {
              if (border_labels[nvox] >= 0)
              {
                if (!mseg)
                {
                  label = border_labels[nvox] ;
                  mseg = &mriseg->segments[label] ;
                  mseg->found = 1 ;
                }
                else
                {
                  mseg2 = &mriseg->segments[border_labels[nvox]] ;
                  if (mseg2->found == 0)
                  {
                    mseg2->found = 1 ;  /* prevent merging more than once */
                    if (mriSegmentMerge(mriseg, label, border_labels[nvox],
                                        mri_labeled)  != NO_ERROR)
                    {
                      MRIsegmentFree(&mriseg) ;
                      return(NULL) ;
                    }
                  }
                }
                  
              }
            }
            for (nvox = 0 ; nvox < NBR_VOX ; nvox++)
              if (border_labels[nvox] >= 0)
                mriseg->segments[border_labels[nvox]].found = 0 ;
            break ;
          }
          /* add it to the existing list */
          if (mseg->nvoxels >= mseg->max_voxels)
          {
            if (mriSegmentReallocateVoxels(mriseg, label, 
                                           nint(mseg->max_voxels*VOX_INCREASE))
                != NO_ERROR)
            {
              MRIsegmentFree(&mriseg) ;
              return(NULL) ;
            }
          }
          mseg->voxels[mseg->nvoxels].x = x ;
          mseg->voxels[mseg->nvoxels].y = y ;
          mseg->voxels[mseg->nvoxels].z = z ;
          mseg->nvoxels++ ;
          mseg->area += voxel_size ;
          if (mseg->nvoxels != (int)mseg->area)
            DiagBreak() ;
          MRISvox(mri_labeled, x, y, z) = label ;
        }
      }
    }
  }
  MRIcompactSegments(mriseg) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_labeled, "labeled.mgh") ;
  MRIfree(&mri_labeled) ;
  mriComputeSegmentStatistics(mriseg) ;
  return(mriseg) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsegmentFree(MRI_SEGMENTATION **pmriseg)
{
  MRI_SEGMENTATION *mriseg ;
  int              i ;

  mriseg = *pmriseg ; *pmriseg = NULL ;
  for (i = 0 ; i < mriseg->max_segments ; i++)
    if (mriseg->segments[i].voxels)
      free(mriseg->segments[i].voxels) ;

  free(mriseg->segments) ;
  free(mriseg) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_SEGMENTATION *
MRIsegmentAlloc(int max_segments, int max_voxels)
{
  MRI_SEGMENTATION *mriseg ;
  int              i ;

  mriseg = (MRI_SEGMENTATION *)calloc(1, sizeof(MRI_SEGMENTATION)) ;
  if (!mriseg)
    ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mriseg") ;
  mriseg->segments = (MRI_SEGMENT *)calloc(max_segments, sizeof(MRI_SEGMENT)) ;
  if (!mriseg->segments)
    ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mrisegments");
  mriseg->max_segments = max_segments ;

  for (i = 0 ; i < max_segments ; i++)
  {
    mriseg->segments[i].area = 0.0f ;
    mriseg->segments[i].nvoxels = 0 ;
    mriseg->segments[i].max_voxels = max_voxels ;
    mriseg->segments[i].voxels = (MSV *)calloc(max_voxels, sizeof(MSV)) ;
    if (!mriseg->segments[i].voxels)
      ErrorExit(ERROR_NOMEMORY, 
                "MRIsegmentAlloc: could not alloc %dth voxel list", i);
  }

  return(mriseg) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Merge segment 1 into segment 0.
------------------------------------------------------*/
static int
mriSegmentMerge(MRI_SEGMENTATION *mriseg, int s0, int s1, MRI *mri_labeled)
{
  MRI_SEGMENT  *mseg0, *mseg1 ;
  int          v, total_voxels, x, y, z ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "merging segments %d and %d\n", s0, s1) ;
  mseg0 = &mriseg->segments[s0] ; mseg1 = &mriseg->segments[s1] ; 
  total_voxels = mseg0->nvoxels+mseg1->nvoxels ;
  if (total_voxels >= mseg0->max_voxels)
  {
    if (mriSegmentReallocateVoxels(mriseg, s0, total_voxels+10) != NO_ERROR)
      return(Gerror) ;
  }

  for (v = mseg0->nvoxels ; v < total_voxels ; v++)
  {
    x = mseg1->voxels[v-mseg0->nvoxels].x ;
    y = mseg1->voxels[v-mseg0->nvoxels].y ;
    z = mseg1->voxels[v-mseg0->nvoxels].z ;
    mseg0->voxels[v].x = x ; mseg0->voxels[v].y = y ; mseg0->voxels[v].z = z ;
    if (mri_labeled)
      MRISvox(mri_labeled, x, y, z) = s0 ;
  }
  mseg0->nvoxels = total_voxels ;
  mseg0->area += mseg1->area ;
  if (mseg0->nvoxels != (int)mseg0->area)
    DiagBreak() ;
  mriseg->nsegments-- ;
  mseg1->nvoxels = 0 ;
  mseg1->area = 0 ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mriSegmentReallocateVoxels(MRI_SEGMENTATION *mriseg, int sno,
                                      int max_voxels)
{
  MRI_SEGMENT *mseg ;
  // MSV         *old_voxels ;
  MSV         *new_voxels;
  int         v ;

  if (max_voxels <= 0)
    max_voxels = MAX_VOXELS ;

  mseg = &mriseg->segments[sno] ;
  // old_voxels = mseg->voxels ;
  
  // fprintf(stderr, "ReallocVoxels:%8d -> %8d\n", mseg->max_voxels, max_voxels);
  // mseg->voxels = (MSV *)calloc(max_voxels, sizeof(MSV)) ;
  new_voxels = (MSV *)realloc(mseg->voxels, max_voxels*sizeof(MSV)) ;
  if (!new_voxels)
  // if (!mseg->voxels)
  {
    // if (old_voxels)
    //   free(old_voxels) ;
    // when fail, mseg->voxels not modified
     ErrorReturn(ERROR_NOMEMORY, 
                (ERROR_NOMEMORY,
                 "mriSegmentReallocateVoxels: could not alloc %d voxels for sno %d",
                 max_voxels, sno)) ;
  }
  // realloc
  mseg->voxels = new_voxels;
  mseg->max_voxels = max_voxels ;
  // realloc lack of initialization
  for (v=mseg->nvoxels; v < mseg->max_voxels; ++v)
  {
    mseg->voxels[v].x = 0;
    mseg->voxels[v].y = 0;
    mseg->voxels[v].z = 0;
  }
  // this is just copying.  realloc makes this unnecessary
  // for (v = 0 ; v < mseg->nvoxels ; v++)
  // {
  //   mseg->voxels[v].x = old_voxels[v].x ;
  //   mseg->voxels[v].y = old_voxels[v].y ;
  //   mseg->voxels[v].z = old_voxels[v].z ;
  // }
  // if (old_voxels)
  //   free(old_voxels) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mriSegmentReallocateSegments(MRI_SEGMENTATION *mriseg, int max_segments)
{
  MRI_SEGMENT *old_segments ;
  int         s ;

  old_segments = mriseg->segments ;
  mriseg->segments = (MRI_SEGMENT *)calloc(max_segments, sizeof(MRI_SEGMENT)) ;
  if (!mriseg->segments)
    ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mrisegments");
  mriseg->max_segments = max_segments ;

  for (s = 0 ; s < mriseg->nsegments ; s++)
  {
    if (s == 1)
      DiagBreak() ;
    mriseg->segments[s].area = old_segments[s].area ;
    mriseg->segments[s].nvoxels = old_segments[s].nvoxels ;
    mriseg->segments[s].max_voxels = old_segments[s].max_voxels ;
    mriseg->segments[s].voxels = old_segments[s].voxels ;
  }
  free(old_segments) ;
  for (s = mriseg->nsegments ; s < mriseg->max_segments ; s++)
  {
    mriseg->segments[s].voxels = (MSV *)calloc(MAX_VOXELS, sizeof(MSV)) ;
    if (!mriseg->segments[s].voxels)
      ErrorExit(ERROR_NOMEMORY, "mriSegmentRealloc: could not alloc voxels");
    mriseg->segments[s].max_voxels = MAX_VOXELS ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIsegmentToImage(MRI *mri_src, MRI *mri_dst, MRI_SEGMENTATION *mriseg, int s)
{
  int          v, x, y, z ;
  MRI_SEGMENT  *mseg ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  mseg = &mriseg->segments[s] ;

  for (v = 0 ; v < mseg->nvoxels ; v++)
  {
    x = mseg->voxels[v].x ; y = mseg->voxels[v].y ; z = mseg->voxels[v].z ; 
    MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z) ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mriSegmentNew(MRI_SEGMENTATION *mriseg)
{
  int s ;

  if (mriseg->nsegments >= mriseg->max_segments)
    mriSegmentReallocateSegments(mriseg, nint(mriseg->max_segments*1.5)) ;
  mriseg->nsegments++ ;
  for (s = 0 ; s < mriseg->nsegments ; s++)
    if (mriseg->segments[s].nvoxels == 0)
      return(s) ;
  return(s) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIcompactSegments(MRI_SEGMENTATION *mriseg)
{
  MRI_SEGMENT *src_segment, *dst_segment ;
  MRI_SEGMENT *newseg;
  int         s, s2 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "compacting segments...\n") ;

  for (s = 0 ; s < mriseg->max_segments ; s++)
  {
    // if voxel count = 0, then
    if (mriseg->segments[s].nvoxels == 0)
    {
      // assgin dst
      dst_segment = &mriseg->segments[s] ;
      // find the next non-zero voxel segments
      for (s2 = s+1 ; s2 < mriseg->max_segments ; s2++)
      {
        if (mriseg->segments[s2].nvoxels > 0)
          break ;
      }
      // if within the boundary
      if (s2 < mriseg->max_segments)
      {
	// assin src
        src_segment = &mriseg->segments[s2] ;
        if (dst_segment->voxels)
          free(dst_segment->voxels) ;
	// copy src to dst
        dst_segment->area = src_segment->area ;
        dst_segment->nvoxels = src_segment->nvoxels ;
        dst_segment->voxels = src_segment->voxels ;
        dst_segment->max_voxels = src_segment->max_voxels ;
        // make the src to have nvoxels = 0
        src_segment->nvoxels = src_segment->max_voxels = 0 ;
        src_segment->area = 0.0 ; src_segment->voxels = NULL ;
	// now we moved the non-zero segment to the front 
      }
    }
  }
  //////////////////////////////////
  // now we have segments zero voxels ones stack at the end
  // now we keep only the non-zero voxel segments
  // find the location of the first non-zero voxel segments
  for (s = 0; s < mriseg->max_segments; s++)
  {
    if (mriseg->segments[s].nvoxels == 0)
      break;
  }
  newseg = (MRI_SEGMENT *) realloc(mriseg->segments, s*sizeof(MRI_SEGMENT));
  if (newseg)
  {
    mriseg->segments = newseg;
    mriseg->max_segments = s;
    mriseg->nsegments = s;
    fprintf(stderr, "segments reduced to %d\n", s);
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsegmentDilate(MRI_SEGMENTATION *mriseg, MRI *mri)
{
  int         x, y, z, xk, yk, zk, xi, yi, zi, segno, v, nvox ;
  MRI_SEGMENT *mseg ;
  MRI         *mri_segments ;
  float       voxel_size ;

  voxel_size = mri->xsize * mri->ysize * mri->zsize ;

  mri_segments = MRIclone(mri, NULL) ;

  /* build image of all other segments to prevent dilation from
     expanding into other segments */
  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
    MRIsegmentToImage(mri, mri_segments, mriseg, segno) ;

  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
  {
    mseg = &mriseg->segments[segno] ;

    nvox = mseg->nvoxels ;
    for (v = 0 ; v < nvox ; v++)
    {
      x = mseg->voxels[v].x ; y = mseg->voxels[v].y ; z = mseg->voxels[v].z ;

      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        for (yk = -1 ; yk <= 1 ; yk++)
        {
          yi = mri->yi[y+yk] ;
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1)
              continue ;
            zi = mri->zi[z+zk] ;
            if (MRIvox(mri, xi, yi, zi) && !MRIvox(mri_segments,xi,yi,zi))
            {
              if (mseg->nvoxels >= mseg->max_voxels)
              {
                if (mriSegmentReallocateVoxels(mriseg, segno, 
                                               mseg->max_voxels*VOX_INCREASE)
                    != NO_ERROR)
                {
                  /*                  MRIsegmentFree(&mseg) ;*/
                  return(Gerror) ;
                }
              }
              mseg->voxels[mseg->nvoxels].x = xi ;
              mseg->voxels[mseg->nvoxels].y = yi ;
              mseg->voxels[mseg->nvoxels].z = zi ;
              MRIvox(mri_segments,xi,yi,zi) = MRIvox(mri, xi, yi, zi) ;
              mseg->nvoxels++ ;
              mseg->area += voxel_size ;
            }
          }
        }
      }
    }
  }
  MRIfree(&mri_segments) ;
  return(NO_ERROR) ;
}

int
MRIremoveSmallSegments(MRI_SEGMENTATION *mriseg, int min_voxels)
{
  int         s ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "compacting segments...\n") ;

  for (s = 0 ; s < mriseg->max_segments ; s++)
    if (mriseg->segments[s].nvoxels < min_voxels)
      mriseg->segments[s].nvoxels = 0 ;
  return(MRIcompactSegments(mriseg)) ;
}
static int
mriComputeSegmentStatistics(MRI_SEGMENTATION *mriseg)
{
  int         segno ;
  MRI_SEGMENT *mseg ;
  int         v, x, y, z ;

  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
  {
    mseg = &mriseg->segments[segno] ;
    mseg->x0 = mseg->y0 = mseg->z0 = 10000 ;
    mseg->x1 = mseg->y1 = mseg->z1 = -10000 ;
    mseg->cx = mseg->cy = mseg->cz = 0.0 ;
    if (!mseg->nvoxels)
      continue ;
    for (v = 0 ; v < mseg->nvoxels ; v++)
    {
      x = mseg->voxels[v].x ; y = mseg->voxels[v].y ; z = mseg->voxels[v].z ;
      mseg->cx += x ; mseg->cy += y ; mseg->cz += z ; 
      if (x < mseg->x0)
        mseg->x0 = x ;
      if (y < mseg->y0)
        mseg->y0 = y ;
      if (z < mseg->z0)
        mseg->z0 = z ;
      if (x > mseg->x1)
        mseg->x1 = x ;
      if (y > mseg->y1)
        mseg->y1 = y ;
      if (z > mseg->z1)
        mseg->z1 = z ;
    }
    mseg->cx /= (float)mseg->nvoxels ;
    mseg->cy /= (float)mseg->nvoxels ;
    mseg->cz /= (float)mseg->nvoxels ;
  }

  return(NO_ERROR) ;
}

int
MRIsegmentDilateThreshold(MRI_SEGMENTATION *mriseg, MRI *mri_binary, 
                          MRI *mri_thresh, int low_thresh, int hi_thresh)
{
  int         x, y, z, xk, yk, zk, xi, yi, zi, segno, v, nvox, val ;
  MRI_SEGMENT *mseg ;
  MRI         *mri_segments ;
  float       voxel_size ;

  voxel_size = mri_binary->xsize * mri_binary->ysize * mri_binary->zsize ;

  mri_segments = MRIclone(mri_binary, NULL) ;

  /* build image of all other segments to prevent dilation from
     expanding into other segments */
  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
    MRIsegmentToImage(mri_binary, mri_segments, mriseg, segno) ;

  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
  {
    mseg = &mriseg->segments[segno] ;

    nvox = mseg->nvoxels ;
    for (v = 0 ; v < nvox ; v++)
    {
      x = mseg->voxels[v].x ; y = mseg->voxels[v].y ; z = mseg->voxels[v].z ;

      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri_binary->xi[x+xk] ;
        for (yk = -1 ; yk <= 1 ; yk++)
        {
          yi = mri_binary->yi[y+yk] ;
          for (zk = -1 ; zk <= 1 ; zk++)
          {
            if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1)
              continue ;
            zi = mri_binary->zi[z+zk] ;
            val = MRIvox(mri_thresh, xi, yi, zi) ;
            if (MRIvox(mri_binary, xi, yi, zi) && 
                !MRIvox(mri_segments,xi,yi,zi) &&
                (val >= low_thresh) && (val <= hi_thresh))
            {
              if (mseg->nvoxels >= mseg->max_voxels)
              {
                if (mriSegmentReallocateVoxels(mriseg, segno, 
                                               nint(mseg->max_voxels
                                                    *VOX_INCREASE)) != 
                    NO_ERROR)
                {
                  /*                  MRIsegmentFree(&mseg) ;*/
                  return(Gerror) ;
                }
              }
              mseg->voxels[mseg->nvoxels].x = xi ;
              mseg->voxels[mseg->nvoxels].y = yi ;
              mseg->voxels[mseg->nvoxels].z = zi ;
              MRIvox(mri_segments,xi,yi,zi) = MRIvox(mri_binary, xi, yi, zi) ;
              mseg->nvoxels++ ;
              mseg->area += voxel_size ;
            }
          }
        }
      }
    }
  }
  MRIfree(&mri_segments) ;
  return(NO_ERROR) ;
}

int
MRIsegmentMax(MRI_SEGMENTATION *mriseg)
{
  int         segno, max_voxels, max_segno, nvox ;
  MRI_SEGMENT *mseg ;

  max_segno = -1 ; max_voxels = 0 ;
  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
  {
    mseg = &mriseg->segments[segno] ;
    if (mseg->ignore)
      continue ;

    nvox = mseg->nvoxels ;
    if (nvox > max_voxels)
    {
      max_segno = segno ;
      max_voxels = nvox ;
    }
  }
  return(max_segno) ;
}

int
MRIsegmentClearIgnoreFlags(MRI_SEGMENTATION *mriseg)
{
  int         segno ;
  MRI_SEGMENT *mseg ;

  for (segno = 0 ; segno < mriseg->nsegments ; segno++)
  {
    mseg = &mriseg->segments[segno] ;
    mseg->ignore = 0 ;
  }
  return(NO_ERROR) ;
}

