/**
 * @brief calc and filter MRI data based on planes of least variance
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "mrisegment.h"
#include "proto.h"

/*-----------------------------------------------------
  MACROS AND CONSTANTS
  -------------------------------------------------------*/

#define VOX_INCREASE 1.1  // 1.25

/*-----------------------------------------------------
  STATIC PROTOTYPES
  -------------------------------------------------------*/

static int mriSegmentReallocateVoxels(MRI_SEGMENTATION *mriseg, int sno, int max_voxels);
static int mriSegmentReallocateSegments(MRI_SEGMENTATION *mriseg, int max_segments);
static int mriSegmentNew(MRI_SEGMENTATION *mriseg);
static int mriSegmentMerge(MRI_SEGMENTATION *mriseg, int s0, int s1, MRI *mri_labeled);
static int mriComputeSegmentStatistics(MRI_SEGMENTATION *mriseg);

/*-----------------------------------------------------
  GLOBAL FUNCTIONS
  -------------------------------------------------------*/
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_SEGMENTS 8000
#define MAX_VOXELS 5
#define NBR_VOX (3 * 3 * 3)

int MRImaxSegmentArea(MRI_SEGMENTATION *mriseg)
{
  int maxarea = -1;
  // int max = -1;
  int k;

  for (k = 0; k < mriseg->max_segments; k++)
    if (mriseg->segments[k].nvoxels > maxarea) {
      // max = k;
      maxarea = mriseg->segments[k].nvoxels;
    }
  return maxarea;
}

// the following only cares about finding the max segments
// ignoring the rest
MRI_SEGMENTATION *MRImaxsegment(MRI *mri, float low_val, float hi_val)
{
  MRI_SEGMENTATION *mriseg;
  MRI_SEGMENT *mseg, *mseg2;
  int x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk, val, border_labels[NBR_VOX], nlabels, label, nvox;
  MRI *mri_labeled;
  float voxel_size;
  int max;
  int count = 0;
  int totcount = 0;
  int mem = 0;
  int maxarea = 0;

  voxel_size = mri->xsize * mri->ysize * mri->zsize;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  mriseg = MRIsegmentAlloc(MAX_SEGMENTS, MAX_VOXELS);
  mriseg->mri = mri;

  for (z = 0; z < depth; z++)
    for (y = 0; y < height; y++)
      for (x = 0; x < width; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (val >= low_val && val <= hi_val) totcount++;
      }

  /* mri_labeled will contain the label number for each voxel (if it
     has been assigned to a segment).
  */
  mri_labeled = MRIalloc(width, height, depth, MRI_SHORT);
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) MRISvox(mri_labeled, x, y, z) = -1;
    }
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        // if within the range
        if (val >= low_val && val <= hi_val) {
          count++;

          // initializer for border_labels
          memset(border_labels, -1, NBR_VOX * sizeof(int));
          for (nvox = 0, zk = -1; zk <= 1; zk++) {
            zi = z + zk;
            if ((zi < 0) || (zi >= depth)) continue;
            for (yk = -1; yk <= 1; yk++) {
              yi = y + yk;
              if ((yi < 0) || (yi >= height)) continue;
              // increment nvox count here
              for (xk = -1; xk <= 1; xk++, nvox++) {
#if 1
                if ((abs(xk) + abs(yk) + abs(zk)) > 1) continue; /* only allow 4 (6 in 3-d) connectivity */
#endif

                xi = x + xk;
                if ((xi < 0) || (xi >= width)) continue;
                // get the neighbor label
                label = MRISvox(mri_labeled, xi, yi, zi);
                if (label >= 0) border_labels[nvox] = label;
              }
            }
          }
          // count nonzero labels in nbrs
          for (nlabels = nvox = 0; nvox < NBR_VOX; nvox++) {
            label = border_labels[nvox];
            if ((label >= 0) && (!mriseg->segments[label].found)) {
              mriseg->segments[label].found = 1;
              nlabels++;
            }
          }
          // reset found
          for (nvox = 0; nvox < NBR_VOX; nvox++) {
            label = border_labels[nvox];
            if (label >= 0) mriseg->segments[label].found = 0; /* for next time */
          }
          //
          label = 0;
          // create labels for those points which are not connected
          switch (nlabels) {
            case 0: /* allocate a new segment */
              if (mriseg->nsegments >= mriseg->max_segments) {
                mem = getMemoryUsed();
                // only Linux gives the correct value
                // (the rest returns -1
                // mem > 800*1024) // 800 Mbytes virtual memory usage
                if (mriseg->max_segments > MAX_SEGMENTS * VOX_INCREASE) {
                  if (mem > 0)  // only Linux can do this
                    fprintf(stdout,
                            "\n             "
                            "heap usage = %d Kbytes.",
                            mem);
                  // find the region with the largest area
                  maxarea = MRImaxSegmentArea(mriseg);
                  fprintf(stdout,
                          "\n             current max segment has "
                          "%d voxels",
                          maxarea);
                  // the second max area can have area up to
                  // (current - maxarea) + (total - maxarea)
                  //  = total + current - 2*maxarea
                  // if this value is less than maxarea,
                  // then the second candidate never
                  // becomes the top. Thus I can remove
                  // small segments
                  // if (count+totcount < 3*maxarea)
                  //
                  // For those which has < 100*N % of maxarea,
                  // the possible max count is
                  //    = maxarea*N + remaining voxels
                  //    = maxarea*N + (total - count)
                  // Thus I can remove those when
                  // maxarea*N + (total-count) < maxarea
                  // or total - count < maxarea*(1 - N)
                  //
                  // Note that if you remove too much,
                  // then possbile merging voxels will be lost
                  //
                  if (totcount - count < maxarea * 0.99) {
                    int i, j, k;
                    fprintf(stdout,
                            "\n             removing "
                            "small segments (less than 1 "
                            "percent of maxarea).");
                    MRIremoveSmallSegments(mriseg, maxarea * 0.01);
                    // this does compactify
                    // go through mri_labeled to
                    // remove this label value s
                    for (k = 0; k < mri_labeled->depth; ++k)
                      for (j = 0; j < mri_labeled->height; ++j)
                        for (i = 0; i < mri_labeled->width; ++i) {
                          if (MRISvox(mri_labeled, i, j, k) >= mriseg->nsegments) MRISvox(mri_labeled, i, j, k) = -1;
                        }
                  }
                }
              }
              label = mriSegmentNew(mriseg);
              // returns this added segment position
              mseg = &mriseg->segments[label];
              if (0 && DIAG_VERBOSE_ON)
                fprintf(stdout, "allocating new label %d (%d total)\n", label, mriseg->nsegments);
              break;
            case 1: /* assign this voxel to the one that it borders */
              for (nvox = 0; nvox < NBR_VOX; nvox++)
                if (border_labels[nvox] >= 0) {
                  label = border_labels[nvox];
                  break;
                }
              // points to the same label position
              mseg = &mriseg->segments[label];
              break;
            default: /* merge segments and assign to lowest number */
              mseg = NULL;
              for (nvox = 0; nvox < NBR_VOX; nvox++) {
                if (border_labels[nvox] >= 0) {
                  // the 1st encountered label case
                  if (!mseg) {
                    // set the first index position
                    label = border_labels[nvox];
                    mseg = &mriseg->segments[label];
                    mseg->found = 1;
                  }
                  // the rest
                  else {
                    mseg2 = &mriseg->segments[border_labels[nvox]];
                    if (mseg2->found == 0) {
                      mseg2->found = 1; /* prevent merging more than once */
                      // merge to the label position
                      if (mriSegmentMerge(mriseg, label, border_labels[nvox], mri_labeled) != NO_ERROR) {
                        // nsegments decreased by one
                        MRIsegmentFree(&mriseg);
                        return (NULL);
                      }
                    }
                  }
                }
              }
              // reset found
              for (nvox = 0; nvox < NBR_VOX; nvox++)
                if (border_labels[nvox] >= 0) mriseg->segments[border_labels[nvox]].found = 0;
              break;
          }
          /* add it to the existing list */
          if (mseg->nvoxels >= mseg->max_voxels) {
            // this max could be the same as mseg->max_voxels
            // this is taken care in mriSegmentReallocateVoxels()
            max = nint(mseg->max_voxels * VOX_INCREASE);
            // if (mriSegmentReallocateVoxels(mriseg, label,
            //            nint(mseg->max_voxels*VOX_INCREASE))
            //    != NO_ERROR)
            if (mriSegmentReallocateVoxels(mriseg, label, max) != NO_ERROR) {
              MRIsegmentFree(&mriseg);
              return (NULL);
            }
          }
          mseg->voxels[mseg->nvoxels].x = x;
          mseg->voxels[mseg->nvoxels].y = y;
          mseg->voxels[mseg->nvoxels].z = z;
          mseg->nvoxels++;
          mseg->area += voxel_size;  // voxel_size = volume
#if 0
          // this is for only 1mm voxel case
          if (mseg->nvoxels != (int)mseg->area)
            DiagBreak() ;
#endif
          MRISvox(mri_labeled, x, y, z) = label;
        }
      }
    }
  }
  mem = getMemoryUsed();
  if (mem > 0)  // only Linux can do this
    fprintf(stdout, "\n             heap usage = %d Kbytes.", mem);
  maxarea = MRImaxSegmentArea(mriseg);
  fprintf(stdout,
          "\n             removing small segments (less "
          "than 1 percent of maxarea).");
  MRIremoveSmallSegments(mriseg, 0.01 * maxarea);
  // MRIcompactSegments(mriseg) ;

  if (0 && Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_labeled, "labeled.mgh");
  MRIfree(&mri_labeled);
  mriComputeSegmentStatistics(mriseg);
  return (mriseg);
}

MRI_SEGMENTATION *MRIsegment(MRI *mri, float low_val, float hi_val)
{
  MRI_SEGMENTATION *mriseg;
  MRI_SEGMENT *mseg, *mseg2;
  int x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk, border_labels[NBR_VOX], nlabels, label, nvox;
  MRI *mri_labeled;
  float voxel_size, val;
  int max;
  int count = 0;

  voxel_size = mri->xsize * mri->ysize * mri->zsize;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  mriseg = MRIsegmentAlloc(MAX_SEGMENTS, MAX_VOXELS);
  mriseg->mri = mri;

  /* mri_labeled will contain the label number for each voxel (if it
     has been assigned to a segment).
  */
  mri_labeled = MRIalloc(width, height, depth, MRI_SHORT);
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) MRISvox(mri_labeled, x, y, z) = -1;
    }
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        // if within the range
        if (val >= low_val && val <= hi_val) {
          count++;

          // initializer for border_labels
          memset(border_labels, -1, NBR_VOX * sizeof(int));
          for (nvox = 0, zk = -1; zk <= 1; zk++) {
            zi = z + zk;
            if ((zi < 0) || (zi >= depth)) continue;
            for (yk = -1; yk <= 1; yk++) {
              yi = y + yk;
              if ((yi < 0) || (yi >= height)) continue;
              // increment nvox count here
              for (xk = -1; xk <= 1; xk++, nvox++) {
#if 1
                if ((abs(xk) + abs(yk) + abs(zk)) > 1) continue; /* only allow 4 (6 in 3-d) connectivity */
#endif

                xi = x + xk;
                if ((xi < 0) || (xi >= width)) continue;
                // get the neighbor label
                label = MRISvox(mri_labeled, xi, yi, zi);
                if (label >= 0) border_labels[nvox] = label;
              }
            }
          }
          // count nonzero labels in nbrs
          for (nlabels = nvox = 0; nvox < NBR_VOX; nvox++) {
            label = border_labels[nvox];
            if ((label >= 0) && (!mriseg->segments[label].found)) {
              mriseg->segments[label].found = 1;
              nlabels++;
            }
          }
          // reset found
          for (nvox = 0; nvox < NBR_VOX; nvox++) {
            label = border_labels[nvox];
            if (label >= 0) mriseg->segments[label].found = 0; /* for next time */
          }
          //
          label = 0;
          // create labels for those points which are not connected
          switch (nlabels) {
            case 0: /* allocate a new segment */
              label = mriSegmentNew(mriseg);
              // returns this added segment position
              mseg = &mriseg->segments[label];
              if (0 && DIAG_VERBOSE_ON)
                fprintf(stdout, "allocating new label %d (%d total)\n", label, mriseg->nsegments);
              break;
            case 1: /* assign this voxel to the one that it borders */
              for (nvox = 0; nvox < NBR_VOX; nvox++)
                if (border_labels[nvox] >= 0) {
                  label = border_labels[nvox];
                  break;
                }
              // points to the same label position
              mseg = &mriseg->segments[label];
              break;
            default: /* merge segments and assign to lowest number */
              mseg = NULL;
              for (nvox = 0; nvox < NBR_VOX; nvox++) {
                if (border_labels[nvox] >= 0) {
                  // the 1st encountered label case
                  if (!mseg) {
                    // set the first index position
                    label = border_labels[nvox];
                    mseg = &mriseg->segments[label];
                    mseg->found = 1;
                  }
                  // the rest
                  else {
                    mseg2 = &mriseg->segments[border_labels[nvox]];
                    if (mseg2->found == 0) {
                      mseg2->found = 1; /* prevent merging more than once */
                      // merge to the label position
                      if (mriSegmentMerge(mriseg, label, border_labels[nvox], mri_labeled) != NO_ERROR) {
                        // nsegments decreased by one
                        MRIsegmentFree(&mriseg);
                        return (NULL);
                      }
                    }
                  }
                }
              }
              // reset found
              for (nvox = 0; nvox < NBR_VOX; nvox++)
                if (border_labels[nvox] >= 0) mriseg->segments[border_labels[nvox]].found = 0;
              break;
          }
          /* add it to the existing list */
          if (mseg->nvoxels >= mseg->max_voxels) {
            // this max could be the same as mseg->max_voxels
            // this is taken care in mriSegmentReallocateVoxels()
            max = nint(mseg->max_voxels * VOX_INCREASE);
            // if (mriSegmentReallocateVoxels(mriseg, label,
            //            nint(mseg->max_voxels*VOX_INCREASE))
            //    != NO_ERROR)
            if (mriSegmentReallocateVoxels(mriseg, label, max) != NO_ERROR) {
              MRIsegmentFree(&mriseg);
              return (NULL);
            }
          }
          mseg->voxels[mseg->nvoxels].x = x;
          mseg->voxels[mseg->nvoxels].y = y;
          mseg->voxels[mseg->nvoxels].z = z;
          mseg->nvoxels++;
          mseg->area += voxel_size;  // voxel_size = volume
#if 0
          // this is for only 1mm voxel case
          if (mseg->nvoxels != (int)mseg->area)
            DiagBreak() ;
#endif
          MRISvox(mri_labeled, x, y, z) = label;
        }
      }
    }
  }
  MRIcompactSegments(mriseg);

  if (0 && Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_labeled, "labeled.mgh");
  MRIfree(&mri_labeled);
  mriComputeSegmentStatistics(mriseg);
  return (mriseg);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIsegmentFree(MRI_SEGMENTATION **pmriseg)
{
  MRI_SEGMENTATION *mriseg;
  int i;

  mriseg = *pmriseg;
  *pmriseg = NULL;
  for (i = 0; i < mriseg->max_segments; i++)
    if (mriseg->segments[i].voxels) free(mriseg->segments[i].voxels);

  free(mriseg->segments);
  free(mriseg);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SEGMENTATION *MRIsegmentAlloc(int max_segments, int max_voxels)
{
  MRI_SEGMENTATION *mriseg;
  int i;

  mriseg = (MRI_SEGMENTATION *)calloc(1, sizeof(MRI_SEGMENTATION));
  if (!mriseg) ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mriseg");
  mriseg->segments = (MRI_SEGMENT *)calloc(max_segments, sizeof(MRI_SEGMENT));
  if (!mriseg->segments) ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mrisegments");
  mriseg->max_segments = max_segments;

  for (i = 0; i < max_segments; i++) {
    mriseg->segments[i].area = 0.0f;
    mriseg->segments[i].nvoxels = 0;
    mriseg->segments[i].max_voxels = max_voxels;
    mriseg->segments[i].voxels = (MSV *)calloc(max_voxels, sizeof(MSV));
    if (!mriseg->segments[i].voxels) ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc %dth voxel list", i);
  }

  return (mriseg);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Merge segment 1 into segment 0.
  ------------------------------------------------------*/
static int mriSegmentMerge(MRI_SEGMENTATION *mriseg, int s0, int s1, MRI *mri_labeled)
{
  MRI_SEGMENT *mseg0, *mseg1;
  int v, total_voxels, x, y, z;

  if (0 && DIAG_VERBOSE_ON) fprintf(stdout, "merging segments %d and %d\n", s0, s1);
  mseg0 = &mriseg->segments[s0];
  mseg1 = &mriseg->segments[s1];
  total_voxels = mseg0->nvoxels + mseg1->nvoxels;
  if (total_voxels >= mseg0->max_voxels) {
    if (mriSegmentReallocateVoxels(mriseg, s0, total_voxels + 10) != NO_ERROR) return (Gerror);
  }

  for (v = mseg0->nvoxels; v < total_voxels; v++) {
    x = mseg1->voxels[v - mseg0->nvoxels].x;
    y = mseg1->voxels[v - mseg0->nvoxels].y;
    z = mseg1->voxels[v - mseg0->nvoxels].z;
    mseg0->voxels[v].x = x;
    mseg0->voxels[v].y = y;
    mseg0->voxels[v].z = z;
    if (mri_labeled) MRISvox(mri_labeled, x, y, z) = s0;
  }
  mseg0->nvoxels = total_voxels;
  mseg0->area += mseg1->area;
#if 0
  // this check is needed only for 1mm voxel size
  if (mseg0->nvoxels != (int)mseg0->area)
    DiagBreak() ;
#endif
  // decrease nsegments by one
  mriseg->nsegments--;
  mseg1->nvoxels = 0;
  mseg1->area = 0;

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mriSegmentReallocateVoxels(MRI_SEGMENTATION *mriseg, int sno, int max_voxels)
{
  MRI_SEGMENT *mseg;
  // MSV         *old_voxels ;
  MSV *new_voxels;
  int v;
  size_t memsize;

  if (max_voxels <= 0) max_voxels = MAX_VOXELS;

  mseg = &mriseg->segments[sno];
  // old_voxels = mseg->voxels ;

  // fprintf(stdout, "ReallocVoxels:%8d -> %8d\n",
  // mseg->max_voxels, max_voxels);
  // mseg->voxels = (MSV *)calloc(max_voxels, sizeof(MSV)) ;

  // if it is the same, then increase by 10
  if (mseg->max_voxels == max_voxels) max_voxels = mseg->max_voxels + 10;

  memsize = max_voxels * sizeof(MSV);
  new_voxels = (MSV *)realloc(mseg->voxels, memsize);
  if (!new_voxels)
  // if (!mseg->voxels)
  {
    // if (old_voxels)
    //   free(old_voxels) ;
    // when fail, mseg->voxels not modified
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "mriSegmentReallocateVoxels: "
                 "could not alloc %d voxels (%d KB) for sno %d",
                 max_voxels,
                 memsize / 1024,
                 sno));
  }
  // realloc
  mseg->voxels = new_voxels;
  mseg->max_voxels = max_voxels;
  // realloc lack of initialization
  for (v = mseg->nvoxels; v < mseg->max_voxels; ++v) {
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

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mriSegmentReallocateSegments(MRI_SEGMENTATION *mriseg, int max_segments)
{
  MRI_SEGMENT *old_segments;
  int s;

  old_segments = mriseg->segments;
  mriseg->segments = (MRI_SEGMENT *)calloc(max_segments, sizeof(MRI_SEGMENT));
  if (!mriseg->segments) ErrorExit(ERROR_NOMEMORY, "MRIsegmentAlloc: could not alloc mrisegments");
  mriseg->max_segments = max_segments;

  for (s = 0; s < mriseg->nsegments; s++) {
    if (s == 1) DiagBreak();
    mriseg->segments[s].area = old_segments[s].area;
    mriseg->segments[s].nvoxels = old_segments[s].nvoxels;
    mriseg->segments[s].max_voxels = old_segments[s].max_voxels;
    mriseg->segments[s].voxels = old_segments[s].voxels;
  }
  free(old_segments);
  for (s = mriseg->nsegments; s < mriseg->max_segments; s++) {
    mriseg->segments[s].voxels = (MSV *)calloc(MAX_VOXELS, sizeof(MSV));
    if (!mriseg->segments[s].voxels) ErrorExit(ERROR_NOMEMORY, "mriSegmentRealloc: could not alloc voxels");
    mriseg->segments[s].max_voxels = MAX_VOXELS;
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIsegmentToImage(MRI *mri_src, MRI *mri_dst, MRI_SEGMENTATION *mriseg, int s)
{
  int v, x, y, z, smin, smax;
  MRI_SEGMENT *mseg;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  if (s >= mriseg->nsegments) ErrorReturn(mri_dst, (ERROR_BADPARM, "MRIsegmentToImage: invalid segment #%d", s));

  if (s < 0)  // do all segments
  {
    smin = 0;
    smax = mriseg->nsegments - 1;
  }
  else  // just do one segment
    smin = smax = s;

  for (s = smin; s <= smax; s++) {
    mseg = &mriseg->segments[s];

    for (v = 0; v < mseg->nvoxels; v++) {
      x = mseg->voxels[v].x;
      y = mseg->voxels[v].y;
      z = mseg->voxels[v].z;
      MRIsetVoxVal(mri_dst, x, y, z, 0, MRIgetVoxVal(mri_src, x, y, z, 0));
    }
  }

  return (mri_dst);
}

int MRIsetSegmentValue(MRI *mri, MRI_SEGMENTATION *mriseg, int s, float val)
{
  int v, x, y, z;
  MRI_SEGMENT *mseg;

  mseg = &mriseg->segments[s];

  for (v = 0; v < mseg->nvoxels; v++) {
    x = mseg->voxels[v].x;
    y = mseg->voxels[v].y;
    z = mseg->voxels[v].z;
    MRIsetVoxVal(mri, x, y, z, 0, val);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mriSegmentNew(MRI_SEGMENTATION *mriseg)
{
  int s;

  if (mriseg->nsegments >= mriseg->max_segments) {
    int newseglength = mriseg->max_segments * VOX_INCREASE;
    // adjust the length when the increase is too small
    // this happens only when (VOX_INCREASE-1)*max_segment < 1
    if (newseglength == mriseg->max_segments) {
      if (mriseg->max_segments < MAX_SEGMENTS / 2)
        newseglength = MAX_SEGMENTS / 2;
      else
        newseglength = mriseg->max_segments + 1024;
    }
    // now allocate
    mriSegmentReallocateSegments(mriseg, newseglength);
  }
  mriseg->nsegments++;
  for (s = 0; s < mriseg->nsegments; s++)
    if (mriseg->segments[s].nvoxels == 0) return (s);
  return (s);  // return mriseg->nsegments-1
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIcompactSegments(MRI_SEGMENTATION *mriseg)
{
  MRI_SEGMENT *src_segment, *dst_segment;
  MRI_SEGMENT *newseg;
  int s, s2;

  if (mriseg->nsegments == 0) return (NO_ERROR);
  if (0 && DIAG_VERBOSE_ON) fprintf(stdout, "compacting segments...\n");

  for (s = 0; s < mriseg->max_segments; s++) {
    // if voxel count = 0, then
    if (mriseg->segments[s].nvoxels == 0) {
      // assgin dst
      dst_segment = &mriseg->segments[s];
      // find the next non-zero voxel segments
      for (s2 = s + 1; s2 < mriseg->max_segments; s2++) {
        if (mriseg->segments[s2].nvoxels > 0) break;
      }
      // if within the boundary
      if (s2 < mriseg->max_segments) {
        // assin src
        src_segment = &mriseg->segments[s2];
        if (dst_segment->voxels) free(dst_segment->voxels);
        // copy src to dst
        dst_segment->area = src_segment->area;
        dst_segment->nvoxels = src_segment->nvoxels;
        dst_segment->voxels = src_segment->voxels;
        dst_segment->max_voxels = src_segment->max_voxels;
        // make the src to have nvoxels = 0
        src_segment->nvoxels = src_segment->max_voxels = 0;
        src_segment->area = 0.0;
        src_segment->voxels = NULL;
        // now we moved the non-zero segment to the front
      }
    }
  }
  //////////////////////////////////
  // now we have segments zero voxels ones stack at the end
  // now we keep only the non-zero voxel segments
  // find the location of the first non-zero voxel segments
  for (s = 0; s < mriseg->max_segments; s++) {
    if (mriseg->segments[s].nvoxels == 0) break;
  }
  mriseg->max_segments = s;
  mriseg->nsegments = s;
  if (mriseg->nsegments == 0) return (NO_ERROR);  // nothing left after compacting

  newseg = (MRI_SEGMENT *)realloc(mriseg->segments, s * sizeof(MRI_SEGMENT));
  if (newseg) {
    mriseg->segments = newseg;
    if (0 && DIAG_VERBOSE_ON) fprintf(stdout, "\n        segments reduced to %d\n", s);
  }
  else {
    MRI_SEGMENT *oldseg = mriseg->segments;

    newseg = (MRI_SEGMENT *)calloc(s, sizeof(MRI_SEGMENT));
    if (newseg == NULL)
      ErrorExit(ERROR_NOMEMORY, "MRIcompactSegments(): could not alloc mriseg with %d segments", mriseg->nsegments);
    memmove(newseg, oldseg, s * sizeof(MRI_SEGMENT));
    free(oldseg);
    mriseg->segments = newseg;
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIsegmentDilate(MRI_SEGMENTATION *mriseg, MRI *mri)
{
  int x, y, z, xk, yk, zk, xi, yi, zi, segno, v, nvox;
  MRI_SEGMENT *mseg;
  MRI *mri_segments;
  float voxel_size;

  voxel_size = mri->xsize * mri->ysize * mri->zsize;

  mri_segments = MRIclone(mri, NULL);

  /* build image of all other segments to prevent dilation from
     expanding into other segments */
  for (segno = 0; segno < mriseg->nsegments; segno++) MRIsegmentToImage(mri, mri_segments, mriseg, segno);

  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];

    nvox = mseg->nvoxels;
    for (v = 0; v < nvox; v++) {
      x = mseg->voxels[v].x;
      y = mseg->voxels[v].y;
      z = mseg->voxels[v].z;

      for (xk = -1; xk <= 1; xk++) {
        xi = mri->xi[x + xk];
        for (yk = -1; yk <= 1; yk++) {
          yi = mri->yi[y + yk];
          for (zk = -1; zk <= 1; zk++) {
            if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1) continue;
            zi = mri->zi[z + zk];
            if (MRIgetVoxVal(mri, xi, yi, zi, 0) && !MRIgetVoxVal(mri_segments, xi, yi, zi, 0)) {
              if (mseg->nvoxels >= mseg->max_voxels) {
                if (mriSegmentReallocateVoxels(mriseg, segno, mseg->max_voxels * VOX_INCREASE) != NO_ERROR) {
                  /*              MRIsegmentFree(&mseg) ;*/
                  return (Gerror);
                }
              }
              mseg->voxels[mseg->nvoxels].x = xi;
              mseg->voxels[mseg->nvoxels].y = yi;
              mseg->voxels[mseg->nvoxels].z = zi;
              MRIsetVoxVal(mri_segments, xi, yi, zi, 0, MRIgetVoxVal(mri, xi, yi, zi, 0));
              mseg->nvoxels++;
              mseg->area += voxel_size;
            }
          }
        }
      }
    }
  }
  MRIfree(&mri_segments);
  return (NO_ERROR);
}

int MRIremoveSmallSegments(MRI_SEGMENTATION *mriseg, int min_voxels)
{
  int s;

  if (0 && DIAG_VERBOSE_ON) fprintf(stdout, "compacting segments...\n");

  for (s = 0; s < mriseg->max_segments; s++)
    if (mriseg->segments[s].nvoxels < min_voxels) {
      if (mriseg->segments[s].voxels) {
        free(mriseg->segments[s].voxels);
        mriseg->segments[s].voxels = NULL;
      }
      mriseg->segments[s].nvoxels = 0;
    }
  return (MRIcompactSegments(mriseg));
}

int MRIeraseSmallSegments(MRI_SEGMENTATION *mriseg, MRI *mri_seg, int min_voxels)
{
  int s, v, nerased = 0, x, y, z, f;
  MRI_SEGMENT *mseg;

  if (0 && DIAG_VERBOSE_ON) fprintf(stdout, "compacting segments...\n");

  for (s = 0; s < mriseg->max_segments; s++) {
    mseg = &mriseg->segments[s];
    if (mseg->nvoxels < min_voxels) {
      if (mseg->voxels)
        for (v = 0; v < mseg->nvoxels; v++) {
          x = mseg->voxels[v].x;
          y = mseg->voxels[v].y;
          z = mseg->voxels[v].z;
          if (x == Gx && y == Gy && z == Gz) printf("MRIeraseSmallSegments: erasing (%d, %d %d)\n", x, y, z);
          for (f = 0; f < mri_seg->nframes; f++) MRIsetVoxVal(mri_seg, x, y, z, f, 0);
          nerased++;
        }
    }
  }
  return (nerased);
}

static int mriComputeSegmentStatistics(MRI_SEGMENTATION *mriseg)
{
  int segno;
  MRI_SEGMENT *mseg;
  int v, x, y, z;

  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];
    mseg->x0 = mseg->y0 = mseg->z0 = 10000;
    mseg->x1 = mseg->y1 = mseg->z1 = -10000;
    mseg->cx = mseg->cy = mseg->cz = 0.0;
    if (!mseg->nvoxels) continue;
    for (v = 0; v < mseg->nvoxels; v++) {
      x = mseg->voxels[v].x;
      y = mseg->voxels[v].y;
      z = mseg->voxels[v].z;
      mseg->cx += x;
      mseg->cy += y;
      mseg->cz += z;
      if (x < mseg->x0) mseg->x0 = x;
      if (y < mseg->y0) mseg->y0 = y;
      if (z < mseg->z0) mseg->z0 = z;
      if (x > mseg->x1) mseg->x1 = x;
      if (y > mseg->y1) mseg->y1 = y;
      if (z > mseg->z1) mseg->z1 = z;
    }
    mseg->cx /= (float)mseg->nvoxels;
    mseg->cy /= (float)mseg->nvoxels;
    mseg->cz /= (float)mseg->nvoxels;
  }

  return (NO_ERROR);
}

int MRIsegmentDilateThreshold(MRI_SEGMENTATION *mriseg, MRI *mri_binary, MRI *mri_thresh, int low_thresh, int hi_thresh)
{
  int x, y, z, xk, yk, zk, xi, yi, zi, segno, v, nvox, val;
  MRI_SEGMENT *mseg;
  MRI *mri_segments;
  float voxel_size;

  voxel_size = mri_binary->xsize * mri_binary->ysize * mri_binary->zsize;

  mri_segments = MRIclone(mri_binary, NULL);

  /* build image of all other segments to prevent dilation from
     expanding into other segments */
  for (segno = 0; segno < mriseg->nsegments; segno++) MRIsegmentToImage(mri_binary, mri_segments, mriseg, segno);

  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];

    nvox = mseg->nvoxels;
    for (v = 0; v < nvox; v++) {
      x = mseg->voxels[v].x;
      y = mseg->voxels[v].y;
      z = mseg->voxels[v].z;

      for (xk = -1; xk <= 1; xk++) {
        xi = mri_binary->xi[x + xk];
        for (yk = -1; yk <= 1; yk++) {
          yi = mri_binary->yi[y + yk];
          for (zk = -1; zk <= 1; zk++) {
            if ((fabs(xk) + fabs(yk) + fabs(zk)) != 1) continue;
            zi = mri_binary->zi[z + zk];
            val = MRIgetVoxVal(mri_thresh, xi, yi, zi, 0);
            if (MRIgetVoxVal(mri_binary, xi, yi, zi, 0) && !MRIgetVoxVal(mri_segments, xi, yi, zi, 0) &&
                (val >= low_thresh) && (val <= hi_thresh)) {
              if (mseg->nvoxels >= mseg->max_voxels) {
                if (mriSegmentReallocateVoxels(mriseg, segno, nint(mseg->max_voxels * VOX_INCREASE)) != NO_ERROR) {
                  /*                  MRIsegmentFree(&mseg) ;*/
                  return (Gerror);
                }
              }
              mseg->voxels[mseg->nvoxels].x = xi;
              mseg->voxels[mseg->nvoxels].y = yi;
              mseg->voxels[mseg->nvoxels].z = zi;
              MRIsetVoxVal(mri_segments, xi, yi, zi, 0, MRIgetVoxVal(mri_binary, xi, yi, zi, 0));
              mseg->nvoxels++;
              mseg->area += voxel_size;
            }
          }
        }
      }
    }
  }
  MRIfree(&mri_segments);
  return (NO_ERROR);
}

int MRIsegmentMax(MRI_SEGMENTATION *mriseg)
{
  int segno, max_voxels, max_segno, nvox;
  MRI_SEGMENT *mseg;

  max_segno = -1;
  max_voxels = 0;
  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];
    if (mseg->ignore) continue;

    nvox = mseg->nvoxels;
    if (nvox > max_voxels) {
      max_segno = segno;
      max_voxels = nvox;
    }
  }
  return (max_segno);
}

int MRIsegmentClearIgnoreFlags(MRI_SEGMENTATION *mriseg)
{
  int segno;
  MRI_SEGMENT *mseg;

  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];
    mseg->ignore = 0;
  }
  return (NO_ERROR);
}

int MRIfindMaxSegmentNumber(MRI_SEGMENTATION *mriseg)
{
  MRI_SEGMENT *mseg;
  int max_segno = 0, segno, max_voxels = 0;

  for (segno = 0; segno < mriseg->nsegments; segno++) {
    mseg = &mriseg->segments[segno];
    if (mseg->nvoxels >= max_voxels) {
      max_voxels = mseg->nvoxels;
      max_segno = segno;
    }
  }
  return (max_segno);
}

MRI *MRIsegmentFill(MRI_SEGMENTATION *mriseg, int s, MRI *mri, float fillval)
{
  int v, x, y, z;
  MRI_SEGMENT *mseg;

  if (s < 0 || s >= mriseg->nsegments) ErrorReturn(NULL, (ERROR_BADPARM, "MRIsegmentFill: invalid segment #%d", s));
  mseg = &mriseg->segments[s];

  if (mri == NULL) mri = MRIclone(mriseg->mri, NULL);

  for (v = 0; v < mseg->nvoxels; v++) {
    x = mseg->voxels[v].x;
    y = mseg->voxels[v].y;
    z = mseg->voxels[v].z;
    MRIsetVoxVal(mri, x, y, z, 0, fillval);
  }

  return (mri);
}

#include "label.h"
LABEL *MRIsegmentToLabel(MRI_SEGMENTATION *mriseg, MRI *mri, int ind)
{
  LABEL *area;
  int i;
  MRI_SEGMENT *mseg;
  double xs, ys, zs;

  mseg = &mriseg->segments[ind];
  area = LabelAlloc(mseg->nvoxels, NULL, NULL);

  area->coords = LABEL_COORDS_TKREG_RAS;
  area->avg_stat = 0;
  for (i = 0; i < mseg->nvoxels; i++) {
    MRIvoxelToSurfaceRAS(mri, mseg->voxels[i].x, mseg->voxels[i].y, mseg->voxels[i].z, &xs, &ys, &zs);
    area->lv[i].x = xs;
    ;
    area->lv[i].y = ys;
    ;
    area->lv[i].z = zs;
    ;
    area->lv[i].vno = -1;
    area->lv[i].stat = MRIgetVoxVal(mri, mseg->voxels[i].x, mseg->voxels[i].y, mseg->voxels[i].z, 0);
    area->avg_stat += area->lv[i].stat;
  }
  area->n_points = mseg->nvoxels;
  if (mseg->nvoxels > 0) area->avg_stat /= mseg->nvoxels;
  strcpy(area->space, "TkReg");

  return (area);
}
