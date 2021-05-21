/*
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


#ifndef MRI_SEGMENT_H
#define MRI_SEGMENT_H

#include "mri.h"

// in order to conserve memory, we change int to short
typedef struct
{
  short x, y, z; // int   x, y, z ;
}
MRI_SEGMENT_VOXEL, MSV ;

typedef struct
{
  float    area ;
  MSV      *voxels ;
  int      nvoxels ;
  int      max_voxels ;
  int      label ;
  int      found ;      /* for use during segmentation */
  float    cx, cy, cz ; /* centroid */
  short x0, x1; // int      x0, x1 ;
  short y0, y1; // int      y0, y1 ;
  short z0, z1; // int      z0, z1 ;
  short ignore; // int      ignore;
}
MRI_SEGMENT ;

typedef struct
{
  MRI_SEGMENT  *segments ;
  int          nsegments ;
  int          max_segments ;
  MRI          *mri ;
}
MRI_SEGMENTATION ;

MRI_SEGMENTATION *MRIsegment(MRI *mri, float low_val, float hi_val) ;
// the next one cares about only finding the maxsegment region
MRI_SEGMENTATION *MRImaxsegment(MRI *mri, float low_val, float hi_val) ;

int              MRIsegmentFree(MRI_SEGMENTATION **pmriseg) ;
MRI_SEGMENTATION *MRIsegmentAlloc(int max_segments, int max_voxels) ;
MRI              *MRIsegmentToImage(MRI *mri_src, MRI *mri_dst,
                                    MRI_SEGMENTATION *mriseg, int s) ;
MRI              *MRIsegmentFill(MRI_SEGMENTATION *mriseg, int s, MRI *mri, float fillval) ;
int              MRIsegmentDilate(MRI_SEGMENTATION *mriseg, MRI *mri) ;
int              MRIsegmentDilateThreshold(MRI_SEGMENTATION *mriseg,
    MRI *mri_binary, MRI *mri_thresh,
    int low_thresh, int hi_thresh) ;
int              MRIcompactSegments(MRI_SEGMENTATION *mriseg) ;
int MRIeraseSmallSegments(MRI_SEGMENTATION *mriseg, MRI *mri_seg, int min_voxels) ;
int              MRIremoveSmallSegments(MRI_SEGMENTATION *mriseg,
                                        int min_voxels) ;
int              MRIsegmentMax(MRI_SEGMENTATION *mriseg) ;
int              MRIsegmentClearIgnoreFlags(MRI_SEGMENTATION *mriseg) ;
int              MRIfindMaxSegmentNumber(MRI_SEGMENTATION *mriseg) ;
int              MRIsetSegmentValue(MRI *mri, MRI_SEGMENTATION *mriseg, int s, float val) ;
#include "label.h"
LABEL            *MRIsegmentToLabel(MRI_SEGMENTATION *mriseg, MRI *mri, int ind) ;

#endif
