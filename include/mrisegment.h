#ifndef MRI_SEGMENT_H
#define MRI_SEGMENT_H

#include "mri.h"

typedef struct
{
  int   x, y, z ;
} MRI_SEGMENT_VOXEL, MSV ;

typedef struct
{
  float    area ;
  MSV      *voxels ;
  int      nvoxels ;
  int      max_voxels ;
  int      label ;
  int      found ;      /* for use during segmentation */
  float    cx, cy, cz ; /* centroid */
  int      x0, x1 ;
  int      y0, y1 ;
  int      z0, z1 ;
} MRI_SEGMENT ;

typedef struct
{
  MRI_SEGMENT  *segments ;
  int          nsegments ;
  int          max_segments ;
  MRI          *mri ;
} MRI_SEGMENTATION ;

MRI_SEGMENTATION *MRIsegment(MRI *mri, float low_val, float hi_val) ;
int              MRIsegmentFree(MRI_SEGMENTATION **pmriseg) ;
MRI_SEGMENTATION *MRIsegmentAlloc(int max_segments, int max_voxels) ;
MRI              *MRIsegmentToImage(MRI *mri_src, MRI *mri_dst, 
                                    MRI_SEGMENTATION *mriseg, int s) ;
int              MRIsegmentDilate(MRI_SEGMENTATION *mriseg, MRI *mri) ;
int              MRIsegmentDilateThreshold(MRI_SEGMENTATION *mriseg, 
                                           MRI *mri_binary, MRI *mri_thresh,
                                           int low_thresh, int hi_thresh) ;
int              MRIcompactSegments(MRI_SEGMENTATION *mriseg) ;
int              MRIremoveSmallSegments(MRI_SEGMENTATION *mriseg,
                                        int min_voxels) ;
int              MRIsegmentMax(MRI_SEGMENTATION *mriseg) ;

#endif
