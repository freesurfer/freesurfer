#ifndef MRI_SEGMENT_H
#define MRI_SEGMENT_H


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


#endif
