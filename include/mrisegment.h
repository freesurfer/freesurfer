/**
 * @file  mrisegment.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2010/11/10 01:45:17 $
 *    $Revision: 1.11 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
int              MRIremoveSmallSegments(MRI_SEGMENTATION *mriseg,
                                        int min_voxels) ;
int              MRIsegmentMax(MRI_SEGMENTATION *mriseg) ;
int              MRIsegmentClearIgnoreFlags(MRI_SEGMENTATION *mriseg) ;
int              MRIfindMaxSegmentNumber(MRI_SEGMENTATION *mriseg) ;
int              MRIsetSegmentValue(MRI *mri, MRI_SEGMENTATION *mriseg, int s, float val) ;

#endif
