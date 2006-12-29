/**
 * @file  voxlist.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.5 $
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


#ifndef VOXLIST_H
#define VOXLIST_H

#include "mri.h"

typedef struct
{
  int     *xi ;
  int     *yi ;
  int     *zi ;
  float   *xd ;  // transformed parameters
  float   *yd ;
  float   *zd ;
  float   *vsrc ;
  float   *vdst ;
  MRI     *mri ;
  MRI     *mri2;
  int     nvox ;
  double  mean;
  double  std;
}
VOXEL_LIST ;

MRI         *VLSTtoMri(VOXEL_LIST *vl, MRI *mri) ;
int         VLSTfree(VOXEL_LIST **pvoxel_list) ;
VOXEL_LIST  *VLSTcreate(MRI *mri, float low_val, float hi_val ,
                        VOXEL_LIST *vl, int skip, int border_only) ;
VOXEL_LIST  *VLSTcreateInRegion(MRI *mri, float low_val, float hi_val ,
                                VOXEL_LIST *vl, int skip, int border_only, MRI_REGION *box) ;
int         VLSTtransform(VOXEL_LIST *vl, MATRIX *m, MRI *mri, int sample_type) ;

MRI         *VLSTcreateMri(VOXEL_LIST *vl, int val) ;
MRI         *VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val) ;
VOXEL_LIST  *VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude) ;
void VLSTcomputeStats(VOXEL_LIST *vl);

#define VL_DILATE_ADD      0
#define VL_DILATE_REPLACE  1

#endif
