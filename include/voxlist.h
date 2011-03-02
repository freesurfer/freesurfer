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
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.12 $
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


#ifndef VOXLIST_H
#define VOXLIST_H

#include "mri.h"


#if defined(__cplusplus)
extern "C" {
#endif


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
  MRI     *mri_grad;
  int     nvox ;
  double  mean;
  double  std;
}
VOXEL_LIST, VOXLIST ;

VOXEL_LIST *VLSTalloc(int nvox) ;
VOXEL_LIST *VLSTcopy(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_index, int num) ;
MRI         *VLSTtoMri(VOXEL_LIST *vl, MRI *mri) ;
int         VLSTfree(VOXEL_LIST **pvoxel_list) ;
VOXEL_LIST  *VLSTcreateFromDifference(MRI *mri1, MRI *mri2, VOXEL_LIST *vl, int target_label) ;
VOXEL_LIST  *VLSTcreate(MRI *mri, float low_val, float hi_val ,
                        VOXEL_LIST *vl, int skip, int border_only) ;
VOXEL_LIST  *VLSTcreateInRegion(MRI *mri, float low_val, float hi_val ,
                                VOXEL_LIST *vl, int skip, int border_only, MRI_REGION *box) ;
int         VLSTtransformCoords(VOXEL_LIST *vl, MATRIX *m, int skip) ;
int         VLSTtransform(VOXEL_LIST *vl, MATRIX *m, MRI *mri, int sample_type) ;

MRI         *VLSTcreateMri(VOXEL_LIST *vl, int val) ;
MRI         *VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val) ;
VOXEL_LIST  *VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude) ;
void VLSTcomputeStats(VOXEL_LIST *vl);
VOXEL_LIST *VLSTsort(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst) ;

#define VL_DILATE_ADD      0
#define VL_DILATE_REPLACE  1



#if defined(__cplusplus)
};
#endif


#endif
