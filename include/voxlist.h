#ifndef VOXLIST_H
#define VOXLIST_H

#include "mri.h"

typedef struct
{
  int *xi ;
  int *yi ;
  int *zi ;
  MRI *mri ;
  MRI *mri2;
  int nvox ;
  double mean;
  double std;
} VOXEL_LIST ;

MRI         *VLSTtoMri(VOXEL_LIST *vl, MRI *mri) ;
int         VLSTfree(VOXEL_LIST **pvoxel_list) ;
VOXEL_LIST  *VLSTcreate(MRI *mri, float low_val, float hi_val , 
			VOXEL_LIST *vl, int skip, int border_only) ;
VOXEL_LIST  *VLSTcreateInRegion(MRI *mri, float low_val, float hi_val , 
			VOXEL_LIST *vl, int skip, int border_only, MRI_REGION *box) ;

MRI         *VLSTcreateMri(VOXEL_LIST *vl, int val) ;
MRI         *VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val) ;
VOXEL_LIST  *VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude) ;
void VLSTcomputeStats(VOXEL_LIST *vl);

#define VL_DILATE_ADD      0
#define VL_DILATE_REPLACE  1

#endif
