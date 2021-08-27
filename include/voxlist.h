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


#ifndef VOXLIST_H
#define VOXLIST_H

#include "mri.h"
#include "label.h"

#define VOXLIST_NORMAL   0
#define VOXLIST_SPLINE   1

typedef struct
{
  int     type ;  // one of the types above. 
  int     *xi ;
  int     *yi ;
  int     *zi ;
  int     *fi ;  // frame index
  float   *xd ;  // transformed parameters
  float   *yd ;
  float   *zd ;
  float   *vsrc ;
  float   *vdst ;
  MRI     *mri ;
  MRI     *mri2;
  MRI     *mri_grad;
  int     nvox ;
  int     max_vox ;
  double  mean;
  double  std;
  float   *mx ;
  float   *my ;
  float   *mz ;  // slopes for splines
  float   *t ;   // parameterization along total curve
}
VOXEL_LIST, VOXLIST ;

VOXEL_LIST *VLSTfromMRI(MRI *mri, int vno) ;
VOXEL_LIST *VLSTalloc(int nvox) ;
VOXEL_LIST *VLSTcopy(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_index, int num) ;
MRI         *VLSTvsrcToMri(VOXEL_LIST *vl, MRI *mri) ;
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
int         VLSTaddUnique(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd);
int         VLSTinList(VOXEL_LIST *vl, int x, int y, int z);
int         VLSTadd(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd) ;
int         VLSTaddWithValue(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd, float vsrc, float vdst) ;

int         VLmostCommonLabel(VOXEL_LIST *vl) ;
int         VLSTwriteLabel(VOXEL_LIST *vl, const char *fname, MRI_SURFACE *mris, MRI *mri) ;
LABEL       *VLSTtoLabel(VOXEL_LIST *vl, MRI_SURFACE *mris, MRI *mri) ;
MRI         *VLSTwriteOrderToMRI(VOXEL_LIST *vl, MRI *mri) ;


#define VL_DILATE_ADD      0
#define VL_DILATE_REPLACE  1


// spline stuff see http://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull.E2.80.93Rom_spline
#define h00(t)   ((1+2*t)*SQR(1-t))
#define h10(t)   (t*SQR(1-t))
#define h01(t)   (t*t*(3-2*t))
#define h11(t)   (t*t*(t-1))
#define X_to_t(xk, xkp1, x)   (((x)-(xk)) / ((xkp1)-(xk)))
#define t_to_X(xk, xkp1, x)   ((t) * ((xkp1)-(xk)) + (xk))


double      VLSTmean(VOXEL_LIST *vl, MRI *mri, double *pvar) ;
int         VLSTsampleFloat(VOXEL_LIST *vl, MRI *mri)  ;
int         VLSTsample(VOXEL_LIST *vl, MRI *mri) ;
VOXEL_LIST  *VLSTsplineFit(VOXEL_LIST *vl, int num_control) ;
VOXEL_LIST  *VLSTinterpolate(VOXEL_LIST *vl, float spacing) ;
int         VLSTinterpolateIntoVolume(VOXEL_LIST *vl, MRI *mri, float val) ;
double      VLSTcomputeEntropy(VOXEL_LIST *vl, MRI *mri, int num) ;
int         VLSTinterpolateSplineIntoVolume(VOXEL_LIST *vl, MRI *mri, double spacing, VOXEL_LIST *vl_total, float val) ;
VOXEL_LIST  *VLSTcopyInto(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_dst_index, int num);

double VLSTcomputeSplineMean(VOXEL_LIST *vl_spline, MRI *mri, double step_size)  ;
double VLSTcomputeSplineSegmentMean(VOXEL_LIST *vl_spline, MRI *mri, double step_size, int start, int stop)  ;
float VLSTcomputeSplineMedian(VOXEL_LIST *vl_spline, MRI *mri, double step_size) ;
  double VLSThausdorffDistance(VOXEL_LIST *vl1, VOXEL_LIST *vl2, double max_dist, MRI **pmri_dist) ;
  double VLSTrmsDistance(VOXEL_LIST *vl1, VOXEL_LIST *vl2, double max_dist, MRI **pmri_dist) ;

#endif
