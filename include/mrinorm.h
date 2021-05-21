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


#ifndef MRINORM_H
#define MRINORM_H

#include "mri.h"
#include "region.h"
#include "histo.h"
#include "label.h"

/* defaults which can be used in MRInormInit call (using 0 will install these)
 */
#define CONTROL_NONE       0
#define CONTROL_MARKED     1
#define CONTROL_TMP        2
#define CONTROL_NBR        3

#define DEFAULT_DESIRED_WHITE_MATTER_VALUE  110
#define DEFAULT_SMOOTH_SIGMA                2.0f
#define DEFAULT_WINDOWS_BELOW_T0            9 /* 4*/
#define DEFAULT_WINDOWS_ABOVE_T0            14 /* 9*/
#define DEFAULT_WINDOW_SIZE                 10

#define WM_MEAN   DEFAULT_DESIRED_WHITE_MATTER_VALUE
#define WHITE_MATTER_MEAN  WM_MEAN


/* minimum size of peak, relative to total # of voxels in histogram
   (excluding background voxels).
*/
#define MIN_HISTO_PCT   0.15f

/* minimum separation of adjacent peaks */
#define HISTO_WINDOW_SIZE  7

/* amount be which window size decreases in superior direction
   via SIZE_MOD ^ (window # above t0)
*/
#define SIZE_MOD                            1.0f /* 0.95*/

/* 20% more of brain behind (0,0,0) than in front */
#define Z_OFFSET_SCALE                      .2

/* 50% overlap in adjacent windows */
#define OVERLAP                             0.5f

/* maximum gradient measured in (intensity change/voxel) */
/* allow about 3 voxel change between adjacent windows */
#define MAX_GRADIENT                        (5.0f/5.0f)
#define MAX_SKIPPED                         2

/* relative to slice at Talairach origin (i.e. 4 slices above it) */
#define STARTING_SLICE                      4
#define SLICE_OFFSET                        3

#define HISTO_BINS                          0
#define BACKGROUND_INTENSITY                30

#define MAX_SPLINE_POINTS                   80

typedef struct
{
  int         windows_above_t0 ;  /* # of windows above Talairach origin */
  int         windows_below_t0 ;  /* # of windows below Talairach origin */
  int         desired_wm_value ;  /* desired intensity value for white matter*/
  float       smooth_sigma ;      /* sigma of gaussian for smoothing histo */
  HISTOGRAM   histograms[MAX_SPLINE_POINTS] ;
  MRI_REGION  regions[MAX_SPLINE_POINTS] ;
  float       max_gradient ;
}
MRI_NORM_INFO, MNI ;


MRI *MRIhistoNormalize(MRI *mri_src,MRI *mri_norm, MRI *mri_template, int low,
                       int high);
MRI *MRIsplineNormalize(MRI *mri_src, MRI *mri_dst, MRI **pmri_bias,
                        float *inputs, float *outputs, int npoints) ;
MRI *MRIadaptiveHistoNormalize(MRI *mri_src, MRI *mri_norm, MRI *mri_template,
                               int wsize, int hsize, int low) ;
MRI *MRIhistoNormalizeRegion(MRI *mri_src, MRI *mri_norm, MRI *mri_template,
                             int low, MRI_REGION *wreg, MRI_REGION *h_src_reg,
                             MRI_REGION *h_tmp_reg) ;
int MRInormInit(MRI *mri, MNI *mni, int windows_above_t0, int windows_below_t0,
                int wsize, int desired_wm_value, float smooth_sigma) ;
int MRInormFillHistograms(MRI *mri, MNI *mni) ;
int MRInormFindPeaks(MNI *mni, float *inputs, float *outputs) ;
MRI *MRInormalize(MRI *mri_src, MRI *mri_dst, MNI *mni) ;
int MRInormCheckPeaks(MNI *mni, float *inputs, float *outputs, int npeaks) ;


MRI *MRInormFindControlPoints(MRI *mri_src, float wm_target,
                              float intensity_above, float intensity_below,
                              MRI *mri_ctrl, int which, int scan_type, MRI *mri_not_control) ;
MRI *MRInormalizeHighSignalLowStd(MRI *mri_src, MRI *mri_dst, float bias_sigma, float wm_target) ;
MRI *MRInormFindHighSignalLowStdControlPoints(MRI *mri_src, MRI *mri_ctrl) ;
MRI *MRInormGentlyFindControlPoints(MRI *mri_src, float wm_target,
                                    float intensity_above,
                                    float intensity_below, MRI *mri_ctrl, MRI *mri_not_control) ;
MRI *MRIbuildBiasImage(MRI *mri_src,MRI *mri_ctrl, MRI *mri_bias, float sigma);
MRI *MRI3dNormalize(MRI *mri_orig, MRI *mri_src, float wm_target, MRI *mri_norm,
                    float intensity_above, float intensity_below,
                    int only_file, int prune, float sigma, int scan_type, MRI *mri_not_control);
MRI *MRI3dGentleNormalize(MRI *mri_src, MRI *mri_bias, float wm_target,
                          MRI *mri_norm, float intensity_above,
                          float intensity_below, int only_file, float bias_sigma, MRI *mri_not_control);
MRI *MRIbuildVoronoiDiagram(MRI *mri_src, MRI *mri_ctrl, MRI *mri_dst);
MRI *MRIsoapBubble(MRI *mri_src, MRI *mri_ctrl, MRI *mri_dst,int niter, float min_change);
MRI *MRIsoapBubbleExpand(MRI *mri_src, MRI *mri_ctrl, MRI *mri_dst,int niter);
int MRI3dUseFileControlPoints(MRI *mri,const char *fname) ;
int MRI3dUseLabelControlPoints(MRI *mri, LABEL *area) ;
int MRI3dWriteControlPoints(char *control_volume_fname) ;
int MRI3dWriteBias(char *bias_volume_fname) ;
MRI *MRIaverageFixedPoints(MRI *mri_src, MRI *mri_ctrl, MRI *mri_dst,int niter) ;
int MRInormAddFileControlPoints(MRI *mri_ctrl, int value, MRI *mri) ;
MRI *MRInormFindControlPointsInWindow(MRI *mri_src, float wm_target,
                                      float intensity_above,
                                      float intensity_below, MRI *mri_ctrl,
                                      float whalf_mm,const char *debug_str,
                                      int *pnctrl, int scan_type, MRI *mri_not_control);
MRI *MRIapplyBiasCorrection(MRI *mri_in, MRI *mri_bias, MRI *mri_out) ;

#endif
