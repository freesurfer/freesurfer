#ifndef MRINORM_H
#define MRINORM_H

#include "mri.h"
#include "region.h"
#include "histo.h"

/* defaults which can be used in MRInormInit call (using 0 will install these)
 */
#define DEFAULT_DESIRED_WHITE_MATTER_VALUE  110
#define DEFAULT_SMOOTH_SIGMA                2.0f
#define DEFAULT_WINDOWS_BELOW_T0            9 /* 4*/
#define DEFAULT_WINDOWS_ABOVE_T0            14 /* 9*/
#define DEFAULT_WINDOW_SIZE                 10

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
#define MAX_GRADIENT                        (3.0f/5.0f)

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
} MRI_NORM_INFO, MNI ;
  

MRI *MRIhistoNormalize(MRI *mri_src,MRI *mri_norm, MRI *mri_template, int low);
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


#endif
