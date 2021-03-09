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


#ifndef MRI_MORPH_H
#define MRI_MORPH_H

#include "histo.h"
#include "mri.h"
#include "matrix.h"
#include "transform.h"
#include "mrisurf.h"
#include "gca.h"

typedef struct
{
  float      neck_x0 ;       /* position of neck in original coordinates */
  float      neck_y0 ;
  float      neck_z0 ;
  float      neck_dx ;       /* vector along neck in 'don't care' direction */
  float      neck_dy ;
  float      neck_dz ;
}
NECK_PARMS ;

#define MAX_LEVELS 10
typedef struct
{
  double     l_intensity ;   /* coefficient of intensity term */
  double     l_area ;        /* coefficient of area term */
  double     l_nlarea ;      /* coefficient of nonlinear area term */
  double     l_nldist ;      /* coefficient of nonlinear distance term */
  double     l_dist ;        /* coefficient of distance term */
  double     l_compression ; /* coefficient of distance compression term */
  double     exp_k ;         /* nonlinear area exponent coefficient */
  LTA        *lta ;          /* octree transform */
  double     dt ;
  double     momentum ;
  int        niterations ;
  int        write_iterations ;
  char       base_name[100] ;
  FILE       *log_fp ;
  int        start_t ;
  float      sigma ;         /* upper bound on amount of blurring */
  double     tol ;           /* for terminating integration */
  int        levels ;
  MRI        *mri_in, *mri_ref ;  /* for diagnostics to get at originals */
  int        navgs ;         /* # of iterations of gradient averaging */
  NECK_PARMS ref_np ;     /* position and orientation of reference neck */
  NECK_PARMS in_np ;      /* position and orientation of input neck */
  int        rigid ;      /* 1 if doing rigid alignment */
  float      trans_mul ;  /* scaling of translation gradient */
  int        morph_skull ;
  int        disable_neck ;
  MRI        *mri_red_in, *mri_red_ref ;  /* current (reduced) volumes */
  GCA        *gca_red ;
  MATRIX     *m_xform_mean ;        /* cross-subject mean of xform parms */
  MATRIX     *m_xform_covariance ;  /* covariance matrix of xform parms */
  MATRIX     *m_inv_cov ;           /* inverse of above */
  double     l_priors ;             /* weighting given to priors term */
  MATRIX     **m_xforms ;           /* one per subject in  priors */
  int        nxforms ;              /* size of previous array */
  int        max_levels ;
  MRI        *mri_crop ;            /* boolean image 1=cropped region */
  int        scout_flag ;
  MRI        *mri_classified ;
  float      factor ;               /* for stabilizing integration */
  GCAS       *gcas ;
  int        nsamples ;
  TRANSFORM  *transform ;
  void       *vgca ;
  double     clamp ;                // saturation threshold for log likelihood
}
MORPH_PARMS, MP ;


#define INTEGRATION_TOL    1e-4  /*5e-5*/
#define NEIGHBORS          6

typedef struct
{
  float      x, y, z ;       /* current coordinates */
}
MORPH_NODE ;

/* this data structure is separate so that it doesn't have to
   be held in memory all the time.
*/
typedef struct
{
  float      dx, dy, dz ;            /* movement delta */
  float      ox, oy, oz ;    /* original coordinates */
  float      orig_dist[NEIGHBORS] ;  /* original distances to 6 neighbors */
  float      orig_area ;             /* original area */
  float      area ;                  /* current area */
  float      tdx, tdy, tdz ;         /* temp. quantities for gradient avging */
}
MORPH_NODE_PROPERTIES, MNP ;

typedef struct
{
  float      node_spacing ;  /* mm between nodes (isotropic) */
  int        width ;         /* number of nodes wide */
  int        height ;
  int        depth ;
  MRI        *mri_in ;       /* source image */
  MRI        *mri_ref ;      /* target of morph */
  LTA        *lta ;          /* octree of linear trasforms */
  MORPH_NODE ***nodes ;
  MORPH_NODE_PROPERTIES ***pnodes ;
  int        neg ;
}
MORPH_3D, M3D ;


#define LTA_TYPE       0
#define MORPH3D_TYPE   1



HISTOGRAM *MRIhorizontalHistogram(MRI *mri, int thresh_low, int thresh_hi) ;
HISTOGRAM *MRIhorizontalBoundingBoxHistogram(MRI *mri, int thresh) ;

int       MRIcountAboveThreshold(MRI *mri, double thresh) ;
MRI       *MRIlabel(MRI *mri_src, MRI *mri_dst, int *nlabels) ;
int       MRIlabelBoundingBoxes(MRI *mri_label,MRI_REGION *bboxes,int nlabels);
int       MRIeraseOtherLabels(MRI *mri_src, MRI *mri_dst, int label) ;
int       MRIeraseLabel(MRI *mri_src, MRI *mri_dst, int label) ;
int       MRIfindHorizontalLabelLimits(MRI *mri, int label,
                                       int *xmins, int *xmaxs) ;
MRI       *MRIfindNeck(MRI *mri_src, MRI *mri_dst, int thresh_low,
                       int thresh_hi, MORPH_PARMS *parms, int dir,
                       NECK_PARMS *np) ;
int       MRIlabelAreas(MRI *mri_label, float *areas, int nlabels) ;
int       MRIlabelCentroid(MRI *mri_label,int l,float *px,float *py,float *pz);
int       MRIlinearAlign(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms);
int       MRIrigidAlign(MRI *mri_in,MRI *mri_ref, MORPH_PARMS *parms,
                        MATRIX *m_L);
int       MRIquasiNewtonAlignVolumes(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MATRIX *m_L);
int       MRIemAlign(MRI *mri_in, GCA *gca, MORPH_PARMS *parms, MATRIX *m_L);
MATRIX    *MRIpowellAlignImages(MRI *mri_in, MRI *mri_target, MATRIX *m_L, 
                                float *pscale_factor, MATRIX *m_constraint, 
                                MRI *mri_source_mask, MRI *mri_target_mask,
                                int map_both_ways) ;
MATRIX *
MRIpowellAlignLabels(MRI *mri_source, MRI *mri_target, MATRIX *m_L);

int       MRIinitTranslation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) ;
int       MRIinitScaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) ;
int       MRIfindMeans(MRI *mri, float *means) ;
int       MRIfindCenterOfBrain(MRI *mri, float *px0, float *py0, float *pz0) ;
MRI       *MRIapply3DMorph(MRI *mri_in, MORPH_3D *m3d, MRI *mri_morphed) ;
MRI       *MRIapplyInverse3DMorph(MRI *mri_ref,MORPH_3D *m3d,MRI *mri_morphed);
MORPH_3D  *MRI3Dmorph(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms) ;
int       MRI3DmorphFree(MORPH_3D **pm3d) ;
int       MRI3DmorphFreeProperties(MORPH_3D *m3d) ;
int       MRI3Dwrite(MORPH_3D *m3d, char *fname) ;
MORPH_3D  *MRI3Dread(char *fname) ;
MORPH_3D  *MRI3DreadSmall(char *fname) ;
int       MRIsample3DmorphOrig(MORPH_3D *m3d, float x, float y, float z,
                               float *pxd, float *pyd, float *pzd);

MRI_SURFACE   *MRISshrinkWrapSkull(MRI *mri, MORPH_PARMS *parms) ;
int           MRIeraseNeck(MRI *mri, NECK_PARMS *np) ;

MATRIX *
MRIfaridAlignImages(MRI *mri_source, MRI *mri_target, MATRIX *m_L);
double MRIcomputeOptimalLinearXform(
  MRI *mri_src,
  MRI *mri_dst,
  MATRIX *m_L,
  float min_angle, float max_angle,
  float min_scale, float max_scale,
  float min_trans, float max_trans,
  float angle_steps, float scale_steps,
  float trans_steps,
  int nreductions, char *base_name, int map_both_ways) ;


#define M3D_MAGIC  0xabcdef42

extern int IMAGE_SIZE ;

#endif
