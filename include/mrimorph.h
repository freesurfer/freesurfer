#ifndef MRI_MORPH_H
#define MRI_MORPH_H

#include "histo.h"
#include "mri.h"
#include "matrix.h"

typedef struct
{
  double     l_intensity ;   /* coefficient of intensity term */
  double     l_area ;        /* coefficient of area term */
  double     l_narea ;       /* coefficient of negative area term */
  double     l_dist ;        /* coefficient of distance term */
  MATRIX     *m_L ;          /* 4x4 linear homogenous tranform */
  double     dt ;
  double     momentum ;
  int        niterations ;
  int        write_iterations ;
  char       base_name[100] ;
  FILE       *log_fp ;
  int        start_t ;
  float      sigma ;         /* upper bound on amount of blurring */
  int        levels ;
  MRI        *mri_in, *mri_ref ;  /* for diagnostics to get at originals */
} MORPH_PARMS, MP ;

#define NEIGHBORS  6
typedef struct
{
  float      ox, oy, oz ;    /* original coordinates */
  float      x, y, z ;       /* current coordinates */
  float      orig_dist[NEIGHBORS] ;  /* original distances to 6 neighbors */
  float      orig_area ;             /* original area */
  float      area ;                  /* current area */
  float      dx, dy, dz ;            /* movement delta */
} MORPH_NODE ;

typedef struct
{
  float      node_spacing ;  /* mm between nodes (isotropic) */
  int        width ;         /* number of nodes wide */
  int        height ;
  int        depth ;
  MRI        *mri_in ;       /* source image */
  MRI        *mri_ref ;      /* target of morph */
  MATRIX     *m_L ;          /* linear trasform */
  MORPH_NODE ***nodes ;
  int        neg ;
} MORPH_3D, M3D ;

HISTOGRAM *MRIhorizontalHistogram(MRI *mri, int thresh_low, int thresh_hi) ;
HISTOGRAM *MRIhorizontalBoundingBoxHistogram(MRI *mri, int thresh) ;

int       MRIcountAboveThreshold(MRI *mri, int thresh) ;
MRI       *MRIlabel(MRI *mri_src, MRI *mri_dst, int *nlabels) ;
int       MRIlabelBoundingBoxes(MRI *mri_label,MRI_REGION *bboxes,int nlabels);
int       MRIeraseOtherLabels(MRI *mri_src, MRI *mri_dst, int label) ;
int       MRIfindHorizontalLabelLimits(MRI *mri, int label, 
                                       int *xmins, int *xmaxs) ;
MRI       *MRIremoveNeck(MRI *mri_src, MRI *mri_dst, int thresh_low, 
                         int thresh_hi, MORPH_PARMS *parms, int dir) ;
int       MRIlabelAreas(MRI *mri_label, float *areas, int nlabels) ;
int       MRIlinearAlign(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms);
int       MRIfindMeans(MRI *mri, float *means) ;
int       MRIfindCenterOfBrain(MRI *mri, float *px0, float *py0, float *pz0) ;
MRI       *MRIapply3DMorph(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d, 
                           MRI *mri_morphed) ;
MORPH_3D  *MRI3Dmorph(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms) ;
int       MRI3DmorphFree(MORPH_3D **pm3d) ;
int       MRI3Dwrite(MORPH_3D *m3d, char *fname) ;
MORPH_3D  *MRI3Dread(char *fname) ;
int      MRIsample3DmorphOrig(MORPH_3D *m3d, float x, float y, float z, 
                                 float *pxd, float *pyd, float *pzd);

#endif
