#ifndef MRI_MORPH_H
#define MRI_MORPH_H

#include "histo.h"
#include "mri.h"
#include "matrix.h"
#include "transform.h"


typedef struct
{
  double     l_intensity ;   /* coefficient of intensity term */
  double     l_area ;        /* coefficient of area term */
  double     l_nlarea ;      /* coefficient of nonlinear area term */
  double     l_dist ;        /* coefficient of distance term */
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
} MORPH_PARMS, MP ;

#define INTEGRATION_TOL    1e-4  /*5e-5*/
#define NEIGHBORS          6

typedef struct
{
  float      x, y, z ;       /* current coordinates */
} MORPH_NODE ;

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
} MORPH_NODE_PROPERTIES, MNP ;

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
} MORPH_3D, M3D ;


#define LTA_TYPE       0
#define MORPH3D_TYPE   1



HISTOGRAM *MRIhorizontalHistogram(MRI *mri, int thresh_low, int thresh_hi) ;
HISTOGRAM *MRIhorizontalBoundingBoxHistogram(MRI *mri, int thresh) ;

int       MRIcountAboveThreshold(MRI *mri, int thresh) ;
MRI       *MRIlabel(MRI *mri_src, MRI *mri_dst, int *nlabels) ;
int       MRIlabelBoundingBoxes(MRI *mri_label,MRI_REGION *bboxes,int nlabels);
int       MRIeraseOtherLabels(MRI *mri_src, MRI *mri_dst, int label) ;
int       MRIeraseLabel(MRI *mri_src, MRI *mri_dst, int label) ;
int       MRIfindHorizontalLabelLimits(MRI *mri, int label, 
                                       int *xmins, int *xmaxs) ;
MRI       *MRIremoveNeck(MRI *mri_src, MRI *mri_dst, int thresh_low, 
                         int thresh_hi, MORPH_PARMS *parms, int dir) ;
int       MRIlabelAreas(MRI *mri_label, float *areas, int nlabels) ;
int       MRIlabelCentroid(MRI *mri_label,int l,float *px,float *py,float *pz);
int       MRIlinearAlign(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms);
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


#define M3D_MAGIC  0xabcdef42

#endif
