#ifndef GCA_MORPH_H
#define GCA_MORPH_H

#include "mri.h"
#include "gca.h"
#include "transform.h"


#define EXP_K            20.0

typedef struct
{
  Real   origx ;      //  (origx,origy,origz) saved
  Real   origy ;
  Real   origz ;
  Real   x ;          //  (x,y,z) saved
  Real   y ;
  Real   z ;
  Real   xs ;         //  not saved
  Real   ys ;
  Real   zs ;
  int    xn ;         /* node coordinates */
  int    yn ;         // (xn, yn, zn) saved
  int    zn ;
  int    label ;      
  int    n ;          /* index in gcan structure */
  float  prior ;
  GC1D   *gc ;
  float  log_p ;         /* current log probability of this sample */
  float  dx, dy, dz;     /* current gradient */
  float  odx, ody, odz ; /* previous gradient */
  float  area ;
  float  orig_area ;
  int    status ;       /* ignore likelihood term */
} GCA_MORPH_NODE, GMN ;

typedef struct
{
  int  width, height ,depth ;
  GCA  *gca ;          // using a separate GCA data (not saved)
  GMN  ***nodes ;
  int  neg ;
  double exp_k ;
  int  spacing ;
  MRI  *mri_xind ;    /* MRI->gca transform */
  MRI  *mri_yind ;
  MRI  *mri_zind ;
  VOL_GEOM   src;            /* src for the transform       */
  VOL_GEOM   dst;            /* dst for the transform       */
} GCA_MORPH, GCAM ;

typedef struct
{
  int    write_iterations ;
  float  dt ;
  float  momentum ;
  int    niterations ;
  char   base_name[STRLEN] ;
  double l_log_likelihood ;
  double l_likelihood ;
  double l_area ;
  double l_jacobian ;
  double l_smoothness ;
  double l_distance ;
  double l_label ;
  double l_map ;
  double tol ;
  int    levels ;
  FILE   *log_fp ;
  int    start_t ;
  MRI    *mri ;
  float  max_grad ;
  double exp_k ;
  double sigma ;
  int    navgs ;
  double label_dist ;
  int    noneg ;
  double ratio_thresh ;
  int    integration_type ;
  int    nsmall ;
  int    relabel ;    /* are we relabeling (i.e. using MAP label, or just max prior label) */
  int    relabel_avgs ; /* what level to start relabeling at */
  int    reset_avgs ;   /* what level to reset metric properties at */
} GCA_MORPH_PARMS, GMP ;

GCA_MORPH *GCAMalloc(int width, int height, int depth) ;
int       GCAMinit(GCA_MORPH *gcam, MRI *mri, GCA *gca, TRANSFORM *transform, int relabel) ;
int       GCAMinitLookupTables(GCA_MORPH *gcam) ;
int       GCAMwrite(GCA_MORPH *gcam, char *fname) ;
GCA_MORPH *GCAMread(char *fname) ;
int       GCAMfree(GCA_MORPH **pgcam) ;
MRI       *GCAMmorphFromAtlas(MRI *mri_src, GCA_MORPH *gcam, MRI *mri_dst) ;
MRI       *GCAMmorphToAtlas(MRI *mri_src, GCA_MORPH *gcam, MRI *mri_dst) ;
int       GCAMregister(GCA_MORPH *gcam, MRI *mri, GCA_MORPH_PARMS *parms) ;
int       GCAMregisterLevel(GCA_MORPH *gcam, MRI *mri, MRI *mri_smooth, 
                            GCA_MORPH_PARMS *parms) ;
int       GCAMsampleMorph(GCA_MORPH *gcam, float x, float y, float z, 
                          float *pxd, float *pyd, float *pzd) ;
int       GCAMcomputeLabels(MRI *mri, GCA_MORPH *gcam) ;
MRI       *GCAMbuildMostLikelyVolume(GCA_MORPH *gcam, MRI *mri) ;
MRI       *GCAMbuildVolume(GCA_MORPH *gcam, MRI *mri) ;
int       GCAMinvert(GCA_MORPH *gcam, MRI *mri) ;
int       GCAMfreeInverse(GCA_MORPH *gcam) ;
int       GCAMcomputeMaxPriorLabels(GCA_MORPH *gcam) ;
int       GCAMcomputeOriginalProperties(GCA_MORPH *gcam) ;
int       GCAMstoreMetricProperties(GCA_MORPH *gcam) ;
int       GCAMcopyNodePositions(GCA_MORPH *gcam, int from, int to) ;

#define GCAM_IGNORE_LIKELIHOOD 0x0001
#define GCAM_USE_LIKELIHOOD    0x0000
#define GCAM_LABEL_NODE        0x0002

int GCAMsetLabelStatus(GCA_MORPH *gcam, int label, int status) ;
int GCAMsetStatus(GCA_MORPH *gcam, int status) ;

#define ORIGINAL_POSITIONS  0
#define ORIG_POSITIONS      ORIGINAL_POSITIONS
#define SAVED_POSITIONS     1
#define CURRENT_POSITIONS   2

#define GCAM_INTEGRATE_OPTIMAL 0
#define GCAM_INTEGRATE_FIXED   1
#define GCAM_INTEGRATE_BOTH    2


#endif
