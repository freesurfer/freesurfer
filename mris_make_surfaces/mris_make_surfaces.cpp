/**
 * @FILE  mris_make_surfaces.c
 * @brief creates white matter and grey matter surface files.
 *
 * This program positions the tessellation of the cortical surface
 * at the white matter surface, then the gray matter surface
 * and generate surface files for these surfaces as well as a
 * 'curvature' file for the cortical thickness, and a surface file
 * which approximates layer IV of the cortical sheet.
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "mri.h"
#include "cma.h"
#include "macros.h"
#include "mrimorph.h"
#include "tags.h"
#include "mrinorm.h"
#include "version.h"
#include "label.h"
#include "voxlist.h"
#include "fsinit.h"
#include "mris_multimodal_refinement.h"
#include "surfgrad.h"
#include "utils.h"
#include "mrisurf_compute_dxyz.h"

#define  MAX_HISTO_BINS 1000


int main(int argc, char *argv[]) ;

#define MIN_NONCORTEX_VERTICES 10
//#define BRIGHT_LABEL         130
//#define BRIGHT_BORDER_LABEL  100

#define MIN_PEAK_PCT 0.1
static double Ghisto_left_outside_peak_pct = MIN_PEAK_PCT ;
static double Ghisto_right_outside_peak_pct = MIN_PEAK_PCT ;
static double Ghisto_left_inside_peak_pct = MIN_PEAK_PCT ;
static double Ghisto_right_inside_peak_pct = MIN_PEAK_PCT ;
static int contrast_type = CONTRAST_T2 ;
static int unpinch = 0 ;
static int find_and_mark_pinched_regions(MRI_SURFACE *mris,
    MRI *mri_T2,
    float nstd_below,
    float nstd_above) ;

static int T2_min_set = 0 ;
static int T2_max_set = 0 ;
static int T2_max_inside_set = 0 ;
static float T2_max_inside;
static float T2_min_inside;
static double T2_min_outside;
static double T2_max_outside;
static double max_outward_dist = 1 ; // for normal cases, for pathology might be much larger

// 0 disables use of max_gray in setting outside_hi in white surfacee deformation
// the larger this is, the closer outside_hi will be to max_gray. Set it high if
// the white surface is settling too far into the white matter
static double max_gray_scale = 0.0 ;  
static int above_set = 0 ;
static int below_set = 0 ;

static int nowmsa = 0 ;

static double l_repulse = 0.0 ;

int MRISremoveSelfIntersections(MRI_SURFACE *mris) ;
//static int  near_bright_structure(MRI *mri_aseg, int xi, int yi, int zi, int whalf) ;
// compute_pial_target_locations() is now MRIScomputePialTargetLocationsMultiModal() in mrisurf_mri.cpp

static int compute_white_target_locations(MRI_SURFACE *mris,
					  MRI *mri_T2,
					  MRI *mri_aseg,
					  float nstd_below,
					  float nstd_above,
					  int contrast_type,
					  float T2_min_inside,
					  float T2_max_inside,
					  int below_set,
					  int above_set,
					  double sigma,
					  double max_out) ;
static int compute_label_normal(MRI *mri_aseg, int x0, int y0, int z0,
                                int label, int whalf,
                                double *pnx, double *pny,
                                double *pnz, int use_abs);

static int edit_aseg_with_surfaces(MRI_SURFACE *mris, MRI *mri_aseg) ;
static double  compute_brain_thresh(MRI_SURFACE *mris,
                                    MRI *mri_ratio,
                                    float nstd) ;
static int fix_midline(MRI_SURFACE *mris,
                       MRI *mri_aseg,
                       MRI *mri_brain,
                       char *hemi, int which, int fix_mtl) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static MRI  *MRIsmoothMasking(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                              int mask_val, int wsize) ;
MRI *MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                      int out_label, MRI *mri_dst);

//int MRISfindExpansionRegions(MRI_SURFACE *mris) ;
MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;

static int fix_mtl = 0 ;
static char *read_pinch_fname = NULL ;
static LABEL *highres_label = NULL ;
static char T1_name[STRLEN] = "brain" ;

static int erase_cerebellum = 0 ;
static int erase_brainstem = 0 ;
static float nsigma = 2.0 ;
static float pial_blur_sigma = 0.0 ;
static float dura_thresh = -1 ;
static float nsigma_above = 6.0 ;
static double wm_weight = 3 ;    // threshold for finding gm outliers: thresh = (wt*wm + gm )/(wt+1)
static float nsigma_below = 6.0 ;
static int remove_contra = 1 ;
static char *write_aseg_fname = NULL ;
static char *white_fname = NULL ;
static int use_mode = 1 ;

const char *Progname ;

static double std_scale = 1.0;

static int label_cortex = 1 ;
static int graymid = 0 ;
static int curvature_avgs = 10 ;
static int create = 1 ;
static int smoothwm = 0 ;
static int overlay = 0 ;
static int inverted_contrast = 0 ;
static const char *filled_name = "filled" ;
static const char *wm_name = "wm" ;
static int auto_detect_stats = 1 ;
static char *dura_echo_name = NULL ;
static char *T2_name = NULL ;
static char *flair_or_T2_name = NULL ;
static int nechos = 0 ;
static bool followGradients = false ;

static int in_out_in_flag = 0 ;  /* for Arthur (as are most things) */

static int nbhd_size = 20 ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;
static float pial_target_offset = 0 ;
static float white_target_offset = 0 ;

static COLOR_TABLE *ctab  = NULL;
static MRI *mri_cover_seg = NULL ;
static const char *aseg_name = "aseg" ;
static const char *aparc_name = "aparc" ;  // for midline and cortex label
static MRI *mri_aseg = NULL;
static int add = 0 ;

static double l_tsmooth = 0.0 ;
static double l_surf_repulse = 5.0 ;

static int smooth_pial = 0 ;
static int smooth = 5 ;
static int vavgs = 5 ;
static int nwhite = 100 /* 25 */ /*5*/ ;
static int ngray = 100  /* 30*/ /*45*/ ;

// orig_name is the surface that will be read first. It will be
// smoothed. If orig_white is specified, the smoothed vertex locations
// will be overwritten by the orig_white surface (which will not be
// smoothed).
static const char *orig_name = ORIG_NAME ; // "orig"

static char *AutoDetSurfName=NULL; // surf to use when autodet stats

// white_only=1, then only place the white, do not place pial
static int white_only = 0 ;

// If nowhite=1, do not place white, forces reading previously placed
// gray/white surface.  This will be the value of orig_white or
// "white" if unspecified. Could have been called pialonly instead.
static int nowhite = 0 ;

// orig_white, if specified, will be used as the initial surface 
// to place the white surface. Whatever surface is specified, it
// will not be smoothed prior to placing the white surface. 
static char *orig_white = NULL ;
static char *orig_sphere = NULL ;

// orig_pial, if specified, will be used as the initial surface 
// to place the white surface
static char *orig_pial = NULL ;

int SaveTarget=0;
int SaveResidual=0;

static int flairwhite = 0 ;
static int nbrs = 2 ;
static int write_vals = 0 ;

static const char *suffix = "" ;
static const char *output_suffix = "" ;
static char *xform_fname = NULL ;
char *ripfilename=NULL, *whitesigmafilename=NULL, *pialsigmafilename=NULL;

static char pial_name[STRLEN] = "pial" ;
static char white_matter_name[STRLEN] = WHITE_MATTER_NAME ; // "white"

static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;

static int max_pial_averages = 16 ;
static int min_pial_averages = 2 ;
static int max_white_averages = 4 ;
static int min_white_averages = 0 ;
static float pial_sigma = 2.0f ;
static float white_sigma = 2.0f ;
static float max_thickness = 5.0 ;

static float variablesigma = 3.0;

#define MAX_WHITE             120
#define MAX_BORDER_WHITE      105
#define MIN_BORDER_WHITE       85
#define MIN_GRAY_AT_WHITE_BORDER  70

#define MAX_GRAY               95
#define MIN_GRAY_AT_CSF_BORDER    40
#define MID_GRAY               ((max_gray + min_gray_at_csf_border) / 2)
#define MAX_GRAY_AT_CSF_BORDER    75
#define MIN_CSF                10
#undef MAX_CSF
#define MAX_CSF                40

static  int   max_border_white_set = 0,
              min_border_white_set = 0,
              min_gray_at_white_border_set = 0,
              max_gray_set = 0,
              max_gray_at_csf_border_set = 0,
              min_gray_at_csf_border_set = 0,
              min_csf_set = 0,
              max_csf_set = 0 ;

static  float   max_border_white = MAX_BORDER_WHITE,
                min_border_white = MIN_BORDER_WHITE,
                min_gray_at_white_border = MIN_GRAY_AT_WHITE_BORDER,
                max_gray = MAX_GRAY,
                max_gray_at_csf_border = MAX_GRAY_AT_CSF_BORDER,
                min_gray_at_csf_border = MIN_GRAY_AT_CSF_BORDER,
                min_csf = MIN_CSF,
                max_csf = MAX_CSF ;
static char sdir[STRLEN] = "" ;

static int MGZ = 1; // for use with MGZ format

static int longitudinal = 0;

static int pial_num = 0 ;
static int pial_nbrs = 0 ;
static int white_num = 0 ;
#define MAX_VERTICES 1000
static float pial_vals[MAX_VERTICES] ;
static int pial_vnos[MAX_VERTICES] ;
static float white_vals[MAX_VERTICES] ;
static int white_vnos[MAX_VERTICES] ;

static int fill_interior = 0 ;
static float check_contrast_direction(MRI_SURFACE *mris,MRI *mri_T1) ;
int UseV6CBV = 0;

int main(int argc, char *argv[])
{
  char          *hemi, *sname, *cp, fname[STRLEN], mdir[STRLEN];
 int           nargs, i, msec, n_averages, j ;
  MRI_SURFACE   *mris, *mrisAutoDet ;
  MRI           *mri_wm, *mri_kernel = NULL;
  MRI *mri_smooth = NULL, *mri_mask = NULL;
  MRI *mri_filled, *mri_T1, *mri_labeled, *mri_T1_white = NULL, *mri_T1_pial ;
  float         max_len ;
  float         white_mean, white_std, gray_mean, gray_std ;
  double        current_sigma, thresh = 0, spring_scale = 1;
  Timer then ;
  M3D           *m3d ;
  INTEGRATION_PARMS old_parms ;
  int memusage[5];
  char *cmdline2, cwd[2000];
  MRIS *mristarget = NULL;
  int vno;
  VERTEX *vgdiag;

  FSinit() ;
  std::string cmdline = getAllInfo(argc, argv, "mris_make_surfaces");

  nargs = handleVersionOption(argc, argv, "mris_make_surfaces");
  getcwd(cwd,2000);
  cmdline2 = argv2cmdline(argc,argv);

  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  // don't let gradient use exterior information (slows things down)
  parms.fill_interior = 0 ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;

  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.2 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 0.5f ;
  parms.l_spring = 0.0f ;
  parms.l_surf_repulse = 0.0 ;
  parms.l_repulse = 5.0 ;
  parms.l_spring_nzr = 0.0 ;
  parms.l_spring_nzr_len = 0.0 ;
  parms.l_hinge = 0;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit() ;
  }

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = nwhite ;
  if (parms.momentum < 0.0)
  {
    parms.momentum = 0.0 /*0.75*/ ;
  }

  then.reset() ;
  sname = argv[1] ;
  hemi = argv[2] ;
  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment.\n", Progname) ;
  strcpy(mdir, cp) ;

  // print out version of this program and mrisurf.c
  printf("%s\n",getVersion().c_str());
  printf("%s\n",getVersion().c_str());
  printf("\n");
  printf("cd %s\n",cwd);
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("%s\n",cmdline2);
  printf("\n");

  fflush(stdout);
  int req = snprintf(fname, STRLEN, "%s/%s/surf/mris_make_surfaces.%s.mrisurf.c.version", 
		     sdir, sname, hemi) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, filled_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ)
  {
    strcat(fname, ".mgz");
  }
  fprintf(stdout, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  ////////////////////////////// we can handle only conformed volumes
//  setMRIforSurface(mri_filled);

  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, T1_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ)
  {
    strcat(fname, ".mgz");
  }
  fprintf(stdout, "reading volume %s...\n", fname) ;
  mri_T1 = mri_T1_pial = MRIread(fname) ;

  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  /////////////////////////////////////////
//  setMRIforSurface(mri_T1);


  // make sure filled and wm are in same voxel coords as T1
  if (MRIcompareHeaders(mri_T1, mri_filled) != 0)
  {
    MRI *mri_tmp ;
    printf("resampling filled volume to be in voxel register with T1\n") ;
    mri_tmp = MRIresample(mri_filled, mri_T1, SAMPLE_NEAREST) ;
    MRIfree(&mri_filled) ; mri_filled = mri_tmp ;
  }
  if (white_fname != NULL)
  {
    int req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, white_fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (MGZ) strcat(fname, ".mgz");
    fprintf(stdout, "reading volume %s...\n", fname) ;
    mri_T1_white = MRIread(fname) ;
    if (!mri_T1_white)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                Progname, fname) ;
    /////////////////////////////////////////
    if (mri_T1_white->type != MRI_UCHAR)
    {
      MRI *mri_tmp ;

      MRIeraseNegative(mri_T1_white, mri_T1_white) ;
      mri_tmp = MRIchangeType(mri_T1_white, MRI_UCHAR, 0, 255, 1) ;
      MRIfree(&mri_T1_white) ;
      mri_T1_white = mri_tmp ; // swap
    }

//    setMRIforSurface(mri_T1_white);
  }

  if (xform_fname)
  {
    char fname[STRLEN], ventricle_fname[STRLEN] ;
    MRI  *mri_lv, *mri_inv_lv ;

    int req = snprintf(fname, STRLEN, "%s/%s/mri/transforms/%s", 
		       sdir, sname,xform_fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "reading transform %s...\n", fname) ;
    m3d = MRI3DreadSmall(fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                Progname, fname) ;
    req = snprintf(ventricle_fname, STRLEN, "%s/average/%s_ventricle.mgz#0@mgh",
		   mdir, !stricmp(hemi, "lh") ? "left" : "right") ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout,"reading ventricle representation %s...\n",
            ventricle_fname);
    mri_lv = MRIread(ventricle_fname) ;
    if (!mri_lv)
      ErrorExit(ERROR_NOFILE,"%s: could not read %s",
                Progname,ventricle_fname);
    fprintf(stdout, "applying inverse morph to ventricle...\n") ;
    mri_inv_lv = MRIapplyInverse3DMorph(mri_lv,m3d,NULL);
    MRI3DmorphFree(&m3d) ;
    MRIfree(&mri_lv) ;
    fprintf(stdout, "filling in ventricle...\n") ;
    mri_lv = MRIfillVentricle(mri_inv_lv,mri_T1,100,
                              DEFAULT_DESIRED_WHITE_MATTER_VALUE, NULL);
    MRIfree(&mri_inv_lv) ;
    MRIunion(mri_lv, mri_T1, mri_T1) ;
    MRIfree(&mri_lv) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      int req = snprintf(fname, STRLEN, "%s/%s/mri/T1_filled", sdir, sname) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      MRIwrite(mri_T1, fname) ;
    }
  }
  if (aseg_name)
  {
    char fname[STRLEN] ;
    int req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, aseg_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (MGZ)
    {
      strcat(fname, ".mgz");
    }

    fprintf(stdout, "reading volume %s...\n", fname) ;
    mri_aseg = MRIread(fname) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation volume %s",
                Progname, fname) ;
    if (nowmsa)
    {
      MRIreplaceValuesOnly(mri_aseg, mri_aseg, Left_WM_hypointensities, Left_Cerebral_White_Matter) ;
      MRIreplaceValuesOnly(mri_aseg, mri_aseg, Right_WM_hypointensities, Right_Cerebral_White_Matter) ;
      MRIreplaceValuesOnly(mri_aseg, mri_aseg, WM_hypointensities, Left_Cerebral_White_Matter) ;
    }
    if (erase_brainstem)
    {
      MRImaskLabel(mri_T1, mri_T1, mri_aseg, Brain_Stem, 0) ;
      MRImaskLabel(mri_T1, mri_T1, mri_aseg,  Left_VentralDC, 0) ;
      MRImaskLabel(mri_T1, mri_T1, mri_aseg,  Right_VentralDC, 0) ;
      if (mri_T1_white && mri_T1_white != mri_T1)
      {
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg, Brain_Stem, 0) ;
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg,  Left_VentralDC, 0) ;
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg,  Right_VentralDC, 0) ;
      }

      if (mri_T1_pial && mri_T1_pial != mri_T1)
      {
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg, Brain_Stem, 0) ;
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg,  Left_VentralDC, 0) ;
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg,  Right_VentralDC, 0) ;
      }

      MRImaskLabel(mri_filled, mri_filled, mri_aseg, Brain_Stem, 0) ;
      MRImaskLabel(mri_filled, mri_filled, mri_aseg,  Left_VentralDC, 0) ;
      MRImaskLabel(mri_filled, mri_filled, mri_aseg,  Right_VentralDC, 0) ;
    }
    if (erase_cerebellum)
    {
      MRImaskLabel(mri_T1, mri_T1, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
      MRImaskLabel(mri_T1, mri_T1, mri_aseg, Left_Cerebellum_Cortex, 0) ;
      MRImaskLabel(mri_T1, mri_T1, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
      MRImaskLabel(mri_T1, mri_T1, mri_aseg, Right_Cerebellum_Cortex, 0) ;

      if (mri_T1_white && mri_T1_white != mri_T1)
      {
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg, Left_Cerebellum_Cortex, 0) ;
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
	MRImaskLabel(mri_T1_white, mri_T1_white, mri_aseg, Right_Cerebellum_Cortex, 0) ;
      }

      if (mri_T1_pial && mri_T1_pial != mri_T1)
      {
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg, Left_Cerebellum_Cortex, 0) ;
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
	MRImaskLabel(mri_T1_pial, mri_T1_pial, mri_aseg, Right_Cerebellum_Cortex, 0) ;
      }

      MRImaskLabel(mri_filled, mri_filled, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
      MRImaskLabel(mri_filled, mri_filled, mri_aseg, Left_Cerebellum_Cortex, 0) ;
      MRImaskLabel(mri_filled, mri_filled, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
      MRImaskLabel(mri_filled, mri_filled, mri_aseg, Right_Cerebellum_Cortex, 0) ;
    }
  }
  else
  {
    mri_aseg = NULL ;
  }

  // Load in wm.mgz (or equivalent)
  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, wm_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stdout, "reading volume %s...\n", fname) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  if (mri_wm->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;
    printf("changing type of input wm volume to UCHAR...\n") ;
    mri_tmp = MRIchangeType(mri_wm, MRI_UCHAR, 0, 255, 1) ;
    MRIfree(&mri_wm) ;
    mri_wm = mri_tmp ;
  }
  if (MRIcompareHeaders(mri_T1, mri_wm) != 0)
  {
    MRI *mri_tmp ;
    printf("resampling wm volume to be in voxel register with T1\n") ;
    mri_tmp = MRIresample(mri_wm, mri_T1, SAMPLE_NEAREST) ;
    MRIfree(&mri_wm) ; mri_wm = mri_tmp ;
  }
  //////////////////////////////////////////
  //  setMRIforSurface(mri_wm);

  // This does not smooth. It clips the maximum WM value
  MRIclipBrightWM(mri_T1, mri_wm) ;

  if(fill_interior == 0)
    mri_labeled = MRIfindBrightNonWM(mri_T1, mri_wm) ;
  else
    mri_labeled = MRIclone(mri_T1, NULL) ;

  if(mri_T1_white)
    MRIclipBrightWM(mri_T1_white, mri_wm) ;


  if (overlay)
  {
    fprintf(stdout, "overlaying editing into T1 volume...\n") ;
    MRImask(mri_T1, mri_wm, mri_T1,
            WM_EDITED_ON_VAL, DEFAULT_DESIRED_WHITE_MATTER_VALUE);
    MRIsmoothMasking(mri_T1, mri_wm, mri_T1, WM_EDITED_ON_VAL, 15) ;
    int req = snprintf(fname, STRLEN, "%s/%s/mri/T1_overlay", sdir, sname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(mri_T1, fname) ;
  }

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, orig_name, suffix) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  printf("reading original surface position from %s...\n", fname) ;
  if(orig_white) printf("  .. but with overwrite the positions with %s...\n",orig_white) ;
  mris = MRISreadOverAlloc(fname, 1.1) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;

  if(AutoDetSurfName){
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s", 
		       sdir, sname, hemi, AutoDetSurfName, suffix) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("reading surface for AutoDet from %s...\n", fname) ;
    mrisAutoDet = MRISreadOverAlloc(fname, 1.1) ;
    if (!mrisAutoDet)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
  }
  else mrisAutoDet = mris;

  MRISfaceMetric(mris,0); // just computes f->norm, does not replace MRISmetricProperties
  MRISedges(mris);
  MRIScorners(mris);
  MRISedgeMetric(mris,0);
  MRIScornerMetric(mris,0);
  printf("Input Surface Quality Stats\n");
  MRISprettyPrintSurfQualityStats(stdout, mris);

  if(parms.l_hinge > 0 || parms.l_spring_nzr > 0){
    if(parms.l_spring_nzr){
      double  *edgestats = MRISedgeStats(mris, 0, NULL, NULL);
      parms.l_spring_nzr_len = edgestats[1];
      free(edgestats);
    }
  }

  if (mris->vg.valid && !FZERO(mris->vg.xsize))
    spring_scale = 3/(mris->vg.xsize + mris->vg.ysize + mris->vg.zsize) ;
  else
  {
    spring_scale = 3/(mri_T1->xsize + mri_T1->ysize + mri_T1->zsize) ;
    mris->vg.xsize = mri_T1->xsize ;
    mris->vg.ysize = mri_T1->ysize ;
    mris->vg.zsize = mri_T1->zsize ;
  }

  if (!FEQUAL(spring_scale, 1.0))
    printf("spring scale = %2.2f\n", spring_scale) ;

  MRISaddCommandLine(mris, cmdline) ;

  if(mris->ct == NULL && ctab)  
    mris->ct = ctab ;  // add user-specified color table to structure

  MRISremoveIntersections(mris,0) ;
  
  if (pial_nbrs > 2)
    MRISsetNeighborhoodSizeAndDist(mris, pial_nbrs) ;
  
  if (auto_detect_stats)
  {
    printf("Auto detecting stats\n");
    MRI *mri_tmp ;
    float white_mode, gray_mode ;

    // Binarize wm.mgz by thresholding at WM_MIN_VAL. Voxels below threshold will 
    // take a value of MRI_NOT_WHITE; those above will get MRI_WHITE.
    printf("Binarizing %s thresholding at %d\n",wm_name,WM_MIN_VAL);
    mri_tmp = MRIbinarize(mri_wm, NULL, WM_MIN_VAL, MRI_NOT_WHITE, MRI_WHITE) ;
    MRISsaveVertexPositions(mrisAutoDet, WHITE_VERTICES) ;
    // WHITE_MATTER_MEAN = 110
    printf("computing class statistics... low=30, hi=%d\n",WHITE_MATTER_MEAN);
    // This computes means and stddevs of voxels near the border of
    // wm.mgz with inside being WM and outside being GM. Seems like
    // the aseg would be better for this than the wm.mgz
    MRIcomputeClassStatistics(mri_T1, mri_tmp, 30, WHITE_MATTER_MEAN,
                              &white_mean, &white_std, &gray_mean, &gray_std) ;
    printf("white_mean = %g +/- %g, gray_mean = %g +/- %g\n",white_mean, white_std, gray_mean,gray_std) ;

    if(use_mode){
      printf("using class modes intead of means, discounting robust sigmas....\n") ;
      //MRIScomputeClassModes(mris, mri_T1, &white_mode, &gray_mode, NULL, &white_std, &gray_std, NULL);
      // This gets stats based on sampling the MRI at 1mm inside (WM) and 1mm outside (GM) of the surface.
      // This makes the identity of mris very important! It will be orig_name by default but will
      // become white_name if white_name specified.
      MRIScomputeClassModes(mrisAutoDet, mri_T1, &white_mode, &gray_mode, NULL, NULL, NULL, NULL);
      white_mean = white_mode ;
      gray_mean = gray_mode ;
      printf("white_mode = %g, gray_mode = %g\n",white_mode, gray_mode);
    }
    printf("std_scale = %g\n",std_scale);

    white_std /= std_scale;
    gray_std /= std_scale;

    //these may be set on the cmd
    if(!min_gray_at_white_border_set)
      min_gray_at_white_border = gray_mean-gray_std ;
    if(!max_border_white_set)
      max_border_white = white_mean+white_std ;
    if(!max_csf_set)
      max_csf = gray_mean - MAX(.5, (variablesigma-1.0))*gray_std ;
    if (!min_border_white_set)
      min_border_white = gray_mean ;

    // apply some sanity checks
    double MAX_SCALE_DOWN = .2;
    printf("Applying sanity checks, MAX_SCALE_DOWN = %g\n",MAX_SCALE_DOWN);

    if (min_gray_at_white_border < MAX_SCALE_DOWN*MIN_GRAY_AT_WHITE_BORDER)
      min_gray_at_white_border = nint(MAX_SCALE_DOWN*MIN_GRAY_AT_WHITE_BORDER) ;
    if (max_border_white < MAX_SCALE_DOWN*MAX_BORDER_WHITE)
      max_border_white = nint(MAX_SCALE_DOWN*MAX_BORDER_WHITE) ;
    if (min_border_white < MAX_SCALE_DOWN*MIN_BORDER_WHITE)
      min_border_white = nint(MAX_SCALE_DOWN*MIN_BORDER_WHITE) ;
    if (max_csf < MAX_SCALE_DOWN*MAX_CSF)
      max_csf = MAX_SCALE_DOWN*MAX_CSF ;

    printf("setting MIN_GRAY_AT_WHITE_BORDER to %2.1f (was %d)\n",
            min_gray_at_white_border, MIN_GRAY_AT_WHITE_BORDER) ;
    printf("setting MAX_BORDER_WHITE to %2.1f (was %d)\n",
            max_border_white, MAX_BORDER_WHITE) ;
    printf("setting MIN_BORDER_WHITE to %2.1f (was %d)\n",
            min_border_white, MIN_BORDER_WHITE) ;
    printf("setting MAX_CSF to %2.1f (was %d)\n",
            max_csf, MAX_CSF) ;

    //these may be set on the cmd
    if (!max_gray_set)
      max_gray = white_mean-white_std ;
    if (!max_gray_at_csf_border_set){
      //max_gray_at_csf_border = gray_mean-0.5*gray_std ;
      max_gray_at_csf_border = gray_mean-1.0*gray_std ;   // changed to push pial surfaces further out BRF 12/10/2015
    }
    if (!min_gray_at_csf_border_set)
      min_gray_at_csf_border = gray_mean - variablesigma*gray_std ;

    if (max_gray < MAX_SCALE_DOWN*MAX_GRAY)
      max_gray = nint(MAX_SCALE_DOWN*MAX_GRAY) ;
    if (max_gray_at_csf_border < MAX_SCALE_DOWN*MAX_GRAY_AT_CSF_BORDER)
      max_gray_at_csf_border = nint(MAX_SCALE_DOWN*MAX_GRAY_AT_CSF_BORDER) ;
    if (min_gray_at_csf_border < MAX_SCALE_DOWN*MIN_GRAY_AT_CSF_BORDER)
      min_gray_at_csf_border = nint(MAX_SCALE_DOWN*MIN_GRAY_AT_CSF_BORDER) ;

    printf("setting MAX_GRAY to %2.1f (was %d)\n",max_gray, MAX_GRAY) ;
    printf("setting MAX_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",max_gray_at_csf_border, MAX_GRAY_AT_CSF_BORDER) ;
    printf("setting MIN_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",min_gray_at_csf_border, MIN_GRAY_AT_CSF_BORDER) ;
    MRIfree(&mri_tmp) ;
  }

  if (dura_echo_name == NULL)
  {
    MRIfree(&mri_wm) ;
  }
  inverted_contrast = (check_contrast_direction(mris,mri_T1) < 0) ;
  if (inverted_contrast)
  {
    printf("inverted contrast detected....\n") ;
  }
  if (highres_label)
  {
    LabelRipRestOfSurface(highres_label, mris) ;
  }
  if (smooth && !nowhite && !dura_echo_name)
  {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    if(orig_white) printf("  .. but with overwrite the positions with %s...\n",orig_white) ;
    /* The effect of smoothing the vertices here is to pull the
       surface away from the boundary toward the center of the
       hemisphere.  For the orig, it has a bigger effect at the crowns
       because that is where the surface bends the sharpest. The
       surface placement then must work to undo this. When 5
       iterations are used (smooth=5), the smoothed surface moves
       quite a bit.  This does not happen when -orig_white is
       specified. Basically, the smoothing is to handle the orig
       surface which may be jagged. */
    MRISaverageVertexPositions(mris, smooth) ;
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  sprintf(parms.base_name, "%s%s%s", white_matter_name, output_suffix, suffix) ;
  if(orig_white){
    printf("reading initial white vertex positions from %s...\n",orig_white) ;
    if (MRISreadVertexPositions(mris, orig_white) != NO_ERROR)
      ErrorExit(Gerror, "reading of orig white failed...");
    MRISremoveIntersections(mris,0) ;
  }
  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
  MRISstoreMetricProperties(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;

  if (add)  {
    fprintf(stdout, "adding vertices to initial tessellation...\n") ;
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris, max_len) > 0) {}
  }
  MRISsetVals(mris, -1) ;  /* clear white matter intensities */

  if (aparc_name) {
    printf("Reading in aparc %s\n",aparc_name);
    if (MRISreadAnnotation(mris, aparc_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation",aparc_name) ;
  }

  if(!nowhite)  {
    printf("repositioning cortical surface to gray/white boundary\n");
    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_BORDER_LABEL, 0) ;
    if (mri_T1_white) {
      MRImask(mri_T1_white, mri_labeled, mri_T1_white, BRIGHT_LABEL, 0) ;
      MRImask(mri_T1_white, mri_labeled, mri_T1_white, BRIGHT_BORDER_LABEL, 0) ;
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_T1, "white_masked.mgz") ;
  }

  if (mri_T1_white)  {
    if (mri_T1 != mri_T1_pial)
      MRIfree(&mri_T1);
    mri_T1 = mri_T1_white ; // T1 and T1_white is swapped
  }

  // ==========================================================================
  // Place white surface. This is a loop where the number of averages
  // descreases by a factor of two until it reaches min_white_averages (0).
  // The sigma also decreases by a factor of two for each iter.
  if(nowhite==0)
    printf("Placing white surface white_sigma=%g, max_white_averages=%d\n",white_sigma, max_white_averages);
  current_sigma = white_sigma ;
  n_averages = max_white_averages;

  for (i = 0 ;  n_averages >= min_white_averages ; n_averages /= 2, current_sigma /= 2, i++)
  {
    if(nowhite) break ; // skip if not placing the white surface

    printf("Iteration %d white =========================================\n",i);
    printf("n_averages=%d, current_sigma=%g\n",n_averages,current_sigma); fflush(stdout);

    // This does not look like it actually smooths anything. It creates the kernel and frees 
    // it without ever applying it. mri_smooth is just a copy of mri_T1. current sigma is used below
    printf("Smoothing T1 volume with sigma = %2.3f\n",current_sigma) ; fflush(stdout);
    parms.sigma = current_sigma ;
    mri_kernel = MRIgaussian1d(current_sigma, 100) ;
    if (!mri_smooth)
      mri_smooth = MRIcopy(mri_T1, NULL) ;
    MRIfree(&mri_kernel) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON){
      char fname[STRLEN] ;
      sprintf(fname, "sigma%.0f.mgz", current_sigma) ;
      fprintf(stdout, "writing smoothed volume to %s...\n", fname) ;
      MRIwrite(mri_smooth, fname) ;
    }

    parms.n_averages = n_averages ;

    MRISprintTessellationStats(mris, stdout) ;
    if(mri_aseg){
      printf("Freezing midline and others\n");  fflush(stdout);
      // This rips the midline vertices so that they do not move (thus
      // "fix"). It may also rip some vertices near the putamen and
      // maybe near lesions. It may set v->marked2 which will
      // influence the cortex.label. At onset, it will set all marked
      // and marked2 to 0. At then end it will set all marked=0.
      // It does not unrip any vertices, so, unless they are unripped
      // at some other point, the number of ripped vertices will 
      // increase.
      if(Gdiag_no > 0){
	vgdiag = &mris->vertices[Gdiag_no];
	printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
      }
      MRISripMidline(mris, mri_aseg, mri_T1, hemi, GRAY_WHITE, 0) ;
      if(Gdiag_no > 0){
	vgdiag = &mris->vertices[Gdiag_no];
	printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
      }
    }

    if (mri_cover_seg) {
      if(i == 0){
	MRI *mri_bin, *mri_tmp;
	int itmp=0;
	if(mri_T1 == mri_T1_pial) itmp=1;
	mri_bin = MRIclone(mri_T1, NULL) ;
	mri_tmp = MRIScoverSeg(mris, mri_bin, mri_cover_seg, GRAY_WHITE);
	MRIfree(&mri_bin) ; 
	MRIfree(&mri_T1) ;
	mri_T1 = mri_tmp ;
	if(itmp) mri_T1_pial = mri_tmp ;
      }
    }
    else if (flair_or_T2_name == NULL){ // otherwise already done
      // This is where most of the recon-all commands end up when placing white

      //outside_hi = max_border_white when max_gray_scale=0 (default)
      double outside_hi = (max_border_white + max_gray_scale*max_gray) / (max_gray_scale+1.0) ;

      if (!FZERO(max_gray_scale))
	printf("setting outside hi = %2.1f (was %2.1f) in white surface deformation\n", outside_hi, max_border_white) ;

      printf("Computing border values for white\n");
      printf("  MAX_WHITE  %d, max_border_white %g, min_border_white %g, min_gray_at_white_border %g\n",
	     MAX_WHITE, max_border_white, min_border_white, min_gray_at_white_border);
      printf("  outside_hi  %g, max_thickness %g, max_gray_scale %g,  max_gray %g\n",
	     outside_hi, max_thickness, max_gray_scale, max_gray); fflush(stdout);
      // MRIScomputeBorderValues() will effectively set all v->marked=1 for all unripped vertices
      if(UseV6CBV == 0){
	MRIScomputeBorderValues(mris, mri_T1, mri_smooth,
				MAX_WHITE, max_border_white, min_border_white,
				min_gray_at_white_border,
				outside_hi, current_sigma,
				2*max_thickness, parms.fp, GRAY_WHITE, NULL, 0, parms.flags,mri_aseg,-1,-1);
      } 
      else {
	MRIScomputeBorderValuesV6(mris, mri_T1, mri_smooth,
				MAX_WHITE, max_border_white, min_border_white,
				min_gray_at_white_border,
				outside_hi, current_sigma,
				2*max_thickness, parms.fp, GRAY_WHITE, NULL, 0, parms.flags,mri_aseg,-1,-1);
      }


      printf("Finding expansion regions\n"); fflush(stdout);
      MRISfindExpansionRegions(mris) ;
      if(Gdiag_no > 0){
	vgdiag = &mris->vertices[Gdiag_no];
	printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
      }
    }

    else if (flairwhite)  {
      MRI  *mri_flair = NULL ;
      char fname[STRLEN] ;
      printf("Code location: flairwhite\n");
      
      if (orig_pial) strcpy(fname, orig_pial) ;
      else           strcpy(fname, "pial") ;

      printf("reading initial pial vertex positions from %s...\n", fname) ;

      if (MRISreadPialCoordinates(mris, fname) != NO_ERROR)
	ErrorExit(Gerror, "reading orig pial positions failed") ;
      MRISremoveIntersections(mris,0) ;
      
      strcpy(fname, flair_or_T2_name) ;
      if (MGZ)
	strcat(fname, ".mgz");
      
      printf("repositioning white surface locations using  %s\n", fname) ;

      mri_flair = MRIread(fname) ;
      if (mri_flair == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load flair volume %s", Progname, fname) ;

      if (orig_white){
	if (i > 0)
	  MRISsaveVertexPositions(mris, WHITE_VERTICES) ; // update estimate of white
	else   // first time - read it in
	{
	  printf("reading initial white vertex positions from %s...\n",orig_white) ;
	  if (MRISreadVertexPositions(mris, orig_white) != NO_ERROR)
	    ErrorExit(Gerror, "reading of orig white failed...");
	  MRISremoveIntersections(mris,0) ;
	}
      }

      if (MRImatchDimensions(mri_flair, mri_aseg) == 0)  {
	MRI *mri_tmp ;
	
	printf("resampling aseg to be in voxel space of T2/FLAIR\n") ;
	mri_tmp = MRIresample(mri_aseg, mri_flair, SAMPLE_NEAREST) ;
	MRIfree(&mri_aseg) ;
	mri_aseg = mri_tmp ;
      }

      // Compute white locations when using FLAIR
      compute_white_target_locations(mris, mri_flair, mri_aseg, 3, 3, CONTRAST_FLAIR, 
				     T2_min_inside, T2_max_inside, below_set, above_set, white_sigma,
				     max_outward_dist) ;
      if (Gdiag & DIAG_WRITE){
	char fname[STRLEN] ;
	static int callno = 0 ;
	sprintf(fname, "%ch.dist.%3.3d.mgz", mris->hemisphere == LEFT_HEMISPHERE ? 'l' : 'r', callno) ;
	printf("writing distances to %s\n", fname) ;
	MRISwriteD(mris, fname) ;
	MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
	MRISrestoreVertexPositions(mris, TARGET_VERTICES) ;
	sprintf(fname, "%ch.targets.%3.3d", mris->hemisphere == LEFT_HEMISPHERE ? 'l' : 'r', callno) ;
	printf("writing target locations to %s\n", fname) ;
	MRISwrite(mris, fname) ;
	MRISrestoreVertexPositions(mris, TMP2_VERTICES) ;
	callno++ ;
      }
    } // end flairwhite

    if(vavgs) {
      fprintf(stdout,"averaging target values for %d iterations...\n",vavgs) ;
      // MRIScomputeBorderValues() sets v->marked=1 for all unripped
      MRISaverageMarkedVals(mris, vavgs) ;
      if (Gdiag_no > 0){
        VERTEX * const v = &mris->vertices[Gdiag_no] ;
        printf("v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f, ripflag=%d\n",
	       Gdiag_no, v->val, v->mean, v->d, v->ripflag) ;
      }
    }

    /*
      there are frequently regions of gray whose intensity is fairly
      flat. We want to make sure the surface settles at the innermost
      edge of this region, so on the first pass, set the target
      intensities artificially high so that the surface will move
      all the way to white matter before moving outwards to seek the
      border (I know it's a hack, but it improves the surface in
      a few areas. The alternative is to explicitly put a gradient-seeking
      term in the cost functional instead of just using one to find
      the target intensities).
    */
    if (write_vals) {
      sprintf(fname, "./%s-white%2.2f.mgz", hemi, current_sigma) ;
      MRISwriteValues(mris, fname);
    }

    // Create a surface of the target surface at the point where the
    // maximum gradient should be. This will not necessarily be a topo
    // correct surface, but it can be used for evaluation.  Have to
    // create this surface here because v->d gets set to 0 after pos surf
    if(mristarget != NULL) MRISfree(&mristarget);
    mristarget = MRISclone(mris);

    MRISfreeDistsButNotOrig(mristarget);

    for(vno=0; vno < mris->nvertices; vno++){
      VERTEX* const v = &(mristarget->vertices[vno]);
      v->d   = mris->vertices[vno].d; // clone does not copy this
      v->val = mris->vertices[vno].val; // clone does not copy this
      MRISsetXYZ(mristarget,vno,
        v->x + (v->d*v->nx), // d is the distance to the max gradient
        v->y + (v->d*v->ny),
        v->z + (v->d*v->nz));
    }

    //sprintf(fname, "./%s.white.target.i%02d", hemi, i);
    //printf("writing white target surface to %s...\n", fname) ;
    //MRISwrite(mristarget, fname);

    // Everything up to now has been leading up to this. This is where
    // the surfaces get placed.
    INTEGRATION_PARMS_copy(&old_parms, &parms) ;

    // This appears to adjust the cost weights based on the iteration but in
    // practice, they never change because of the PARMS_copy above and below
    parms.l_nspring *= spring_scale ; 
    parms.l_spring *= spring_scale ; 
    // This line with tspring being ajusted twice was probably originally a typo
    // but it has existed this way for a long time. It was changed after 
    // version 6, but I just changed it back for consistency. 
    parms.l_tspring *= spring_scale ;  parms.l_tspring *= spring_scale ;

    //parms.l_tspring = MIN(1.0,parms.l_tspring) ; // This had a bad effect on highres and no effect on 1mm
    parms.l_nspring = MIN(1.0, parms.l_nspring) ;
    parms.l_spring = MIN(1.0, parms.l_spring) ;
    printf("Positioning White Surface: tspring = %g, nspring = %g, spring = %g, niters = %d ",
	   parms.l_tspring,parms.l_nspring,parms.l_spring,parms.niterations); 
    printf("l_repulse = %g, checktol = %d\n",parms.l_repulse,parms.check_tol);fflush(stdout);
    if(Gdiag_no > 0){
      vgdiag = &mris->vertices[Gdiag_no];
      printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d, xyz=[%g,%g,%g]; nxyz=[%g,%g,%g];\n",
	     Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag,
	     vgdiag->x,vgdiag->y,vgdiag->z,vgdiag->nx,vgdiag->ny,vgdiag->nz);
    }
    MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
    if(Gdiag_no > 0){
      vgdiag = &mris->vertices[Gdiag_no];
      printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d, xyz=[%g,%g,%g]; nxyz=[%g,%g,%g];\n",
	     Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag,
	     vgdiag->x,vgdiag->y,vgdiag->z,vgdiag->nx,vgdiag->ny,vgdiag->nz);
    }

    old_parms.start_t = parms.start_t ;
    INTEGRATION_PARMS_copy(&parms, &old_parms) ;

    if (add) {
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
        while (MRISdivideLongEdges(mris, max_len) > 0) {}
    }

    if (!n_averages)
      break ; 

  } // end major loop placing the white surface using
  // ==========================================================================

  if(!nowhite) {
    printf("Done placing white\n");
    if(ripfilename){
      printf("Saving ripflag overlay to %s\n",ripfilename);
      MRI *mritmp = MRIcopyMRIS(NULL, mris, 0, "ripflag");
      MRIwrite(mritmp,ripfilename);
      MRIfree(&mritmp);
    }
    if(whitesigmafilename){
      printf("Saving white sigma overlay to %s\n",whitesigmafilename);
      MRI *mritmp = MRIcopyMRIS(NULL, mris, 0, "val2");
      MRIwrite(mritmp,whitesigmafilename);
      MRIfree(&mritmp);
    }
    printf("Unripping vertices\n");
    MRISunrip(mris) ;
  }

  if(nowhite){ /* read in previously generated white matter surface */
    if (orig_white)    {
      sprintf(fname, "%s%s", orig_white, suffix) ;
      printf("reading white vertex positions from %s...\n",orig_white) ;
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read white matter surface.",Progname) ;
      MRISremoveIntersections(mris,0) ;
    }
    else // read default white (something needs to be read if nowhite was created)
    {
      // if you don't like the default, give an error message here and exit,
      // to force passing the -orig_white white
      sprintf(fname, "%s%s", white_matter_name, suffix) ;
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read white matter surface.",
                  Progname) ;
    }
    MRIScomputeMetricProperties(mris) ;
  } // end nowhite (read in previously generated white matter surface)


  if (mri_aseg) //update aseg using either generated or orig_white
  {
    //    MRISripMidline(mris, mri_aseg, mri_T1, hemi, GRAY_CSF, fix_mtl) ;   //moved to later
    if (write_aseg_fname) {
      edit_aseg_with_surfaces(mris, mri_aseg) ;
      printf("writing corrected aseg to %s\n", write_aseg_fname) ;
      MRIwrite(mri_aseg, write_aseg_fname) ;
    }
  }

  // NJS HACK: if filename passed to -white is "NOWRITE", then dont write
  // the white, curv, area, and cortex.label files.  this is in lieu of
  // -nowhite not creating pial surfaces that match those created
  // w/o the -nowhite option.
  if (!nowhite && strcmp(white_matter_name,"NOWRITE"))  {
    // white_matter_name != "NOWRITE"

    if(SaveResidual){
      RmsValErrorRecord = 1;
      mrisRmsValError(mris, mri_T1);
      RmsValErrorRecord = 0;
      MRI *ValResid;
      ValResid = MRIcopyMRIS(NULL, mristarget, 3, "d"); // target distance
      MRIcopyMRIS(ValResid, mristarget, 2, "val"); // target val
      MRIcopyMRIS(ValResid, mris, 1, "valbak"); // value sampled at vertex
      MRIcopyMRIS(ValResid, mris, 0, "val2bak"); // residual = sample-target
      if (getenv("FS_POSIX")) {
        int req = snprintf(fname,STRLEN,"./%s.%s.res%s%s.mgz", 
			   hemi, white_matter_name, output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      } else {
        int req = snprintf(fname,STRLEN,"%s/%s/surf/%s.%s.res%s%s.mgz", 
			   sdir, sname, hemi, white_matter_name, output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      printf("Saving white value residual to %s\n",fname);
      MRIwrite(ValResid,fname);
      MRIfree(&ValResid);
    }

    if(Gdiag_no > 0){
      vgdiag = &mris->vertices[Gdiag_no];
      printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",
	     Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
    }
    if(SaveTarget) {
      if (getenv("FS_POSIX")) {
        int req = snprintf(fname, STRLEN, "./%s.%s.target%s%s",
			   hemi, white_matter_name, output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      } else {
        int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s.target%s%s",
			   sdir, sname, hemi, white_matter_name, output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      printf("writing white target surface to %s...\n", fname) ;
      MRISwrite(mristarget, fname) ;
    }

    MRISremoveIntersections(mris,0) ;
    if(smoothwm > 0){
      printf("Averaging white surface by %d iterations\n",smoothwm);
      MRISaverageVertexPositions(mris, smoothwm); // "smoothwm" is a bad name
    }
    if (getenv("FS_POSIX")) {
      int req = snprintf(fname,STRLEN,"./%s.%s%s%s", hemi, white_matter_name, output_suffix, suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname,STRLEN,"%s/%s/surf/%s.%s%s%s",
			 sdir, sname, hemi, white_matter_name, output_suffix, suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("writing white surface to %s...\n", fname) ;
    MRISwrite(mris, fname) ;

    printf("Output Surface Quality Stats\n");
    MRISprettyPrintSurfQualityStats(stdout, mris);

    if(mri_aseg && label_cortex) {
      // Label cortex based on aseg (4=ndilate,4=nerode)
      LABEL *lcortex = MRIScortexLabelDECC(mris, mri_aseg, 4, 4, -1, 0) ;
      if (getenv("FS_POSIX")) {
        int req = snprintf(fname,STRLEN,"./%s.%s%s%s.label",
			   hemi, "cortex", output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      } else {
        int req = snprintf(fname,STRLEN,"%s/%s/label/%s.%s%s%s.label",
			   sdir, sname, hemi, "cortex", output_suffix, suffix);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      printf("Writing cortical label to %s\n",fname);
      LabelWrite(lcortex, fname) ;
      LabelFree(&lcortex) ;
    }// done writing out white surface

    MRISfree(&mristarget);

    if (create)   /* write out curvature and area files */
    {
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, curvature_avgs) ;
      int req = snprintf(fname, STRLEN, "%s.curv%s%s",
			 mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
			 suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "writing smoothed curvature to %s\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      req = snprintf(fname, STRLEN, "%s.area%s%s",
		     mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
		     suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "writing smoothed area to %s\n", fname) ;
      MRISwriteArea(mris, fname) ;
      MRISprintTessellationStats(mris, stderr) ;
    }
  } // done with white surface placement

  if (white_only)  {
    msec = then.milliseconds() ;
    printf("refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
    MRIfree(&mri_T1);
    printf("#VMPC# mris_make_surfaces VmPeak  %d\n",GetVmPeak());
    exit(0) ;
  }

  // ==========================================================================
  // Place pial surface. #pial
  printf("\n\n\nPlacing pial surface (pialpial)\n");

  MRISsetVal2(mris, 0) ;   // will be marked for vertices near lesions
  MRISunrip(mris) ;

  // Rip vertices in the midline and other places
  if(mri_aseg) {
    if(Gdiag_no > 0){
      vgdiag = &mris->vertices[Gdiag_no];
      printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
    }
    printf("Freezing midline and others\n");  fflush(stdout);
    MRISripMidline(mris, mri_aseg, mri_T1, hemi, GRAY_CSF, fix_mtl) ;
    if(Gdiag_no > 0){
      vgdiag = &mris->vertices[Gdiag_no];
      printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d\n",Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag);
    }
  }

  parms.t = parms.start_t = 0 ;
  sprintf(parms.base_name, "%s%s%s", pial_name, output_suffix, suffix) ;
  parms.niterations = ngray ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
  parms.l_surf_repulse = l_surf_repulse ;

  MRISsetVals(mris, -1) ;  /* clear target intensities */

  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    if(orig_pial)  printf(" .. but this will be overwritten using vertex positions from %s...\n", orig_pial) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  if (fill_interior)  { // off by default
    MRI *mri_interior ;
    if (fill_interior== 255 && mri_T1_pial->type == MRI_UCHAR)
    {
      // see MRIcomputeBorderValues which changes all 255 to 0
      printf("!!!!!!!! WARNING - 255 is a protected value for uchar volumes. Changing to 254... !!!!!!!!!!!!!!!!!!!!!!\n") ;
      fill_interior = 254 ;
    }
    mri_interior = MRIcopy(mri_T1_pial, NULL) ;
    MRIScopyVolGeomFromMRI(mris, mri_T1_pial) ;
    MRISfillInterior(mris, mri_T1_pial->xsize, mri_interior) ;
    MRImask(mri_T1_pial, mri_interior, mri_T1_pial, 1, fill_interior) ;
    printf("writing filled_interior.mgz\n") ;
    MRIwrite(mri_T1_pial, "filled_interior.mgz") ; 
    MRIfree(&mri_interior) ;
  }

  if(orig_pial)  {
    // often white.preaparc is passed as orig_pial
    printf("reading initial pial vertex positions from %s...\n", orig_pial) ;

    if (longitudinal)    {
      //save final white location into TMP_VERTICES (v->tx, v->ty, v->tz)
      MRISsaveVertexPositions(mris, TMP_VERTICES);
    }

    if (MRISreadVertexPositions(mris, orig_pial) != NO_ERROR)
      ErrorExit(Gerror, "reading orig pial positions failed") ;

    if(smooth_pial){ // off by default
      printf("smoothing pial surface for %d iterations before deformation\n", smooth_pial) ;
      MRISaverageVertexPositions(mris, smooth_pial) ;
    }

    MRISremoveIntersections(mris,0) ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;

    if (longitudinal) {
      //reset the starting position to be
      //slightly inside the orig_pial in the longitudinal case
      //between final white and orig pial
      MRISblendXYZandTXYZ(mris, 0.75f, 0.25f);
    }

    MRIScomputeMetricProperties(mris) ; //shouldn't this be done whenever
    // orig_pial is used??? Maybe that's why the cross-intersection
    // was caused
  }

  /*    parms.l_convex = 1000 ;*/
  mri_T1 = mri_T1_pial ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_T1, "p.mgz") ;

  if (dura_echo_name)
  {
    #define MAX_VOLUMES 100
    char fname[STRLEN], fmt[STRLEN] ;
    MRI *mri_ratio, *mri_echos[MAX_VOLUMES] ;
    int  e ;

    printf("masking dura from pial surface locations\n") ;
    sprintf(parms.base_name, "%s%s%s", pial_name, output_suffix, suffix) ;
    for (e = 0 ; e < nechos ; e++)
    {
      if (e != 0 && e != nechos-1)
      {
        mri_echos[e] = NULL ;
        continue ;
      }
      int req = snprintf(fmt, STRLEN, "%s/%s/mri/%s", sdir, sname, dura_echo_name) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      req = snprintf(fname, STRLEN, fmt, e) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      mri_echos[e] = MRIread(fname) ;
      if (mri_echos[e] == NULL)
        ErrorExit
        (ERROR_BADPARM,
         "%s: could not read %dth echo for dura localization from %s",
         Progname, e, fname) ;
      if (mri_echos[e]->type != MRI_FLOAT)
      {
	MRI *mri_tmp = MRIchangeType(mri_echos[e], MRI_FLOAT, 0, 1, 1) ;
	printf("changing echo voxel type to FLOAT\n") ;
	MRIfree(&mri_echos[e]) ;
	mri_echos[e] = mri_tmp ;
      }
      if (MRImatchDimensions(mri_T1, mri_echos[e]) == 0)
      {
	MRI *mri_tmp = MRIresample(mri_echos[e], mri_T1, SAMPLE_TRILINEAR) ;
	printf("echo not conformed - reslicing to match conformed volumes\n") ;
	MRIfree(&mri_echos[e]) ;
	mri_echos[e] = mri_tmp ;
      }
    }

    mri_ratio = MRIdivide(mri_echos[0], mri_echos[nechos-1], NULL) ;
    if (dura_thresh < 0)
      thresh = compute_brain_thresh(mris, mri_ratio, nsigma) ;
    else
      thresh = dura_thresh ;
    mri_mask = mri_ratio ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      int req = snprintf(fname, STRLEN, "%s_ratio.mgz", parms.base_name) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing dura ratio image to %s...\n", fname) ;
      MRIwrite(mri_ratio, fname) ;
    }
    parms.mri_dura = mri_ratio ;
    parms.dura_thresh = thresh ;
    parms.l_dura = 10 ;
    printf("setting dura threshold to %2.3f\n", thresh) ;

    for (e = 0 ; e < nechos ; e++)
      if (mri_echos[e])
      {
        MRIfree(&mri_echos[e]) ;
      }
  } // end dura_echo_name
  else
    mri_mask = NULL ;

  sprintf(parms.base_name, "%s%s%s", pial_name, output_suffix, suffix) ;
  fprintf(stdout, "repositioning cortical surface to gray/csf boundary.\n") ;
  parms.l_repulse = l_repulse ;
  MRISclearMark2s(mris) ;  // will be used to keep track of which vertices found pial surface in previous cycle
  //  for (j = 0 ; j <= 0 ; parms.l_intensity *= 2, j++)  /* only once for now */

  for (j = 0 ; j < 1 ; j++)  /* only once for now */
  {
    current_sigma = pial_sigma ;

    for (n_averages = max_pial_averages, i = 0 ;
         n_averages >= min_pial_averages ;
         n_averages /= 2, current_sigma /= 2, i++)
    {

      printf("Iteration %d pial =========================================\n",i);
      printf("j=%d, i=%d, sigma=%g\n",j,i,current_sigma);

      if(flair_or_T2_name) {
        // Note: l_intensity is turned off at the end
        MRI  *mri_flair = nullptr;
        LABEL **labels = nullptr;
        int nlabels = 0;
        char fname[STRLEN] ;

        strcpy(fname, flair_or_T2_name) ;
        if (MGZ) strcat(fname, ".mgz");

        printf("flair_or_T2_name: repositioning pial surface locations using  %s\n", fname) ;
        if (mri_flair)  MRIfree(&mri_flair) ;

	if (n_averages <= min_pial_averages)
	  parms.l_repulse = 0 ;  // for final round turn off the 'unpinching' term

        mri_flair = MRIread(fname) ;
        if (mri_flair == NULL)
          ErrorExit(ERROR_NOFILE, "%s: could not load flair volume %s", Progname, fname) ;

	if (MRImatchDimensions(mri_flair, mri_aseg) == 0){
	  MRI *mri_tmp ;
	  printf("resampling to be in voxel space of T2/FLAIR\n") ;
	  mri_tmp = MRIresample(mri_aseg, mri_flair, SAMPLE_NEAREST) ;
	  MRIfree(&mri_aseg) ;
	  mri_aseg = mri_tmp ;
	}

	if(mri_aseg && erase_cerebellum){ // off by default
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg,  Left_Cerebellum_Cortex, 0) ;
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg,  Right_Cerebellum_White_Matter, 0) ;
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg, Right_Cerebellum_Cortex, 0) ;
	}

	if (mri_aseg && erase_brainstem) {  // off by default
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg, Brain_Stem, 0) ;
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg,  Left_VentralDC, 0) ;
	  MRImaskLabel(mri_flair, mri_flair, mri_aseg,  Right_VentralDC, 0) ;
	}

        if (read_pinch_fname) {  // off by default
          char marked_fname[STRLEN] ;
          MRISreadVertexPositions(mris, read_pinch_fname) ;
          sprintf(marked_fname, "%s.marked.mgz", read_pinch_fname) ;
          MRISreadMarked(mris, marked_fname) ;
          MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
        }
        else if (unpinch){  // off by default
          if (mri_aseg && erase_cerebellum) {
            MRImaskLabel(mri_flair, mri_flair, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
            MRImaskLabel(mri_flair, mri_flair, mri_aseg, Left_Cerebellum_Cortex, 0) ;
            MRImaskLabel(mri_flair, mri_flair, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
            MRImaskLabel(mri_flair, mri_flair, mri_aseg, Right_Cerebellum_Cortex, 0) ;
          }
          find_and_mark_pinched_regions(mris, mri_flair,
                                        nsigma_below, nsigma_above) ;
          if (MRIScountMarked(mris, 1) > 0)
          {
            INTEGRATION_PARMS saved_parms ;

            MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
            INTEGRATION_PARMS_copy(&saved_parms, &parms) ;

            printf("%d vertices found in imminent self-intersecting "
                   "regions, deforming to remove...\n",
                   MRIScountMarked(mris,1));
            MRISinvertMarks(mris) ;
            MRISerodeMarked(mris, 2) ;
            MRISripMarked(mris) ;

            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISremoveCompressedRegions(mris, .5) ;
            MRISwrite(mris, "after_uncompress") ;
            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISsoapBubbleVertexPositions(mris, 25) ;
            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISwrite(mris, "after_soap") ;

            MRISunrip(mris) ;
            MRISremoveIntersections(mris,0) ;
          }

          MRISwriteMarked(mris, "distant.mgz") ;
          MRIScomputeMetricProperties(mris) ;
        } //end if(unpinch)
        else
          nlabels = 0 ;

	MRISremoveIntersections(mris,0) ;
	
	GetMemUsage(memusage);  printf("Pre Pial Targ Loc VmPeak %d\n",memusage[1]);fflush(stdout);

	printf("CPTL: nlabels=%d, T2_min_inside=%g, T2_max_inside=%g, T2_min_outside=%g, T2_max_outside=%g, max_out_d=%g\n",
           nlabels,T2_min_inside, T2_max_inside, T2_min_outside, T2_max_outside, max_outward_dist);
	printf("CPTL: inside_peak_pct=%g,%g outside_peak_pct=%g,%g , wm_weight=%g\n",
               Ghisto_left_inside_peak_pct, Ghisto_right_inside_peak_pct, Ghisto_left_outside_peak_pct, Ghisto_right_outside_peak_pct, wm_weight);
	if ( followGradients)
	{
		std::cout << "T2/FLAIR compute target based on gradients "<< std::endl;
		int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s",
				   sdir, sname, hemi, orig_sphere, suffix) ;
		if( req >= STRLEN ) {
		  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
		}
		printf("reading sphere position from %s...\n", fname) ;
		//MRIS* sph = MRISread(fname);
		//MRIStoParameterization(sph,NULL, 1,0);

		MRIS_MultimodalRefinement* refine = new MRIS_MultimodalRefinement();
		MRI* whiteMR= MRIcopy(mri_T1_pial,NULL);
		MRI* vesselMR= MRIcopy(mri_T1_pial,NULL);
		refine->SegmentWM(mri_T1_pial,mri_flair, whiteMR, contrast_type);
		refine->SegmentVessel(mri_T1_pial,mri_flair, vesselMR, contrast_type);

		refine->SetStep(.4);
		refine->SetNumberOfSteps(8);
		refine->SetGradientSigma(.3);
		refine->SetSegmentation(mri_aseg);
		refine->FindMaximumGradient(contrast_type== CONTRAST_T2);
		refine->addImage(mri_T1_pial);
		refine->addImage(mri_flair);
		refine->SetWhiteMR(whiteMR);
		refine->SetVesselMR(vesselMR);
		//refine->SetSphere(sph);
		refine->getTarget(mris); //, debugVertex);
		
		MRIfree(&whiteMR);
		MRIfree(&vesselMR);
		delete refine;
	}
	else
	{
	  MRIScomputePialTargetLocationsMultiModal(mris, mri_flair, labels, nlabels,
				      contrast_type, mri_aseg, 
				      T2_min_inside, T2_max_inside, 
				      T2_min_outside, T2_max_outside, 
				      max_outward_dist,
				      Ghisto_left_inside_peak_pct, Ghisto_right_inside_peak_pct, 
				      Ghisto_left_outside_peak_pct, Ghisto_right_outside_peak_pct, 
 				      wm_weight, pial_sigma, mri_T1_pial) ;
	}
	//	parms.l_location /= spring_scale ; 
	GetMemUsage(memusage);  printf("Post Pial Targ Loc VmPeak %d\n",memusage[1]);fflush(stdout);

        if(Gdiag & DIAG_WRITE) {
          char fname[STRLEN] ;
	  static int n = 0 ;
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrestoreVertexPositions(mris, TARGET_VERTICES) ;
          int req = snprintf(fname, STRLEN, "%s.flair.target.%3.3d", parms.base_name, n++) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          printf("writing surface targets to %s\n", fname) ;
          MRISwrite(mris, fname) ;
          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        }
	//  parms.l_histo = 1 ;
        parms.l_intensity = 0 ;
	//  parms.l_max_spring = .1 ;
	//  parms.l_parea = .05 ;
      } // end if(FLAIR_or_T2)

      // T2=========================================
      if(T2_name) {
	// Note: when T2 is specified, it is assigned to flair_or_T2_name, and this
	// branch is not activated. It does not appear that there is a way to 
	// activate this branch.
        // Note: l_location is set to 0.5 at the end and l_intensity is turned off
        MRI  *mri_T2 = NULL ;
        int n = 0 ;
        LABEL             **labels ;
        int               nlabels ;

        printf("removing non-brain from pial surface locations using T2 volume %s\n", T2_name) ;
        if (mri_T2)  // first time
          MRIfree(&mri_T2) ;

        mri_T2 = MRIread(T2_name) ;
        if (mri_T2 == NULL)
          ErrorExit(ERROR_NOFILE,"%s: could not load T2 volume %s", Progname, T2_name) ;

        if (read_pinch_fname) {
          char marked_fname[STRLEN] ;

          MRISreadVertexPositions(mris, read_pinch_fname) ;
          int req = snprintf(marked_fname, STRLEN, "%s.marked.mgz", read_pinch_fname) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          MRISreadMarked(mris, marked_fname) ;
          MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
        }
        else if (unpinch)  {
          if (mri_aseg && erase_cerebellum)
          {
            MRImaskLabel(mri_T2, mri_T2, mri_aseg, Left_Cerebellum_White_Matter, 0) ;
            MRImaskLabel(mri_T2, mri_T2, mri_aseg, Left_Cerebellum_Cortex, 0) ;
            MRImaskLabel(mri_T2, mri_T2, mri_aseg, Right_Cerebellum_White_Matter, 0) ;
            MRImaskLabel(mri_T2, mri_T2, mri_aseg, Right_Cerebellum_Cortex, 0) ;
          }
          find_and_mark_pinched_regions(mris, mri_T2,
                                        nsigma_below, nsigma_above) ;
          if (MRIScountMarked(mris, 1) > 0)
          {
            INTEGRATION_PARMS saved_parms ;

            MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
            INTEGRATION_PARMS_copy(&saved_parms, &parms) ;

            printf("%d vertices found in imminent self-intersecting regions, deforming to remove...\n", 
                   MRIScountMarked(mris,1));
            MRISinvertMarks(mris) ;
            MRISerodeMarked(mris, 2) ;
            MRISripMarked(mris) ;

            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISremoveCompressedRegions(mris, .5) ;
            MRISwrite(mris, "after_uncompress") ;
            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISsoapBubbleVertexPositions(mris, 25) ;
            MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
            MRISwrite(mris, "after_soap") ;

            MRISunrip(mris) ;
            MRISremoveIntersections(mris,0) ;
          }

          MRISwriteMarked(mris, "distant.mgz") ;
          MRIScomputeMetricProperties(mris) ;
        }// if unpinch
        else
        {
          nlabels = 0 ;
        }

        if (0)
        {
		std::cout << " nooooooooooooooooooooooooooooooooooooooooooooooooooo" <<std::endl;
          MRIScomputePialTargetLocationsMultiModal(mris, mri_T2,
                                        labels, nlabels,
                                        contrast_type, mri_aseg, T2_min_inside, T2_max_inside, 
					T2_min_outside, T2_max_outside,max_outward_dist,
					Ghisto_left_inside_peak_pct,Ghisto_right_inside_peak_pct, 
					Ghisto_left_outside_peak_pct,Ghisto_right_outside_peak_pct, 
					wm_weight, pial_sigma, mri_T1_pial) ;
        }

        if(Gdiag & DIAG_WRITE)
        {
          char fname[STRLEN] ;
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrestoreVertexPositions(mris, TARGET_VERTICES) ;
          int req = snprintf(fname, STRLEN, "%s.T2.target.%3.3d", parms.base_name, n++) ;
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          printf("writing surface targets to %s\n", fname) ;
          MRISwrite(mris, fname) ;
          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        }
	//        parms.l_histo = 1 ;
        parms.l_location = .5 ;
        parms.l_intensity = 0 ;
      } // if T2_name ==================================================

      parms.sigma = current_sigma ;
      // This mri_kernel does not appear to do anything
      mri_kernel = MRIgaussian1d(current_sigma, 100) ;
      printf("smoothing T1 volume with sigma = %2.3f\n",current_sigma) ;

      parms.n_averages = n_averages ;
      parms.l_tsmooth = l_tsmooth ;

      if (fill_interior == 0){ // this happens by default default
        /* Replace bright stuff such as eye sockets with 255.  Simply
        zeroing it out would make the border always go through the
        sockets, and ignore subtle local minima in intensity at the
        border of the sockets.  Will set to 0 after border values have
        been computed so that it doesn't mess up gradients.      */
	printf("removing bright stuff\n");
	MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 255) ;
	MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_BORDER_LABEL, MID_GRAY) ;
	printf("done removing bright stuff\n");
      }

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_T1, "pial_masked.mgz") ; 

      if (mri_cover_seg) {// off by default
	if(i == 0) {
	  MRI *mri_bin, *mri_tmp;
	  int itmp=0;
	  if(mri_T1 == mri_T1_pial) itmp=1;
	  mri_bin = MRIclone(mri_T1, NULL) ;
	  mri_tmp = MRIScoverSeg(mris, mri_bin, mri_cover_seg, GRAY_CSF);
	  MRIfree(&mri_bin) ; 
	  MRIfree(&mri_T1) ;
	  mri_T1 = mri_tmp ;
	  if(itmp) mri_T1_pial = mri_tmp ;
	}
      }
      else if (flair_or_T2_name == NULL) 
      {
        // This is the main CBV for pial
        // Intensity is ignored with T2 or FLAIR
        // inside_hi = max_gray (eg, 99.05)
        // border_hi = max_gray_at_csf_border = meanGM-1stdGM (eg, 65.89)
        // border_low = min_gray_at_csf_border = meanGM-V*stdGM (V=3) (eg, 45.7)
        // outside_low = min_csf = meanGM - MAX(0.5,(V-1)*stdGM) (eg, 10)
        // outside_hi  = (max_csf+max_gray_at_csf_border)/2 (eg, 60.8)
	printf("Starting ComputeBorderVals for pial\n");
        MRIScomputeBorderValues(mris, mri_T1, mri_smooth, max_gray,
         max_gray_at_csf_border, min_gray_at_csf_border,
         min_csf,(max_csf+max_gray_at_csf_border)/2,
         current_sigma, 2*max_thickness, parms.fp,
 	 GRAY_CSF, mri_mask, thresh, parms.flags,mri_aseg,-1,-1) ;
        MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
      }

      if (!FZERO(pial_target_offset)){
	MRISaddToValues(mris, pial_target_offset) ;
	printf("adding %2.2f to pial targets\n", pial_target_offset);
      }

      {// open scope -- this section generally has no effect because pial_num=0
        int ii, vno, n, vtotal ;
        for (ii = 0; ii < pial_num ; ii++)
        {
          vno = pial_vnos[ii] ;
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
          VERTEX                * const v  = &mris->vertices         [vno];
          v->val = pial_vals[ii] ;
          v->marked = 1 ;
	  if (v->val2 > 0)
	    v->val2 = 3*current_sigma ;
	  else
	    v->val2 = current_sigma ;

          vtotal = 0 ;
          switch (pial_nbrs)
          {
          case 1:
            vtotal = vt->vnum ;
            break ;
          case 2:
            vtotal = vt->v2num ;
            break ;
          case 3:
            vtotal = vt->v3num ;
            break ;
          default:
            break ;
          }
          for (n = 0 ; n < vtotal ; n++)
          {
            mris->vertices[vt->v[n]].val = pial_vals[ii] ;
            mris->vertices[vt->v[n]].marked = 1 ;
            mris->vertices[vt->v[n]].val2 = current_sigma ;
          }
        }
       } // just an end of scope

      if (pial_blur_sigma>0) { // off by default
	MRI *mri_tmp ;
	printf("blurring input image with sigma = %2.3f\n", pial_blur_sigma) ;
	mri_kernel = MRIgaussian1d(pial_blur_sigma, 100) ;
	mri_tmp = MRIconvolveGaussian(mri_T1_pial, NULL, mri_kernel) ;
	MRIfree(&mri_kernel) ; MRIfree(&mri_T1_pial) ;
	mri_T1_pial = mri_T1 = mri_tmp ;
	MRIwrite(mri_T1_pial, "pblur.mgz") ;
      }

      if(vavgs){ // 5 by default, should do nothing when T2/FLAIR
        printf("averaging target values for %d iterations...\n",vavgs) ;
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no >= 0){
          VERTEX const * const v = &mris->vertices[Gdiag_no] ;
          printf("v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f, ripflag=%d\n",
           Gdiag_no, v->val, v->mean, v->d, v->ripflag) ;
        }
      }

      if (write_vals) {
        int req = snprintf(fname, STRLEN, "./%s-gray%2.2f.mgz", hemi, current_sigma) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        MRISwriteValues(mris, fname) ;
      }

      // probably don't need this, does not appear that anything is smoothed
      if (!mri_smooth)
        mri_smooth = MRIcopy(mri_T1, NULL) ;

      // First time increase smoothness to prevent pinching when pushing out from white surface
      // Generally has no effect with T2/FLAIR because orig_pial usually != NULL
      if(j == 0 && i == 0 && orig_pial == NULL)  {
        // Problem: this stage is not run in recon-all because orig_pial is specified, but
        // orig_pial is actually white.preaparc, which is not a real pial surface. 
	INTEGRATION_PARMS  saved_parms ;
	int k, start_t ;

	printf("perforing initial smooth deformation to move away from white surface\n") ;

	// Copy the params to restore after this stage
	INTEGRATION_PARMS_copy(&saved_parms, &parms) ;

	// parms.l_intensity /= 5 ;
	// parms.l_spring = 1 ;
	parms.dt /= 10 ;   // take small steps to unkink pinched areas
	parms.niterations = 10 ;
	parms.l_surf_repulse /= 5 ;

	for (k = 0 ; k < 3 ; k++){
	  INTEGRATION_PARMS_copy(&old_parms, &parms) ;
	  parms.l_tspring *= spring_scale ; parms.l_nspring *= spring_scale ; parms.l_spring *= spring_scale ;
	  parms.l_tspring = MIN(1.0,parms.l_tspring) ;
	  parms.l_nspring = MIN(1.0, parms.l_nspring) ;
	  parms.l_spring = MIN(1.0, parms.l_spring) ;
	  MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
	  old_parms.start_t = parms.start_t ;
	  INTEGRATION_PARMS_copy(&parms, &old_parms) ;
	  if (!FZERO(parms.l_intensity)){ // l_intensity should be 0 for T2/FLAIR

	    MRIScomputeBorderValues
	      (mris, mri_T1, mri_smooth, max_gray,
	       max_gray_at_csf_border, min_gray_at_csf_border,
	       min_csf,(max_csf+max_gray_at_csf_border)/2,
	       current_sigma, 2*max_thickness, parms.fp,
	       GRAY_CSF, mri_mask, thresh, parms.flags,mri_aseg,-1,-1) ;
	  }
        }// loop over k
	start_t = parms.start_t ;
	INTEGRATION_PARMS_copy(&parms, &saved_parms) ;
	parms.start_t = start_t ;
	parms.l_surf_repulse = l_surf_repulse ;
      }// end if (j == 0 && i == 0 && orig_pial == NULL)  

      if(mristarget != NULL) MRISfree(&mristarget);
      
      mristarget = MRISclone(mris);

      float *px,*py,*pz;
      MRISexportXYZ(mristarget, &px, &py, &pz);
      
      for(vno=0; vno < mris->nvertices; vno++){
        VERTEX * const v = &(mristarget->vertices[vno]);
        if(flair_or_T2_name == NULL) {
          v->d   = mris->vertices[vno].d;   // clone does not copy this
          v->val = mris->vertices[vno].val; // clone does not copy this
          px[vno] += (v->d*v->nx);          // d is the distance to the max gradient
          py[vno] += (v->d*v->ny);
          pz[vno] += (v->d*v->nz);
        }
	else{
	  VERTEX const * const v2 = &(mris->vertices[vno]);
          px[vno] = v2->targx; // clone does not copy this
          py[vno] = v2->targy; 
          pz[vno] = v2->targz; 
        }
      }

      MRISimportXYZ(mristarget,  px,  py,  pz);
      freeAndNULL(px);
      freeAndNULL(py);
      freeAndNULL(pz);

      INTEGRATION_PARMS_copy(&old_parms, &parms) ;
      parms.l_tspring *= spring_scale ; parms.l_nspring *= spring_scale ; parms.l_spring *= spring_scale ;
      parms.l_tspring = MIN(1.0,parms.l_tspring) ;
      parms.l_nspring = MIN(1.0, parms.l_nspring) ;
      parms.l_spring = MIN(1.0, parms.l_spring) ;
      MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
      old_parms.start_t = parms.start_t ;
      INTEGRATION_PARMS_copy(&parms, &old_parms) ;

      if (parms.l_location > 0)  { // 0 by default, but turned on with T2
	int vno ;
	double dist ;

	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
	  VERTEX * const v = &mris->vertices[vno] ;
	  dist = sqrt(SQR(v->targx - v->x) + SQR(v->targy-v->y) + SQR(v->targz-v->z)) ;
	  v->cropped = nint((float)v->cropped * dist) ;
	}
	if (Gdiag_no > 0)
	  MRISwriteCropped(mris, "cropped.mgz") ;
	DiagBreak() ;
	MRISsetCroppedToZero(mris) ;
      }

      /*  parms.l_nspring = 0 ;*/

      if (!n_averages)
        break ;

    } // end loop over number of averages 
  } // end loop over j (only 1 for now)

  MRISremoveIntersections(mris,0) ;
  if (getenv("FS_POSIX")) {
    int req = snprintf(fname, STRLEN, "./%s.%s%s%s", hemi, pial_name, output_suffix, suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } else {
    int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", 
		       sdir, sname, hemi, pial_name, output_suffix, suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
  fprintf(stdout, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;

  if(pialsigmafilename){
    printf("Saving pial sigma overlay to %s\n",pialsigmafilename);
    MRI *mritmp = MRIcopyMRIS(NULL, mris, 0, "val2");
    MRIwrite(mritmp,pialsigmafilename);
    MRIfree(&mritmp);
  }

  if(create)   /* write out curvature and area files */
  {
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    MRISuseMeanCurvature(mris) ;
    MRISaverageCurvatures(mris, curvature_avgs) ;
    int req = snprintf(fname, STRLEN, "%s.curv.pial%s%s",
		       mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
		       suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing smoothed curvature to %s\n", fname) ;
    MRISwriteCurvature(mris, fname) ;
    req = snprintf(fname, STRLEN, "%s.area.pial%s%s",
		   mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
		   suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing smoothed area to %s\n", fname) ;
    MRISwriteArea(mris, fname) ;
    MRISprintTessellationStats(mris, stdout) ;
  }

  if(SaveResidual && flair_or_T2_name==NULL){
    RmsValErrorRecord = 1;
    mrisRmsValError(mris, mri_T1);
    RmsValErrorRecord = 0;
    MRI *ValResid;
    ValResid = MRIcopyMRIS(NULL, mristarget, 3, "d"); // target distance
    MRIcopyMRIS(ValResid, mristarget, 2, "val"); // target val
    MRIcopyMRIS(ValResid, mris, 1, "valbak"); // value sampled at vertex
    MRIcopyMRIS(ValResid, mris, 0, "val2bak"); // residual = sample-target
    if (getenv("FS_POSIX")) {
      int req = snprintf(fname, STRLEN, "./%s.%s.res%s%s.mgz",
			 hemi, pial_name, output_suffix, suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s.res%s%s.mgz",
			 sdir, sname, hemi, pial_name, output_suffix, suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("Saving pial value residual to %s\n",fname);
    MRIwrite(ValResid,fname);
    MRIfree(&ValResid);
  }

  if(SaveTarget){
    if (getenv("FS_POSIX")) {
      int req = snprintf(fname, STRLEN, "./%s.%s.target%s%s",
			 hemi, pial_name, output_suffix, suffix);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s.target%s%s",
			 sdir, sname, hemi, pial_name,output_suffix, suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    printf("writing pial target surface to %s...\n", fname) ;
    MRISwrite(mristarget, fname) ;
  }

  if (in_out_in_flag)// off by default
  {
    int req = snprintf(parms.base_name, STRLEN, "%s%s%s",
		       white_matter_name, output_suffix, suffix) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    fprintf(stdout,
            "repositioning cortical surface to gray/white boundary\n");

    MRISsetVals(mris, -1) ;  /* clear white matter intensities */

    current_sigma = white_sigma ;
    MRImask(mri_T1, mri_labeled, mri_T1, BRIGHT_LABEL, 0) ;
    for (n_averages = max_white_averages, i = 0 ;
         n_averages >= min_white_averages ;
         n_averages /= 2, current_sigma /= 2, i++)
    {
      if (nowhite)
      {
        break ;
      }

      parms.sigma = current_sigma ;
      mri_kernel = MRIgaussian1d(current_sigma, 100) ;
      fprintf(stdout, "smoothing T1 volume with sigma = %2.3f\n",
              current_sigma) ;
      if (!mri_smooth)
      {
        mri_smooth = MRIclone(mri_T1, NULL) ;
      }
      MRIfree(&mri_kernel) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        char fname[STRLEN] ;
        sprintf(fname, "sigma%.0f.mgz", current_sigma) ;
        fprintf(stdout, "writing smoothed volume to %s...\n", fname) ;
        MRIwrite(mri_smooth, fname) ;
      }

      parms.n_averages = n_averages ;
      MRISprintTessellationStats(mris, stdout) ;
      if (flair_or_T2_name == NULL) // otherwise already done
	// within if(in_out_in_flag)
        MRIScomputeBorderValues
        (mris, mri_T1, mri_smooth, MAX_WHITE,
         max_border_white, min_border_white,
         min_gray_at_white_border, max_border_white /*max_gray*/,
         current_sigma, 2*max_thickness, parms.fp,
         GRAY_WHITE, NULL, 0, parms.flags, mri_aseg,-1,-1) ;
      MRISfindExpansionRegions(mris) ;
      if (vavgs)
      {
        fprintf
        (stdout,
         "averaging target values for %d iterations...\n",vavgs);
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no > 0)
        {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf
          (stdout,
           "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f, ripflag=%d\n",
           Gdiag_no, v->val, v->mean, v->d, v->ripflag) ;
        }
      }

      if (write_vals)
      {
        sprintf(fname, "./%s-white%2.2f.mgz", hemi, current_sigma) ;
        MRISwriteValues(mris, fname) ;
      }
      if(Gdiag_no > 0){
	vgdiag = &mris->vertices[Gdiag_no];
	printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d(%g,%g,%g)\n",
	       Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag,vgdiag->x,vgdiag->y,vgdiag->z);
      }
      MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);
      if(Gdiag_no > 0){
	vgdiag = &mris->vertices[Gdiag_no];
	printf("vno=%d  v->val=%g v->d=%g v->marked=%d, v->ripflag=%d(%g,%g,%g)\n",
	       Gdiag_no,vgdiag->val,vgdiag->d,vgdiag->marked,vgdiag->ripflag,vgdiag->x,vgdiag->y,vgdiag->z);
      }
      if (parms.l_location > 0)
      {
	int vno ;
	VERTEX *v ;
	double dist ;

	for (vno = 0 ; vno < mris->nvertices ; vno++)
	{
	  v = &mris->vertices[vno] ;
	  dist = sqrt(SQR(v->targx - v->x) + SQR(v->targy-v->y) + SQR(v->targz-v->z)) ;
	  v->cropped = nint((float)v->cropped * dist) ;
	}
	if (Gdiag_no > 0)
	  MRISwriteCropped(mris, "cropped.mgz") ;
	DiagBreak() ;
	MRISsetCroppedToZero(mris) ;
      }
      if (!n_averages)
      {
        break ;
      }
    }
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* gray/white surface */
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* pial surface */
    MRISremoveIntersections(mris,0) ;
    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s2%s", 
		   sdir, sname,hemi,white_matter_name,
		   suffix);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing gray/white surface to %s...\n", fname) ;
    MRISwrite(mris, fname) ;
  } // end if(in_and_out)
  MRIfree(&mri_T1);

  /*  if (!(parms.flags & IPFLAG_NO_SELF_INT_TEST))*/
  {
    fprintf(stdout, "measuring cortical thickness...\n") ;
    if (longitudinal)
      MRISmeasureCorticalThickness(mris, nbhd_size, 5.0) ;
    else
      MRISmeasureCorticalThickness(mris, nbhd_size, max_thickness) ;

    if(create){
      printf("writing cortical thickness estimate to 'thickness' file.\n") ;
      int req = snprintf(fname, STRLEN, "thickness%s%s", output_suffix, suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      MRISwriteCurvature(mris, fname) ;
    }

    /* at this point, the v->curv slots contain the cortical surface. Now
       move the white matter surface out by 1/2 the thickness as an estimate
       of layer IV.
    */
    if (graymid)
    {
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      mrisFindMiddleOfGray(mris) ;
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s",
			 sdir, sname, hemi, GRAYMID_NAME,
			 suffix) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "writing layer IV surface to %s...\n", fname) ;
      MRISwrite(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
  }

  msec = then.milliseconds() ;
  fprintf(stdout,"positioning took %2.1f minutes\n", (float)msec/(60*1000.0f));

  printf("#VMPC# mris_make_surfaces VmPeak  %d\n",GetVmPeak());

  printf("mris_make_surface done\n");

  exit(0) ;
  return(0) ;  /* for ansi */
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "nocerebellum") || !stricmp(option, "-erasecerebellum") || !stricmp(option, "erase_cerebellum"))
  {
    erase_cerebellum = 1 ;
    printf("erasing cerebellum in aseg before deformation\n") ;
  }
  else if (!stricmp(option, "nobrainstem") || !stricmp(option, "-erasebrainstem") || !stricmp(option, "erase_brainstem"))
  {
    erase_brainstem = 1 ;
    printf("erasing brainstem in aseg before deformation\n") ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nowmsa"))
  {
    nowmsa = 1 ;
    fprintf(stderr,  "removing WMSA from input data\n") ;
  }
  else if (!stricmp(option, "min_peak_pct"))
  {
    Ghisto_left_inside_peak_pct = atof(argv[2]) ;
    Ghisto_right_inside_peak_pct = atof(argv[3]) ;
    Ghisto_left_outside_peak_pct = atof(argv[4]) ;
    Ghisto_right_outside_peak_pct = atof(argv[5]) ;
    fprintf(stderr,  "using peak pct %2.2f,%2.2f and %2.2f,%2.2f (default=%2.2f) to find bounds of GM range\n",
	    Ghisto_left_inside_peak_pct, Ghisto_right_inside_peak_pct,
	    Ghisto_left_outside_peak_pct, Ghisto_right_outside_peak_pct, MIN_PEAK_PCT) ;
    nargs = 4 ;
  }
  else if (!stricmp(option, "repulse"))
  {
    parms.l_repulse = atof(argv[2]) ;
    l_repulse = parms.l_repulse;
    fprintf(stderr,  "setting l_repulse = %2.2f\n", parms.l_repulse) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "no-unitize"))
  {
    UnitizeNormalFace = 0;
    printf("Turning off face normal unitization\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "first-peak-d1"))
  {
    CBVfindFirstPeakD1 = 1;
    printf("Turning on FindFirstPeakD1 option for MRIScomputeBorderValues_new()\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "first-peak-d2"))
  {
    CBVfindFirstPeakD2 = 1;
    printf("Turning on FindFirstPeakD2 option for MRIScomputeBorderValues_new()\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "max_gray_scale"))
  {
    max_gray_scale = atof(argv[2]) ;
    fprintf(stderr,  "setting max_gray_scale (mgs) = %2.1f in calculation of outside_hi = (max_border_white + mgs*max_gray) / (mgs+1.0)\n", max_gray_scale) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_min_inside"))
  {
    T2_min_inside = atof(argv[2]) ; T2_min_set = 1 ;
    fprintf(stderr,  "using min interior T2 gray matter threshold %2.1f\n", T2_min_inside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_min"))
  {
    T2_min_outside = T2_min_inside = atof(argv[2]) ; T2_min_set = 1 ;
    fprintf(stderr,  "using min T2 gray matter (interior and exterior) threshold %2.1f\n", T2_min_inside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_out"))
  {
    max_outward_dist = atof(argv[2]) ;
    fprintf(stderr,  "using max outward dist %2.1f\n", max_outward_dist) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_max_inside"))
  {
    T2_max_inside_set = 1;
    T2_max_inside = atof(argv[2]) ; T2_max_set = 1 ;
    fprintf(stderr,  "using max T2 gray matter (interior and exterior) threshold %2.1f\n", T2_max_inside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_max_outside"))
  {
    T2_max_outside = atof(argv[2]) ; T2_max_set = 1 ;
    fprintf(stderr,  "using max T2 gray matter exterior threshold %2.1f\n", T2_max_outside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_min_outside"))
  {
    T2_min_outside = atof(argv[2]) ;
    fprintf(stderr,  "using min T2 exterior gray matter threshold %2.1f\n", T2_min_outside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T2_min_outside"))
  {
    T2_min_outside = atof(argv[2]) ;
    fprintf(stderr,  "using max T2 exterior gray matter threshold %2.1f\n", T2_min_outside) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "grad_dir"))
  {
    int i = atoi(argv[2]) ;
    parms.grad_dir =  i < 0 ? -1 : i > 0 ? 1 : 0 ;
    fprintf(stderr,  "setting grad dir to %d\n", parms.grad_dir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ct"))
  {
    printf("reading color table from %s\n", argv[2]) ;
    ctab = CTABreadASCII(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "soap"))
  {
    parms.smooth_intersections = 1 ;
    printf("using soap bubble smoothing to remove vertex intersections\n") ;
  }
  else if (!stricmp(option, "read_pinch"))
  {
    read_pinch_fname = argv[2] ;
    printf("reading pinch initialization from %s\n", read_pinch_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "pial_offset"))
  {
    pial_target_offset = atof(argv[2]) ;
    fprintf(stderr,  "offseting pial target vals by %2.2f\n",
            pial_target_offset) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "white_offset"))
  {
    white_target_offset = atof(argv[2]) ;
    fprintf(stderr,  "offseting white target vals by %2.0f\n",
            white_target_offset) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "both"))
  {
    remove_contra = 0 ;
    fprintf(stderr,  "not removing contralateral hemi\n") ;
  }
  else if (!stricmp(option, "nsigma") || !stricmp(option, "nsigmas"))
  {
    nsigma = atof(argv[2]) ;
    fprintf(stderr,
            "using dura threshold of %2.2f sigmas from mean (default=2)\n",
            nsigma) ;
    nargs = 1;
  }
  else if (!stricmp(option, "dura_thresh") || !stricmp(option, "dura+thresh"))
  {
    dura_thresh = atof(argv[2]) ;
    fprintf(stderr, "setting dura threshold to %2.2f instead of estimating it\n", dura_thresh) ;
    nargs = 1;
  }
  else if (!stricmp(option, "nsigma_above") ||
           !stricmp(option, "nsigmas_above"))
  {
    above_set = 1 ;
    nsigma_above = atof(argv[2]) ;
    fprintf(stderr,
            "using T2 threshold of %2.2f sigmas above the mean (default=2)\n",
            nsigma_above) ;
    nargs = 1;
  }
  else if (!stricmp(option, "wm_weight"))
  {
    wm_weight = atof(argv[2]) ;
    fprintf(stderr,
            "using wm weight to compute T2 threshold of %2.2f (default=3)\n", wm_weight) ;
    nargs = 1;
  }
  else if (!stricmp(option, "nsigma_below") ||
           !stricmp(option, "nsigmas_below"))
  {
    below_set = 1 ;
    nsigma_below = atof(argv[2]) ;
    fprintf(stderr,
            "using T2 threshold of %2.2f sigmas below the mean (default=2)\n",
            nsigma_below) ;
    nargs = 1;
  }
  else if (!stricmp(option, "dura"))
  {
    dura_echo_name = argv[2] ;
    nechos = atoi(argv[3]) ;
    fprintf(stderr,
            "detecting dura using %d echos from %s\n",
            nechos, dura_echo_name) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "T2dura") || !stricmp(option, "T2"))
  {
    flair_or_T2_name = argv[2] ;
    contrast_type = CONTRAST_T2 ;
    parms.l_location = 0.5 ;
    parms.l_nspring = 1.0 ;
    parms.l_curv = 0.5 ;
    parms.check_tol = 1 ;
    wm_weight = 3 ;
    // See mris_place_surface.cpp for a discussion on these limits. If
    // you change these defaults, then change the values in
    // mris_place_surface for completeness
    if (T2_max_inside_set == 0) T2_max_inside = 300 ;
    T2_min_inside = 110 ;
    T2_min_outside = 130 ;  
    T2_max_outside = 300 ;
    max_pial_averages = 4 ;
    max_outward_dist = 3 ;
    Ghisto_left_inside_peak_pct = 0.1 ;
    Ghisto_right_inside_peak_pct = 0.01 ;
    Ghisto_left_outside_peak_pct = 0.5 ;
    Ghisto_right_outside_peak_pct = 0.5 ;
    l_repulse = .025 ;

    fprintf(stderr,
            "refining pial surfaces placement using T2 volume %s\n", 
            flair_or_T2_name) ;
    parms.check_tol = 1 ;
    nargs = 1 ;
    printf("setting default T2 deformation parameters to: \n") ;
    printf("\tT2_max_inside = %2.1f\n", T2_max_inside) ;
    printf("\tT2_min_inside = %2.1f\n", T2_min_inside) ;
    printf("\tT2_min_outside = %2.1f\n", T2_min_outside) ;
    printf("\tT2_max_outside = %2.1f\n", T2_max_outside) ;
    printf("\tmax_pial_averages = %d\n", max_pial_averages) ;
    printf("\tmax_outward_dist = %2.1f\n", max_outward_dist) ;
    printf("\tGhisto_inside_peak_pct = %2.2f, %2.2f\n", Ghisto_left_inside_peak_pct,Ghisto_right_inside_peak_pct) ;
    printf("\tGhisto_outside_peak_pct = %2.2f, %2.2f\n", Ghisto_left_outside_peak_pct,Ghisto_right_outside_peak_pct) ;
    printf("\tparms.l_repulse  = %2.3f\n", l_repulse) ;
    printf("\tparms.l_curv     = %2.1f\n", parms.l_curv) ;
    printf("\tparms.l_nspring  = %2.1f\n", parms.l_nspring) ;
    printf("\tparms.l_location = %2.1f\n", parms.l_location) ;
  }
  else if (!stricmp(option, "norms") || !stricmp(option, "no_rms"))
  {
    parms.check_tol = 0 ;
    printf("not checking rms error decrease\n") ;
  }
  else if (!stricmp(option, "flair"))
  {
    contrast_type = CONTRAST_FLAIR ;
    parms.l_location = 0.5 ;
    parms.l_nspring = 1.0 ;
    parms.l_curv = 0.5 ;
    parms.check_tol = 1 ;
    flair_or_T2_name = argv[2] ;
    fprintf(stderr,
            "deforming surfaces based on FLAIR volume %s\n", flair_or_T2_name) ;
    nargs = 1 ;

    // See mris_place_surface.cpp for a discussion on these limits. If
    // you change these defaults, then change the values in
    // mris_place_surface for completeness
    if(T2_max_inside_set == 0) T2_max_inside = 200;
    T2_min_inside =  50;
    T2_min_outside = 10; // outside not so important for FLAIR
    T2_max_outside = 50;

    max_pial_averages = 4 ;
    max_outward_dist = 3 ;
    Ghisto_left_inside_peak_pct = 0.01 ;
    Ghisto_right_inside_peak_pct = 0.01 ;
    Ghisto_left_outside_peak_pct = 0.5 ;
    Ghisto_right_outside_peak_pct = 0.5 ;
    l_repulse = .025 ;
    printf("setting default FLAIR deformation parameters to: \n") ;
    printf("\tT2_max_inside = %2.1f\n", T2_max_inside) ;
    printf("\tT2_min_inside = %2.1f\n", T2_min_inside) ;
    printf("\tT2_min_outside = %2.1f\n", T2_min_outside) ;
    printf("\tT2_max_outside = %2.1f\n", T2_max_outside) ;
    printf("\tmax_pial_averages = %d\n", max_pial_averages) ;
    printf("\tmax_outward_dist = %2.1f\n", max_outward_dist) ;
    printf("\tGhisto_min_inside_peak_pct = %2.2f, %2.2f\n", Ghisto_left_inside_peak_pct,Ghisto_right_inside_peak_pct) ;
    printf("\tGhisto_min_outside_peak_pct = %2.2f, %2.2f\n", Ghisto_left_outside_peak_pct, Ghisto_right_outside_peak_pct) ;
    printf("\tparms.l_repulse  = %2.3f\n", l_repulse) ;
    printf("\tparms.l_curv     = %2.1f\n", parms.l_curv) ;
    printf("\tparms.l_nspring  = %2.1f\n", parms.l_nspring) ;
    printf("\tparms.l_location = %2.1f\n", parms.l_location) ;
  }
  else if (!stricmp(option, "location"))
  {
    parms.l_location = atof(argv[2]) ;
    printf("setting location term coefficent to l_location = %2.2f\n", parms.l_location) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "cortex"))
  {
    // -cortex 0 to turn off
    label_cortex = atoi(argv[2]) ;
    if((label_cortex != 0 && label_cortex != 1 ) || isalpha(argv[2][0]) || argv[2][0] == '-'){
      printf("ERROR: -cortex flag takes one argument that must be 0 or 1; you passed %s\n",argv[2]);
      exit(1);
    }
    printf("%sgenerating cortex label to subject's "
           "label/?h.cortex.label file\n", label_cortex ? "" : "not ") ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "fix_mtl"))
  {
    fix_mtl = 1 ;
    printf("not allowing deformations in hippocampus or amygdala "
           "when estimating pial surface\n") ;
  }
  else if (!stricmp(option, "mode"))
  {
    use_mode = atoi(argv[2]) ;
    printf("%susing class modes instead of means...\n",
           use_mode ? "" : "NOT ") ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "aseg"))
  {
    aseg_name = argv[2] ;
    printf("using aseg volume %s to prevent surfaces crossing the midline\n",
           aseg_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "cover_seg"))
  {
    auto_detect_stats = 0 ;
    printf("creating surfaces to cover  segmented volume %s\n", argv[2]) ;
    mri_cover_seg = MRIread(argv[2]) ;
    if (mri_cover_seg == NULL)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not read segmentation volume %s", argv[2]) ;
    }
    parms.grad_dir = 1 ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "write_aseg"))
  {
    write_aseg_fname = argv[2] ;
    printf("writing corrected  aseg volume to %s\n", write_aseg_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "noaseg"))
  {
    aseg_name = NULL ;
    printf("not using aseg volume to prevent surfaces "
           "crossing the midline\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "noaparc"))
  {
    aparc_name = NULL ;
    printf("not using aparc to prevent surfaces "
           "crossing the midline\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "wval"))
  {
    if (white_num >= MAX_VERTICES)
    {
      ErrorExit(ERROR_NOMEMORY,
                "%s: too many white vertex vals specified", Progname) ;
    }
    white_vnos[white_num] = atoi(argv[2]) ;
    white_vals[white_num] = atof(argv[3]) ;
    printf("constraining white surface val for vno %d to be %2.0f\n",
           white_vnos[white_num], white_vals[white_num]) ;
    white_num++ ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "pval"))
  {
    if (pial_num >= MAX_VERTICES)
    {
      ErrorExit(ERROR_NOMEMORY,
                "%s: too many pial vertex vals specified", Progname) ;
    }
    pial_vnos[pial_num] = atoi(argv[2]) ;
    pial_vals[pial_num] = atof(argv[3]) ;
    printf("constraining pial surface val for vno %d to be %2.0f\n",
           pial_vnos[pial_num], pial_vals[pial_num]) ;
    pial_num++ ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "pnbrs"))
  {
    pial_nbrs = atoi(argv[2]) ;
    printf("setting pvals out to %d nbrs\n", pial_nbrs) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T1") || !stricmp(option, "gvol"))
  {
    strcpy(T1_name, argv[2]) ;
    printf("using %s as T1 volume...\n", T1_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "wvol"))
  {
    white_fname = argv[2] ;
    printf("using %s as volume for white matter deformation...\n",
           white_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "hires") || !stricmp(option, "highres"))
  {
    highres_label = LabelRead(NULL, argv[2]) ;
    if (!highres_label)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read highres label %s", Progname, argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "long"))
  {
    longitudinal = 1;
    printf("Using longitudinal scheme\n");
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    setenv("SUBJECTS_DIR",sdir,1);
    nargs = 1 ;
  }
  else if (!stricmp(option, "orig_white"))
  {
    orig_white = argv[2] ;
    printf("using %s starting white location...\n", orig_white) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "orig_pial"))
  {
    orig_pial = argv[2] ;
    printf("using %s starting pial locations...\n", orig_pial) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "orig_sphere"))
  {
    orig_sphere = argv[2] ;
    printf("using %s sphere locations...\n", orig_sphere) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_border_white"))
  {
    max_border_white_set = 1 ;
    max_border_white = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "min_border_white"))
  {
    min_border_white_set = 1 ;
    min_border_white = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "wlo"))     // same flag name as mri_segment
  {
    min_border_white_set = 1 ;
    min_border_white = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "scale_std"))
  {
    std_scale = atof(argv[2]);
    printf("scale the estimated WM and GM std by %g \n", std_scale) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "min_gray_at_white_border"))
  {
    min_gray_at_white_border_set = 1 ;
    min_gray_at_white_border = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_gray"))
  {
    max_gray_set = 1 ;
    max_gray = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ghi"))     // same flag name as mri_segment
  {
    max_gray_set = 1 ;
    max_gray = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_gray_at_csf_border"))
  {
    max_gray_at_csf_border_set = 1 ;
    max_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "min_gray_at_csf_border"))
  {
    min_gray_at_csf_border_set = 1 ;
    min_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "variablesigma"))
  {
    variablesigma = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "min_csf"))
  {
    min_csf_set = 1 ;
    min_csf = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_csf"))
  {
    max_csf_set = 1 ;
    max_csf = atof(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "noauto"))
  {
    auto_detect_stats = 0 ;
    fprintf(stderr, "disabling auto-detection of border ranges...\n") ;
  }
  else if (!stricmp(option, "autodetsurf")){
    AutoDetSurfName = argv[2];
    nargs = 1 ;
  }

  else if (!stricmp(option, "inoutin"))
  {
    in_out_in_flag = 1 ;
    fprintf(stderr, "applying final white matter deformation after pial\n") ;
  }
  else if (!stricmp(option, "graymid"))
  {
    graymid = 1 ;
    fprintf(stderr, "generating graymid surface...\n") ;
  }
  else if (!strcmp(option, "rval"))
  {
    rh_label = atoi(argv[2]) ;
    nargs = 1 ;
   
    fprintf(stderr,"using %d as fill val for right hemisphere.\n", rh_label);
  }
  else if (!stricmp(option, "fill_interior"))
  {
    fill_interior = atoi(argv[2]) ;
    nargs = 1 ;
    printf("filling interior of surface with %d\n", fill_interior) ;
  }
  else if (!strcmp(option, "nbhd_size"))
  {
    nbhd_size = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d size nbhd for thickness calculation.\n",
            nbhd_size);
  }
  else if (!strcmp(option, "lval"))
  {
    lh_label = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d as fill val for left hemisphere.\n", lh_label);
  }

  else if (!stricmp(option, "white"))
  {
    strcpy(white_matter_name, argv[2]) ;
    nargs = 1 ;
    // NJS HACK: if filename passed to -white is "NOWRITE", then dont write
    // the white, curv, area, and cortex.label files.
    // this is in lieu of -nowhite not creating pial surfaces that
    // match those created w/o the -nowhite option.
    if (!strcmp(white_matter_name,"NOWRITE"))
    {
      fprintf(stderr, "-white NOWRITE indicates that white, curv, area, "
              "and cortex.label files will not be written...\n") ;
    }
    else
    {
      fprintf(stderr, "using %s as white matter name...\n",
              white_matter_name) ;
    }
  }
  else if (!stricmp(option, "whiteonly"))
  {
    white_only = 1 ;
    fprintf(stderr,  "only generating white matter surface\n") ;
  }
  else if (!stricmp(option, "nowhite"))
  {
    // Do not place the white surface. Only place the pial surface. It
    // reads in a white surface based on whatever orig_white is set to.
    // if that is null, then uses whatever white_matter_name is set to.
    // The default for this is white.
    nowhite = 1 ;
    printf("reading previously compute gray/white surface\n") ;
  }

  else if (!stricmp(option, "overlay"))
  {
    overlay = !overlay ;
    fprintf(stderr,  "%soverlaying T1 volume with edited white matter\n",
            overlay ? "" : "not") ;
  }
  else if (!stricmp(option, "pial"))
  {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "writing pial surface to file named %s\n", pial_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "write_vals"))
  {
    write_vals = 1 ;
    fprintf(stderr,  "writing gray and white surface targets to .mgz files\n") ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "spring_nzr"))
  {
    parms.l_spring_nzr = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_spring_nzr = %2.3f\n", parms.l_spring_nzr) ;
  }
  else if (!stricmp(option, "hinge"))
  {
    parms.l_hinge = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_hinge = %2.3f\n", parms.l_hinge);
  }
  else if (!stricmp(option, "spring"))
  {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "tsmooth"))
  {
    l_tsmooth = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tsmooth = %2.3f\n", l_tsmooth) ;
  }
  else if (!stricmp(option, "grad"))
  {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_grad = %2.3f\n", parms.l_grad) ;
  }
  else if (!stricmp(option, "tspring"))
  {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  }
  else if (!stricmp(option, "unpinch"))
  {
    unpinch = 1 ;
    fprintf(stderr, "removing pinches from surface before deforming\n") ;
  }
  else if (!stricmp(option, "nltspring"))
  {
    parms.l_nltspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nltspring = %2.3f\n", parms.l_nltspring) ;
  }
  else if (!stricmp(option, "nspring"))
  {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "smooth"))
  {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing for %d iterations\n", smooth) ;
  }
  else if (!stricmp(option, "smooth_pial"))
  {
    smooth_pial = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing pial surface for %d iterations before deforming\n", smooth_pial) ;
  }
  else if (!stricmp(option, "output"))
  {
    output_suffix = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "appending %s to output names...\n", output_suffix) ;
  }
  else if (!stricmp(option, "rip"))
  {
    // This is the full file name, so you need to add lh/rh, suffix, etc
    ripfilename = argv[2] ;
    printf("saving ripflag to %s\n",ripfilename);
    nargs = 1 ;
  }
  else if (!stricmp(option, "location-mov-len")){
    double locationmovlen;
    sscanf(argv[2],"%lf",&locationmovlen);
    mrisDxyzSetLocationMoveLen(locationmovlen);
    printf("Setting LOCATION_MOVE_LEN to %g\n",locationmovlen);
    // Used in mrisComputeTargetLocationTerm()
    nargs = 1 ;
  }
  else if (!stricmp(option, "sigma-white"))
  {
    // This is the full file name, so you need to add lh/rh, suffix, etc
    whitesigmafilename = argv[2] ;
    printf("saving white surf sigma to %s\n",whitesigmafilename);
    nargs = 1 ;
  }
  else if (!stricmp(option, "sigma-pial"))
  {
    // This is the full file name, so you need to add lh/rh, suffix, etc
    pialsigmafilename = argv[2] ;
    printf("saving white surf sigma to %s\n",pialsigmafilename);
    nargs = 1 ;
  }
  else if (!stricmp(option, "vavgs"))
  {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing values for %d iterations\n", vavgs) ;
  }
  else if (!stricmp(option, "intensity"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "nwhite"))
  {
    nwhite = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,
            "integrating gray/white surface positioning for %d time steps\n",
            nwhite) ;
  }
  else if (!stricmp(option, "flair_white") || !stricmp(option, "flairwhite"))
  {
    contrast_type = CONTRAST_FLAIR ;
    flairwhite = 1 ;
    nbrs = 3 ;
    parms.l_location = .5 ;
    parms.l_intensity = 0 ;
    fprintf(stderr, "deforming white matter surfaces to match FLAIR volume\n") ;
  }
  else if (!stricmp(option, "smoothwm"))
  {
    smoothwm = atoi(argv[2]) ;
    fprintf(stderr, "writing smoothed (%d iterations) wm surface\n",
            smoothwm) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "filled"))
  {
    filled_name = (argv[2]) ;
    fprintf(stderr, "using %s as filled name\n", filled_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "wm"))
  {
    wm_name = (argv[2]) ;
    fprintf(stderr, "using %s as filled name\n", wm_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ngray"))
  {
    ngray = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,
            "integrating pial surface positioning for %d time steps\n",
            ngray) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, 
            "setting tol = %2.1f%% change in sse/rms\n", parms.tol) ;
  }
  else if (!stricmp(option, "wsigma"))
  {
    white_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            white_sigma) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "psigma"))
  {
    pial_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            pial_sigma) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "pblur"))
  {
    pial_blur_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing pial volume with Gaussian sigma = %2.1f\n",
            pial_blur_sigma) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "pa"))
  {
    max_pial_averages = atoi(argv[2]) ;
    fprintf(stderr, "using max pial averages = %d\n", max_pial_averages) ;
    nargs = 1 ;
    if (isdigit(*argv[3]))
    {
      min_pial_averages = atoi(argv[3]) ;
      fprintf
      (stderr, "using min pial averages = %d\n", min_pial_averages) ;
      nargs++ ;
    }
  }
  else if (!stricmp(option, "wa"))
  {
    max_white_averages = atoi(argv[2]) ;
    fprintf(stderr, "using max white averages = %d\n", max_white_averages) ;
    nargs = 1 ;
    if (isdigit(*argv[3]))
    {
      min_white_averages = atoi(argv[3]) ;
      fprintf
      (stderr, "using min white averages = %d\n", min_white_averages) ;
      nargs++ ;
    }
  }
  else if (!stricmp(option, "add"))
  {
    add = 1 ;
    fprintf(stderr, "adding vertices to tessellation during deformation.\n");
  }
  else if (!stricmp(option, "max") || !stricmp(option, "max_thickness"))
  {
    max_thickness = atof(argv[2]) ;
    nargs = 1 ;
    printf("using max_thickness = %2.1f\n", max_thickness) ;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
      Gx = atoi(argv[2]) ;
      Gy = atoi(argv[3]) ;
      Gz = atoi(argv[4]) ;
      nargs = 3 ;
      printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "v6-cbv")){
    UseV6CBV = 1;
    printf("Using v6 CBV\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "save-res")){
    SaveResidual = 1;
    printf("Saving residual\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "no-save-res")){
    SaveResidual = 0;
    printf("Not saving residual\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "save-target")){
    SaveTarget = 1;
    printf("Saving target surface\n");
    nargs = 0 ;
  }
  else if (!stricmp(option, "no-save-target")){
    SaveTarget = 0;
    printf("Not saving target surface\n");
    nargs = 0 ;
  }
  else if(!stricmp(option, "debug-vertex")){
    // 1=vno 2=surf 3=mri 4=aseg 5=which 6=insidehi 7=blow 8=bhi 9=olow 10=ohi 11=sigma
    Gdiag_no = atoi(argv[2]) ;
    MRI_SURFACE *mris = MRISread(argv[3]);
    if(mris==NULL) exit(1);
    MRIScomputeMetricProperties(mris); /* surface normals */
    MRISsaveVertexPositions(mris, ORIG_VERTICES) ;
    MRI *mri = MRIread(argv[4]);
    if(mri==NULL) exit(1);
    MRI *mriaseg = NULL;
    if(strcmp(argv[5],"noaseg")!=0){
      mriaseg = MRIread(argv[5]);
      if(mriaseg==NULL) exit(1);
    }
    int which          = atoi(argv[6]); //1=white, 2=pial
    double inside_hi   = atof(argv[7]);
    double border_low  = atof(argv[8]);
    double border_hi   = atof(argv[9]);
    double outside_low = atof(argv[10]);
    double outside_hi  = atof(argv[11]);
    double sigma       = atof(argv[12]);
    MRIScomputeBorderValues(mris,mri,NULL,inside_hi,border_hi,border_low,outside_low,outside_hi,sigma,
			    2*max_thickness,NULL,which,NULL,0,0,mriaseg,Gdiag_no,Gdiag_no+1);
    exit(1);
  }
  else if (!stricmp(option, "openmp")) 
  {
    char str[STRLEN] ;
    sprintf(str, "OMP_NUM_THREADS=%d", atoi(argv[2]));
    putenv(str) ;
#ifdef HAVE_OPENMP
    omp_set_num_threads(atoi(argv[2]));
#else
    fprintf(stderr, "Warning - built without openmp support\n");
#endif
    nargs = 1 ;
    fprintf(stderr, "Setting %s\n", str) ;
  }
  else if (!stricmp(option, "mgz"))
  {
    MGZ = 1;
    printf("INFO: assuming MGZ format for volumes.\n");
 } else if (!stricmp(option, "followgradients"))
  {
   	 followGradients= true;
         printf("INFO: Following gradients to compute refine pial surface.\n");
   }
  else switch (toupper(*option))
    {
    case 'D':
      Gx = atoi(argv[2]) ;
      Gy = atoi(argv[3]) ;
      Gz = atoi(argv[4]) ;
      nargs = 3 ;
      printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
      break ;
    case 'S':
      suffix = argv[2] ;
      fprintf(stderr, "using %s as suffix\n", suffix) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      xform_fname = argv[2] ;
      nargs = 1;
      fprintf(stderr, "applying ventricular xform %s\n", xform_fname);
      break ;
    case 'O':
      orig_name = argv[2] ;
      nargs = 1 ;
      fprintf
      (stderr, "reading original vertex positions from %s\n", orig_name);
      break ;
    case 'Q':
      parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
      fprintf(stderr,
              "doing quick (no self-intersection) surface positioning.\n") ;
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      if (parms.momentum >= 1 || parms.momentum < 0)
	ErrorExit(ERROR_BADPARM, "invalid parms.momentum  %f\n",parms.momentum);
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
      break ;
    case 'R':
      l_surf_repulse = atof(argv[2]) ;
      fprintf(stderr, "l_surf_repulse = %2.3f\n", l_surf_repulse) ;
      nargs = 1 ;
      break ;
    case 'B':
      base_dt_scale = atof(argv[2]) ;
      parms.base_dt = base_dt_scale*parms.dt ;
      nargs = 1;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'C':
      create = 0;
      fprintf(stderr,
              "%screating area and curvature files for wm surface...\n",
              create ? "" : "not ") ;
      break ;
    case 'W':
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      printf("write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "niterations = %d\n", parms.niterations) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
      exit(1) ;
      break ;
    }
  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

#include "mris_make_surfaces.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_make_surfaces_help_xml,mris_make_surfaces_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

MRI *
MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_T1, float thresh,
                 int out_label, MRI *mri_dst)
{
  BUFTYPE   *pdst, *pinv_lv, out_val, inv_lv_val ;
  int       width, height, depth, x, y, z,
            ventricle_voxels;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_T1, NULL) ;
  }

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
  */


  ventricle_voxels = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pinv_lv = &MRIvox(mri_inv_lv, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        inv_lv_val = *pinv_lv++ ;
        out_val = 0 ;
        if (inv_lv_val >= thresh)
        {
          ventricle_voxels++ ;
          out_val = out_label ;
        }
        *pdst++ = out_val ;
      }
    }
  }

  return(mri_dst) ;
}

static MRI *
MRIsmoothMasking(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int mask_val,
                 int wsize)
{
  int      width, height, depth, x, y, z, xi, yi, zi, xk, yk, zk, whalf,
           nvox, mean, avg ;
  BUFTYPE  *psrc, *pdst ;

  whalf = (wsize-1) / 2 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
  }

  for ( z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        mean = *psrc++ ;
        nvox = 1 ;

        /* this is a hack to prevent smoothing of non-white values */
        if (MRIvox(mri_mask, x, y, z) > WM_MIN_VAL)
        {
          avg = 0 ;  /* only average if a masked
                        value is close to this one */
          for (zk = -whalf ; zk <= whalf ; zk++)
          {
            zi = mri_mask->zi[z+zk] ;
            for (yk = -whalf ; yk <= whalf ; yk++)
            {
              yi = mri_mask->yi[y+yk] ;
              for (xk = -whalf ; xk <= whalf ; xk++)
              {
                xi = mri_mask->xi[x+xk] ;
                if (MRIvox(mri_mask, xi, yi, zi) == mask_val)
                {
                  avg = 1 ;
                  break ;
                }
              }
              if (avg)
              {
                break ;
              }
            }
            if (avg)
            {
              break ;
            }
          }
          if (avg)
          {
            for (zk = -whalf ; zk <= whalf ; zk++)
            {
              zi = mri_mask->zi[z+zk] ;
              for (yk = -whalf ; yk <= whalf ; yk++)
              {
                yi = mri_mask->yi[y+yk] ;
                for (xk = -whalf ; xk <= whalf ; xk++)
                {
                  xi = mri_mask->xi[x+xk] ;
                  if (MRIvox(mri_mask, xi, yi, zi) >=
                      WM_MIN_VAL)
                  {
                    mean += MRIvox(mri_src, xi, yi, zi) ;
                    nvox++ ;
                  }
                }
              }
            }
          }
        }
        *pdst++ = (BUFTYPE)nint((float)mean/(float)nvox) ;
      }
    }
  }
  return(mri_dst) ;
}

/*!
  \fn int MRISfindExpansionRegions(MRI_SURFACE *mris)
  \brief Sets v->curv=0 unless the given vertex has more than 25% of
  the neighbors whose v->d value is greater than mean+2*std (mean and
  std are the global mean and stddev distances). The dist, eg, is the
  distance along the normal to point of the max gradient. The v->val
  is the max gradient.
  Hidden parameters: 0.25 and mean+2*std
 */
int DeleteMRISfindExpansionRegions(MRI_SURFACE *mris)
{
  int    vno, num, n, num_long, total ;
  float  d, dsq, mean, std, dist ;

  // Compute the mean and stddev of the distance to max gradient
  d = dsq = 0.0f ;
  for (total = num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->ripflag || v->val <= 0)
    {
      continue ;
    }
    num++ ;
    dist = fabs(v->d) ;
    d += dist ;
    dsq += (dist*dist) ;
  }
  mean = d / num ;
  std = sqrt(dsq/num - mean*mean) ;
  printf("mean absolute distance = %2.2f +- %2.2f\n", mean, std); fflush(stdout);

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    v->curv = 0 ;

    if (v->ripflag || v->val <= 0) continue ;

    if (fabs(v->d) < mean+2*std) continue ;

    // Only gets here if distance is not too big

    // Count number of neighbors with big distances
    for (num_long = num = 1, n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->val <= 0 || v->ripflag) continue ;
      if (fabs(vn->d) >= mean+2*std)
        num_long++ ;
      num++ ;
    }

    // If the number of big dist neighbors is greater than 25%
    // Set v->curv = fabs(v->d) and increment total
    if ((float)num_long / (float)num > 0.25)
    {
      v->curv = fabs(v->d) ;
      total++ ; // not used for anything except diagnostic
    }

  } // end loop over vertices

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d vertices more than 2 sigmas from mean.\n", total) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteCurvature(mris, "long") ;

  return(NO_ERROR) ;
}


MRI *
DeletemeMRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm)
{
  int     width, height, depth, x, y, z, nlabeled, nwhite,
          xk, yk, zk, xi, yi, zi;
  BUFTYPE *pwm, val, wm ;
  MRI     *mri_labeled, *mri_tmp ;

  mri_labeled = MRIclone(mri_T1, NULL) ;
  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pwm = &MRIvox(mri_wm, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        wm = *pwm++ ;

        if (x == Gx && y == Gy && z == Gz)  /* T1=127 */
        {
          DiagBreak() ;
        }
        /* not white matter and bright (e.g. eye sockets) */
        if ((wm < WM_MIN_VAL) && (val > 125))
        {
          nwhite = 0 ;
          for (xk = -1 ; xk <= 1 ; xk++)
          {
            xi = mri_T1->xi[x+xk] ;
            for (yk = -1 ; yk <= 1 ; yk++)
            {
              yi = mri_T1->yi[y+yk] ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_T1->zi[z+zk] ;
                if (MRIvox(mri_wm, xi, yi, zi) >= WM_MIN_VAL)
                {
                  nwhite++ ;
                }
              }
            }
          }
#define MIN_WHITE  ((3*3*3-1)/2)
          if (nwhite < MIN_WHITE)
          {
            MRIvox(mri_labeled, x, y, z) = BRIGHT_LABEL ;
          }
        }
      }
    }
  }

  /* find all connected voxels that are above 115 */
  MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL, 10,115);
  MRIclose(mri_labeled, mri_labeled) ;

  /* expand once more to all neighboring voxels that are bright. At
     worst we will erase one voxel of white matter.
  */
  mri_tmp =
    MRIdilateThreshLabel(mri_labeled, mri_T1, NULL, BRIGHT_LABEL,1,100);
  MRIxor(mri_labeled, mri_tmp, mri_tmp, 1, 255) ;
  MRIreplaceValues(mri_tmp, mri_tmp, 1, BRIGHT_BORDER_LABEL) ;
  MRIunion(mri_tmp, mri_labeled, mri_labeled) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_labeled, "label.mgz") ;
  }
  /*    MRIwrite(mri_tmp, "tmp.mgz") ;*/
  nlabeled = MRIvoxelsInLabel(mri_labeled, BRIGHT_LABEL) ;
  fprintf(stdout, "%d bright non-wm voxels segmented.\n", nlabeled) ;

  /* dilate outwards if exactly 0 */
  MRIdilateInvThreshLabel
  (mri_labeled, mri_T1, mri_labeled, BRIGHT_LABEL, 3, 0) ;

  MRIfree(&mri_tmp) ;
  return(mri_labeled) ;
}
static float
check_contrast_direction(MRI_SURFACE *mris,MRI *mri_T1)
{
  int     vno, n ;
  VERTEX  *v ;
  double  x, y, z, xw, yw, zw, val, mean_inside, mean_outside ;

  mean_inside = mean_outside = 0.0 ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag != 0)
    {
      continue ;
    }
    x = v->x+0.5*v->nx ;
    y = v->y+0.5*v->ny ;
    z = v->z+0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_T1, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_T1, xw, yw, zw, &val) ;
    mean_outside += val ;

    x = v->x-0.5*v->nx ;
    y = v->y-0.5*v->ny ;
    z = v->z-0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_T1, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_T1, xw, yw, zw, &val) ;
    mean_inside += val ;
    n++ ;
  }
  mean_inside /= (float)n ;
  mean_outside /= (float)n ;
  printf("mean inside = %2.1f, mean outside = %2.1f\n",
         mean_inside, mean_outside) ;
  return(mean_inside - mean_outside) ;
}


/*!
  \fn int fix_midline()

  \brief Note: "fix" here means to "hold in place", not to "repair
  something that is broken".  It also does more than midline.  

  Sets the v->ripped flag at certain vertices (which then holds them
  in place later on). Does not unrip vertices; if a vertex is ripped,
  then it is not procsssed. All unripped v->marked are set to 0 at the
  beginning; all v->marked (ripped or otherwise) are set to 0 at end
  of the function. All unripped v->marked2 are set to 0 at the
  begining of the function. The ripped vertices will have
  v->val=intensity at the vertex and v->d=0.

  It finds the midline mostly by finding subcortical GM, ventricular
  CSF, or CC in the aseg.presurf within a few mm of the vertex. These
  vertices are marked (v->marked=1).  It will also mark some vertices
  that are near putamen or a lesion. All ripped vertices are marked.

  The marks are spatially filtered by eliminating clusters that have
  fewer than 5 vertices or, if annotation is present, that are less
  than 60% "unknown". This last step unmarks the putamen marks done
  before (bug).

  If there is an annotation, it will only make a difference in two
  places: (1) for entorhinal (where vertices are never marked/ripped),
  and (2) in some areas of "unknown" cortex (see above). #1 probably
  does not have much of an impact because because there is usually
  some WM between entorhinal and hippo. #2 probably removes a few
  vertices (and causes the putamen bug). 

  v->marked2 - affects MRIScortexLabel(), forces labeling as cortex. Also
  influences this program the second time around in that vertices with
  v->marked2=1 are not reprocessed. 

  v->val2=1 for Lesions when fitting GRAY_CSF
  #FIX
*/
static int fix_midline(MRI_SURFACE *mris, MRI *mri_aseg, MRI *mri_brain, char *hemi,
            int which, int fix_mtl)
{
  int vno, label, contra_wm_label, nvox=0, total_vox=0, adjacent=0;
  int wm_label, gm_label, nlabels, n, index, annotation, entorhinal_index ;
  VERTEX   *v ;
  double   xv, yv, zv, val, xs, ys, zs, d, nx, ny, nz ;
  LABEL    **labels ;
  int nmarked, nmarked2, nripped;

  printf("Entering: fix_midline(): inhibiting deformation at non-cortical midline structures...\n") ;
  printf("  which=%d, fix_mtl=%d\n",which, fix_mtl);
  printf("#FML0# fix_midline(): nripped=%d\n",MRIScountRipped(mris));
  fflush(stdout);

  if (stricmp(hemi, "lh") == 0){
    contra_wm_label = Right_Cerebral_White_Matter ;
    wm_label = Left_Cerebral_White_Matter ;
    gm_label = Left_Cerebral_Cortex ;
  }
  else{
    contra_wm_label = Left_Cerebral_White_Matter ;
    wm_label = Right_Cerebral_White_Matter ;
    gm_label = Right_Cerebral_Cortex ;
  }

  // Clear the deck
  MRISclearMarks(mris) ;  // Sets v->marked=0  for all unripped vertices
  MRISclearMark2s(mris) ; // Sets v->marked2=0 for all unripped vertices

  if (mris->ct)
    CTABfindName(mris->ct, "entorhinal", &entorhinal_index);
  else
    entorhinal_index = -1 ;

  // Loop over vertices ==================================
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

    if(v->ripflag || v->marked2 > 0){
      // It was ripped previously - it should still be excluded, but
      // still mark it.  Note: all unripped marked2 are set to zero
      // above, so if(v->marked2>0) is meaningless.
      v->marked = 1 ;
      continue ;
    }

    if (mris->ct)
      CTABfindAnnotation(mris->ct, v->annotation, &index);

    if (vno == Gdiag_no ) {
      printf("vno %d: annotation %x, index %d, EC %d\n", vno, v->annotation, index, entorhinal_index) ;
      DiagBreak() ;
    }

    // don't freeze vertices that are in EC and very close to hippocampus
    // How do you know from this that the vertex is close to hippo?
    if (mris->ct && index == entorhinal_index)  
      continue ;

    // search outwards. Hidden parameters 2mm and 0.5mm
    for (d = 0 ; d <= 2 ; d += 0.5) {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;

      // Sample the aseg at this distance
      MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
      if (vno == Gdiag_no){
	printf("vno %d: dist %2.2f - %s (%d)\n", vno, d, cma_label_to_name(label), label) ;
	DiagBreak() ;
      }

      if (d > 0 && label == gm_label)
	break ;  // cortical gray matter of this hemisphere - doesn't matter what is outside it

      // If the vertex lands on a certain label, then v->marked=1 and
      // v->d=0 and v->val=val at this location. The val may be
      // changed at later distances but marked and d will not even if
      // the point lands in an acceptable label. marked2 may be set to
      // 1 if the label is a lesion or WMSA.  Labels excluded from
      // thist list include ipsilateral cortex and WM.
      if(label == contra_wm_label ||
          label == Left_Lateral_Ventricle ||
          label == Left_vessel ||
          label == Right_vessel ||
          label == Optic_Chiasm ||
          label == Left_choroid_plexus ||
          label == Right_choroid_plexus ||
          label == Third_Ventricle ||
          label == Right_Lateral_Ventricle ||
          ((label == Left_Accumbens_area ||
            label == Right_Accumbens_area) &&  // only for gray/white
           which == GRAY_WHITE) ||
          ((label == Left_Lesion ||
            label == Right_Lesion ||
            label == WM_hypointensities ||
            label == Left_WM_hypointensities ||
            label == Right_non_WM_hypointensities ||
            label == Left_non_WM_hypointensities ||
            label == Right_WM_hypointensities) &&  // only for gray/white
           which == GRAY_WHITE) ||
          label == Left_Caudate ||
          label == Right_Caudate ||
          label == Left_Pallidum ||
          label == Right_Pallidum ||
          IS_CC(label) ||
          ((IS_HIPPO(label)  || IS_AMYGDALA(label)) && fix_mtl)  ||
          label == Right_Thalamus_Proper ||
          label == Left_Thalamus_Proper ||
          label == Brain_Stem ||
          label == Left_VentralDC ||
          label == Right_VentralDC)
      {
	// below are labels where the intensities aren't useful, so just freeze surface there
	if (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label))
	  v->marked2 = 1 ; // afects the cortex.label
        if (label == Left_Putamen || label == Right_Putamen)
          DiagBreak() ;
        if (vno == Gdiag_no)
          DiagBreak() ;
	// Sample brain intensity at vertex (not at distance from vertex)
        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ; // not redundant
        MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
        v->val = val ;
        v->d = 0 ;
        v->marked = 1 ;
      }
    } // end loop over distance

    if (vno == Gdiag_no)
      DiagBreak() ;

    // Compute the normal to the edge of the label if this is putamen
    MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
    label = nint(val) ;
    if (IS_PUTAMEN(label)) {
      // 3=whalf in voxels, hidden parameter
      // 1=use_abs, hidden parameter
      // nx, ny, nz are in vox coords making orientation a hidden parameter
      // Are nx,ny,nz ever used?
      if(vno == Gdiag_no) printf("vno=%d, putamen label: running compute_label_normal\n",vno);
      compute_label_normal(mri_aseg, xv, yv, zv, label, 3, &nx, &ny, &nz, 1) ;
    }
    else
    {
      nx = ny = nz = 0 ;
    }

    /*
      for gray/white surface, if we are in insula, don't want to let the
      surface diverge into the putamen.
    */
    if (which == GRAY_WHITE) {

      // search inwards. Hidden parameters 2mm and 0.5mm
      for (d = 0 ; d <= 2 ; d += 0.5) {

	// Sample the aseg at this distance
        xs = v->x - d*v->nx ;
        ys = v->y - d*v->ny ;
        zs = v->z - d*v->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
	if (vno == Gdiag_no){
	  printf("vno %d: dist %2.2f - %s (%d)\n", vno, -d, cma_label_to_name(label), label) ;
	  DiagBreak() ;
	}

	// these are labels where the intensities aren't useful, so just freeze surface there
	// hidden parameter d<1mm. v->marked2 affects the cortical label
	if (d < 1 && (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label)))
	  v->marked = v->marked2 = 1 ;

        if (IS_PUTAMEN(label) || IS_ACCUMBENS(label) || IS_CLAUSTRUM(label)){
	  // 3=whalf in voxels, hidden parameter. Voxel resolution is
	  //   a hidden parameter since whalf is in voxels
	  // 1=use_abs, hidden parameter, but abs used regardless below
	  // nx, ny, nz are in vox coords making orientation a hidden parameter
          compute_label_normal(mri_aseg, xv, yv, zv, label, 3,&nx, &ny, &nz, 1) ;
          if(fabs(nx) > fabs(ny) && fabs(nx) > fabs(nz))  {
	    // edge is mostly oriented in column direction
            if(vno == Gdiag_no)
              DiagBreak() ;
	    if(IS_ACCUMBENS(label)){
	      // Same as put and claust but val not set
	      v->d = 0 ;
	      v->marked = 1 ;
	      if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
		printf("marking vertex %d as adjacent to accumbens on midline\n", vno);
	    }
	    else {
	      // Sample brain intensity at the vertex (not the distance)
	      MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ; // not redundant
	      MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
	      v->val = val ;
	      v->d = 0 ;
	      v->marked = 1 ;
	      if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
		printf("marking vertex %d as adjacent to putamen/claustrum in insula\n", vno);
	    }
	    break ;
          }
        } // if putamen, accumbens, or claustrum

      } // loop over distance
    }// if GRAY_WHITE


    // search inwards 2mm step 0.5mm (hidden parameters)
    for (d = 0 ; d <= 2 ; d += 0.5) {

      // Sample the aseg at this distance
      xs = v->x - d*v->nx ;
      ys = v->y - d*v->ny ;
      zs = v->z - d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
 
      // Found ispilateral GM or WM next to surface (1.1mm hidden
      // parameter), so don't go any deeper.  This is probably why the
      // surface near entorhinal is not affected (or not much) even
      // when no annotation is passed.
      if(d < 1.1 && (label == wm_label || label == gm_label))
        break ;  

      if ((which == GRAY_CSF) && (d < 1) && (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label)))
	v->val2 = 1 ;

      if ((label == contra_wm_label ||
           label == Left_vessel ||
           label == Right_vessel ||
           label == Optic_Chiasm ||
           label == Left_choroid_plexus ||
           label == Right_choroid_plexus ||
           label == Left_Lateral_Ventricle ||
           label == Third_Ventricle ||
           label == Right_Lateral_Ventricle ||
           ((label == Left_Accumbens_area ||
             label == Right_Accumbens_area) &&
            which == GRAY_WHITE)||
           label == Left_Caudate ||
           label == Right_Caudate ||
           label == Left_Pallidum ||
           IS_CC(label) ||
           label == Right_Thalamus_Proper ||
           label == Left_Thalamus_Proper ||
           label == Right_Pallidum ||
           label == Brain_Stem ||
           label == Left_VentralDC ||
           label == Right_VentralDC) ||
          // putamen can be adjacent to insula in aseg for pial
          (which == GRAY_WHITE && (d < 1.1) &&
           (label == Left_Putamen || label == Right_Putamen)))

      {
        if (label == Left_Putamen || label == Right_Putamen)
          DiagBreak() ;

        if((label == Left_Lateral_Ventricle || label == Right_Lateral_Ventricle) && d > 1)  
        {
	  // In calcarine ventricle can be pretty close to wm surface
	  // but could affect any vertex that is near L/R Lat Vent.
	  // d>1mm hidden parameter
          break ;
        }
        if (vno == Gdiag_no)
          DiagBreak() ;

        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
        MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
        v->val = val ;
        v->d = 0 ;
        v->marked = 1 ;
      }

    } // end loop over inward distance

    /* Now check for putamen superior to this vertex. If there's a lot
       of it there, then we are in basal forebrain and not cortex. */
    if (which == GRAY_WHITE) {
      // Project in the row direction (not the normal direction) 10mm step 0.5 mm (hidden par)
      // For a conformed volume, row direction = superior/inferior
      adjacent = total_vox = nvox = 0;
      for (d = 0 ; d <= 10 ; d += 0.5, total_vox++){
	// Sample superiorly. This leads to a dependence on how the
	// head is oriented with respect to the voxel coordinates
	// (hidden parameter) since the surface coordinates are
	// aligned with the voxel coords.  It also leads to a
	// dependence on voxel axis orientation since the row
	// direction is not always superior/inf.
        xs = v->x ;
        ys = v->y ;
        zs = v->z + d ;  // z is along the row direction
        MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
        if (label == Left_Putamen || label == Right_Putamen) {
          nvox++ ;
          if (d < 1.5) // d<1.5mm hidden parameter
            adjacent = 1 ;  // right next to putamen
        }
      } // loop over distance

      // if more than 50% of the samples are putamen and vertex is within 1.5mm of putamen
      // and edge of putamen is mostly in the row direction
      if (adjacent && (double)nvox/(double)total_vox > 0.5)
      {
        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
	// 3=whalf in voxels, hidden parameter. Voxel resolution is
	//   a hidden parameter since whalf is in voxels
	// 1=use_abs, hidden parameter, but abs used regardless below
	// nx, ny, nz are in vox coords making orientation a hidden parameter
        compute_label_normal(mri_aseg, xv, yv, zv, label, 3, &nx, &ny, &nz, 1) ;
        if (ny > 0 && fabs(ny) > fabs(nx) &&  fabs(ny) > fabs(nz))  {
	  // Normal is mostly in the row direction. Note that ny will always be >=0 since use_abs=1
          if (vno == Gdiag_no)
            DiagBreak() ;
          MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
          MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
          v->val = val ;
          v->d = 0 ;
          v->marked = 1 ;
        }
      }

    } // end if GRAY_WHITE

  } // end loop over vertices =================================================

  if (Gdiag_no >= 0) {
    int index ;
    v = &mris->vertices[Gdiag_no] ;
    if (mris->ct)
      CTABfindAnnotation(mris->ct, v->annotation, &index);
    else
      index = -1 ;
    printf("v %d: ripflag = %d before connected components, annot %d (%d)\n",
           Gdiag_no, v->marked, v->annotation, index) ;
    if (v->marked == 0)
      DiagBreak() ;
    else
      DiagBreak() ;
  }

  // Dilate and erode the marked vertices
  MRISdilateMarked(mris, 3) ;
  MRISerodeMarked(mris, 3) ;

  // Get a list of marked clusters (ie, connected components) whose
  // area is greater than 1mm2 (hidden parameter)
  MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
  
  if (Gdiag_no > 0)
    printf("v %d: ripflag = %d after morphology\n", Gdiag_no, mris->vertices[Gdiag_no].marked) ;

  // Go through all the clusters
  for (n = 0 ; n < nlabels ; n++)  {

    // If there are fewer than 5 points in the cluster, then discard it
    if(labels[n]->n_points < 5) {
      int i, threadno = 0 ;
      #ifdef HAVE_OPENMP 
      threadno = omp_get_thread_num();
      #endif
      printf("removing %d vertices from ripped group in thread:%d\n",labels[n]->n_points,threadno); 
      for (i = 0 ; i < labels[n]->n_points ; i++) 
        mris->vertices[labels[n]->lv[i].vno].marked = 0 ;
    }

    // If there is an annoation that includes "unknown", then ...
    if (mris->ct && CTABfindName(mris->ct, "unknown", &index) == NO_ERROR) {
      double pct_unknown;
      int    i ;
      VERTEX *v ;
      CTABannotationAtIndex(mris->ct, index, &annotation) ;

      // Compute the fraction of this cluster that is in the unknown annotation
      for (pct_unknown = 0.0, i = 0 ; i < labels[n]->n_points ; i++) {
	v = &mris->vertices[labels[n]->lv[i].vno] ;
        if ((v->annotation == annotation || v->annotation == 0) &&
	    (v->marked2 == 0))  // will be 1 if a lesion or WMSA (keep frozen if so)
        {
          pct_unknown = pct_unknown + 1 ;
        }
      }
      pct_unknown /= (double)labels[n]->n_points ;

      // If this fraction is < 0.6 (hidden parameter), unmark all vertices in the cluster
      // Won't this undo the marks near putamen or lesion? Answer: yes, this is a bug
      // that makes the surface worse when the annotation is loaded.
      if (pct_unknown < .6) {
        printf("deleting segment %d with %d points - only %2.2f%% unknown, v=%d\n",
               n,labels[n]->n_points,100*pct_unknown, labels[n]->lv[0].vno) ;
        for (i = 0 ; i < labels[n]->n_points ; i++) {
	  v = &mris->vertices[labels[n]->lv[i].vno] ;
	  if (v->marked2 == 0){
	    mris->vertices[labels[n]->lv[i].vno].marked = 0 ;
	    if (labels[n]->lv[i].vno  == Gdiag_no)
	      printf("removing mark from v %d due to non-unknown aparc\n",Gdiag_no) ;
	  }
        }
      }
    }

    LabelFree(&labels[n]) ;

  } // end loop over clusters
  free(labels) ;

  if (Gdiag_no > 0)
    printf("v %d: ripflag = %d after connected components\n",
           Gdiag_no, mris->vertices[Gdiag_no].marked) ;

  // This is where the effects of this routine are instantiated
  // Sets v->ripflag=1 if v->marked==1  Note: does not unrip any vertices.
  MRISripMarked(mris) ; 

  // Count up the final number of marked vertices
  nmarked=0, nmarked2=0, nripped=0;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    if(mris->vertices[vno].marked)  nmarked++;
    if(mris->vertices[vno].marked2) nmarked2++;
    if(mris->vertices[vno].ripflag) nripped++;
  }
  printf("#FML# fix_midline(): nmarked=%d, nmarked2=%d, nripped=%d\n",nmarked,nmarked2,nripped);
  fflush(stdout);

  // Sets all marked vertices to 0, even the ripped ones
  MRISsetAllMarks(mris, 0) ;

  return(NO_ERROR) ;
}


static double
compute_brain_thresh(MRI_SURFACE *mris, MRI *mri_ratio, float nstd)
{
  double    val, xs, ys, zs, xv, yv, zv, d, mean, std ;
  int       vno, num ;
  VERTEX    *v ;
  double    thresh ;
  FILE      *logfp = NULL ;
  MRI       *mri_filled ;

  mri_filled = MRIclone(mri_ratio, NULL) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    logfp = fopen("gm.plt", "w") ;
  }
  mean = std = 0.0 ;
  num = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    MRISvertexToVoxel(mris, v, mri_ratio, &xv, &yv, &zv) ;
    MRIsampleVolume(mri_ratio, xv, yv, zv, &val) ;
    for (d = .5 ; d <= 1.0 ; d += 0.5)
    {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_ratio, xs, ys, zs, &xv, &yv, &zv);
      if (MRIgetVoxVal(mri_filled, nint(xv), nint(yv), nint(zv), 0) > 0)
	continue ;  // already visited this voxel
      MRIsetVoxVal(mri_filled, nint(xv), nint(yv), nint(zv), 0, 1) ;
      val = MRIgetVoxVal(mri_ratio, nint(xv), nint(yv), nint(zv), 0) ;
      if (val < 0)
        continue ;

      if (logfp)
        fprintf(logfp, "%f\n", val) ;

      mean += val ;
      std += val*val ;
      num++ ;
    }
  }

  if (logfp)
  {
    fclose(logfp) ;
  }
  mean /= num ;
  std = sqrt(std/num - mean*mean) ;
  thresh = mean+nstd*std ;
  MRIfree(&mri_filled) ;
  return(thresh) ;
}
#undef SAMPLE_DIST
#define SAMPLE_DIST .1
#define PERCENTILE   0.9

#define PIAL_SURFACE_SMOOTHING_MM 0.5


#define WM_BORDER    1
#define WM_INTERIOR  2
#define WM_NBR       3
#define GRAD_BINS 256


static int
compute_white_target_locations(MRI_SURFACE *mris,
			       MRI *mri_T2,
			       MRI *mri_aseg,
			       float nstd_below,
			       float nstd_above,
			       int contrast_type, float T2_min_inside, float T2_max_inside, 
			       int below_set, int above_set, double wsigma, double max_out)
{
  int       vno, vnum, outer_nbhd_size = 11, inner_nbhd_size = 3, n ;
  int       num_in, num_out ;
  VERTEX    *v, *vn ;
  double    min_grad, max_grad, min_inside, max_inside, min_outside, max_outside, sample_dist ;
  double    nx, ny, nz, d, thickness, grad_scale = 2 ;
  double    grad_val, inside_val, outside_val, xv, yv, zv, xs, ys, zs  ;
  HISTOGRAM *h_inside, *h_outside, *h_grad, *hs ;

  int const vlistCapacity = mris->nvertices;
  int * vlist = (int *)malloc(vlistCapacity*sizeof(vlist[0])) ;
  int * hops  = (int *)malloc(vlistCapacity*sizeof(hops [0])) ;
  if (!vlist|| !hops)
    ErrorExit(ERROR_NOFILE, "compute_white_target_locations: could not allocate vlist") ;

  sample_dist = MIN(SAMPLE_DIST, mri_T2->xsize/2) ;

  MRISsaveVertexPositions(mris, TMP2_VERTICES) ;

  // first compute bounds on gradient and inner and outer values
  min_outside = min_inside = min_grad = 1e10 ;
  max_outside = max_inside = max_grad = -1e-10 ;
  d = 2*sample_dist ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = v->pialx - v->whitex ; ny = v->pialy - v->whitey ; nz = v->pialz - v->whitez ;
    thickness = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
    if (FZERO(thickness))  // pial and white in same place - no cortex here
      continue ;

    MRISvertexToVoxel(mris, v, mri_T2, &xv, &yv, &zv) ;
    MRISvertexNormalToVoxelScaled(mris, v, mri_T2, &nx, &ny, &nz) ;
    MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, 2) ;
    xs = v->x - d*v->nx ; ys = v->y - d*v->ny ; zs = v->z - d*v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
    MRIsampleVolume(mri_T2, xv, yv, zv, &inside_val) ;
    xs = v->x + d*v->nx ; ys = v->y + d*v->ny ; zs = v->z + d*v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
    MRIsampleVolume(mri_T2, xv, yv, zv, &outside_val) ;

    if (inside_val > 250)
      DiagBreak() ;

    if (outside_val < min_outside)
      min_outside = outside_val ;
    if (outside_val > max_outside)
      max_outside = outside_val ;

    if (grad_val < min_grad)
      min_grad = grad_val ;
    if (grad_val > max_grad)
      max_grad = grad_val ;

    if (inside_val < min_inside)
      min_inside = inside_val ;
    if (inside_val > max_inside)
      max_inside = inside_val ;
  }

  h_inside = HISTOalloc((int)ceil(max_inside)+1) ;
  HISTOinit(h_inside, h_inside->nbins, 0, ceil(max_inside)) ;
  h_outside = HISTOalloc((int)ceil(max_outside)+1) ;
  HISTOinit(h_outside, h_outside->nbins, 0, ceil(max_outside)) ;
  h_grad = HISTOalloc(GRAD_BINS) ;
  HISTOinit(h_grad, h_grad->nbins, min_grad, ceil(max_grad)) ;

  num_in = num_out = 0 ;
  MRISsaveVertexPositions(mris, TARGET_VERTICES) ;
  
  MRISclearMarks(mris) ;    // used by the old implementation of MRISfindNeighborsAtVertex
  MRISclearMark2s(mris) ;
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    double dist_from_current_white, pin, pout, best_dist, best_dist_prelim ;
    double pin0, pout0, mean_gm, mean_wm, sigma_wm, sigma_gm, grad_outside, grad_inside ;
    int  found ;

    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      v->targx = v->x ; v->targy = v->y ; v->targz = v->z ;
      continue ;
    }

    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = v->pialx - v->whitex ; ny = v->pialy - v->whitey ; nz = v->pialz - v->whitez ;
    thickness = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
    if (FZERO(thickness))  // pial and white in same place - no cortex here
      continue ;

    d = 2*sample_dist ;
    vnum = MRISfindNeighborsAtVertex(mris, vno, outer_nbhd_size, vlistCapacity, vlist, hops);
    for (n = 0 ; n < vnum ; n++)
    {
      vn = &mris->vertices[vlist[n]] ;
      if (vn->ripflag || hops[n] <= inner_nbhd_size)
	continue ;

      MRISvertexToVoxel(mris, vn, mri_T2, &xv, &yv, &zv) ;
      MRISvertexNormalToVoxelScaled(mris, vn, mri_T2, &nx, &ny, &nz) ;
      MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, grad_scale) ;
      if (grad_val < -2)
	DiagBreak() ;

      // sample further into interior of WM to avoid biasing distribution by misplaced boundaries
      xs = vn->x - 4*d*vn->nx ; ys = vn->y - 4*d*vn->ny ; zs = vn->z - 4*d*vn->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolume(mri_T2, xv, yv, zv, &inside_val) ;

      xs = vn->x + 2*d*vn->nx ; ys = vn->y + 2*d*vn->ny ; zs = vn->z + 2*d*vn->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolume(mri_T2, xv, yv, zv, &outside_val) ;

      // sample inside and outside and make sure we are at a local max of the gradient. If not, don't trust
      xs = vn->x - 2*d*vn->nx ; ys = vn->y - 2*d*vn->ny ; zs = vn->z - 2*d*vn->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_inside, grad_scale) ;
      xs = vn->x + 2*d*vn->nx ; ys = vn->y + 2*d*vn->ny ; zs = vn->z + 2*d*vn->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_outside, grad_scale) ;
      if ((grad_val > 1) && (grad_val > grad_inside && grad_val > grad_outside)) // only regions we are confident of
      {
	double dist ;

	HISTOaddSample(h_outside, outside_val, 0, 0) ;
	HISTOaddSample(h_grad, grad_val, 0, 0) ;

	// sample a number of points in the interior to get a good histogram
	for (dist = -2*sample_dist ; dist >= -8*d ; dist -= sample_dist)
	{
	  xs = vn->x + dist*vn->nx ; ys = vn->y + dist*vn->ny ; zs = vn->z + dist*vn->nz ;
	  MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
	  MRIsampleVolume(mri_T2, xv, yv, zv, &inside_val) ;
	  MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, grad_scale) ;
	  if (grad_val < 0)
	    break ;
	  HISTOaddSample(h_inside, inside_val, 0, 0) ;
	}
      }
    }

    HISTOfillHoles(h_grad) ; HISTOfillHoles(h_inside) ; HISTOfillHoles(h_outside) ;
    hs = HISTOsmooth(h_inside, NULL, 2) ; HISTOfree(&h_inside) ; h_inside = hs ;
    hs = HISTOsmooth(h_outside, NULL, 2) ; HISTOfree(&h_outside) ; h_outside = hs ;
    hs = HISTOsmooth(h_grad, NULL, 2) ; HISTOfree(&h_grad) ; h_grad = hs ;
    HISTOmakePDF(h_grad, h_grad) ; HISTOmakePDF(h_inside, h_inside) ; HISTOmakePDF(h_outside, h_outside) ;

    d = 2*sample_dist ;   
    MRISvertexNormalToVoxelScaled(mris, v, mri_T2, &nx, &ny, &nz) ;
    MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, grad_scale) ;
    xs = v->x - d*v->nx ; ys = v->y - d*v->ny ; zs = v->z - d*v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
    MRIsampleVolume(mri_T2, xv, yv, zv, &inside_val) ;
    xs = v->x + d*v->nx ; ys = v->y + d*v->ny ; zs = v->z + d*v->nz ;
    MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
    MRIsampleVolume(mri_T2, xv, yv, zv, &outside_val) ;
    HISTOrobustGaussianFit(h_inside, .9, &mean_wm, &sigma_wm) ;
    HISTOrobustGaussianFit(h_outside, .9, &mean_gm, &sigma_gm) ;
    sigma_gm = sigma_wm = 1 ; // don't trust variances
#define MIN_GM_DIST_FROM_WM_PEAK 15
    if (abs(mean_wm - DEFAULT_DESIRED_WHITE_MATTER_VALUE) > 10)
      mean_wm = DEFAULT_DESIRED_WHITE_MATTER_VALUE ;
    if (mean_gm - mean_wm < MIN_GM_DIST_FROM_WM_PEAK)
      mean_gm = DEFAULT_DESIRED_WHITE_MATTER_VALUE + MIN_GM_DIST_FROM_WM_PEAK ;

    if (FZERO(inside_val) || FZERO(outside_val))
      continue ;  // ventricles are near edge of volume - don't deform

    pin0 = HISTOgetCount(h_inside, inside_val) ;
    pout0 = HISTOgetCount(h_outside, outside_val) ;
    if (vno == Gdiag_no)
    {
      printf("vno %d: vertex in=%2.1f (%2.3f), out=%2.1f (%2.3f), grad=%2.1f\n", vno, inside_val, pin0, outside_val, pout0, grad_val) ;
      HISTOplot(h_inside, "hi.plt") ;
      HISTOplot(h_outside, "ho.plt") ;
      printf(" wm = %2.1f +- %2.1f, gm = %2.1f +- %2.1f\n", mean_wm, sigma_wm, mean_gm, sigma_gm) ;
    }

    best_dist = 0 ; max_grad = 0 ;
    for (dist_from_current_white = -max_out ; dist_from_current_white <= max_out ; dist_from_current_white += sample_dist)
    {
      double xs0, ys0, zs0, gm_dist_in, gm_dist_out, wm_dist_in, wm_dist_out ;

      xs0 = v->x + dist_from_current_white*v->nx ; 
      ys0 = v->y + dist_from_current_white*v->ny ; 
      zs0 = v->z + dist_from_current_white*v->nz ;

      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs0, ys0, zs0, &xv, &yv, &zv);
      MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, grad_scale) ;
      if (grad_val < 0)
	continue ; // don't consider locations with gradient pointing the wrong way

      xs = xs0 - d*v->nx ; ys = ys0 - d*v->ny ; zs = zs0 - d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolume(mri_T2, xv, yv, zv, &inside_val) ;
      xs = xs0 + d*v->nx ; ys = ys0 + d*v->ny ; zs = zs0 + d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolume(mri_T2, xv, yv, zv, &outside_val) ;

      pin = HISTOgetCount(h_inside, inside_val) ;
      pout = HISTOgetCount(h_outside, outside_val) ;
      wm_dist_in = abs(mean_wm-inside_val) / sigma_wm ;
      gm_dist_in = abs(mean_gm-inside_val) / sigma_gm ;  
      wm_dist_out = abs(mean_wm-outside_val) / sigma_wm ;
      gm_dist_out = abs(mean_gm-outside_val) / sigma_gm ;
      if (vno == Gdiag_no)
      {
	MRISsurfaceRASToVoxelCached(mris, mri_T2, xs0, ys0, zs0, &xv, &yv, &zv);
	printf("vno %d: dist %2.2f, vox=(%2.1f %2.1f %2.1f),  in=%2.1f (%2.3f, WM/GM=%2.1f/%2.1f), out=%2.1f (%2.3f, WM/GM=%2.1f/%2.1f), grad=%2.1f\n", vno, dist_from_current_white, xv, yv, zv, inside_val, pin, wm_dist_in, gm_dist_in, outside_val, pout, wm_dist_out, gm_dist_out, grad_val) ;
      }
      
      if (dist_from_current_white > 0 && (wm_dist_in > 6 || gm_dist_in < wm_dist_in))
      {
	if (vno == Gdiag_no)
	  printf("terminating search due to unlikely WM found (interior wm_dist=%2.1f, gm_dist=%2.1f)\n", wm_dist_in,gm_dist_in) ;
	break ;
      }
      wm_dist_out = abs(mean_wm-outside_val) / sigma_wm ;  // now compute dist to exterior val
      gm_dist_out = abs(mean_gm-outside_val) / sigma_gm ;  // now compute dist to exterior val
      if (dist_from_current_white < 0 && gm_dist_out > 4 && !FZERO(best_dist))  // found illegal interior val
      {
	if (vno == Gdiag_no)
	  printf("resetting search due to non-allowable gm val %2.1f (dist=%2.1f)\n", outside_val, gm_dist_out);
	max_grad = 0 ; best_dist = 0 ;
      }
      if (grad_val > max_grad)    // this is the largest gradient found yet, so if vals are reasonable
      {
	if ((wm_dist_in < 1 && gm_dist_out < 1) ||   // either it really looks correct
	    (
	      (wm_dist_in < 2 && gm_dist_out < 2) && // or it  almost looks correct and is better than current
	      (pin0*pout0 < pin*pout)
	      )
	    ||
	    ((wm_dist_in < gm_dist_in) && (wm_dist_out > gm_dist_out)) // inside more likely wm and out more likely gm
	  )
	{
	  found = 1 ;
	  max_grad = grad_val ;
	  best_dist = dist_from_current_white ;
	}
      }
    }

    v->marked2 = found ;
    best_dist_prelim = best_dist;
    if (vno == Gdiag_no)
    {
      printf("vno %d: found %d, prelim best distance %2.2f, searching for local grad max\n", vno, found, best_dist_prelim) ;
      DiagBreak() ;
    }
    // tweak location to be at peak of gradient mag
    for (dist_from_current_white = best_dist_prelim-2*sample_dist ; dist_from_current_white <= best_dist_prelim + 2*sample_dist ; dist_from_current_white += sample_dist)
    {
      double xs0, ys0, zs0 ;

      xs0 = v->x + dist_from_current_white*v->nx ; 
      ys0 = v->y + dist_from_current_white*v->ny ; 
      zs0 = v->z + dist_from_current_white*v->nz ;

      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs0, ys0, zs0, &xv, &yv, &zv);
      MRISvertexNormalToVoxelScaled(mris, v, mri_T2, &nx, &ny, &nz) ;
      MRIsampleVolumeDerivativeScale(mri_T2, xv, yv, zv, nx, ny, nz, &grad_val, grad_scale) ;
      if (grad_val > max_grad)
      {
	if (vno == Gdiag_no)
	{
	  printf("vno %d: updating max grad dist to %2.2f, vox=(%2.1f %2.1f %2.1f),  grad=%2.1f\n", 
		 vno, dist_from_current_white, xv, yv, zv, grad_val) ;
	  DiagBreak() ;
	}
	best_dist = dist_from_current_white ;
	max_grad = grad_val ;
      }
    }

    if (vno == Gdiag_no)
    {
      printf("vno %d: best distance %2.2f\n", vno, best_dist) ;
      DiagBreak() ;
    }
    v->targx = v->x+best_dist*v->nx ; v->targy = v->y+best_dist*v->ny ; v->targz = v->z+best_dist*v->nz ;
    v->d = best_dist ;

    // reset things for next iteration
    HISTOclearCounts(h_inside, h_inside) ;
    HISTOclearCounts(h_outside, h_outside) ;
    HISTOclearCounts(h_grad, h_grad) ;
    
    if (best_dist > 0)
      num_out++ ;
    else if (best_dist < 0)
      num_in++ ;
  }

  MRIScopyMarked2ToMarked(mris) ;
  MRISwriteMarked(mris, "found.mgz") ;
  HISTOfree(&h_inside) ; HISTOfree(&h_outside) ; HISTOfree(&h_grad) ;

  printf("%d surface locations found to contain inconsistent values (%d in, %d out)\n",
         num_in+num_out, num_in, num_out) ;

  freeAndNULL(hops);
  freeAndNULL(vlist);
  
  return(NO_ERROR) ;
}


#define HISTO_NBINS 256

static int
find_and_mark_pinched_regions(MRI_SURFACE *mris,
                              MRI *mri_T2,
                              float nstd_below,
                              float nstd_above)
{
  double    val, xs, ys, zs, xv, yv, zv, d, mean, std ;
  int       vno, num_in, num_out, found_bad_intensity, done, bin, niter ;
  VERTEX    *v ;
  double    min_gray, max_gray, thickness, nx, ny, nz, mn, mx, sig, sample_dist ;
  HISTOGRAM *h, *hcdf ;
  MRI       *mri_filled ;

  sample_dist = MIN(SAMPLE_DIST, mri_T2->xsize/2) ;
  h = HISTOalloc(HISTO_NBINS) ;
  mx = HISTO_NBINS-1 ;
  HISTOinit(h, HISTO_NBINS, 0, mx) ;

  mean = std = 0.0 ;
  num_in = 0 ;
  niter = 0 ;
  do
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      nx = v->x - v->whitex ;
      ny = v->y - v->whitey ;
      nz = v->z - v->whitez ;
      thickness = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
      if (FZERO(thickness))  // pial and white in same place - no cortex here
      {
        continue ;
      }
      MRISvertexToVoxel(mris, v, mri_T2, &xv, &yv, &zv) ;
      MRIsampleVolume(mri_T2, xv, yv, zv, &val) ;
      for (d = thickness/2 ; d <= thickness ; d += sample_dist)
      {
        xs = v->whitex + d*v->nx ;
        ys = v->whitey + d*v->ny ;
        zs = v->whitez + d*v->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
        MRIsampleVolumeType(mri_T2, xv, yv, zv, &val, SAMPLE_TRILINEAR) ;
        if (val <= 0)
        {
          continue ;
        }

        mean += val ;
        std += val*val ;
        num_in++ ;
        HISTOaddSample(h, val, 0, mx) ;
      }
    }
    hcdf = HISTOmakeCDF(h, NULL) ;
    bin = HISTOfindBinWithCount(hcdf, PERCENTILE);
    if (bin < ceil(PERCENTILE*HISTO_NBINS)/4)  // data range is too compressed for histogram to represent
    {
      done = (niter > 10) ;
      if (niter++ > 0)
      {
        mx /= 2 ;
      }
      else
      {
        mx = 2*bin/(.9) ;  // first time - take an educated guess
      }
      HISTOinit(h, HISTO_NBINS, 0, mx) ;
      printf("compressed histogram detected, changing bin size to %f\n",
             h->bin_size) ;
    }
    else
    {
      done = 1 ;
    }
  }
  while (!done) ;
  mean /= num_in ;
  std = sqrt(std/num_in - mean*mean) ;
  max_gray = mean+nstd_above*std ;
  min_gray = mean-nstd_below*std ;

  HISTOrobustGaussianFit(h, .9, &mn, &sig) ;
  if (Gdiag & DIAG_WRITE)
    HISTOplot(h, "h.plt") ;
  if (T2_max_inside < 0)
    max_gray = mn+nstd_above*sig ;
  else
    max_gray = T2_max_inside ;
  if (T2_min_inside < 0)
    min_gray = mn-nstd_below*sig ;
  else
    min_gray = T2_min_inside ;
  printf("locating cortical regions not in the range [%2.2f %2.2f], "
         "gm=%2.2f+-%2.2f, and vertices in regions > %2.1f\n",
         min_gray, max_gray, mn, sig, mn-.5*sig) ;

  MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
  MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  mri_filled = MRIclone(mri_T2, NULL) ;
  MRISfillInterior(mris, mri_T2->xsize, mri_filled) ;
  MRISrestoreVertexPositions(mris, TMP2_VERTICES) ;

  for (num_in = num_out = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->targx = v->x ;
    v->targy = v->y ;
    v->targz = v->z ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    nx = v->x - v->whitex ;
    ny = v->y - v->whitey ;
    nz = v->z - v->whitez ;
    thickness = sqrt(SQR(nx)+SQR(ny)+SQR(nz)) ;
    if (FZERO(thickness))
    {
      continue ;
    }
    MRISvertexToVoxel(mris, v, mri_T2, &xv, &yv, &zv) ;
    nx /= thickness ;
    ny /= thickness ;
    nz /= thickness ;
    found_bad_intensity = 0 ;
    for (d = thickness/2 ; d <= thickness ; d += sample_dist)
    {
      if (MRIgetVoxVal(mri_filled, nint(xv), nint(yv), nint(zv), 0) >0) // interior of white surface
      {
        continue ;
      }
      xs = v->whitex + d*nx ;
      ys = v->whitey + d*ny ;
      zs = v->whitez + d*nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_T2, xv, yv, zv, &val, SAMPLE_TRILINEAR) ;
      if (val < 0)
      {
        continue ;
      }
      if (val < min_gray || val > max_gray)
      {
        found_bad_intensity = 1 ;
        break ;
      }
    }
    if (found_bad_intensity)
    {
      num_in++ ;
      // target surface so interior is good value and exterior is bad gm value
      v->targx = xs - (sample_dist/2*nx) ;
      v->targy = ys - (sample_dist/2*ny) ;
      v->targz = zs - (sample_dist/2*nz) ;
      if (vno == Gdiag_no)
      {
        printf("vno %d: resetting target location to be (%2.1f %2.1f %2.1f), "
               "val = %2.0f\n",
               vno, v->targx, v->targy, v->targz, val) ;
        DiagBreak() ;
      }
    }
    else  // no invalid intensities found in the interior, check for valid ones in the exterior?
    {
#define MAX_OUTWARD_DIST2 10
      for (d = 0 ; d <= MAX_OUTWARD_DIST2 ; d += sample_dist)
      {
        xs = v->x + d*v->nx ;
        ys = v->y + d*v->ny ;
        zs = v->z + d*v->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
        if (MRIgetVoxVal(mri_filled, nint(xv), nint(yv), nint(zv), 0) >0) // interior of white surface
        {
          continue ;
        }
        MRIsampleVolumeType(mri_T2, xv, yv, zv, &val, SAMPLE_TRILINEAR) ;
        if (val < 0)
        {
          continue ;
        }
        if (val < mn-sig ||  val > mn+sig)  // only look for a very narrow range of intensities
        {
          break ;
        }
      }
      if (d > MAX_OUTWARD_DIST2)  // couldn't find pial surface
      {
        d = 0 ;
      }

      if (d > 0)
      {
        d -= sample_dist ;
        num_out++ ;
      }
      v->targx = v->x+d*v->nx ;
      v->targy = v->y+d*v->ny ;
      v->targz = v->z+d*v->nz ;
      {
        int xk, yk, zk, whalf, found ;
        float xi, yi, zi ;

        MRISsurfaceRASToVoxelCached(mris, mri_T2, xs, ys, zs, &xv, &yv, &zv);
        for (whalf = 1 ; whalf <= 6 ; whalf++)
        {
          found = 0 ;
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            xi = (xv + xk) ;
            for (yk = -whalf ; yk <= whalf ; yk++)
            {
              yi = (yv + yk) ;
              for (zk = -whalf ; zk <= whalf ; zk++)
              {
                zi = (zv + zk) ;
                if (MRIgetVoxVal(mri_filled, xi, yi, zi, 0) > 0)
                {
                  continue ;  // don't sample interior to white surface
                }
                MRIsampleVolumeType(mri_T2,
                                    xi, yi, zi,
                                    &val, SAMPLE_TRILINEAR) ;
                if (val < 0)
                {
                  continue ;
                }
                if (val < mn-.5*sig)  // found something that could be csf in FLAIR
                {
                  found = 1 ;
                  break ;
                }
              }
              if (found)
              {
                break ;
              }
            }
            if (found)
            {
              break ;
            }
          }
          if (found)
          {
            break ;
          }
        }
        v->marked = whalf ;
      }
    }
  }
#define OUT_DIST 3
  // find regions that don't have any reasonable CSF-like intensity nearby
  MRISmarkedToCurv(mris) ;
  MRISaverageCurvatures(mris, 10) ;
  MRISthresholdCurvature(mris, OUT_DIST, 1) ;
  MRIScurvToMarked(mris) ;
  MRISdilateMarked(mris, 7) ;      // must smooth a nbhd around each pinch
  MRIScopyMarkedToMarked2(mris) ;  // save mark = bad vertex
  MRISwriteMarked(mris, "bad") ;
  MRISinvertMarks(mris) ;
  if (vno >= 0)
  {
    printf("before soap bubble smoothing:\n") ;
    MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
  }
//#define K_bp   0.1
#define K_bp   0.5
#define Lambda .3
#define Mu     (1.0)/((K_bp)-1.0/Lambda)
  MRISripMarked(mris) ;
  MRISweightedSoapBubbleVertexPositions(mris, 500) ;
  MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
  MRISunrip(mris) ;
  //  MRISweightedSoapBubbleVertexPositions(mris, 500) ;
  if (vno >= 0)
  {
    printf("after soap bubble smoothing:\n") ;
    MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES) ;
    MRISprintVertexStats(mris, Gdiag_no, Gstdout, ORIGINAL_VERTICES) ;
    MRISprintVertexStats(mris, Gdiag_no, Gstdout, WHITE_VERTICES) ;
  }
  MRIScomputeMetricProperties(mris) ;
  MRIScopyMarked2ToMarked(mris) ;  // restore mark = bad vertex
  for (num_in = num_out = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == 0)
    {
      continue ;
    }
    v->targx = v->x+v->nx*OUT_DIST ;
    v->targy = v->y+v->ny*OUT_DIST ;
    v->targz = v->z+v->nz*OUT_DIST ;
    v->val = mn-sig ;   // a mildly CSF intensity in a flair image to seek
    if (vno == Gdiag_no)
    {
      double xv, yv, zv, xv0, yv0, zv0 ;
      MRISsurfaceRASToVoxelCached(mris, mri_T2,
                                  v->x, v->y, v->z,
                                  &xv0, &yv0, &zv0);
      MRISsurfaceRASToVoxelCached(mris, mri_T2,
                                  v->targx, v->targy, v->targz,
                                  &xv, &yv, &zv);
      printf("vno %d: target = (%2.1f, %2.1f, %2.1f) --> "
             "vox (%2.1f, %2.1f, %2.1f), current=(%2.1f, %2.1f, %2.1f)\n",
             vno, v->targx, v->targy, v->targz, xv, yv, zv, xv0, yv0, zv0) ;
      DiagBreak() ;
    }
  }
  printf("%d surface locations found to contain inconsistent "
         "values (%d in, %d out)\n",
         num_in+num_out, num_in, num_out) ;
  MRIfree(&mri_filled) ;
  return(NO_ERROR) ;
}


#include "mrisegment.h"

static int labels_to_correct[] = { Left_Hippocampus,
                                   Right_Hippocampus,
                                   Left_Amygdala,
                                   Right_Amygdala
                                 } ;
#define NLABELS (sizeof(labels_to_correct)/sizeof(labels_to_correct[0]))

static int
edit_aseg_with_surfaces(MRI_SURFACE *mris, MRI *mri_aseg)
{
  MRI              *mri_filled, *mri_hires_aseg ;
  MRI_SEGMENTATION *mseg1, *mseg2 ;
  MRI_SEGMENT      *mseg ;
  int              label, *counts, x, y, z, max_seg_no, sno,vno, alabel, l ;
  MATRIX           *m_vox2vox ;
  VECTOR           *v1, *v2 ;

  printf("correcting aseg with surfaces...\n");
  mri_filled = MRISfillInterior(mris, mri_aseg->xsize/4, NULL) ;
  mri_filled->c_r += mri_aseg->c_r ;
  mri_filled->c_a += mri_aseg->c_a ;
  mri_filled->c_s += mri_aseg->c_s ;
  mri_hires_aseg = MRIresample(mri_aseg, mri_filled, SAMPLE_NEAREST);
  /*  MRIdilate(mri_filled, mri_filled) ; */// fill small breaks
  MRIcopyLabel(mri_filled, mri_hires_aseg, 1) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_hires_aseg, "ha.mgz") ;
    MRIwrite(mri_filled, "hs.mgz") ;
  }

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_hires_aseg, mri_aseg) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  counts = MRIhistogramLabels(mri_aseg,  NULL, MAX_CMA_LABEL+1) ;

  for (l = 0 ; l < NLABELS ; l++)
  {
    label = labels_to_correct[l] ;
    if (counts[label] == 0)
    {
      continue ;
    }
    mseg1 = MRIsegment(mri_aseg, label, label) ;
    if (mseg1->nsegments != 1) // wasn't topologically correct
    {
      MRIsegmentFree(&mseg1) ;
      continue ;
    }
    mseg2 = MRIsegment(mri_hires_aseg, label, label) ;
    if (mseg2->nsegments == 1)  // topology already correct
    {
      MRIsegmentFree(&mseg2) ;
      continue ;
    }

    // turn off the other pieces of the label
    max_seg_no = MRIfindMaxSegmentNumber(mseg2) ;
    for (sno = 0 ; sno < mseg2->nsegments ; sno++)
    {
      if (sno == max_seg_no)
      {
        continue ;
      }
      mseg = &mseg2->segments[sno] ;
      printf("label %s: removing %d voxels in segment %d\n",
             cma_label_to_name(label), mseg->nvoxels, sno) ;
      for (vno = 0 ; vno < mseg->nvoxels ; vno++)
      {
        V3_X(v1) = mseg->voxels[vno].x ;
        V3_Y(v1) = mseg->voxels[vno].y ;
        V3_Z(v1) = mseg->voxels[vno].z ;
        MatrixMultiply(m_vox2vox, v1, v2) ; // to lowres coords
        x = nint(V3_X(v2)) ;
        y = nint(V3_Y(v2)) ;
        z = nint(V3_Z(v2)) ;
        alabel = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (alabel == label)
        {
          MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_undetermined) ;
        }
      }
    }

    MRIsegmentFree(&mseg1) ;
    MRIsegmentFree(&mseg2) ;
  }
  free(counts) ;
  MRIfree(&mri_hires_aseg) ;
  MRIfree(&mri_filled) ;

  return(NO_ERROR) ;
}

/*!
  \fn compute_label_normal()
  \brief This function searches whalf voxels around the given voxel
  for voxels that are on the edge of the given label. It then computes
  the average normal of the edge over these voxels. If use_abs=1, then
  the average of the absolute value of the normals is computed. The
  normal is in voxel units, and if the orientation changes, then the
  interpretation of the normal changes making orientation effectively
  a hidden parameter.
*/
static int
compute_label_normal(MRI *mri_aseg, int x0, int y0, int z0,
                     int label, int whalf, double *pnx, double *pny,
                     double *pnz, int use_abs)
{
  int xi, yi, zi, xk, yk, zk, nvox = 0, val, dx, dy, dz, xn, yn, zn ;
  double  nx, ny, nz, mag ;

  nx = ny = nz = 0.0 ;

  // Search an area of 3x3x3 around xyz0
  for (xk = -whalf ; xk <= whalf ; xk++)  {
    xi = mri_aseg->xi[x0+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)    {
      yi = mri_aseg->yi[y0+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)      {
        zi = mri_aseg->zi[z0+zk] ;

        val = (int)MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
        if (val != label) continue ;

	// If this point is within the label, then look at only the 6 face neighbors
        for (dx = -1 ; dx <= 1 ; dx++)  {
          for (dy = -1 ; dy <= 1 ; dy++)  {
            for (dz = -1 ; dz <= 1 ; dz++)   {
              if (fabs(dx) + fabs(dy) + fabs(dz) != 1) continue ;  // only 8-connected nbrs (??)
                
              xn = mri_aseg->xi[xi+dx] ;
              yn = mri_aseg->yi[yi+dy] ;
              zn = mri_aseg->zi[zi+dz] ;
              val = (int)MRIgetVoxVal(mri_aseg, xn, yn, zn, 0) ;
              if (val != label) {
		// This voxel not the target label but is at the edge of the target label
		// "surface" of label - interface between label and non-label
                nvox++ ;
                if(use_abs){
                  nx += fabs(dx) ;
                  ny += fabs(dy) ;
                  nz += fabs(dz) ;
                }
                else
                {
                  nx += dx ;
                  ny += dy ;
                  nz += dz ;
                }
              }
            }
          }
        }
      }
    }
  }

  if (nvox > 0)
  {
    nx /= nvox ;
    ny /= nvox ;
    nz /= nvox ;
  }

  mag = sqrt(nx*nx + ny*ny + nz*nz) ;
  if (mag > 0)
  {
    nx /= mag ;
    ny /= mag ;
    nz /= mag ;
  }

  *pnx = nx ;
  *pny = ny ;
  *pnz = nz ;
  return(NO_ERROR) ;
}


int
MRISremoveSelfIntersections(MRI_SURFACE *mris)
{
  int              n, vno, fno ;
  VERTEX           *v ;

//  return(NO_ERROR) ;
  MRISremoveIntersections(mris,0);

  return(NO_ERROR) ;

  MRIS_HASH_TABLE* mht = MHTcreateFaceTable(mris) ;
  
  MRISsetMarks(mris, 1) ;  // fixed vertices
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (MHTdoesFaceIntersect(mht, mris, fno))
    {
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      {
	vno = mris->faces[fno].v[n] ;
	v = &mris->vertices[vno] ;
	v->marked = 0 ;
	if (vno == Gdiag_no)
	  DiagBreak() ;
      }
    }
  }

  MRISdilateMarked(mris, 1) ;
  MRISwriteMarked(mris, "marked.mgz") ;
  MRISsoapBubbleVertexPositions(mris, 500) ;
  MRISwrite(mris, "after_soap") ;

  MHTfree(&mht) ;
  return(NO_ERROR) ; ;
}



