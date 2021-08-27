/**
 * @brief refines the surfaces lh(rh).pial and lh(rh).white
 *
 * This program refines the surfaces lh(rh).pial and lh(rh).white
 * around the region specified by the label file, producing lh(rh).pialhires
 * and lh(rh).whitehires. The subject must be processed beforehand to have
 * mri/filled, mri/wm, surf/lh(rh).orig, surf/lh(rh).white, surf/lh(rh).pial 
 * in the subject directory.  The lowtohires.xfm is an optional argument 
 * to give.
 */
/*
 * Original Author: Bruce Fischl (June 16, 1998)
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

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "version.h"
#include "label.h"
#include "registerio.h"


int debug__ = 0; /// tosa debug


int main(int argc, char *argv[]) ;

//#define BRIGHT_LABEL         130
//#define BRIGHT_BORDER_LABEL  100

//static int  MRIScomputeClassStatistics(MRI_SURFACE *mris, MRI *mri, 
//float *pwhite_mean, float *pwhite_std, float *pgray_mean, float *pgray_std) ;
//static int compute_border_gradients(MRI_SURFACE *mris, MRI *mri, double thresh) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

MRI *MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_hires, float thresh,
                      int out_label, MRI *mri_dst);

int LocalMRISfindExpansionRegions(MRI_SURFACE *mris) ;

static double pial_errfunc_gradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double pial_errfunc_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double pial_errfunc_rms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;

static MRI *make_pial_location_mask(MRI_SURFACE *mris, MRI *mri_hires, float max_thickness, MRI *mri_dst) ;
static LABEL *hires_label = NULL ;
static TRANSFORM *hires_xform = NULL ;
static LTA *hires_lta = 0 ;

static const char *orig_white = "white";
static char *orig_pial = NULL ;

const char *Progname ;

static char *reg_fname = NULL ;
static int graymid = 0 ;
static int curvature_avgs = 10 ;
static int create = 1 ;
static int smoothwm = 0 ;
static int white_only = 0 ;
static int inverted_contrast = 0 ;
static int auto_detect_stats = 0 ;
static int apply_median_filter = 0 ;
static int nbhd_size = 20 ;

#define MAX_OTHER_LABELS 100
static int num_other_labels = 0 ;
static char *other_label_names[MAX_OTHER_LABELS] ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0

static int add = 0 ;

static double l_tsmooth = 0.0 ;
static double l_surf_repulse = 5.0 ;

static int smooth = 0 ;
static int vavgs = 5 ;
static int nwhite = 50 /*5*/ ;
static int ngray = 50 /*45*/ ;

static int nowhite = 0 ;
static int nbrs = 2 ;
static int write_vals = 0 ;

static const char *orig_name = ORIG_NAME ;

static char output_suffix[STRLEN] = "hires" ;
static char pial_name[STRLEN] = "pial" ;
static char white_matter_name[STRLEN] = WHITE_MATTER_NAME ;

#if 0
static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;
#endif

static int max_pial_averages = 16 ;
static int min_pial_averages = 2 ;
static int max_white_averages = 4 ;
static int min_white_averages = 0 ;
static float pial_sigma = 2.0f ;
static float white_sigma = 1.0 ;
static float max_thickness = 5.0 ;

static int num_dilate = 0 ;

#define MAX_WHITE             120
#define MAX_BORDER_WHITE      105
#define MIN_BORDER_WHITE       85
#define MIN_GRAY_AT_WHITE_BORDER  70

#define MAX_GRAY               95
#define MIN_GRAY_AT_CSF_BORDER    40
#define MID_GRAY               ((max_gray + min_gray_at_csf_border) / 2)
#define MAX_GRAY_AT_CSF_BORDER    75
#define MIN_CSF                10
#define MAX_CSF                40

static MRI *find_wm(MRI_SURFACE *mris, MRI *mri, MRI *mri_wm) ;
static double max_white = MAX_WHITE ;
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
static char suffix[STRLEN] = "" ;

static int MGZ = 1; // for use with MGZ format
static float check_contrast_direction
(MRI_SURFACE *mris,MRI *mri_hires) ; //defined at the end

#define MIN_BORDER_DIST 2.0  // mm from border

int
main(int argc, char *argv[]) {
  char          *hemi, *sname, *cp=0, fname[STRLEN], med_fname[STRLEN], *dot ;
  int           nargs, i/*, label_val, replace_val*/, msec, n_averages, j ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm = NULL, *mri_smooth = NULL, *mri_hires, 
                *mri_hires_pial, *mri_mask = NULL, *mri_hires_orig ;
  /*  MRI           *mri_tran = NULL;*/
  float         max_len ;
  float         white_mean, white_std, gray_mean, gray_std ;
  double        current_sigma ;
  Timer then ;
  MATRIX        *m_reg = NULL ;
  /*  LT            *lt =0;*/

  nargs = handleVersionOption(argc, argv, "mris_refine_surfaces");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.projection = NO_PROJECTION ;
  parms.fill_interior = 0 ;  // don't let gradient use exterior information (slows things down)
  parms.tol = 1e-4 ;
  parms.dt = 0.5f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;
  parms.l_spring = 1.0f ;
  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.1 ;
  parms.l_spring = 0.0f ;
  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.2 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 1 ;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.l_surf_repulse = 0.0 ;
  parms.l_repulse = 1 ;
#if 0
  parms.l_grad = 1 ;
  parms.l_intensity = 0 ;
#endif

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  /////////////////////////////////////////////////////////////////////////
  // make sure that subjects_dir is set
  /////////////////////////////////////////////////////////////////////////
  if (!strlen(sdir))  // -sdir option
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  if (argc < 4)
    usage_exit() ;
  else {
    printf("\nArguments given are:\n");
    printf("\tSUBJECTS_DIR = %s\n", sdir);
    printf("\tsubject_name = %s\n", argv[1]);
    printf("\themisphere  = %s\n", argv[2]);
    printf("\thires volume = %s\n", argv[3]);
    printf("\tlabel file   = %s\n", argv[4]);
    if (argc==5)
      printf("\txfm          = assumes lowres to highres "
             "ras-to-ras identity\n");
    else if (argc==6)
      printf("\txfm          = %s\n", argv[5]);
    else if (argc > 6) {
      print_usage();
      ErrorExit(ERROR_BADPARM, "%s: gave too many arguments\n", Progname);
    }
    if (MGZ)
      printf("\t-mgz         = use .mgz extension for reading volumes\n");
    if (strlen(suffix))
      printf("\t-suffix      = add suffix %s to the final surfaces\n", suffix);
    printf("\n");
  }
  /////////////////////////////////////////////////////////////////////////////
  // arg1 = subject_name, arg2 = hemi, arg3 = hires vol, arg4 = label file
  // arg5 = lowtohires.xfm is an option
  /////////////////////////////////////////////////////////////////////////////
  sname = argv[1] ;
  hemi = argv[2] ;

  ///////////////////////////////////////////////////////////////////////
  // read orig surface
  ///////////////////////////////////////////////////////////////////////
  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, sname, hemi, orig_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading original surface position from %s...\n", fname) ;
  mris = MRISreadOverAlloc(fname, 1.1) ;
  // this tries to get src c_(ras) info
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  ///////////////////////////////////////////////////////////////////////
  // read hires volume
  ///////////////////////////////////////////////////////////////////////
  req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, sname, argv[3]) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading hires volume %s...\n", fname) ;

  strcpy(med_fname, fname) ;
  dot = strchr(med_fname, '.') ;
  *dot = 0 ;
  strcat(med_fname, ".median.mgz") ;
  if (0)
  {
    MRI_REGION box ;
    double     dist = MIN_BORDER_DIST-.5 ;
    MRI        *mri_tmp ;
      
    mri_hires = mri_hires_orig = MRIread(fname);
    if (!mri_hires)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                Progname, fname) ;

    box.x = nint(dist/mri_hires->xsize) ;
    box.y = nint(dist/mri_hires->ysize) ;
    box.z = nint(dist/mri_hires->zsize) ;
    box.dx = mri_hires->width - 2*nint(dist/mri_hires->xsize) ;
    box.dy = mri_hires->height - 2*nint(dist/mri_hires->ysize) ;
    box.dz = mri_hires->depth - 2*nint(dist/mri_hires->zsize) ;
      
    mri_tmp = MRIextract(mri_hires, NULL, box.x, box.y, box.z,
                         box.dx, box.dy, box.dz) ;
    //    MRIfree(&mri_hires) ;  kept in mri_hires_orig
    mri_hires = mri_hires_pial = mri_tmp ;

    if (apply_median_filter) /* -median option ... modify 
                                mri_hires using the filter */
    {
      fprintf(stderr, "applying median filter to T1 image...\n") ;
      if (1)
      {
        mri_tmp = MRIcopy(mri_hires, NULL) ;  // copy stuff outside of range
        MRImedian(mri_hires, mri_tmp, 3, NULL) ;
      //      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        printf("writing median to %s\n", med_fname) ;
        MRIwrite(mri_tmp, med_fname) ;
      }
      else
      {
        printf("reading median from %s\n", med_fname) ;
        mri_tmp = MRIread(med_fname) ;
        if (!mri_tmp)
          ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                    Progname, med_fname) ;

      }
      MRIfree(&mri_hires) ;
      mri_hires = mri_hires_pial = mri_tmp ;
    }
  }
  else
  {
    printf("reading median filtered volume from %s\n", med_fname) ;
    mri_hires = mri_hires_orig = mri_hires_pial = MRIread(med_fname) ;
    if (!mri_hires)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname, med_fname) ;
  }


  /////////////////////////////////////////////////////////////////////////////
  // read or setup lowrestohires transform
  /////////////////////////////////////////////////////////////////////////////
  if (argc <= 5) {
    LINEAR_TRANSFORM *lt = 0;
    MATRIX *m_L = 0;
    fprintf(stderr, "allocating identity RAS-to-RAS xform...\n") ;
    // allocate hires_xform->xform
    hires_xform = TransformAlloc(MNI_TRANSFORM_TYPE, NULL);
    if (!hires_xform)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not allocate hires xform %s", Progname, argv[3]) ;
    hires_lta = (LTA *) (hires_xform->xform);
    lt = &hires_lta->xforms[0];
    lt->sigma = 1.0f ;
    lt->x0 = lt->y0 = lt->z0 = 0 ;
    hires_lta->type = LINEAR_RAS_TO_RAS;
    m_L = lt->m_L;
    MatrixIdentity(4, m_L);
    // geometry for src and dst are not specified yet
  } else if (argc ==6) {
    fprintf(stderr, "reading lowrestohires.xfm %s....\n", argv[5]);
    hires_xform = TransformRead(argv[5]) ;
    if (!hires_xform)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read hires xform %s", Progname, argv[3]) ;
    hires_lta = (LTA *)(hires_xform->xform) ;
  }

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = nwhite ;
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  //////////////////////////////////////////////////////////////////////////
  // timer starts here
  //////////////////////////////////////////////////////////////////////////
  then.reset() ;

  ////////////////////////////////////////////////////////////////////////////
  // move the vertex positions to the lh(rh).white positions (low res)
  ////////////////////////////////////////////////////////////////////////////
  printf("reading initial white vertex positions from %s...\n", orig_white) ;
  if (reg_fname && mri_hires_orig)
  {
    char *subject ;
    int  type ;
    LTA  *lta ;

    type = TransformFileNameType(reg_fname) ;
    printf("reading registration from %s (type = %d)\n", reg_fname, type);
    switch (type)
    {
    case REGISTER_DAT:
      m_reg = regio_read_surfacexform_from_register_dat(reg_fname, mris, mri_hires_orig, &subject);
      if (m_reg == NULL)
        exit(Gerror) ;
      break ;
    case LINEAR_VOX_TO_VOX:
    case LINEAR_RAS_TO_RAS:
    case TRANSFORM_ARRAY_TYPE:
      lta = LTAreadEx(reg_fname) ;
      m_reg = MatrixCopy(lta->xforms[0].m_L, NULL) ;
      LTAfree(&lta) ;
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "%s: cannot read registration file %s type %d", Progname, fname, type) ;
    }
  }

  if (MRISreadVertexPositions(mris, orig_white) != NO_ERROR)
    ErrorExit(Gerror, "reading of orig white failed...");

  if (m_reg)
    MRISmatrixMultiply(mris, m_reg) ;

  /////////////////////////////////////////////////////////////////////////////
  // read label file
  /////////////////////////////////////////////////////////////////////////////
  if (argc > 4)  // label specified explicitly
  {
    req = snprintf(fname, STRLEN, "%s/%s/label/%s", sdir, sname, argv[4]) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "reading the label file %s...\n", fname);
    hires_label = LabelRead(NULL, fname) ;
    if (!hires_label)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read hires label %s", Progname, argv[2]) ;
    if (num_dilate > 0)
      LabelDilate(hires_label, mris, num_dilate, CURRENT_VERTICES) ;

    for (i = 0 ; i < num_other_labels ; i++) {
      LABEL *olabel ;
      req = snprintf(fname, STRLEN, "%s/%s/label/%s", sdir, sname, other_label_names[i]) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "reading the label file %s...\n", fname);
      olabel = LabelRead(NULL, fname) ;
      if (!olabel)
        ErrorExit(ERROR_NOFILE, 
                  "%s: could not read label %s", Progname, fname) ;
      if (num_dilate > 0)
        LabelDilate(olabel, mris, num_dilate, CURRENT_VERTICES) ;
      hires_label = LabelCombine(olabel, hires_label) ;
      LabelFree(&olabel) ;
    }
  }
  else   // build a label that contains the whole field of view
  {
    hires_label = LabelInFOV(mris, mri_hires, 2*mri_hires->xsize) ;
  }

  printf("deforming %d vertices\n", hires_label->n_points) ;

  ////////////////////////////////////////////////////////////////////////////
  // mark vertices near label positions as ripped=0 but the rest rippped = 1
  ////////////////////////////////////////////////////////////////////////////
  LabelToFlat(hires_label, mris) ;
  LabelRipRestOfSurface(hires_label, mris) ;

  ////////////////////////////////////////////////////////////////////////
  // check contrast direction
  ////////////////////////////////////////////////////////////////////////
  inverted_contrast = (check_contrast_direction(mris, mri_hires) < 0) ;
  if (inverted_contrast) {
    fprintf(stderr, "inverted contrast detected....\n") ;
  }

  mri_wm = find_wm(mris, mri_hires, NULL) ;


  /////////////////////////////////////////////////////////////////////
  // detect grey/white stats
  /////////////////////////////////////////////////////////////////////
  if (auto_detect_stats) {
    MRI *mri_tmp=0 ;

    // now mri_wm and mri_hires are in the same resolution space
    mri_tmp = MRIbinarize(mri_wm, NULL, WM_MIN_VAL, MRI_NOT_WHITE, MRI_WHITE) ;
    if (debug__ ==1)
      MRIwrite(mri_tmp, "./wmtranbinarized.mgz");

    fprintf(stderr, "computing class statistics...\n");
#if 0
    MRIcomputeClassStatistics(mri_hires, mri_tmp, 30, WHITE_MATTER_MEAN,
                              &white_mean, &white_std, &gray_mean,
                              &gray_std) ;
#else
    MRIScomputeClassStatistics(mris, mri_hires,
                               &white_mean, 
                               &white_std, 
                               &gray_mean, 
                               &gray_std) ;
#endif

    //    if (abs(white_mean-110) > 50)
    //   parms.l_nspring *= sqrt(110/white_mean) ;

    if (!min_gray_at_white_border_set)
      min_gray_at_white_border = gray_mean-gray_std ;
    if (!max_border_white_set)
      max_border_white = white_mean+white_std ;
    if (!max_csf_set)
      max_csf = gray_mean-2*gray_std ;
    if (!min_border_white_set)
      min_border_white = white_mean - white_std ;
    fprintf(stderr, "setting MIN_GRAY_AT_WHITE_BORDER to %2.1f (was %d)\n",
            min_gray_at_white_border, MIN_GRAY_AT_WHITE_BORDER) ;
    fprintf(stderr, "setting MAX_BORDER_WHITE to %2.1f (was %d)\n",
            max_border_white, MAX_BORDER_WHITE) ;
    fprintf(stderr, "setting MIN_BORDER_WHITE to %2.1f (was %d)\n",
            min_border_white, MIN_BORDER_WHITE) ;
    fprintf(stderr, "setting MAX_CSF to %2.1f (was %d)\n",
            max_csf, MAX_CSF) ;
    if (min_border_white > max_white) {
      max_white = white_mean+white_std ;
      printf("setting MAX_WHITE to %2.1f (was %d)\n", max_white, MAX_WHITE);
    }

    if (!max_gray_set)
      max_gray = white_mean-white_std ;
    if (!max_gray_at_csf_border_set)
      max_gray_at_csf_border = gray_mean-0.5*gray_std ;
    if (!min_gray_at_csf_border_set)
      min_gray_at_csf_border = gray_mean - 3*gray_std ;
    fprintf(stderr, "setting MAX_GRAY to %2.1f (was %d)\n",
            max_gray, MAX_GRAY) ;
    fprintf(stderr, "setting MAX_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",
            max_gray_at_csf_border, MAX_GRAY_AT_CSF_BORDER) ;
    fprintf(stderr, "setting MIN_GRAY_AT_CSF_BORDER to %2.1f (was %d)\n",
            min_gray_at_csf_border, MIN_GRAY_AT_CSF_BORDER) ;
    MRIfree(&mri_tmp) ;
  }
  /*  MRIfree(&mri_wm);*/
  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  req = snprintf(parms.base_name, STRLEN, "%s%s%s", white_matter_name, output_suffix, suffix);
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////////
  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
  MRISstoreMetricProperties(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

  if (add) {
    fprintf(stderr, "adding vertices to initial tessellation...\n") ;
    for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
    while (MRISdivideLongEdges(mris, max_len) > 0) {}
  }

  MRISsetVals(mris, -1) ;  /* clear white matter intensities */

  if (!nowhite) {
    fprintf(stderr, "repositioning cortical surface to gray/white boundary\n");

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_hires, "white_masked.mgh") ;
  }
  //////////////////////////////////////////////////////////////////////
  // vertices are lh(rh).white positions now
  ///////////////////////////////////////////////////////////////////////

#if 0
  //////////////////////////////////////////////////////////////////////
  // convert surface into hires volume surface
  //////////////////////////////////////////////////////////////////////
  MRISsurf2surf(mris, mri_hires, hires_lta);
#endif

  ///////////////////////////////////////////////////////////////////////
  // loop here
  ///////////////////////////////////////////////////////////////////////
  current_sigma = white_sigma ;  // in mm not voxels
  for (n_averages = max_white_averages, i = 0 ;
       n_averages >= min_white_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    if (nowhite)     // if nowhite set, don't do the loop
      break ;

    parms.sigma = current_sigma ;
    fprintf(stderr, "smoothing hires volume with sigma = %2.3f\n",
            current_sigma) ;
    if (!mri_smooth)
      mri_smooth = MRIclone(mri_hires, NULL) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN] ;
      sprintf(fname, "sigma%.0f.mgh", current_sigma) ;
      fprintf(stderr, "writing smoothed volume to %s...\n", fname) ;
      MRIwrite(mri_smooth, fname) ;
    }
    parms.n_averages = n_averages ;
    MRISprintTessellationStats(mris, stderr) ;
#if 0
    if (inverted_contrast)
      MRIScomputeInvertedGrayWhiteBorderValues(mris, mri_hires, mri_smooth,
          max_white, max_border_white, min_border_white,
          min_gray_at_white_border,
          max_border_white /*max_gray*/, current_sigma,
          2*max_thickness, parms.fp) ;
    else
      MRIScomputeBorderValues(mris, mri_hires, mri_smooth,
                              max_white, max_border_white, min_border_white,
                              min_gray_at_white_border,
                              max_border_white /*max_gray*/, current_sigma,
                              2*max_thickness, parms.fp, GRAY_WHITE, NULL, 0, parms.flags) ;
#else
    MRIScomputeMaxGradBorderValues(mris,mri_hires, mri_smooth,
                                   current_sigma, max_thickness, 1, parms.fp,
                                   mri_wm, i) ;
    //    compute_border_gradients(mris, mri_hires, 20) ;
#endif    
    LocalMRISfindExpansionRegions(mris) ;
    if (vavgs) {
      fprintf(stderr, "averaging target values for %d iterations...\n",vavgs) ;
      MRISaverageMarkedVals(mris, vavgs) ;
      if (Gdiag_no > 0) {
        VERTEX *v ;
        v = &mris->vertices[Gdiag_no] ;
        fprintf(stderr,"v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                Gdiag_no, v->val, v->mean, v->d) ;
      }
    }


    if (write_vals) {
      sprintf(fname, "./%s-white%2.2f.mgz", hemi, current_sigma) ;
      MRISwriteValues(mris, fname) ;
    }
    if (i == 0)
    {
      parms.niterations = 15 ;
      parms.grad_dir = 1 ;
      MRISpositionSurface(mris, mri_hires, mri_smooth,&parms);
      MRIScomputeMaxGradBorderValues(mris,mri_hires, mri_smooth,
                                     current_sigma, max_thickness, 1, parms.fp,
                                     mri_wm, i) ;
      //    compute_border_gradients(mris, mri_hires, 20) ;
      LocalMRISfindExpansionRegions(mris) ;
      if (vavgs) {
        fprintf(stderr, "averaging target values for %d iterations...\n",vavgs) ;
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no > 0) {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf(stderr,"v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                  Gdiag_no, v->val, v->mean, v->d) ;
        }
      }

      printf("enabling interior filling and gradient direction calculation...\n") ;
      parms.niterations = nwhite ;
      parms.grad_dir = 0 ;
      parms.fill_interior = 1 ;
    }
    MRISpositionSurface(mris, mri_hires, mri_smooth,&parms);
    if (add) {
      for (max_len = 1.5*8 ; max_len > 1 ; max_len /= 2)
      while (MRISdivideLongEdges(mris, max_len) > 0) {}
    }
    if (!n_averages)
      break ;
  }
  // end of loop
  /////////////////////////////////////////////////////////////////////////

  if (!nowhite) {
    req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname,hemi,white_matter_name,
            output_suffix,suffix);
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "writing white matter surface to %s...\n", fname) ;
    /////////////////////////////////////////////////////////////////////////
    // convert surface into lowres volume surface again for more processing
    //////////////////////////////////////////////////////////////////////////

    // default smoothwm = 0 and thus don't do anything
    MRISaverageVertexPositions(mris, smoothwm) ;
    // write out lh(rh).whitehires
    MRISwrite(mris, fname) ;
    if (create)   /* write out curvature and area files */
    {
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, curvature_avgs) ;
      req = snprintf(fname, STRLEN, "%s.curv%s%s",
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
              suffix);
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "writing smoothed curvature to %s\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      req = snprintf(fname, STRLEN, "%s.area%s%s",
              mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh", output_suffix,
              suffix);
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      MRISprintTessellationStats(mris, stderr) ;

      //  restore to hires for further processing
      MRISsurf2surfAll(mris, mri_hires, hires_lta);
    }
  } else   /* read in previously generated white matter surface */
  {
    req = snprintf(fname, STRLEN, "%s", white_matter_name) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: could not read white matter surfaces.",
                Progname) ;
    //  restore to hires for further processing
    //    MRISsurf2surfAll(mris, mri_hires, hires_lta);
    MRIScomputeMetricProperties(mris) ;
    if (m_reg)
      MRISmatrixMultiply(mris, m_reg) ;
  }

  if (white_only) {
    msec = then.milliseconds() ;
    fprintf(stderr,
            "refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
    MRIfree(&mri_hires);
    exit(0) ;
  }

  ///////////////////////////////////////////////////////////////////////
  // now pial surface
  ///////////////////////////////////////////////////////////////////////
  parms.t = parms.start_t = 0 ;
  req = snprintf(parms.base_name, STRLEN, "%s%s%s", pial_name, output_suffix, suffix);
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  parms.niterations = ngray ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter 
                                                        positions */
  parms.l_surf_repulse = l_surf_repulse ;

  MRISsetVals(mris, -1) ;  /* clear target intensities */

  if (smooth && !nowhite) {
    printf("smoothing surface for %d iterations...\n", smooth) ;
    MRISaverageVertexPositions(mris, smooth) ;
  }

  fprintf(stderr, "repositioning cortical surface to gray/csf boundary.\n") ;
  parms.l_repulse = 0 ;

  if (orig_pial != NULL) {
    //////////////////////////////////////////////////////////////////////////
    // set vertices to orig_pial positions
    //////////////////////////////////////////////////////////////////////////
    printf("reading initial pial vertex positions from %s...\n", orig_pial) ;
    if (MRISreadVertexPositions(mris, orig_pial) != NO_ERROR)
      ErrorExit(Gerror, "reading orig pial positions failed") ;
    if (m_reg)
      MRISmatrixMultiply(mris, m_reg) ;
  }

  // we need to retain the mask
  mri_hires = mri_hires_pial ;

  /////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////
  for (j = 0 ; j <= 0 ; parms.l_intensity *= 2, j++)  /* only once for now */
  {
    current_sigma = pial_sigma ;
    for (n_averages = max_pial_averages, i = 0 ;
         n_averages >= min_pial_averages ;
         n_averages /= 2, current_sigma /= 2, i++) {
      parms.sigma = current_sigma ;
      fprintf(stderr, "smoothing hires volume with sigma = %2.3f\n",
              current_sigma) ;
      parms.n_averages = n_averages ;
      parms.l_tsmooth = l_tsmooth ;
      /*
        replace bright stuff such as eye sockets with 255. 
        Simply zeroing it out
        would make the border always go through the sockets, and ignore subtle
        local minima in intensity at the border of the sockets. Will set to 0
        after border values have been computed so that it 
        doesn't mess up gradients.
      */
      if (inverted_contrast == 0) {
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
          MRIwrite(mri_hires, "pial_masked.mgh") ;
      }
      MRISunrip(mris) ;
      if (hires_label)
        LabelFree(&hires_label) ;
      hires_label = LabelInFOV(mris, mri_hires, 2*mri_hires->xsize) ;
      LabelRipRestOfSurface(hires_label, mris) ;
      if (mri_mask == NULL)
        mri_mask = make_pial_location_mask(mris, mri_hires, max_thickness, NULL) ;
      if (inverted_contrast)
        MRIScomputeMaxGradBorderValuesPial(mris,mri_hires, mri_smooth,
                                           current_sigma, max_thickness, -1, parms.fp,i, mri_mask) ;
      else
        MRIScomputeBorderValues
          (mris, mri_hires, mri_smooth, max_gray,
           max_gray_at_csf_border, min_gray_at_csf_border,
           min_csf,(max_csf+max_gray_at_csf_border)/2,
           current_sigma, 2*max_thickness, parms.fp,
           GRAY_CSF, NULL, 0, parms.flags, NULL, -1, -1) ;

      if (i == 0)
      {
        gMRISexternalGradient = pial_errfunc_gradient ;
        gMRISexternalSSE = pial_errfunc_sse ;
        gMRISexternalRMS = pial_errfunc_rms ;
        parms.l_external = parms.l_intensity ;
        parms.l_intensity = 0 ;
      }
      else
      {
        gMRISexternalGradient = gMRISexternalSSE = gMRISexternalRMS = NULL ;
        parms.l_intensity = parms.l_external ;
        parms.l_external = 0 ;
      }

      if (vavgs) {
        fprintf(stderr, 
                "averaging target values for %d iterations...\n",vavgs) ;
        MRISaverageMarkedVals(mris, vavgs) ;
        if (Gdiag_no > 0) {
          VERTEX *v ;
          v = &mris->vertices[Gdiag_no] ;
          fprintf(stderr,
                  "v %d, target value = %2.1f, mag = %2.1f, dist=%2.2f\n",
                  Gdiag_no, v->val, v->mean, v->d) ;
        }
      }

      if (write_vals) {
        sprintf(fname, "./%s-pial%2.2f.mgz", hemi, current_sigma) ;
        MRISwriteValues(mris, fname) ;
      }
      if (!mri_smooth)
        mri_smooth = MRIcopy(mri_hires, NULL) ;
      MRISpositionSurface(mris, mri_hires, mri_smooth,&parms);
      /*    parms.l_nspring = 0 ;*/
      if (!n_averages)
        break ;
    }
  }

  req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s%s%s", sdir, sname, hemi, pial_name,
          output_suffix, suffix) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  /////////////////////////////////////////////////////////////////////////
  // write out pial surface
  /////////////////////////////////////////////////////////////////////////
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;

  ////////////////////////////////////////////////////////////////////////
  // free up memory
  ////////////////////////////////////////////////////////////////////////
  MRIfree(&mri_hires);

  /*  if (!(parms.flags & IPFLAG_NO_SELF_INT_TEST))*/
  ////////////////////////////////////////////////////////////////////////
  // measure surface statistics
  ////////////////////////////////////////////////////////////////////////
  {
    fprintf(stderr, "measuring cortical thickness...\n") ;
    MRISmeasureCorticalThickness(mris, nbhd_size, max_thickness) ;
    fprintf(stderr,
            "writing cortical thickness estimate to 'thickness' file.\n") ;
    req = snprintf(fname, STRLEN, "thickness%s%s", output_suffix, suffix) ;
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwriteCurvature(mris, fname) ;

    /* at this point, the v->curv slots contain the cortical surface. Now
       move the white matter surface out by 1/2 the thickness as an estimate
       of layer IV.
    */
    if (graymid) {
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      mrisFindMiddleOfGray(mris) ;
      sprintf(fname, "%s/%s/surf/%s.%s%s", sdir, sname, hemi, GRAYMID_NAME,
              suffix) ;
      fprintf(stderr, "writing layer IV surface to %s...\n", fname) ;
      MRISwrite(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
  }
  MRISfree(&mris);
  msec = then.milliseconds() ;
  fprintf(stderr,"positioning took %2.1f minutes\n", (float)msec/(60*1000.0f));

  return(0) ;  /* for ansi */
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "fill_interior"))
  {
    parms.fill_interior = 1 ;  // don't let gradient use exterior information (slows things down)
    printf("limiting gradient calculations to interior of surface\n") ;
  }
  else if (!stricmp(option, "reg"))
  {
    reg_fname = argv[2] ;
    printf("using registration file %s to transform surfaces\n", reg_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nowhite") || !stricmp(option, "pialonly")) {
    nowhite = 1 ;
    fprintf(stderr, "reading previously computed gray/white surface\n") ;
  } else if (!stricmp(option, "wa")) {
    max_white_averages = atoi(argv[2]) ;
    fprintf(stderr, "using max white averages = %d\n", max_white_averages) ;
    nargs = 1 ;
    if (isdigit(*argv[3])) {
      min_white_averages = atoi(argv[3]) ;
      fprintf(stderr, "using min white averages = %d\n", min_white_averages) ;
      nargs++ ;
    }
  } else if (!stricmp(option, "orig_pial")) {
    orig_pial = argv[2] ;
    printf("using %s starting pial locations...\n", orig_pial) ;
    nargs = 1 ;
  } else if (!stricmp(option, "orig_white")) {
    orig_white = argv[2] ;
    printf("using %s starting white location...\n", orig_white) ;
    nargs = 1 ;
  } else if (!stricmp(option, "median")) {
    nargs = 0 ;
    printf("applying median filter\n") ;
    apply_median_filter = 1 ;
  } else if (!stricmp(option, "whiteonly") || !stricmp(option, "nopial")) {
    white_only = 1 ;
    fprintf(stderr,  "only generating white matter surface\n") ;
  } else if (!stricmp(option, "dilate")) {
    num_dilate = atoi(argv[2]) ;
    fprintf(stderr,  "dilating labels %d times\n", num_dilate) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_border_white")) {
    min_border_white_set = 1 ;
    min_border_white = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "wsigma")) {
    white_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            white_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "psigma")) {
    pial_sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n",
            pial_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_gray_at_white_border")) {
    min_gray_at_white_border_set = 1 ;
    min_gray_at_white_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_gray")) {
    max_gray_set = 1 ;
    max_gray = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_gray_at_csf_border")) {
    max_gray_at_csf_border_set = 1 ;
    max_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_gray_at_csf_border")) {
    min_gray_at_csf_border_set = 1 ;
    min_gray_at_csf_border = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "min_csf")) {
    min_csf_set = 1 ;
    min_csf = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max_csf")) {
    max_csf_set = 1 ;
    max_csf = atof(argv[2]) ;
    nargs = 1 ;
  } else if (!stricmp(option, "noauto")) {
    auto_detect_stats = 0 ;
    fprintf(stderr, "disabling auto-detection of border ranges...\n") ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "intensity")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "spring")) {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
  } else if (!stricmp(option, "tspring")) {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  } else if (!stricmp(option, "smooth")) {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing for %d iterations\n", smooth) ;
  } else if (!stricmp(option, "tsmooth")) {
    l_tsmooth = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tsmooth = %2.3f\n", l_tsmooth) ;
  } else if (!stricmp(option, "nspring")) {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  } else if (!stricmp(option, "curv")) {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  } else if (!stricmp(option, "MGZ")) {
    MGZ=1;
  } else if (!stricmp(option, "SUFFIX")) {
    strcpy(suffix, argv[2]);
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'O':
      strcpy(output_suffix, argv[2]) ;
      printf("using output suffix %s...\n", output_suffix) ;
      nargs = 1 ;
      break ;
    case 'L':
      other_label_names[num_other_labels++] = argv[2] ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    }
  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr, "usage: %s [options] <subject_name> <hemi> "
          "<hires volume> <label> [<lowtohires.xfm>]\n",
          Progname) ;
  fprintf(stderr, "options:\n");
  fprintf(stderr, "       -sdir $(SUBJECTS_DIR) : specifies SUBJECTS_DIR\n");
  fprintf(stderr, "       -mgz                  : use .mgz volumes\n");
  fprintf(stderr, "       -suffix $(SUFFIX)     : "
          "add $(SUFFIX) to the final surfaces\n");
  fprintf(stderr, "\nsample: mris_refine_surface agt rh "
          "hires rh_hires.label\n");
  fprintf(stderr, "        where hires is located at "
          "$(SUBJECTS_DIR)/agt/mri/hires, \n");
  fprintf(stderr, "        rh_hires.label is located at "
          "$(SUBJECTS_DIR)/agt/label/rh_hires.label.\n");
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr, "\nThis program refines the surfaces lh(rh).pial "
          "and lh(rh).white\n");
  fprintf(stderr, "around the region specified by the label file, "
          "producing lh(rh).pialhires\n");
  fprintf(stderr, "and lh(rh).whitehires. The subject must be "
          "processed beforehand to have\n");
  fprintf(stderr, "mri/filled, mri/wm, surf/lh(rh).orig, "
          "surf/lh(rh).white, surf/lh(rh).pial in the \n");
  fprintf(stderr, "subject directory.  The lowtohires.xfm is "
          "an optional argument to give.\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

MRI *
MRIfillVentricle(MRI *mri_inv_lv, MRI *mri_hires, float thresh,
                 int out_label, MRI *mri_dst) {
  BUFTYPE   *pdst, *pinv_lv, out_val, inv_lv_val, *pT1 ;
  int       width, height, depth, x, y, z,
  ventricle_voxels;

  if (!mri_dst)
    mri_dst = MRIclone(mri_hires, NULL) ;

  width = mri_hires->width ;
  height = mri_hires->height ;
  depth = mri_hires->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
  */


  ventricle_voxels = 0 ;
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      pT1 = &MRIvox(mri_hires, 0, y, z) ;
      pinv_lv = &MRIvox(mri_inv_lv, 0, y, z) ;
      for (x = 0 ; x < width ; x++) {
        pT1++ ;
        inv_lv_val = *pinv_lv++ ;
        out_val = 0 ;
        if (inv_lv_val >= thresh) {
          ventricle_voxels++ ;
          out_val = out_label ;
        }
        *pdst++ = out_val ;
      }
    }
  }

#if 0
  MRIfillRegion(mri_hires, mri_dst, 30, out_label, 2*ventricle_voxels) ;
  MRIdilate(mri_dst, mri_dst) ;
  MRIdilate(mri_dst, mri_dst) ;
#endif
  return(mri_dst) ;
}


int
LocalMRISfindExpansionRegions(MRI_SURFACE *mris) {
  int    vno, num, n, num_long, total ;
  float  d, dsq, mean, std, dist ;

  d = dsq = 0.0f ;
  for (total = num = vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->ripflag || v->val <= 0)
      continue ;
    num++ ;
    dist = fabs(v->d) ;
    d += dist ;
    dsq += (dist*dist) ;
  }

  mean = d / num ;
  std = sqrt(dsq/num - mean*mean) ;
  fprintf(stderr, "mean absolute distance = %2.2f +- %2.2f\n", mean, std) ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    v->curv = 0 ;
    if (v->ripflag || v->val <= 0)
      continue ;
    if (fabs(v->d) < mean+2*std)
      continue ;
    for (num_long = num = 1, n = 0 ; n < vt->vnum ; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->val <= 0 || v->ripflag)
        continue ;
      if (fabs(vn->d) >= mean+2*std)
        num_long++ ;
      num++ ;
    }

    if ((float)num_long / (float)num > 0.25) {
      v->curv = fabs(v->d) ;
      total++ ;
#if 0
      fprintf(stderr, "v %d long: (%d of %d)\n", vno, num_long, num) ;
#endif
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d vertices more than 2 sigmas from mean.\n", total) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteCurvature(mris, "long") ;
  return(NO_ERROR) ;
}

// mris here is the one in mri_hires space
static float
check_contrast_direction(MRI_SURFACE *mris,MRI *mri_hires) {
  int     vno, n ;
  VERTEX  *v ;
  double  x, y, z, xw, yw, zw, val, mean_inside, mean_outside ;

  mean_inside = mean_outside = 0.0 ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag != 0)
      continue ;

    x = v->x+0.5*v->nx ;
    y = v->y+0.5*v->ny ;
    z = v->z+0.5*v->nz ;
    // check new function

#if 0
    if (mris->useRealRAS)
      MRIworldToVoxel(mri_hires, x, y, z, &xw, &yw, &zw);
    else
      MRIsurfaceRASToVoxel(mri_hires, x, y, z, &xw, &yw, &zw);
#else
    MRISsurfaceRASToVoxelCached(mris, mri_hires, x, y, z, &xw, &yw, &zw) ;
#endif
    MRIsampleVolume(mri_hires, xw, yw, zw, &val) ;
    mean_outside += val ;

    x = v->x-0.5*v->nx ;
    y = v->y-0.5*v->ny ;
    z = v->z-0.5*v->nz ;
#if 0
    if (mris->useRealRAS)
      MRIworldToVoxel(mri_hires, x, y, z, &xw, &yw, &zw);
    else
      MRIsurfaceRASToVoxel(mri_hires, x, y, z, &xw, &yw, &zw);
#else
    MRISsurfaceRASToVoxelCached(mris, mri_hires, x, y, z, &xw, &yw, &zw) ;
#endif
    MRIsampleVolume(mri_hires, xw, yw, zw, &val) ;
    mean_inside += val ;
    n++ ;
  }
  mean_inside /= (float)n ;
  mean_outside /= (float)n ;
  printf("mean inside = %2.1f, mean outside = %2.1f\n", 
         mean_inside, mean_outside) ;
  return(mean_inside - mean_outside) ;
}

#if 0
static int
MRIScomputeClassStatistics(MRI_SURFACE *mris, 
                           MRI *mri, 
                           float *pwhite_mean, 
                           float *pwhite_std, 
                           float *pgray_mean, 
                           float *pgray_std) 
{
  double  val, x, y, z, xw, yw, zw ;
  int     total_vertices, vno ;
  VERTEX  *v ;
  double  mean_white, mean_gray, std_white, std_gray, nsigma, gw_thresh  ;
  FILE    *fpwm, *fpgm ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    fpwm = fopen("wm.dat", "w") ;
    fpgm = fopen("gm.dat", "w") ;
  } else
    fpgm = fpwm = NULL ;

  std_white = std_gray = mean_white = mean_gray = 0.0 ;
  for (total_vertices = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    total_vertices++ ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = v->x+1.0*v->nx ;
    y = v->y+1.0*v->ny ;
    z = v->z+1.0*v->nz ;
    MRISsurfaceRASToVoxel(mris, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri, xw, yw, zw, &val) ;
    if (fpgm)
      fprintf(fpgm, "%d %2.1f %2.1f %2.1f %f\n", vno, xw, yw, zw, val) ;
    mean_gray += val ;
    std_gray += (val*val) ;

    x = v->x-0.5*v->nx ;
    y = v->y-0.5*v->ny ;
    z = v->z-0.5*v->nz ;
    MRISsurfaceRASToVoxel(mris, mri, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri, xw, yw, zw, &val) ;
    if (fpwm)
      fprintf(fpwm, "%d %2.1f %2.1f %2.1f %f\n", vno, xw, yw, zw, val) ;
    mean_white += val ;
    std_white += (val*val) ;
  }

  *pwhite_mean = mean_white /= (float)total_vertices ;
  *pwhite_std = std_white = 
    sqrt(std_white / (float)total_vertices - mean_white*mean_white) ;
  *pgray_mean = mean_gray /= (float)total_vertices ;
  *pgray_std = std_gray = 
    sqrt(std_gray / (float)total_vertices - mean_gray*mean_gray) ;
  nsigma = (mean_gray-mean_white) / (std_gray+std_white) ;
  gw_thresh = mean_white + nsigma*std_white ;
  printf("white %2.1f +- %2.1f,    gray %2.1f +- %2.1f, "
         "G/W boundary at %2.1f\n",
         mean_white, std_white, mean_gray, std_gray, gw_thresh) ;

  if (fpwm) {
    fclose(fpgm) ;
    fclose(fpwm) ;
  }
  return(NO_ERROR) ;
}
static int
compute_border_gradients(MRI_SURFACE *mris, MRI *mri, double thresh)
{
  int    vno, n ;
  VERTEX *v, *vn ;
  double  dif, max_dif ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag || v->marked == 0)
      continue ;
    max_dif = 0 ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag)
        continue ;
      dif = fabs(vn->val-v->val) ;
      if (dif > max_dif)
        max_dif = dif ;
    }
    v->val2 = max_dif ;
    if (max_dif > thresh)
      DiagBreak() ;
  }

  MRIScopyValuesToImagValues(mris) ;
  MRIScopyVal2ToVal(mris) ;
  MRISwriteValues(mris, "border_grad.mgz") ;
  MRIScopyImagValuesToValues(mris);
  MRISwriteValues(mris, "border.mgz") ;
  return(NO_ERROR) ;
}
#endif
static MRI *
find_wm(MRI_SURFACE *mris, MRI *mri, MRI *mri_wm)
{
  MRI       *mri_interior, *mri_ctrl, *mri_kernel ;
  int       vno, x, y, z ;
  VERTEX    *v ;
  double    xv, yv, zv, val ;
  
  return(NULL) ;
#define SDIST 0.5
  mri_interior = MRIclone(mri, NULL) ;
  mri_ctrl = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR) ;
  MRIcopyHeader(mri, mri_ctrl) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    MRISsurfaceRASToVoxel(mris, mri, 
                          v->x-SDIST*v->nx, v->y-SDIST*v->ny, 
                          v->z-SDIST*v->nz, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    for (x = floor(xv) ; x <= ceil(xv) ; x++)
      for (y = floor(yv) ; y <= ceil(yv) ; y++)
        for (z = floor(zv) ; z <= ceil(zv) ; z++)
        {
          if (MRIindexNotInVolume(mri, x, y, z) == 0)
          {
            MRIsetVoxVal(mri_ctrl, x, y, z, 0, 1) ;
            MRIsetVoxVal(mri_interior, x, y, z, 0, val) ;
          }
        }
  }

  printf("writing 1...\n") ;
  MRIwrite(mri_interior, "i1.mgz") ;
  MRIbuildVoronoiDiagram(mri, mri_ctrl, mri_interior) ;
  printf("writing 2...\n") ;
  MRIwrite(mri_interior, "i2.mgz") ;
  MRIwrite(mri_ctrl, "c.mgz") ;
  mri_kernel = MRIgaussian1d(.25/mri->xsize, 0) ;
  mri_wm = MRIconvolveGaussian(mri_interior, NULL, mri_kernel) ;
  //  MRIsoapBubble(mri_interior, mri_ctrl, mri_interior, 5) ;
  printf("writing 3...\n") ;
  MRIwrite(mri_wm, "i3.mgz") ;

  MRIfree(&mri_interior) ; MRIfree(&mri_ctrl) ; MRIfree(&mri_kernel) ;
  return(mri_wm) ;
}

static double
pial_errfunc_gradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int    vno ;
  VERTEX *v ;
  double dist, sse, dx, dy, dz, l_external ;

  l_external = parms->l_external ;
  for (sse = 0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked == 0 || v->ripflag)
      continue ;
    dx = v->tx - v->x ;
    dy = v->ty - v->y ;
    dz = v->tz - v->z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (vno == Gdiag_no)
      printf("v %d: pial_grad d = (%2.3f, %2.3f, %2.3f)\n", vno, dx, dy, dz) ;
    v->dx += l_external * dx ;
    v->dy += l_external * dy ;
    v->dz += l_external * dz ;
    sse += dist ;
  }
  return(sse) ;
}

static double
pial_errfunc_sse(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int    vno ;
  VERTEX *v ;
  double dist, sse, dx, dy, dz ;

  for (sse = 0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked == 0 || v->ripflag)
      continue ;
    dx = v->tx - v->x ;
    dy = v->ty - v->y ;
    dz = v->tz - v->z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    sse += dist ;
  }
  
  return(sse) ;
}
static double
pial_errfunc_rms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double sse ;
  int    nv ;

  sse = pial_errfunc_gradient(mris, parms) ;
  nv = MRISmarkedVertices(mris) ;
  if (nv == 0)
    nv = 1 ;
  return(sqrt(sse)/nv) ;
}

static MRI *
make_pial_location_mask(MRI_SURFACE *mris, MRI *mri_hires, float max_thickness, MRI *mri_dst)
{
  mri_dst = MRIScomputeDistanceToSurface(mris, NULL, 0.5) ;

  MRIthresholdRangeInto(mri_dst, mri_dst, 0, max_thickness) ;
  MRIwrite(mri_dst, "dist.mgz") ;
  return(mri_dst) ;
}

