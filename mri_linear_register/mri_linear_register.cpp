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




#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "matrix.h"
#include "version.h"

const char         *Progname ;
static MORPH_PARMS  parms ;

static int invert_flag = 0 ;
static int voxel_coords = 0 ;
static int nopca = 0 ;
static int full_res = 0 ;

static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static MATRIX *initialize_transform(MRI *mri_in, MRI *mri_ref, MP *parms) ;
static MRI *find_cropping(MRI *mri_in, MRI *mri_ref, MP *parms) ;
static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in, MRI *mri_ref, MP *parms, MATRIX *m_L) ;
static int order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors) ;
static float window_size = 0 ;
static unsigned char thresh_low = 40 ;
static int binarize = 1 ;
static int check_crop_flag = 0 ;
static int use_gradient = 1 ;

static char *var_fname = NULL ;
static char *xform_mean_fname = NULL ;
static char *xform_covariance_fname = NULL ;

#if 0
static unsigned char thresh_hi = 120 ;
#endif

static MATRIX *pca_matrix(MATRIX *m_in_evectors, double in_means[3],
                          MATRIX *m_ref_evectors, double ref_means[3]) ;
static MATRIX *compute_pca(MRI *mri_in, MRI *mri_ref) ;
static int init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) ;
#if 0
static int init_translation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L);
#endif

static int nreductions = 1 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;

static double blur_sigma = 2.0f ;
static double l_priors = 1 ;  /* weighting for prior term (if used) */

/*
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

int
main(int argc, char *argv[]) {
  char         *ref_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
  MRI          *mri_ref, *mri_in, *mri_orig, *mri_in_red, *mri_ref_red,
  *mri_in_tmp, *mri_ref_tmp, *mri_ref_orig, *mri_in_orig ;
  int          ac, nargs, i, msec, minutes, seconds ;
  Timer start ;
  MATRIX       *m_L ;

  nargs = handleVersionOption(argc, argv, "mri_linear_register");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  parms.mri_crop = NULL ;
  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.dt = 5e-6 ;  /* was 5e-6 */
  parms.tol = 1e-3 ;
  parms.momentum = 0.8 ;
  parms.max_levels = MAX_LEVELS ;
  parms.factor = 1.0 ;
  parms.niterations = 25 ;
  Progname = argv[0] ;


  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <in brain> <template> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  ref_fname = argv[2] ;
  if (xform_mean_fname) {
    int   sno, nsubjects ;
    FILE  *fp ;

    parms.m_xform_mean = MatrixAsciiRead(xform_mean_fname, NULL) ;
    if (!parms.m_xform_mean)
      ErrorExit(Gerror, "%s: could not read parameter means from %s",
                Progname, xform_mean_fname) ;

    fp = fopen(xform_covariance_fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read covariances from %s",
                Progname, xform_covariance_fname) ;

    fscanf(fp, "nsubjects=%d", &nsubjects) ;
    printf("reading %d transforms...\n", nsubjects) ;

    parms.m_xforms = (MATRIX **)calloc(nsubjects, sizeof(MATRIX *)) ;
    if (!parms.m_xforms)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate array of %d xforms",
                Progname, nsubjects) ;
    for (sno = 0 ; sno < nsubjects ; sno++) {
      parms.m_xforms[sno] = MatrixAsciiReadFrom(fp, NULL) ;
      if (!parms.m_xforms[sno])
        ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %dth xform",
                  Progname, sno) ;

    }
    parms.m_xform_covariance = MatrixAsciiReadFrom(fp, NULL) ;
    if (!parms.m_xform_covariance)
      ErrorExit(Gerror, "%s: could not read parameter covariance from %s",
                Progname, xform_covariance_fname) ;
    fclose(fp) ;
    parms.l_priors = l_priors ;
    parms.nxforms = nsubjects ;
  }
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  start.reset() ;
  fprintf(stderr, "reading '%s'...\n", ref_fname) ;
  fflush(stderr) ;
  mri_ref = MRIread(ref_fname) ;
  if (!mri_ref)
    ErrorExit(ERROR_NOFILE, "%s: could not open reference volume %s.\n",
              Progname, ref_fname) ;
  if (mri_ref->type != MRI_UCHAR) {
    MRI *mri_tmp ;

    mri_tmp = MRIchangeType(mri_ref, MRI_UCHAR, 0.0, 0.999, FALSE) ;
    MRIfree(&mri_ref) ;
    mri_ref = mri_tmp ;
  }

  if (var_fname)  /* read in a volume of standard deviations */
  {
    MRI *mri_var, *mri_tmp ;

    fprintf(stderr, "reading '%s'...\n", var_fname) ;
    mri_var = MRIread(var_fname) ;
    if (!mri_var)
      ErrorExit(ERROR_NOFILE, "%s: could not open variance volume %s.\n",
                Progname, var_fname) ;
    mri_tmp = MRIconcatenateFrames(mri_ref, mri_var, NULL) ;
    MRIfree(&mri_var) ;
    MRIfree(&mri_ref) ;
    mri_ref = mri_tmp ;
  }
  fprintf(stderr, "reading '%s'...\n", in_fname) ;
  fflush(stderr) ;
  mri_orig = mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;
  if (mri_in->type != MRI_UCHAR) {
    MRI *mri_tmp ;

    mri_orig = mri_tmp = MRIchangeType(mri_in, MRI_UCHAR, 0.0, 0.999, FALSE) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  /* make sure they are the same size */
  if (mri_in->width  != mri_ref->width ||
      mri_in->height != mri_ref->height  ||
      mri_in->depth  != mri_ref->depth) {
    int  width, height, depth ;
    MRI  *mri_tmp ;

    width = MAX(mri_in->width, mri_ref->width) ;
    height = MAX(mri_in->height, mri_ref->height) ;
    depth = MAX(mri_in->depth, mri_ref->depth) ;
    mri_tmp = MRIalloc(width, height, depth, MRI_UCHAR) ;
    MRIextractInto(mri_in, mri_tmp, 0, 0, 0,
                   mri_in->width, mri_in->height, mri_in->depth, 0, 0, 0) ;
#if 0
    MRIfree(&mri_in) ;
#else
    parms.mri_in = mri_in ;
#endif
    mri_in = mri_orig = mri_tmp ;

    mri_tmp = MRIallocSequence(width, height,depth,MRI_UCHAR,mri_ref->nframes);
    MRIextractInto(mri_ref, mri_tmp, 0, 0, 0,
                   mri_ref->width, mri_ref->height, mri_ref->depth, 0, 0, 0) ;
#if 0
    MRIfree(&mri_ref) ;
#else
    parms.mri_in = mri_in ;
#endif
    mri_ref = mri_tmp ;
  }


  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz)) {
    MRI *mri_tmp ;

    fprintf(stderr, "translating second volume by (%2.1f, %2.1f, %2.1f)\n",
            tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  if (!FZERO(rzrot)) {
    MRI *mri_tmp ;

    fprintf(stderr,
            "rotating second volume by %2.1f degrees around Z axis\n",
            (float)DEGREES(rzrot)) ;
    mri_tmp = MRIrotateZ_I(mri_in, NULL, rzrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(rxrot)) {
    MRI *mri_tmp ;

    fprintf(stderr,
            "rotating second volume by %2.1f degrees around X axis\n",
            (float)DEGREES(rxrot)) ;
    mri_tmp = MRIrotateX_I(mri_in, NULL, rxrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(ryrot)) {
    MRI *mri_tmp ;

    fprintf(stderr,
            "rotating second volume by %2.1f degrees around Y axis\n",
            (float)DEGREES(ryrot)) ;
    mri_tmp = MRIrotateY_I(mri_in, NULL, ryrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  if (!transform_loaded)   /* wasn't preloaded */
    parms.lta = LTAalloc(1, mri_in) ;

  if (!FZERO(blur_sigma)) {
    MRI *mri_kernel, *mri_tmp ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    mri_in = mri_tmp ;
    MRIfree(&mri_kernel) ;
  }
  MRIscaleMeanIntensities(mri_in, mri_ref, mri_in);

  mri_ref_orig = mri_ref ;
  mri_in_orig = mri_in ;
  if (nreductions > 0) {
    mri_in_red = mri_in_tmp = MRIcopy(mri_in, NULL) ;
    mri_ref_red = mri_ref_tmp = MRIcopy(mri_ref, NULL) ;
    for (i = 0 ; i < nreductions ; i++) {
      mri_in_red = MRIreduceByte(mri_in_tmp, NULL) ;
      mri_ref_red = MRIreduceMeanAndStdByte(mri_ref_tmp,NULL);
      MRIfree(&mri_in_tmp);
      MRIfree(&mri_ref_tmp) ;
      mri_in_tmp = mri_in_red ;
      mri_ref_tmp = mri_ref_red ;
    }
    mri_in = mri_in_red ;
    mri_ref = mri_ref_red ;
  }
  /* for diagnostics */
  if (full_res) {
    parms.mri_ref = mri_ref ;
    parms.mri_in = mri_in ;
  } else {
    parms.mri_ref = mri_ref_orig ;
    parms.mri_in = mri_in_orig ;
  }

  m_L = initialize_transform(mri_in, mri_ref, &parms) ;

  if (use_gradient) {
    MRI  *mri_in_mag, *mri_ref_mag, *mri_grad, *mri_mag ;

    printf("computing gradient magnitude of input image...\n") ;
    mri_mag = MRIalloc(mri_in->width, mri_in->height, mri_in->depth,MRI_FLOAT);
    MRIcopyHeader(mri_in, mri_mag) ;
    mri_grad = MRIsobel(mri_in, NULL, mri_mag) ;
    MRIfree(&mri_grad) ;

    /* convert it to ubytes */
    MRIvalScale(mri_mag, mri_mag, 0.0f, 255.0f) ;
    mri_in_mag = MRIclone(mri_in, NULL) ;
    MRIcopy(mri_mag, mri_in_mag) ;
    MRIfree(&mri_mag) ;

    /* now compute gradient of ref image */
    printf("computing gradient magnitude of reference image...\n") ;
    mri_mag = MRIalloc(mri_ref->width, mri_ref->height, mri_ref->depth,MRI_FLOAT);
    MRIcopyHeader(mri_ref, mri_mag) ;
    mri_grad = MRIsobel(mri_ref, NULL, mri_mag) ;
    MRIfree(&mri_grad) ;

    /* convert it to ubytes */
    MRIvalScale(mri_mag, mri_mag, 0.0f, 255.0f) ;
    mri_ref_mag = MRIclone(mri_ref, NULL) ;
    MRIcopy(mri_mag, mri_ref_mag) ;
    MRIfree(&mri_mag) ;

    register_mri(mri_in_mag, mri_ref_mag, &parms, m_L) ;
    MRIfree(&mri_in_mag) ;
    MRIfree(&mri_ref_mag) ;
  }
  register_mri(mri_in, mri_ref, &parms, m_L) ;
  if (check_crop_flag)  /* not working yet! */
  {
    printf("searching for cropped regions in the input image...\n") ;
    parms.mri_crop = find_cropping(mri_orig, mri_ref, &parms) ;
    MRIwrite(parms.mri_crop, "crop.mgh") ;
    register_mri(mri_in, mri_ref, &parms, m_L) ;
  }

  if (voxel_coords) {
    printf("transforming xform to voxel coordinates...\n") ;
    MRIrasXformToVoxelXform(mri_in_orig, mri_ref_orig,
                            parms.lta->xforms[0].m_L,
                            parms.lta->xforms[0].m_L);
    if (Gdiag & DIAG_WRITE) {
      MRI *mri_tmp ;

      mri_tmp = MRIlinearTransform(mri_in_orig, NULL,parms.lta->xforms[0].m_L);
      MRIwriteImageViews(mri_tmp, "morphed", IMAGE_SIZE) ;
      MRIfree(&mri_tmp) ;
    }
  }
  // save src and target info in lta
  getVolGeom(mri_in_orig, &parms.lta->xforms[0].src);
  getVolGeom(mri_ref_orig, &parms.lta->xforms[0].dst);
  fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
  if (invert_flag) {
    MATRIX *m_tmp ;

    m_tmp = MatrixInverse(parms.lta->xforms[0].m_L, NULL) ;
    MatrixFree(&parms.lta->xforms[0].m_L) ;
    // change src and dst
    getVolGeom(mri_in_orig, &parms.lta->xforms[0].dst);
    getVolGeom(mri_ref_orig, &parms.lta->xforms[0].src);
    parms.lta->xforms[0].m_L = m_tmp ;
  }
  //
  LTAwriteEx(parms.lta, out_fname) ;
  //
  if (mri_ref)
    MRIfree(&mri_ref) ;
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "registration took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  StrUpper(option) ;
  if (!stricmp(option, "DIST") || !stricmp(option, "DISTANCE")) {
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.2f\n", parms.l_dist) ;
  } else if (!stricmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt = %2.2e\n", parms.dt) ;
  } else if (!stricmp(option, "INVERT")) {
    invert_flag = 1 ;
    fprintf(stderr, "inverting transform before writing...\n") ;
  } else if (!stricmp(option, "crop")) {
    check_crop_flag = 1 ;
    nargs = 1 ;
    fprintf(stderr, "checking for cropping....\n") ;
  } else if (!stricmp(option, "nlevels")) {
    parms.max_levels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "nlevels = %d\n", parms.max_levels) ;
  } else if (!stricmp(option, "image_size")) {
    IMAGE_SIZE = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting default image size to %d\n", IMAGE_SIZE) ;
  } else if (!stricmp(option, "TOL")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  } else if (!stricmp(option, "NUM")) {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "finding a total of %d linear transforms\n", num_xforms) ;
  } else if (!stricmp(option, "SCOUT")) {
    parms.scout_flag = 1 ;
    printf("limitting domain of integration to central slices...\n") ;
  } else if (!stricmp(option, "AREA")) {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.2f\n", parms.l_area) ;
  } else if (!stricmp(option, "WINDOW")) {
    window_size = atof(argv[2]) ;
    fprintf(stderr, "applying Hanning window (R=%2.1f) to images...\n",
            window_size) ;
    nargs = 1 ;
  } else if (!stricmp(option, "NLAREA")) {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nlarea = %2.2f\n", parms.l_nlarea) ;
  } else if (!stricmp(option, "LEVELS")) {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "levels = %d\n", parms.levels) ;
  } else if (!stricmp(option, "INTENSITY") || !stricmp(option, "CORR")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.2f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "thresh")) {
    thresh_low = atoi(argv[2]) ;
#if 1
    fprintf(stderr, "setting threshold to %d\n", thresh_low) ;
    nargs = 1 ;
#else
    thresh_hi = atoi(argv[3]) ;
    fprintf(stderr, "thresholds set to %d --> %d\n", thresh_low, thresh_hi) ;
    nargs = 2 ;
#endif
  } else if (!stricmp(option, "reduce")) {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "reducing input images %d times before aligning...\n",
            nreductions) ;
  } else if (!stricmp(option, "priors")) {
    l_priors = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using %2.2f as weight for prior term\n", l_priors) ;
  } else if (!stricmp(option, "voxel")) {
    voxel_coords = 1 ;
    fprintf(stderr, "outputting transform in voxel coordinates\n") ;
  } else if (!stricmp(option, "full_res")) {
    full_res = 1 ;
    fprintf(stderr, "outputting full resolution images\n") ;
  } else if (!stricmp(option, "nopca")) {
    nopca = 1 ;
    fprintf(stderr, "disabling pca\n") ;
  } else if (!stricmp(option, "factor")) {
    parms.factor = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using time step factor of %2.2f\n",parms.factor) ;
  } else switch (*option) {
    case 'X':
      xform_mean_fname = argv[2] ;
      xform_covariance_fname = argv[3] ;
      printf("reading means (%s) and covariances (%s) of xforms\n",
             xform_mean_fname, xform_covariance_fname) ;
      nargs = 2 ;
      break ;
    case 'D':
      tx = atof(argv[2]) ;
      ty = atof(argv[3]) ;
      tz = atof(argv[4]) ;
      nargs = 3 ;
      break ;
    case 'R':
      rxrot = RADIANS(atof(argv[2])) ;
      ryrot = RADIANS(atof(argv[3])) ;
      rzrot = RADIANS(atof(argv[4])) ;
      nargs = 3 ;
      break ;
    case 'T':
      parms.lta = LTAread(argv[2]) ;
      if (!parms.lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                  Progname, argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using previously computed transform %s\n", argv[2]) ;
      transform_loaded = 1 ;
      break ;
    case 'B':
      blur_sigma = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "blurring input image with sigma=%2.3f\n", blur_sigma);
      break ;
    case 'S':
      parms.sigma = atof(argv[2]) ;
      fprintf(stderr, "using sigma=%2.3f as upper bound on blurring.\n",
              parms.sigma) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      printf("usage: %s <in volume> <template volume> <output transform>\n",
             argv[0]) ;
      exit(1) ;
      break ;
    case 'V':
      var_fname = argv[2] ;
      fprintf(stderr, "reading variance image from %s...\n", var_fname) ;
      nargs = 1 ;
      break ;
    case 'N':
      parms.niterations = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "niterations = %d\n", parms.niterations) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

#if 1
static int
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MATRIX *m_L) {
  MRI     *mri_in_windowed, *mri_ref_windowed ;

  fprintf(stderr, "aligning volume with average...\n") ;

  if (window_size > 0) {
    double in_means[3], ref_means[3] ;

    MRIcenterOfMass(mri_in, in_means, 0) ;
    MRIcenterOfMass(mri_ref, ref_means, 0) ;
    printf("windowing ref around (%d, %d, %d) and input around (%d, %d, %d)\n",
           nint(ref_means[0]), nint(ref_means[1]), nint(ref_means[2]),
           nint(in_means[0]), nint(in_means[1]), nint(in_means[2])) ;
    mri_in_windowed =
      MRIwindow(mri_in, NULL, WINDOW_HANNING,nint(in_means[0]),
                nint(in_means[1]), nint(in_means[2]),window_size);
    mri_ref_windowed =
      MRIwindow(mri_ref,NULL,WINDOW_HANNING,nint(ref_means[0]),
                nint(ref_means[1]), nint(ref_means[2]),window_size);
    mri_in = mri_in_windowed ;
    mri_ref = mri_ref_windowed ;
  }

  MRIrigidAlign(mri_in, mri_ref, parms, m_L) ;

  fprintf(stderr, "final transform:\n") ;
  MatrixPrint(stderr, parms->lta->xforms[0].m_L) ;
  fprintf(stderr, "\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri_aligned ;

    mri_aligned =
      MRIapplyRASlinearTransform(mri_in, NULL, parms->lta->xforms[0].m_L) ;
    MRIwriteImageViews(mri_aligned, "after_alignment", IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }


  MatrixCopy(parms->lta->xforms[0].m_L, m_L) ;
  return(NO_ERROR) ;
}
static MATRIX *
compute_pca(MRI *mri_in, MRI *mri_ref) {
  int    row, col, i ;
  float  dot ;
  MATRIX *m_ref_evectors = NULL, *m_in_evectors = NULL ;
  float  in_evalues[3], ref_evalues[3] ;
  double  ref_means[3], in_means[3] ;

  if (!m_ref_evectors)
    m_ref_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  if (!m_in_evectors)
    m_in_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;

  if (binarize) {
    MRIbinaryPrincipleComponents(mri_ref, m_ref_evectors, ref_evalues,
                                 ref_means, thresh_low);
    MRIbinaryPrincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,
                                 thresh_low);
  } else {
    MRIprincipleComponents(mri_ref, m_ref_evectors, ref_evalues, ref_means,
                           thresh_low);
    MRIprincipleComponents(mri_in,m_in_evectors,in_evalues,in_means,
                           thresh_low);
  }

  order_eigenvectors(m_in_evectors, m_in_evectors) ;
  order_eigenvectors(m_ref_evectors, m_ref_evectors) ;

  /* check to make sure eigenvectors aren't reversed */
  for (col = 1 ; col <= 3 ; col++) {
#if 0
    float theta ;
#endif

    for (dot = 0.0f, row = 1 ; row <= 3 ; row++)
      dot += m_in_evectors->rptr[row][col] * m_ref_evectors->rptr[row][col] ;

    if (dot < 0.0f) {
      fprintf(stderr, "WARNING: mirror image detected in eigenvector #%d\n",
              col) ;
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        m_in_evectors->rptr[row][col] *= -1.0f ;
    }
#if 0
    theta = acos(dot) ;
    fprintf(stderr, "angle[%d] = %2.1f\n", col, DEGREES(theta)) ;
#endif
  }
  fprintf(stderr, "ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  fprintf(stderr, "\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;

  return(pca_matrix(m_in_evectors, in_means,m_ref_evectors, ref_means)) ;
}

static int
order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors) {
  int    row, col, xcol, ycol, zcol ;
  double mx ;

  if (m_src_evectors == m_dst_evectors)
    m_src_evectors = MatrixCopy(m_src_evectors, NULL) ;

  /* find columx with smallest dot product with unit x vector */
  mx = fabs(*MATRIX_RELT(m_src_evectors, 1, 1)) ;
  xcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 1, col)) > mx) {
      xcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 1, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 2, 1)) ;
  ycol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (*MATRIX_RELT(m_src_evectors, 2, col) > mx) {
      ycol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 2, col)) ;
    }

  mx = fabs(*MATRIX_RELT(m_src_evectors, 3, 1)) ;
  zcol = 1 ;
  for (col = 2 ; col <= 3 ; col++)
    if (fabs(*MATRIX_RELT(m_src_evectors, 3, col)) > mx) {
      zcol = col ;
      mx = fabs(*MATRIX_RELT(m_src_evectors, 3, col)) ;
    }

  for (row = 1 ; row <= 3 ; row++) {
    *MATRIX_RELT(m_dst_evectors,row,1) = *MATRIX_RELT(m_src_evectors,row,xcol);
    *MATRIX_RELT(m_dst_evectors,row,2) = *MATRIX_RELT(m_src_evectors,row,ycol);
    *MATRIX_RELT(m_dst_evectors,row,3) = *MATRIX_RELT(m_src_evectors,row,zcol);
  }
  return(NO_ERROR) ;
}


static MATRIX *
pca_matrix(MATRIX *m_in_evectors, double in_means[3],
           MATRIX *m_ref_evectors, double ref_means[3]) {
  float   dx, dy, dz ;
  MATRIX  *mRot, *m_in_T, *mOrigin, *m_L, *m_R, *m_T, *m_tmp ;
  double  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy ;
  int     row, col ;

  m_in_T = MatrixTranspose(m_in_evectors, NULL) ;
  mRot = MatrixMultiply(m_ref_evectors, m_in_T, NULL) ;

  r11 = mRot->rptr[1][1] ;
  r21 = mRot->rptr[2][1] ;
  r31 = mRot->rptr[3][1] ;
  r32 = mRot->rptr[3][2] ;
  r33 = mRot->rptr[3][3] ;
  y_angle = atan2(-r31, sqrt(r11*r11+r21*r21)) ;
  cosy = cos(y_angle) ;
  z_angle = atan2(r21 / cosy, r11 / cosy) ;
  x_angle = atan2(r32 / cosy, r33 / cosy) ;

#define MAX_X_ANGLE  (RADIANS(35))
#define MAX_Y_ANGLE  (RADIANS(15))
#define MAX_Z_ANGLE  (RADIANS(15))
  if (fabs(x_angle) > MAX_X_ANGLE || fabs(y_angle) > MAX_Y_ANGLE ||
      fabs(z_angle) > MAX_Z_ANGLE) {
    MATRIX *m_I ;

    /*    MatrixFree(&m_in_T) ; MatrixFree(&mRot) ;*/
    fprintf(stderr,
            "eigenvector swap detected (%2.0f, %2.0f, %2.0f): ignoring rotational PCA...\n",
            DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

    m_I = MatrixIdentity(3, NULL) ;
    MatrixCopy(m_I, mRot) ;
    MatrixFree(&m_I) ;
    x_angle = y_angle = z_angle = 0.0 ;
  }

  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = ref_means[0] ;
  mOrigin->rptr[2][1] = ref_means[1] ;
  mOrigin->rptr[3][1] = ref_means[2] ;

  fprintf(stderr, "reference volume center of mass at (%2.1f,%2.1f,%2.1f)\n",
          ref_means[0], ref_means[1], ref_means[2]) ;
  fprintf(stderr, "input volume center of mass at     (%2.1f,%2.1f,%2.1f)\n",
          in_means[0], in_means[1], in_means[2]) ;
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;

  fprintf(stderr, "translating volume by %2.1f, %2.1f, %2.1f\n",
          dx, dy, dz) ;
  fprintf(stderr, "rotating volume by (%2.2f, %2.2f, %2.2f)\n",
          DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;

  /* build full rigid transform */
  m_R = MatrixAlloc(4,4,MATRIX_REAL) ;
  m_T = MatrixAlloc(4,4,MATRIX_REAL) ;
  for (row = 1 ; row <= 3 ; row++) {
    for (col = 1 ; col <= 3 ; col++) {
      *MATRIX_RELT(m_R,row,col) = *MATRIX_RELT(mRot, row, col) ;
    }
    *MATRIX_RELT(m_T,row,row) = 1.0 ;
  }
  *MATRIX_RELT(m_R, 4, 4) = 1.0 ;

  /* translation so that origin is at ref eigenvector origin */
  dx = -ref_means[0] ;
  dy = -ref_means[1] ;
  dz = -ref_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ;
  *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ;
  *MATRIX_RELT(m_T, 4, 4) = 1 ;
  m_tmp = MatrixMultiply(m_R, m_T, NULL) ;
  *MATRIX_RELT(m_T, 1, 4) = -dx ;
  *MATRIX_RELT(m_T, 2, 4) = -dy ;
  *MATRIX_RELT(m_T, 3, 4) = -dz ;
  MatrixMultiply(m_T, m_tmp, m_R) ;

  /* now apply translation to take in centroid to ref centroid */
  dx = ref_means[0] - in_means[0] ;
  dy = ref_means[1] - in_means[1] ;
  dz = ref_means[2] - in_means[2] ;
  *MATRIX_RELT(m_T, 1, 4) = dx ;
  *MATRIX_RELT(m_T, 2, 4) = dy ;
  *MATRIX_RELT(m_T, 3, 4) = dz ;
  *MATRIX_RELT(m_T, 4, 4) = 1 ;

  m_L = MatrixMultiply(m_R, m_T, NULL) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("m_T:\n") ;
    MatrixPrint(stdout, m_T) ;
    printf("m_R:\n") ;
    MatrixPrint(stdout, m_R) ;
    printf("m_L:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  MatrixFree(&m_R) ;
  MatrixFree(&m_T) ;

  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(m_L) ;
}
#else
static int
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms) {
  MRI  *mri_in_red, *mri_ref_red ;

  mri_in_red = MRIreduceByte(mri_in, NULL) ;
  mri_ref_red = MRIreduceMeanAndStdByte(mri_ref,NULL);

  /*  parms->write_iterations = 0 ; */
  if (!parms->niterations)
    parms->niterations = 1000 ;
  if (transform_loaded)  /* don't recompute rotation based on neck */
  {
    if (MRIfindNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, NULL,
                    -1,NULL) == NULL)
      ErrorExit(Gerror, "%s: could not find subject neck.\n", Progname) ;
    if (MRIfindNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, NULL,1,
                    NULL) == NULL)
      ErrorExit(Gerror, "%s: could not find neck in reference volume.\n",
                Progname) ;
  } else {
    if (MRIfindNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, parms, -1,
                    &parms->in_np) == NULL)
      ErrorExit(Gerror, "%s: could not find subject neck.\n", Progname) ;
    if (MRIfindNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, parms, 1,
                    &parms->ref_np) == NULL)
      ErrorExit(Gerror, "%s: could not find neck in reference volume.\n",
                Progname) ;
  }

if (full_res) {}
  else {
    parms->mri_ref = mri_ref_red ;
    parms->mri_in = mri_in_red ;  /* for diagnostics */
  }
  while (parms->lta->num_xforms < num_xforms)
    LTAdivide(parms->lta, mri_in_red) ;
  fprintf(stderr,"computing %d linear transformation%s...\n",
          parms->lta->num_xforms, parms->lta->num_xforms>1?"s":"");
  MRIlinearAlign(mri_in_red, mri_ref_red, parms) ;

  MRIfree(&mri_in_red) ;
  MRIfree(&mri_ref_red) ;
  return(NO_ERROR) ;
}
#endif

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_DX   1.2
#define MAX_DY   1.2
#define MAX_DZ   1.2
#define MIN_DX   (1.0/MAX_DX)
#define MIN_DY   (1.0/MAX_DY)
#define MIN_DZ   (1.0/MAX_DZ)
#define MAX_RATIO 1.2

static int
init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) {
  MATRIX      *m_scaling ;
  float       sx, sy, sz, dx, dy, dz ;
  MRI_REGION  in_bbox, ref_bbox ;

  m_scaling = MatrixIdentity(4, NULL) ;

  MRIboundingBox(mri_in, 60, &in_bbox) ;
  MRIboundingBox(mri_ref, 60, &ref_bbox) ;
  sx = (float)ref_bbox.dx / (float)in_bbox.dx ;
  sy = (float)ref_bbox.dy / (float)in_bbox.dy ;
  sz = (float)ref_bbox.dz / (float)in_bbox.dz ;
  dx = (ref_bbox.x+ref_bbox.dx-1)/2 - (in_bbox.x+in_bbox.dx-1)/2 ;
  dy = (ref_bbox.y+ref_bbox.dy-1)/2 - (in_bbox.y+in_bbox.dy-1)/2 ;
  dz = (ref_bbox.z+ref_bbox.dz-1)/2 - (in_bbox.z+in_bbox.dz-1)/2 ;

  if (sx > MAX_DX)
    sx = MAX_DX ;
  if (sx < MIN_DX)
    sx = MIN_DX ;
  if (sy > MAX_DY)
    sy = MAX_DY ;
  if (sy < MIN_DY)
    sy = MIN_DY ;
  if (sz > MAX_DZ)
    sz = MAX_DZ ;
  if (sz < MIN_DZ)
    sz = MIN_DZ ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "initial scaling: (%2.2f, %2.2f, %2.2f) <-- "
            "(%d/%d,%d/%d,%d/%d)\n",
            sx,sy,sz, ref_bbox.dx, in_bbox.dx, ref_bbox.dy, in_bbox.dy,
            ref_bbox.dz, in_bbox.dz) ;
  *MATRIX_RELT(m_scaling, 1, 1) = sx ;
  *MATRIX_RELT(m_scaling, 2, 2) = sy ;
  *MATRIX_RELT(m_scaling, 3, 3) = sz ;

#if 0
  *MATRIX_RELT(m_L, 1, 4) = dx ;
  *MATRIX_RELT(m_L, 2, 4) = dy ;
  *MATRIX_RELT(m_L, 3, 4) = dz ;
#endif
  MatrixMultiply(m_scaling, m_L, m_L) ;
  return(NO_ERROR) ;
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
init_translation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L) {
  MATRIX *m_translation ;
  float  in_means[4], ref_means[4] ;
  double dx, dy, dz ;

  m_translation = MatrixIdentity(4, NULL) ;
  MRIfindCenterOfBrain(mri_in, in_means, in_means+1, in_means+2) ;
  MRIfindCenterOfBrain(mri_ref, ref_means, ref_means+1, ref_means+2) ;
  dx = (double)(ref_means[0] - in_means[0]) * mri_in->thick ;
  dy = (double)(ref_means[1] - in_means[1]) * mri_in->thick ;
  dz = (double)(ref_means[2] - in_means[2]) * mri_in->thick ;
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr, "centering template around (%d,%d,%d) and input around"
            " (%d,%d,%d)\n",
            (int)ref_means[0], (int)ref_means[1], (int)ref_means[2],
            (int)in_means[0], (int)in_means[1], (int)in_means[2]) ;
    fprintf(stderr, "initial translation = (%2.0f, %2.0f, %2.0f).\n",dx,dy,dz);
  }
  *MATRIX_RELT(m_translation, 1, 4) = dx ;
  *MATRIX_RELT(m_translation, 2, 4) = dy ;
  *MATRIX_RELT(m_translation, 3, 4) = dz ;
  MatrixMultiply(m_translation, m_L, m_L) ;
  MatrixFree(&m_translation) ;
  return(NO_ERROR) ;
}
#endif

#define BORDER_LABEL   1
#define CROP_LABEL     2

static MRI *
find_cropping(MRI *mri_in, MRI *mri_ref, MP *parms) {
  MRI    *mri_crop, *mri_tmp ;
  int    x, y, z, xi, yi, zi, xk, yk, zk, width, height, depth, num_on,
  ncropped, nfilled ;

  ncropped = 0 ;
  mri_crop = MRIclone(mri_in, NULL) ;
  mri_tmp =
    MRIapplyRASinverseLinearTransform(mri_ref, NULL,parms->lta->xforms[0].m_L);
  if (!mri_tmp)
    return(NULL) ;
  mri_ref = mri_tmp ;   /* reference volume transformed into input space */

  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth;

  /*
     first find regions in which the source image has 0s next to non-zeros
     and the reference image is non-zero
  */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == 0 && y == 254 && z == 127)
          DiagBreak() ;


        if (x == 126 && y == 124 && z == 217)
          DiagBreak() ;
        if (x == 126 && y == 241 && z == 234)
          DiagBreak() ;
        if ((MRIvox(mri_in, x, y, z) > 0) || (MRIvox(mri_ref,x,y,z) == 0))
          continue ;

        /* find at least one non-zero voxel in input volume */
        num_on = 0 ;
        for (zk = -1 ; zk <= 1 ; zk++) {
          zi = mri_in->zi[z+zk] ;
          for (yk = -1 ; yk <= 1 ; yk++) {
            yi = mri_in->yi[y+yk] ;
            for (xk = -1 ; xk <= 1 ; xk++) {
              xi = mri_in->xi[x+xk] ;
              if (MRIvox(mri_in, xi, yi, zi) > 0) {
                num_on++ ;
                break ;
              }
            }
            if (num_on)
              break ;
          }
          if (num_on)
            break ;
        }
        if (num_on) {
          MRIvox(mri_crop, x, y, z) = BORDER_LABEL ;
          ncropped++ ;
        }
      }
    }
  }

  /* now do a flood fill outward from the seed points to include all
     connected voxels that are 0 in the input volume and non-zero in
     the reference volume.
  */

  do {
    nfilled = 0 ;
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == 0 && y == 254 && z == 127)
            DiagBreak() ;

          if (x == 126 && y == 124 && z == 217)
            DiagBreak() ;
          if (x == 126 && y == 241 && z == 234)
            DiagBreak() ;
          if ((MRIvox(mri_in, x, y, z) > 0) ||
              (MRIvox(mri_ref,x,y,z) <= 10) ||
              (MRIvox(mri_crop,x,y,z) > 0))
            continue ;

          /* find at least one non-zero voxel in input volume */
          num_on = 0 ;
          for (zk = -1 ; zk <= 1 ; zk++) {
            zi = mri_in->zi[z+zk] ;
            for (yk = -1 ; yk <= 1 ; yk++) {
              yi = mri_in->yi[y+yk] ;
              for (xk = -1 ; xk <= 1 ; xk++) {
                xi = mri_in->xi[x+xk] ;
                if (MRIvox(mri_crop, xi, yi, zi) > 0) {
                  num_on++ ;
                  break ;
                }
              }
              if (num_on)
                break ;
            }
            if (num_on)
              break ;
          }
          if (num_on) {
            MRIvox(mri_crop, x, y, z) = CROP_LABEL ;
            nfilled++ ;
          }
        }
      }
    }
    ncropped += nfilled ;
    fprintf(stderr, "%d cropped voxels detected - %d total\n", nfilled,
            ncropped) ;
  } while (nfilled > 0) ;


  /* now remove all border points that don't have at least one neighbor
     that is cropped and non-border.
  */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == 0 && y == 254 && z == 127)
          DiagBreak() ;


        if (x == 126 && y == 124 && z == 217)
          DiagBreak() ;
        if (x == 126 && y == 241 && z == 234)
          DiagBreak() ;
        if (MRIvox(mri_crop, x, y, z) != BORDER_LABEL)
          continue ;

        /* find at least one non-zero voxel in input volume */
        num_on = 0 ;
        for (zk = -1 ; zk <= 1 ; zk++) {
          zi = mri_in->zi[z+zk] ;
          for (yk = -1 ; yk <= 1 ; yk++) {
            yi = mri_in->yi[y+yk] ;
            for (xk = -1 ; xk <= 1 ; xk++) {
              xi = mri_in->xi[x+xk] ;
              if ((MRIvox(mri_crop, xi, yi, zi) > 0) &&
                  MRIvox(mri_crop, xi, yi, zi) != BORDER_LABEL) {
                num_on++ ;
                break ;
              }
            }
            if (num_on)
              break ;
          }
          if (num_on)
            break ;
        }
        if (num_on == 0) {
          MRIvox(mri_crop, x, y, z) = 0 ;
          ncropped-- ;
        }
      }
    }
  }

  printf("%d cropped points detected...\n", ncropped) ;
  MRIfree(&mri_ref) ; /* NOT the one passed in - the inverse transformed one */
  return(mri_crop) ;
}


static MATRIX *
initialize_transform(MRI *mri_in, MRI *mri_ref, MP *parms) {
  MATRIX *m_L ;

  fprintf(stderr, "initializing alignment using PCA...\n") ;
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
    MRIwriteImageViews(parms->mri_ref, "ref", IMAGE_SIZE) ;
    MRIwriteImageViews(parms->mri_in, "before_pca", IMAGE_SIZE) ;
  }

  if (nopca)
    m_L = MatrixIdentity(3, NULL) ;
  else
    m_L = compute_pca(mri_in, mri_ref) ;

  init_scaling(mri_in, mri_ref, m_L) ;
#if 0
  init_translation(mri_in, mri_ref, m_L) ; /* in case PCA failed */
#endif

  /* convert it to RAS mm coordinates */
  MRIvoxelXformToRasXform(mri_in, mri_ref, m_L, m_L) ;

  if (Gdiag & DIAG_SHOW) {
    printf("initial transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }
  if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
    MRI *mri_aligned ;

    mri_aligned = MRIapplyRASlinearTransform(parms->mri_in, NULL, m_L) ;
    MRIwriteImageViews(mri_aligned, "after_pca", IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }
  return(m_L) ;
}
