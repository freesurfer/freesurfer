

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
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "gcamorph.h"
#include "transform.h"

char         *Progname ;
static GCA_MORPH_PARMS  parms ;

static char *mask_fname = NULL ;
static char *norm_fname = NULL ;

static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static int register_wm_flag = 0 ;

static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static char *tl_fname = NULL ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;
static char *normalized_transformed_sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 0 ;
static int relabel = 0 ;

static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static double tol = 1 ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

static FILE *diag_fp = NULL ;

static int translation_only = 0 ;
static int get_option(int argc, char *argv[]) ;
static int write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname) ;

static char *renormalization_fname = NULL ;
static char *tissue_parms_fname = NULL ;
static int center = 1 ;
static int nreductions = 1 ;
static char *xform_name = NULL ;
static int noscale = 0 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;
static TRANSFORM  *transform = NULL ;
static char *vf_fname = NULL ;

static double blur_sigma = 0.0f ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

#define NPARMS           12
#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;

int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
  MRI          *mri_in ;
  GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
  int          ac, nargs ;
  int          msec, minutes, seconds/*, min_left_cbm, min_right_cbm*/ ;
  struct timeb start ;
  GCA_MORPH    *gcam ;

  parms.l_likelihood = 1.0f ;
  parms.niterations = 3 ;
  parms.levels = 3 ;   /* use default */
  parms.dt = 0.1 ;  /* was 5e-6 */
  parms.tol = tol ;  /* at least 1% decrease in sse */
  parms.l_distance = 0.0 ;
  parms.l_jacobian = 40.0 ;
  parms.l_area = 0 ;
	parms.l_label = 1.0 ;
	parms.l_map = 1.0 ;
	parms.label_dist = 3.0 ;
  parms.l_smoothness = 0.0 ;
  parms.start_t = 0 ;
  parms.max_grad = 0.3 ;
  parms.sigma = 1.0f ;
  parms.exp_k = 20 ;
	parms.navgs = 0 ;

  Progname = argv[0] ;
  setRandomSeed(-1L) ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM, 
              "usage: %s <in brain> <template> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  gca_fname = argv[2] ;
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  printf("logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  printf("reading '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;
  printf("freeing gibbs priors...") ;
  GCAfreeGibbs(gca) ;
  printf("done.\n") ;

  if (novar)
    GCAunifyVariance(gca) ;

  if (renormalization_fname)
  {
    FILE   *fp ;
    int    *labels, nlines, i ;
    float  *intensities, f1, f2 ;
    char   *cp, line[STRLEN] ;

    fp = fopen(renormalization_fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                Progname, renormalization_fname) ;

    cp = fgetl(line, 199, fp) ;
    nlines = 0 ;
    while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ; intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
        DiagBreak() ;
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ; free(intensities) ;
  }


  printf("reading '%s'...\n", in_fname) ;
  fflush(stdout) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;

  if (mask_fname)
  {
    MRI *mri_mask ;

    mri_mask = MRIread(mask_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, mask_fname) ;

    MRImask(mri_in, mri_mask, mri_in, 0, 0) ;
    MRIfree(&mri_mask) ;
  }

    
  if (alpha > 0)
    mri_in->flip_angle = alpha ;
  if (TR > 0)
    mri_in->tr = TR ;
  if (TE > 0)
    mri_in->te = TE ;

  if (example_T1)
  {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE,"%s: could not read example segmentation from %s",
                Progname, example_segmentation) ;
    mri_T1 = MRIread(example_T1) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE,"%s: could not read example T1 from %s",
                Progname, example_T1) ;
    printf("scaling atlas intensities using specified examples...\n") ;
    MRIeraseBorderPlanes(mri_seg) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ; MRIfree(&mri_T1) ;
  }

  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;

  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
  {
    MRI *mri_tmp ;
    
    printf("translating second volume by (%2.1f, %2.1f, %2.1f)\n",
            tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  if (!FZERO(rzrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around Z axis\n",
            (float)DEGREES(rzrot)) ;
    mri_tmp = MRIrotateZ_I(mri_in, NULL, rzrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(rxrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around X axis\n",
            (float)DEGREES(rxrot)) ;
    mri_tmp = MRIrotateX_I(mri_in, NULL, rxrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(ryrot))
  {
    MRI *mri_tmp ;
    
    printf(
            "rotating second volume by %2.1f degrees around Y axis\n",
            (float)DEGREES(ryrot)) ;
    mri_tmp = MRIrotateY_I(mri_in, NULL, ryrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }

  if (!transform_loaded)   /* wasn't preloaded */
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;

  if (!FZERO(blur_sigma))
  {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_in) ; mri_in = mri_tmp ;
  }


	if (xform_name)
	{
		gcam = GCAMread(xform_name) ;
		if (!gcam)
			ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_name) ;
	}
	else
		gcam = GCAMalloc(gca->prior_width, gca->prior_height, gca->prior_depth) ;
	if (tl_fname)
	{
		GCA *gca_tl ;

		gca_tl = GCAread(tl_fname) ;
		if (!gca_tl)
			ErrorExit(ERROR_NOFILE, "%s: could not temporal lobe gca %s",
								Progname, tl_fname) ;
		GCAMinit(gcam, mri_in, gca_tl, transform, 0) ;
		if (parms.write_iterations != 0)
		{
			char fname[STRLEN] ;
			MRI  *mri_gca ;
			mri_gca = MRIclone(mri_in, NULL) ;
			GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
			sprintf(fname, "%s_target", parms.base_name) ;
			MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
			sprintf(fname, "%s_target.mgh", parms.base_name) ;
			printf("writing target volume to %s...\n", fname) ;
			MRIwrite(mri_gca, fname) ;
			MRIfree(&mri_gca) ;
		}
		GCAMregister(gcam, mri_in, &parms) ;
		printf("temporal lobe registration complete - registering whole brain...\n") ;
		GCAfree(&gca_tl) ;
	}

	if (!xform_name)  /* only if the transform wasn't previously created */
		GCAMinit(gcam, mri_in, gca, transform, relabel) ;
	else
	{
		gcam->gca = gca ;
		GCAMcomputeOriginalProperties(gcam) ;
		GCAMcomputeMaxPriorLabels(gcam) ;
	}
  if (parms.write_iterations != 0)
  {
    char fname[STRLEN] ;
    MRI  *mri_gca ;
    mri_gca = MRIclone(mri_in, NULL) ;
    GCAMbuildMostLikelyVolume(gcam, mri_gca) ;
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgh", parms.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
    MRIwrite(mri_gca, fname) ;
    MRIfree(&mri_gca) ;
  }

#if 0
  {
    GCA_SAMPLE *gcas ;
    int  nsamples ;

    gcas = GCAfindAllSamples(gca, &nsamples) ;
    GCAtransformAndWriteSamples(gca, mri_in, gcas, nsamples, 
                                "gcas_fsamples.mgh", transform) ;
    free(gcas) ;
  }
#endif

	if (tl_fname == NULL && register_wm_flag)
	{
		GCAMsetStatus(gcam, GCAM_IGNORE_LIKELIHOOD) ; /* disable everything */
		GCAMsetLabelStatus(gcam, Left_Cerebral_White_Matter, GCAM_USE_LIKELIHOOD) ;
		GCAMsetLabelStatus(gcam, Right_Cerebral_White_Matter, GCAM_USE_LIKELIHOOD) ;
		GCAMsetLabelStatus(gcam, Left_Cerebellum_White_Matter, GCAM_USE_LIKELIHOOD) ;
		GCAMsetLabelStatus(gcam, Right_Cerebellum_White_Matter, GCAM_USE_LIKELIHOOD) ;

		printf("initial white matter registration...\n") ;
		GCAMregister(gcam, mri_in, &parms) ;
		GCAMsetStatus(gcam, GCAM_USE_LIKELIHOOD) ; /* disable everything */
		printf("initial white matter registration complete - full registration...\n") ;
	}
  GCAMregister(gcam, mri_in, &parms) ;

  printf("writing output transformation to %s...\n", out_fname) ;
  if (vf_fname)
    write_vector_field(mri_in, gcam, vf_fname) ;
  if (GCAMwrite(gcam, out_fname) != NO_ERROR)
    ErrorExit(Gerror, "%s: GCAMwrite(%s) failed", Progname, out_fname) ;

#if 0
  if (gca)
    GCAfree(&gca) ;
#endif
  GCAMfree(&gcam) ;
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  if (diag_fp)
    fclose(diag_fp) ;
  exit(0) ;
  return(0) ;
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
  StrUpper(option) ;
  if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE"))
  {
    parms.l_distance = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_distance) ;
  }
  else if (!strcmp(option, "SMOOTH") || !strcmp(option, "SMOOTHNESS"))
  {
    parms.l_smoothness = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_smoothness = %2.2f\n", parms.l_smoothness) ;
  }
  else if (!strcmp(option, "SAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "ISIZE") || !strcmp(option, "IMAGE_SIZE"))
  {
    IMAGE_SIZE = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting diagnostic image size to %d\n", IMAGE_SIZE) ;
  }
  else if (!strcmp(option, "WM"))
  {
		register_wm_flag = 1 ;
    printf("registering white matter in initial pass...\n") ;
  }
  else if (!strcmp(option, "TL"))
  {
    tl_fname = argv[2] ;
    nargs = 1 ;
    printf("reading temporal lobe atlas from %s...\n", tl_fname) ;
  }
  else if (!strcmp(option, "RELABEL"))
  {
    relabel = atoi(argv[2]) ;
    nargs = 1 ;
    printf("%srelabeling nodes with MAP estimates\n", 
					 relabel ? "" : "not ") ;
  }
  else if (!strcmp(option, "VF"))
  {
    vf_fname = argv[2] ;
    nargs = 1 ;
    printf("writing vector field to %s...\n", vf_fname) ;
  }
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "DIAG"))
  {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open diag file %s for writing",
                Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!strcmp(option, "EXAMPLE"))
  {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  else if (!strcmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!strcmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!strcmp(option, "FSAMPLES") || !strcmp(option, "ISAMPLES"))
  {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transformed control points to %s...\n", 
            transformed_sample_fname) ;
  }
  else if (!strcmp(option, "NSAMPLES"))
  {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n", 
            normalized_transformed_sample_fname) ;
  }
  else if (!strcmp(option, "CONTRAST"))
  {
    use_contrast = 1 ;
    printf("using contrast to find labels...\n") ;
  }
  else if (!strcmp(option, "RENORM"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!strcmp(option, "FLASH"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    translation_only = 1 ;
    printf("only computing translation parameters...\n") ;
  }
  else if (!strcmp(option, "WRITE_MEAN"))
  {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    printf("writing gca means to %s...\n", gca_mean_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "TOL"))
  {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "CENTER"))
  {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  }
  else if (!strcmp(option, "NOSCALE"))
  {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  }
  else if (!strcmp(option, "LEVELS"))
  {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  }
  else if (!strcmp(option, "LIKELIHOOD"))
  {
    parms.l_likelihood = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_likelihood = %2.2f\n", parms.l_likelihood) ;
  }
  else if (!strcmp(option, "LABEL"))
  {
    parms.l_label = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_label = %2.2f\n", parms.l_label) ;
  }
  else if (!strcmp(option, "MAP"))
  {
    parms.l_map = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_map = %2.2f\n", parms.l_map) ;
  }
  else if (!strcmp(option, "LDIST") || !strcmp(option, "LABEL_DIST"))
  {
    parms.label_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("label_dist = %2.2f\n", parms.label_dist) ;
  }
  else if (!stricmp(option, "REDUCE"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
            nreductions) ;
  }
  else if (!stricmp(option, "DEBUG_NODE"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gvx = atoi(argv[2]) ;
    Gvy = atoi(argv[3]) ;
    Gvz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gvx, Gvy, Gvz) ;
  }
  else if (!stricmp(option, "norm"))
  {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  }
  else if (!stricmp(option, "avgs"))
  {
		parms.navgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing gradient with %d averages...\n", parms.navgs) ;
  }
  else switch (*option)
  {
  case 'J':
    parms.l_jacobian = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_jacobian=%2.3f\n", parms.l_jacobian) ;
    break ;
  case 'A':
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("using l_area=%2.3f\n", parms.l_area) ;
    break ;
  case 'F':
    ctl_point_fname = argv[2] ;
    nargs = 1 ;
    printf("reading manually defined control points from %s\n", ctl_point_fname) ;
    break ;
	case 'X':
		xform_name = argv[2] ;
		nargs = 1 ;
		printf("reading previous transform from %s...\n", xform_name) ;
		break ;
  case 'D':
    tx = atof(argv[2]) ; ty = atof(argv[3]) ; tz = atof(argv[4]) ;
    nargs = 3 ;
    break ;
  case 'K':
    parms.exp_k = atof(argv[2]) ;
    printf("setting exp_k to %2.2f (default=%2.2f)\n",
           parms.exp_k, EXP_K) ;
    nargs = 1 ;
    break ;
  case 'R':
    rxrot = RADIANS(atof(argv[2])) ; 
    ryrot = RADIANS(atof(argv[3])) ;
    rzrot = RADIANS(atof(argv[4])) ; 
    nargs = 3 ;
    break ;
  case 'T':
    transform = TransformRead(argv[2]) ;
    if (!transform)
      ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                Progname, argv[2]) ;
    nargs = 1 ;
    printf("using previously computed transform %s\n", argv[2]) ;
    transform_loaded = 1 ;
    break ;
  case 'B':
    blur_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("blurring input image with sigma=%2.3f\n", blur_sigma);
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'S':
    parms.sigma = atof(argv[2]) ;
    printf("using sigma=%2.3f as upper bound on blurring.\n", 
            parms.sigma) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s <in volume> <template volume> <output transform>\n", 
           argv[0]) ;
    exit(1) ;
    break ;
  case 'N':
    parms.niterations = atoi(argv[2]) ;
    nargs = 1 ;
    printf("niterations = %d\n", parms.niterations) ;
    break ;
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    printf("write iterations = %d\n", parms.write_iterations) ;
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'P':
    ctl_point_pct = atof(argv[2]) ;
    nargs = 1 ;
    printf("using top %2.1f%% wm points as control points....\n",
           100.0*ctl_point_pct) ;
    break ;
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("momentum = %2.2f\n", parms.momentum) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
static int
write_vector_field(MRI *mri, GCA_MORPH *gcam, char *vf_fname)
{
  FILE            *fp ;
  int             x, y, z ;
  GCA_MORPH_NODE  *gcamn ;

  fp = fopen(vf_fname, "w") ;

  for (x = 0 ; x < gcam->width ; x++)
  {
    for (y = 0 ; y < gcam->height ; y++)
    {
      for (z = 0 ; z < gcam->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        gcamn = &gcam->nodes[x][y][z] ;
        fprintf(fp, "%f %f %f %f\n", 
                gcamn->x-gcamn->origx,
                gcamn->y-gcamn->origy,
                gcamn->z-gcamn->origz, 
                gcamn->mean) ;
      }
    }
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

