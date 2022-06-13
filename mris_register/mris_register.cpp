/**
 * @brief cmd line utility for registering a subject's surface with an atlas,
 *
 * Command line utility for registering a subject's surface with an atlas.
 * Can also be used with -1 to register a surface to another individual surfce.
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
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "macros.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"

#include "romp_support.h"

#include "timer.h"
#include "error.h"
#include "diag.h"
#include "tags.h"
#include "proto.h"
#include "macros.h"
#include "version.h"
#include "gcsa.h"

#define PARAM_IMAGES (IMAGES_PER_SURFACE * SURFACES)


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  compute_area_ratios(MRI_SURFACE *mris) ;
static double gcsaSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;

static const char *surface_names[] =
{
  "inflated",
  "smoothwm",
  "smoothwm"
} ;
static const char *curvature_names[] =
{
  "inflated.H",
  "sulc",
  NULL
} ;

#define MAX_SIGMAS 10
static int nsigmas=0 ;
static float sigmas[MAX_SIGMAS] ;

#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])

static char *starting_reg_fname = NULL ;
static int multi_scale = 0 ;
static int which_norm = NORM_MEAN ;
static int navgs = 0 ;
static int single_surf = 0 ;
static double l_ocorr = 1.0 ;
static char *annot_name = NULL ;
static int atlas_size=3;
static int max_passes = 4 ;
static float min_degrees = 0.5 ;
static float max_degrees = 64.0 ;
static int nangles = 8 ;
static int nbrs = 1 ;
static float scale = 1.0f ;

static int reverse_flag = 0 ;

static float dalpha = 0.0f ;
static float dbeta = 0.0f ;
static float dgamma = 0.0f ;

#define MAX_OVERLAYS 1000
static int noverlays = 0 ;
static char *overlays[MAX_OVERLAYS]  ;
const char *Progname ;
static char curvature_fname[STRLEN] = "" ;
static const char *orig_name = "smoothwm" ;
static const char *canon_name = "sphere" ;
static char *jacobian_fname = NULL ;
static char *inflated_name = NULL ;

#define MAX_LABELS 100
static int nlabels = 0 ;
static LABEL *labels[MAX_LABELS] ;
static char  *label_names[MAX_LABELS] ;
static GCSA  *label_gcsa[MAX_LABELS] ;
static int   label_indices[MAX_LABELS] ;
static int   label_annots[MAX_LABELS] ;

/* multiframe registration */
static int multiframes = 0;
static int use_initial_registration = 0;

static void initParms(void);

static int use_defaults = 1 ;

static INTEGRATION_PARMS  parms ;
static int remove_negative = 1 ;
char *rusage_file=NULL;
char *regfile = NULL;

int
main(int argc, char *argv[])
{
  char **av, *surf_fname, *template_fname, *out_fname, fname[STRLEN],*cp;
  int ac, nargs,err, msec ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp_template ;

  char cwd[2000],*cmdline2 ;
  std::string cmdline = getAllInfo(argc, argv, "mris_register");

  nargs = handleVersionOption(argc, argv, "mris_register");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Timer start;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  cmdline2 = argv2cmdline(argc,argv);
  getcwd(cwd,2000);
  printf("\ncwd %s\n",cwd);
  printf("cmdline %s\n\n",cmdline2);

  parms.projection = PROJECT_SPHERE ;
  parms.flags |= IP_USE_CURVATURE ;
  parms.trinarize_thresh = 0.0 ;  // disabled by default
  parms.tol = 0.5 ;    // was 1e-0*2.5
  parms.min_averages = 0 ;
  parms.l_area = 0.0 ;
  parms.l_parea = 0.1f ;  // used to be 0.2
  parms.l_dist = 5.0 ; // used to be 0.5, and before that 0.1
  parms.l_corr = 1.0f ;
  parms.l_nlarea = 1 ;
  parms.l_pcorr = 0.0f ;
  parms.niterations = 25 ;
  parms.n_averages = 1024 ;   // used to be 256
  parms.first_pass_averages = 1024*16 ;  // only used in first pass
  parms.write_iterations = 100 ;
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.99 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.dt_increase = 1.0 ;
  parms.dt_decrease = 1.0 ;
  parms.l_external = 10000 ;   /* in case manual label is specified */
  parms.error_ratio = 1.1 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_ADAPTIVE ;
  parms.integration_type = INTEGRATE_MOMENTUM /*INTEGRATE_LINE_MINIMIZE*/ ;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.dt = 0.9 ;
  parms.momentum = 0.95 ;
  parms.desired_rms_height = -1.0 ;
  parms.nbhd_size = -10 ;
  parms.max_nbrs = 10 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  MRISprintCurvatureNames(stdout);

  if (nsigmas > 0)
  {
    MRISsetRegistrationSigmas(sigmas, nsigmas) ;
  }
  parms.which_norm = which_norm ;
  if (argc < 4)
  {
    usage_exit() ;
  }

  std::cout << getVersion() << std::endl;
  printf("  %s\n",getVersion().c_str());
  fflush(stdout);

  surf_fname = argv[1] ;
  template_fname = argv[2] ;
  out_fname = argv[3] ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
    {
      strcpy(parms.base_name, cp+1) ;
    }
    else
    {
      strcpy(parms.base_name, "sphere") ;
    }
  }

  fprintf(stderr, "reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if(regfile){
    printf("Reading in reg file %s\n",regfile);
    LTA *lta = LTAread(regfile);
    if(lta==NULL) exit(1);
    printf("Extracting rotational components\n");
    LTAmat2RotMat(lta);
    printf("Applying rotation matrix to surface\n");
    //MatrixPrint(stdout,lta->xforms[0].m_L);
    int err = MRISltaMultiply(mris, lta);
    if(err) exit(1);
    LTAfree(&lta);
  }

  if (parms.var_smoothness)
  {
    parms.vsmoothness = (float *)calloc(mris->nvertices, sizeof(float)) ;
    if (parms.vsmoothness == NULL)
    {
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate vsmoothness array",
                Progname) ;
    }
    parms.dist_error = (float *)calloc(mris->nvertices, sizeof(float)) ;
    if (parms.dist_error == NULL)
    {
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate dist_error array",
                Progname) ;
    }
    parms.area_error = (float *)calloc(mris->nvertices, sizeof(float)) ;
    if (parms.area_error == NULL)
    {
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate area_error array",
                Progname) ;
    }
    parms.geometry_error = (float *)calloc(mris->nvertices, sizeof(float)) ;
    if (parms.geometry_error == NULL)
    {
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate geometry_error array",
                Progname) ;
    }
  }

  MRISresetNeighborhoodSize(mris, 1) ;
  if (annot_name)
  {
    if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
      ErrorExit(ERROR_BADPARM,
                "%s: could not read annot file %s",
                Progname, annot_name) ;
    MRISripMedialWall(mris) ;
  }

  MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
  MRISaddCommandLine(mris, cmdline) ;
  if (!FZERO(dalpha) || !FZERO(dbeta) || !FZERO(dgamma))
    MRISrotate(mris, mris, RADIANS(dalpha), RADIANS(dbeta),
               RADIANS(dgamma)) ;

  if (curvature_fname[0])
  {
    fprintf(stderr, "reading source curvature from %s\n",curvature_fname) ;
    MRISreadCurvatureFile(mris, curvature_fname) ;
    if (parms.nonmax)
    {
      printf("suppressing nonmaxima in %s\n", curvature_fname);
      MRISnonmaxSuppress(mris) ;
    }
  }
  if (single_surf)
  {
    char        fname[STRLEN], *cp, surf_dir[STRLEN], hemi[10]  ;
    MRI_SURFACE *mris_template ;
    int         sno, tnbrs=3 ;

    FileNamePath(template_fname, surf_dir) ;
    cp = strrchr(template_fname, '/') ;
    if (cp == NULL) // no path - start from beginning of file name
    {
      cp = template_fname ;
    }
    cp = strchr(cp, '.') ;
    if (cp == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could no scan hemi from %s",
                Progname, template_fname) ;
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
    fprintf(stderr, "reading spherical surface %s...\n", template_fname) ;
    mris_template = MRISread(template_fname) ;
    if (mris_template == NULL)
    {
      ErrorExit(ERROR_NOFILE, "") ;
    }
#if 0
    if (reverse_flag)
    {
      MRISreverse(mris_template, REVERSE_X, 1) ;
    }
#endif
    MRISsaveVertexPositions(mris_template, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris_template) ;
    MRISstoreMetricProperties(mris_template) ;

    if (noverlays > 0)
    {
      mrisp_template = MRISPalloc(scale, IMAGES_PER_SURFACE*noverlays);
      for (sno = 0; sno < noverlays ; sno++)
      {
        int req = snprintf(fname, STRLEN, "%s/../label/%s.%s", surf_dir, hemi, overlays[sno]) ; 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        if (MRISreadValues(mris_template, fname)  != NO_ERROR)
          ErrorExit(ERROR_NOFILE,
                    "%s: could not read overlay from %s",
                    Progname, fname) ;
        MRIScopyValuesToCurvature(mris_template) ;
        MRISaverageCurvatures(mris_template, navgs) ;
        MRISnormalizeCurvature(mris_template, which_norm) ;
        fprintf(stderr,
                "computing parameterization for overlay %s...\n",
                fname);
        MRIStoParameterization(mris_template, mrisp_template, scale, sno*3) ;
        MRISPsetFrameVal(mrisp_template, sno*3+1, 1.0) ;
      }
    }
    else
    {
      mrisp_template = MRISPalloc(scale, PARAM_IMAGES);
      for (sno = 0; sno < SURFACES ; sno++)
      {
        if (curvature_names[sno])  /* read in precomputed curvature file */
        {
          int req = snprintf(fname, STRLEN, "%s/%s.%s", surf_dir, hemi, curvature_names[sno]) ; 
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          if (MRISreadCurvatureFile(mris_template, fname) != NO_ERROR)
            ErrorExit(Gerror,
                      "%s: could not read curvature file '%s'\n",
                      Progname, fname) ;

          /* the two next lines were not in the original code */
          MRISaverageCurvatures(mris_template, navgs) ;
          MRISnormalizeCurvature(mris_template, which_norm) ;
        }
        else                         /* compute curvature of surface */
        {
          int req = snprintf(fname, STRLEN, "%s/%s.%s", surf_dir, hemi, surface_names[sno]) ; 
	  if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	  }
          if (MRISreadVertexPositions(mris_template, fname) != NO_ERROR)
            ErrorExit(ERROR_NOFILE,
                      "%s: could not read surface file %s",
                      Progname, fname) ;

          if (tnbrs > 1)
          {
            MRISresetNeighborhoodSize(mris_template, tnbrs) ;
          }
          MRIScomputeMetricProperties(mris_template) ;
          MRIScomputeSecondFundamentalForm(mris_template) ;
          MRISuseMeanCurvature(mris_template) ;
          MRISaverageCurvatures(mris_template, navgs) ;
          MRISrestoreVertexPositions(mris_template, CANONICAL_VERTICES) ;
          MRISnormalizeCurvature(mris_template, which_norm) ;
        }
	if (parms.nonmax)
	{
	  printf("suppressing nonmaxima in %s\n", fname);
	  MRISnonmaxSuppress(mris) ;
	}
        fprintf(stderr,
                "computing parameterization for surface %s...\n",
                fname);
        MRIStoParameterization(mris_template, mrisp_template, scale, sno*3) ;
        MRISPsetFrameVal(mrisp_template, sno*3+1, 1.0) ;
      }
    }
  }
  else
  {
    fprintf(stderr, "reading template parameterization from %s...\n",
            template_fname) ;
    mrisp_template = MRISPread(template_fname) ;
    if (!mrisp_template)
      ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;
    if (noverlays > 0)
    {
      if (mrisp_template->Ip->num_frame != IMAGES_PER_SURFACE*noverlays)
        ErrorExit(ERROR_BADPARM,
                  "template frames (%d) doesn't match input (%d x %d) = %d\n",
                  mrisp_template->Ip->num_frame, IMAGES_PER_SURFACE,noverlays,
                  IMAGES_PER_SURFACE*noverlays) ;
    }
  }

  if (use_defaults)
  {
    if (*IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2) <= 1.0)  /* 1st time */
    {
      parms.l_dist = 5.0 ;
      parms.l_corr = 1.0 ;
      parms.l_parea = 0.2 ;
    }
    else   /* subsequent alignments */
    {
      parms.l_dist = 5.0 ;
      parms.l_corr = 1.0 ;
      parms.l_parea = 0.2 ;
    }
  }

  if (nbrs > 1)
  {
    MRISresetNeighborhoodSize(mris, nbrs) ;
  }
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  mris->status = MRIS_PARAMETERIZED_SPHERE ;
  MRIScomputeMetricProperties(mris) ;
  if (MRIScountNegativeFaces(mris) > nint(.8*mris->nfaces))
  {
    printf("!!!!!!!!!  everted surface detected - correcting !!!!!!!!!!!!!!\n") ;
    MRISevertSurface(mris) ;
  }

  if (!FZERO(parms.l_dist))
  {
    MRISscaleDistances(mris, scale) ;
  }
#if 0
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISzeroNegativeAreas(mris) ;
  MRISstoreMetricProperties(mris) ;
#endif
  MRISstoreMeanCurvature(mris) ;  /* use curvature from file */
  MRISsetOriginalFileName(orig_name) ;
  if (inflated_name)
  {
    MRISsetInflatedFileName(inflated_name) ;
  }
  err = MRISreadOriginalProperties(mris, orig_name) ;
  if (err != 0)
  {
    printf("ERROR %d from MRISreadOriginalProperties().\n",err);
    exit(1);
  }

  if (MRISreadCanonicalCoordinates(mris, canon_name) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not read canon surface %s",
              Progname, canon_name) ;

  if (reverse_flag)
  {
    MRISreverse(mris, REVERSE_X, 1) ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISreverse(mris, REVERSE_X, 0) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
  }
#if 0
  MRISsaveVertexPositions
  (mris, CANONICAL_VERTICES) ;  // uniform spherical positions
#endif
  if (starting_reg_fname)
    if (MRISreadVertexPositions(mris, starting_reg_fname) != NO_ERROR)
    {
      exit(Gerror) ;
    }

  if (multiframes)
  {
    if (use_initial_registration)
      MRISvectorRegister(mris, mrisp_template, &parms, max_passes,
                         min_degrees, max_degrees, nangles) ;
    parms.l_corr=parms.l_pcorr=0.0f;
#if 0
    parms.l_dist = 0.0 ;
    parms.l_corr = 0.0 ;
    parms.l_parea = 0.0 ;
    parms.l_area = 0.0 ;
    parms.l_parea = 0.0f ;
    parms.l_dist = 0.0 ;
    parms.l_corr = 0.0f ;
    parms.l_nlarea = 0.0f ;
    parms.l_pcorr = 0.0f ;
#endif
    MRISvectorRegister(mris,
                       mrisp_template,
                       &parms,
                       max_passes,
                       min_degrees,
                       max_degrees,
                       nangles) ;
  }
  else
  {
    double l_dist = parms.l_dist ;
    if (multi_scale > 0)
    {
      int i ;

      parms.l_dist = l_dist * pow(5.0, (multi_scale-1.0)) ;
      parms.flags |= IPFLAG_NOSCALE_TOL ;
      parms.flags &= ~IP_USE_CURVATURE ;
      for (i = 0 ; i < multi_scale ; i++)
      {
        printf("*************** round %d, l_dist = %2.3f **************\n", i,
               parms.l_dist) ;
        MRISregister(mris, mrisp_template,
                     &parms, max_passes,
                     min_degrees, max_degrees, nangles) ;
        parms.flags |= IP_NO_RIGID_ALIGN ;
        parms.flags &= ~IP_USE_INFLATED ;
        parms.l_dist /= 5 ;
      }

      if (parms.nbhd_size < 0)
      {
        parms.nbhd_size *= -1 ;
        printf("**** starting 2nd epoch, with long-range distances *****\n");
        parms.l_dist = l_dist * pow(5.0, (multi_scale-2.0)) ;
        for (i = 1 ; i < multi_scale ; i++)
        {
          printf("*********** round %d, l_dist = %2.3f *************\n", i,
                 parms.l_dist) ;
          MRISregister(mris, mrisp_template,
                       &parms, max_passes,
                       min_degrees, max_degrees, nangles) ;
          parms.l_dist /= 5 ;
        }
      }
      printf("****** final curvature registration ***************\n") ;
      if (parms.nbhd_size > 0)
      {
        parms.nbhd_size *= -1 ;  // disable long-range stuff
      }
      parms.l_dist *= 5 ;
      parms.flags |= (IP_USE_CURVATURE | IP_NO_SULC);
      MRISregister(mris, mrisp_template,
                   &parms, max_passes,
                   min_degrees, max_degrees, nangles) ;
    }
    else
      MRISregister(mris, mrisp_template,
                   &parms, max_passes,
                   min_degrees, max_degrees, nangles) ;

  }
#if 0
  parms.l_dist *= 50 ;
  MRISregister(mris, mrisp_template,
                   &parms, max_passes,
                   min_degrees, max_degrees, nangles) ;
  parms.l_dist /= 50 ;
  MRISregister(mris, mrisp_template,
                   &parms, max_passes,
                   min_degrees, max_degrees, nangles) ;
#endif

  if (remove_negative)
  {
    parms.niterations = 1000 ;
    MRISremoveOverlapWithSmoothing(mris,&parms) ;
  }
  fprintf(stderr, "writing registered surface to %s...\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  if (jacobian_fname)
  {
    MRIScomputeMetricProperties(mris) ;
    compute_area_ratios(mris) ;  /* will put results in v->curv */
#if 0
    MRISwriteArea(mris, jacobian_fname) ;
#else
    MRISwriteCurvature(mris, jacobian_fname) ;
#endif
  }

  MRISPfree(&mrisp_template) ;
  MRISfree(&mris) ;

  msec = start.milliseconds() ;
  printf("registration took %2.2f hours\n",(float)msec/(1000.0f*60.0f*60.0f));

  printf("#VMPC# mris_register VmPeak  %d\n",GetVmPeak());

  // Output formatted so it can be easily grepped
#ifdef HAVE_OPENMP
  int n_omp_threads = omp_get_max_threads();
  printf("FSRUNTIME@ mris_register %7.4f hours %d threads\n",msec/(1000.0*60.0*60.0),n_omp_threads);
#else
  printf("FSRUNTIME@ mris_register %7.4f hours %d threads\n",msec/(1000.0*60.0*60.0),1);
#endif


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
  int    n,nargs = 0 ;
  char   *option ;
  float  f ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "trinarize"))
  {
    parms.trinarize_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("binarizing curvatures with threshold = %f\n", parms.trinarize_thresh) ;
  }
  else if (!stricmp(option, "median"))
  {
    which_norm = NORM_MEDIAN ;
    printf("using median normalization\n") ;
  }
  else if (!stricmp(option, "vnum") || !stricmp(option, "distances"))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  }
  else if (!stricmp(option, "nonorm"))
  {
    which_norm = NORM_NONE ;
    printf("disabling normalization\n") ;
  }
  else if (!stricmp(option, "vsmooth"))
  {
    parms.var_smoothness = 1 ;
    printf("using space/time varying smoothness weighting\n") ;
  }
  else if (!stricmp(option, "sigma"))
  {
    if (nsigmas >= MAX_SIGMAS)
    {
      ErrorExit(ERROR_NOMEMORY, "%s: too many sigmas specified (%d)\n",
                Progname, nsigmas) ;
    }
    sigmas[nsigmas] = atof(argv[2]) ;
    nargs = 1 ;
    nsigmas++ ;
  }
  else if (!stricmp(option, "vector"))
  {
    fprintf(stderr,
            "\nMultiframe Mode:\n");
    fprintf(stderr,
            "Use -addframe option to add extra-fields into average atlas\n");
    fprintf(stderr,
            "\t-addframe which_field where_in_atlas l_corr l_pcorr\n");
    fprintf(stderr,
            "\tfield code:\n");
    for (n=0 ; n < NUMBER_OF_VECTORIAL_FIELDS ; n++)
      fprintf(stderr,
              "\t      field %d is '%s' (type = %d)\n",
              n,ReturnFieldName(n),IsDistanceField(n));
    exit(1);
  }
  else if (!stricmp(option, "annot"))
  {
    annot_name = argv[2] ;
    fprintf(stderr,"zeroing medial wall in %s\n", annot_name) ;
    nargs=1;
  }
  else if (!stricmp(option, "init"))
  {
    use_initial_registration=1;
    fprintf(stderr,"use initial registration\n");
  }
  else if (!stricmp(option, "addframe"))
  {
    int which_field,where_in_atlas;
    float l_corr,l_pcorr;

    if (multiframes==0)
    {
      /* activate multiframes mode */
      initParms();
      multiframes = 1;
    }

    which_field=atoi(argv[2]);
    where_in_atlas=atoi(argv[3]);
    l_corr=atof(argv[4]);
    l_pcorr=atof(argv[5]);

    fprintf(stderr,
            "adding field %d (%s) with location %d in the atlas\n",
            which_field,ReturnFieldName(which_field),where_in_atlas) ;
    /* check if this field exist or not */
    for (n = 0 ; n < parms.nfields ; n++)
    {
      if (parms.fields[n].field==which_field)
      {
        fprintf(stderr,"ERROR: field already exists\n");
        exit(1);
      }
    }
    /* adding field into parms */
    n=parms.nfields++;
    SetFieldLabel(&parms.fields[n],
                  which_field,
                  where_in_atlas,
                  l_corr,
                  l_pcorr,
                  0,
                  which_norm);
    nargs = 4 ;
  }
  /* else if (!stricmp(option, "hippocampus")) */
  /*   { */
  /*   if(multiframes==0){ */
  /*    multiframes = 1 ; */
  /*    initParms(); */
  /*   } */
  /*   where=atoi(argv[2]); */
  /*   parms.l_corrs[]=atof(argv[3]); */
  /*   parms.l_pcorrs[]=atof(argv[4]); */
  /*   parms.frames[]=; */
  /*   parms. */
  /*     fprintf(stderr, "using hippocampus distance map\n") ; */
  /*   nargs=3; */
  /*   } */
  else if (!stricmp(option, "topology"))
  {
    parms.flags |= IPFLAG_PRESERVE_SPHERICAL_POSITIVE_AREA;
    fprintf(stderr, "preserving the topology of positive area triangles\n");
  }
  else if (!stricmp(option, "vnum") || !stricmp(option, "distances"))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  }
  else if (!stricmp(option, "rotate"))
  {
    dalpha = atof(argv[2]) ;
    dbeta = atof(argv[3]) ;
    dgamma = atof(argv[4]) ;
    fprintf(stderr, "rotating brain by (%2.2f, %2.2f, %2.2f)\n",
            dalpha, dbeta, dgamma) ;
    nargs = 3 ;
  }
  else if (!stricmp(option, "reg"))
  {
    regfile = argv[2];
    nargs = 1 ;
  }
  else if (!stricmp(option, "reverse"))
  {
    reverse_flag = 1 ;
    fprintf(stderr, "mirror image reversing brain before morphing...\n") ;
  }
  else if (!stricmp(option, "min_degrees"))
  {
    min_degrees = atof(argv[2]) ;
    fprintf(stderr,
            "setting min angle for search to %2.2f degrees\n",
            min_degrees) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "max_degrees"))
  {
    max_degrees = atof(argv[2]) ;
    fprintf(stderr,
            "setting max angle for search to %2.2f degrees\n",
            max_degrees) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nangles"))
  {
    nangles = atoi(argv[2]) ;
    fprintf(stderr, "setting # of angles/search per scale to %d\n", nangles) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "jacobian"))
  {
    jacobian_fname = argv[2] ;
    nargs = 1 ;
    printf("writing out jacobian of mapping to %s\n", jacobian_fname) ;
  }
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    use_defaults = 0 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  }
  else if (!stricmp(option, "norot"))
  {
    fprintf(stderr, "disabling initial rigid alignment...\n") ;
    parms.flags |= IP_NO_RIGID_ALIGN ;
  }
  else if (!stricmp(option, "inflated"))
  {
    fprintf(stderr, "using inflated surface for initial alignment\n") ;
    parms.flags |= IP_USE_INFLATED ;
  }
  else if (!stricmp(option, "noinflated"))
  {
    fprintf(stderr, "using inflated surface for initial alignment\n") ;
    parms.flags &= ~IP_USE_INFLATED ;
  }
  else if (!stricmp(option, "multi_scale"))
  {
    multi_scale = atoi(argv[2]) ;
    fprintf(stderr, "using %d scales for morphing\n", multi_scale) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nsurfaces"))
  {
    parms.nsurfaces = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,
            "using %d surfaces/curvatures for alignment\n",
            parms.nsurfaces) ;
  }
  else if (!stricmp(option, "infname"))
  {
    char fname[STRLEN] ;
    inflated_name = argv[2] ;
    surface_names[0] = argv[2] ;
    nargs = 1 ;
    printf("using %s as inflated surface name, "
           "and using it for initial alignment\n", inflated_name) ;
    sprintf(fname, "%s.H", argv[2]) ;
    curvature_names[0] = (char *)calloc(strlen(fname)+1, sizeof(char)) ;
    strcpy(const_cast<char*>(curvature_names[0]), fname) ; // strcpy _and_ const_cast in a single line....
    parms.flags |= IP_USE_INFLATED ;
  }
  else if (!stricmp(option, "nosulc"))
  {
    fprintf(stderr, "disabling initial sulc alignment...\n") ;
    parms.flags |= IP_NO_SULC ;
  }
  else if (!stricmp(option, "sulc"))
  {
    curvature_names[1] = argv[2] ;
    fprintf(stderr, "using %s to replace 'sulc' alignment\n",
            curvature_names[1]) ;
    nargs = 1 ;
    MRISsetSulcFileName(argv[2]);
  }

  else if (!stricmp(option, "surf0"))
  {
    surface_names[0]       = argv[2];
    fprintf(stderr, "using %s as input surface 0.\n",
            surface_names[0]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "surf1"))
  {
    surface_names[1]       = argv[2];
    fprintf(stderr, "using %s as input surface 1.\n",
            surface_names[1]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "surf2"))
  {
    surface_names[2]       = argv[2];
    fprintf(stderr, "using %s as input surface 2.\n",
            surface_names[2]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "curv0"))
  {
    curvature_names[0]  = argv[2];
    MRISsetCurvatureName(0, curvature_names[0]);
    fprintf(stderr, "using %s as curvature function for surface 0.\n",
            curvature_names[0]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "curv1"))
  {
    curvature_names[1]  = argv[2];
    MRISsetCurvatureName(1, curvature_names[1]);
    fprintf(stderr, "using %s as curvature function for surface 1.\n",
            curvature_names[1]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "curv2"))
  {
    curvature_names[2]  = argv[2];
    MRISsetCurvatureName(2, curvature_names[2]);
    fprintf(stderr, "using %s as curvature function for surface 2.\n",
            curvature_names[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "threads")){
    int nthreads;
    sscanf(argv[2],"%d",&nthreads);
    #ifdef _OPENMP
    omp_set_num_threads(nthreads);
    #endif
    nargs = 1;
  } 
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "search"))
  {
    parms.integration_type = INTEGRATE_LM_SEARCH ;
    fprintf(stderr, "integrating with binary search line minimization\n") ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = .2*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
  }
  else if (!stricmp(option, "area"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "parea"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_parea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_parea = %2.3f\n", parms.l_parea) ;
  }
  else if (!stricmp(option, "nlarea"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_nlarea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_nlarea = %2.3f\n", parms.l_nlarea) ;
  }
  else if (!stricmp(option, "spring"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "corr"))
  {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_corr) ;
    nargs = 1 ;
    fprintf(stderr, "using l_corr = %2.3f\n", parms.l_corr) ;
  }
  else if (!stricmp(option, "remove_negative"))
  {
    remove_negative = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "%sremoving negative triangles with iterative smoothing\n",
            remove_negative ? "" : "not ") ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.flags |= IP_USE_CURVATURE ;
    fprintf(stderr, "using smoothwm curvature for final alignment\n") ;
  }
  else if (!stricmp(option, "nocurv"))
  {
    parms.flags &= ~IP_USE_CURVATURE ;
    fprintf(stderr, "NOT using smoothwm curvature for final alignment\n") ;
  }
  else if (!stricmp(option, "sreg"))
  {
    starting_reg_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "starting registration with coordinates in  %s\n",
            starting_reg_fname) ;
  }
  else if (!stricmp(option, "adaptive"))
  {
    parms.integration_type = INTEGRATE_ADAPTIVE ;
    fprintf(stderr, "using adaptive time step integration\n") ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "tol"))
  {
    if (sscanf(argv[2], "%e", &f) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan tol from %s",
                Progname, argv[2]) ;
    parms.tol = (double)f ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.2e\n", (float)parms.tol) ;
  }
  else if (!stricmp(option, "error_ratio"))
  {
    parms.error_ratio = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "error_ratio=%2.3f\n", parms.error_ratio) ;
  }
  else if (!stricmp(option, "dt_inc"))
  {
    parms.dt_increase = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_increase=%2.3f\n", parms.dt_increase) ;
  }
  else if (!stricmp(option, "lap") || !stricmp(option, "lap"))
  {
    parms.l_lap = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_laplacian = %2.3f\n", parms.l_lap) ;
  }
  else if (!stricmp(option, "vnum"))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  }
  else if (!stricmp(option, "dt_dec"))
  {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else if (!stricmp(option, "ocorr"))
  {
    l_ocorr = atof(argv[2]) ;
    printf("setting overlay correlation coefficient to %2.1f\n", l_ocorr) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else if (!stricmp(option, "rusage"))
  {
    // resource usage
    rusage_file = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "overlay"))
  {
    int navgs ;
    if (noverlays == 0)
    {
      atlas_size = 0 ;
    }
    if (multiframes == 0)
    {
      initParms() ;
      multiframes = 1 ;
    }
    overlays[noverlays++] = argv[2] ;
    navgs = atof(argv[3]) ;
    printf("reading overlay from %s and smoothing it %d times\n",
           argv[2], navgs) ;
    n=parms.nfields++;
    SetFieldLabel(&parms.fields[n],
                  OVERLAY_FRAME,
                  atlas_size,
                  l_ocorr,
                  0.0,
                  navgs,
                  which_norm);
    SetFieldName(&parms.fields[n], argv[2]) ;
    atlas_size++ ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "distance"))
  {
    int navgs ;
    if (noverlays == 0)
    {
      atlas_size = 0 ;
    }
    if (multiframes == 0)
    {
      initParms() ;
      multiframes = 1 ;
    }
    overlays[noverlays++] = argv[2] ;
    navgs = atof(argv[3]) ;
    printf("reading overlay from %s and smoothing it %d times\n",
           argv[2], navgs) ;
    n=parms.nfields++;
    SetFieldLabel(&parms.fields[n],
                  DISTANCE_TRANSFORM_FRAME,
                  atlas_size,
                  l_ocorr,
                  0.0,
                  navgs,
                  NORM_MAX) ;
    SetFieldName(&parms.fields[n], argv[2]) ;
    atlas_size++ ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "canon"))
  {
    canon_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using %s for canonical properties...\n", canon_name) ;
  }
  else if (!stricmp(option, "nonmax"))
  {
    parms.nonmax = 1 ;
    fprintf(stderr, "applying nonmax supporession prior to registration\n") ;
  }
  else if (!stricmp(option, "overlay-dir"))
  {
    parms.overlay_dir = strcpyalloc(argv[2]) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
      break ;
    case 'L':
      if (nlabels >= MAX_LABELS-1)
        ErrorExit(ERROR_NO_MEMORY,
                  "%s: too many labels specified (%d max)",
                  Progname, MAX_LABELS) ;
      nargs = 3 ;
      labels[nlabels] = LabelRead(NULL, argv[2]) ;
      if (labels[nlabels] == NULL)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read label file %s",
                  Progname, argv[2]) ;
      label_gcsa[nlabels] = GCSAread(argv[3]) ;
      if (label_gcsa[nlabels] == NULL)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read GCSA file %s",
                  Progname, argv[3]) ;
      label_names[nlabels] = argv[4] ;
      CTABfindName(label_gcsa[nlabels]->ct, argv[4],&label_indices[nlabels]) ;

      if (label_indices[nlabels] < 0)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not map name %s to index",
                  Progname, argv[3]) ;
      CTABannotationAtIndex(label_gcsa[nlabels]->ct, label_indices[nlabels],
                            &label_annots[nlabels]);
      nlabels++ ;
      gMRISexternalSSE = gcsaSSE ;
      break ;
    case 'E':
      parms.l_external = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting l_external = %2.1f\n", parms.l_external) ;
      break ;
    case 'C':
      strcpy(curvature_fname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'A':
      sscanf(argv[2], "%d", &parms.n_averages) ;
      nargs = 1 ;
      fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
      break ;
    case '1':
      single_surf = True ;
      printf("treating target as a single subject's surface...\n") ;
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      fprintf(stderr, "scaling distances by %2.2f\n", scale) ;
      nargs = 1 ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
      break ;
    case 'W':
      Gdiag |= DIAG_SHOW ;
      Gdiag |= DIAG_WRITE ;
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr,
              "using write iterations = %d\n",
              parms.write_iterations) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'O':
      orig_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "using %s for original properties...\n", orig_name) ;
      break ;
    case 'P':
      max_passes = atoi(argv[2]) ;
      fprintf(stderr, "limiting unfolding to %d passes\n", max_passes) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_help() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_help() ;
  exit(1) ;
}

#include "mris_register.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_register_help_xml,mris_register_help_xml_len);
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

static int
compute_area_ratios(MRI_SURFACE *mris)
{
  VERTEX  *v ;
  int     vno ;
  float   area_scale ;

  area_scale = mris->total_area / mris->orig_area  ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }

    v->curv = v->area / (v->origarea*area_scale) ;
  }

  return(NO_ERROR) ;
}

static double
gcsaSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int       vno, ano, lno, vno_prior, n, found ;
  VERTEX    *v, *v_prior ;
  double    sse ;
  LABEL     *area ;
  CP_NODE   *cpn ;
  CP        *cp ;
  GCSA      *gcsa ;

  for (sse = 0.0, ano = 0 ; ano < nlabels ; ano++)
  {
    area = labels[ano] ;
    gcsa = label_gcsa[ano] ;
    for (lno = 0 ; lno < area->n_points ; lno++)
    {
      vno = area->lv[lno].vno ;
      if (vno < 0)
      {
        continue ;
      }
      v = &mris->vertices[vno] ;
      found = 0 ;
      v_prior = GCSAsourceToPriorVertex(gcsa, v) ;
      vno_prior = v_prior - gcsa->mris_priors->vertices ;
      cpn = &gcsa->cp_nodes[vno_prior] ;
      for (n = 0 ; n < cpn->nlabels ; n++)
      {
        cp = &cpn->cps[n] ;
        if (cpn->labels[n] == label_annots[ano])
        {
          found = 1 ;
          break ;
        }
      }
      if (found == 0)
      {
        sse += parms->l_external ;
      }
    }
  }

  //printf("external SSE: %2.3f\n", sse) ;
  return(sse) ;
}

void initParms(void)
{
  int n;
  parms.flags |= IP_USE_MULTIFRAMES;
  parms.nfields=0;
  for (n = 0 ; n < MNOFIV ; n++)
  {
    InitFieldLabel(&parms.fields[n]);
  }
}
