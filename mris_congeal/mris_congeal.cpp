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

#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "tags.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "gcsa.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
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
static int navgs = 0 ;
static float sigmas[MAX_SIGMAS] ;
static int Gsno = 0 ;
static int Giter = 0 ;

static double max_angle = 16.0 ;
#define IMAGE_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

static char *starting_reg_fname = NULL ;
static int multi_scale = 0 ;
static int which_norm = NORM_MEAN ;
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

static int ico_order = 7 ;
static char sdir[STRLEN] = "" ;
#define MAX_SUBJECTS 200

MRI_SP *MRIScongeal(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, int nsubjects, 
                    MRI_SP *mrisp, INTEGRATION_PARMS *parms) ;
int
main(int argc, char *argv[])
{
  char         **av, *surf_name, *out_fname, fname[STRLEN],*cp, *hemi, *fs_home;
  int          ac, nargs, msec, nsubjects, i ;
  MRI_SURFACE  *mris_array[MAX_SUBJECTS], *mris_ico ;
  MRI_SP       *mrisp_template ;


  std::string cmdline = getAllInfo(argc, argv, "mris_congeal");

  nargs = handleVersionOption(argc, argv, "mris_congeal");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Timer start;
  Gdiag = DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = PROJECT_SPHERE ;
  parms.flags |= IP_USE_CURVATURE ;
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

  fs_home = getenv("FREESURFER_HOME") ;
  if (!fs_home)
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment.\n", Progname) ;
  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  nsubjects = argc-4 ;
  out_fname = argv[argc-1] ;
  printf("processing %d subjects and writing output to %s\n", nsubjects, out_fname) ;

  if (nsigmas > 0)
    MRISsetRegistrationSigmas(sigmas, nsigmas) ;
  parms.which_norm = which_norm ;
  if (argc < 4) usage_exit() ;

  std::cout << getVersion() << std::endl;
  printf("  %s\n",getVersion().c_str());
  fflush(stdout);

  hemi = argv[1] ;
  surf_name = argv[2] ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "sphere") ;
  }

  for (i = 0 ; i < nsubjects ; i++)
  {
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, argv[i+3], hemi, surf_name) ;
    fprintf(stderr, "reading surface from %s...\n", fname) ;
    mris_array[i] = MRISread(fname) ;
    if (!mris_array[i])
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISaddCommandLine(mris_array[i], cmdline) ;
  }

  sprintf(fname, "%s/lib/bem/ic%d.tri", fs_home, ico_order) ;
  mris_ico = NULL;
  mris_ico = MRISread(fname) ;
  if (mris_ico == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read icosahedron from %s", fname) ;
  }
  mrisp_template = MRIScongeal(mris_ico, mris_array, nsubjects, NULL, &parms) ;

  msec = start.milliseconds() ;
  if (Gdiag & DIAG_SHOW)
    printf("registration took %2.2f hours\n",
           (float)msec/(1000.0f*60.0f*60.0f));
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
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
  else if (!stricmp(option, "max_angle"))
  {
    max_angle = atof(argv[2]) ;
    printf("setting max search angle to %2.2f degrees\n", max_angle) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "vsmooth"))
  {
    parms.var_smoothness = 1 ;
    printf("using space/time varying smoothness weighting\n") ;
  }
  else if (!stricmp(option, "sigma"))
  {
    if (nsigmas >= MAX_SIGMAS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many sigmas specified (%d)\n", Progname, nsigmas) ;
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
    { /* activate multiframes mode */
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
    strcpy(const_cast<char*>(curvature_names[0]), fname) ; // const_cast AND strcpy
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
    fprintf(stderr, "starting registration with coordinates in  %s\n", starting_reg_fname) ;
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
  else if (!stricmp(option, "overlay"))
  {
    int navgs ;
    if (noverlays == 0)
      atlas_size = 0 ;
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
      atlas_size = 0 ;
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
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "canon"))
  {
    canon_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using %s for canonical properties...\n", canon_name) ;
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

static void
print_usage(void)
{
  printf("\nUSAGE:\n"
         "%s [options] <input surface name> <hemi> <subject 1>  <subject 2> ... <output surface name>\n",
         Progname) ;
  printf("\n") ;
  printf("Options:\n\n") ;
  printf("  -SDIR SUBJECTS_DIR \n");
  printf("  -norot            : disable initial rigid alignment\n");
  printf("  -nosulc           : disable initial sulc alignment\n");
  printf("  -curv             : use smoothwm curvature for final alignment\n");
  printf("  -jacobian <fname> : write-out jacobian overlay data to fname\n");
  printf("  -dist <num>       : specify distance term\n");
  printf("  -l <label file> <atlas (*.gcs)> <label name>\n"
         "                    : this option will specify a manual label\n"
         "                      to align with atlas label <label name>\n");
  printf("  -addframe <which_field> <where_in_atlas> <l_corr> <l_pcorr>\n");
  printf("  -overlay <surfvals> <navgs> : subject/labels/hemi.surfvals\n");
  printf("  -overlay-dir <dir> : subject/dir/hemi.surfvals\n");
  printf("  -1                 : target specifies a subject's surface,\n"
         "                       not a template file\n") ;
}

static void
print_help(void)
{
  print_usage() ;
  printf("\nThis program registers a set of input surfaces together, and generates an atlas.\n\n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
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
        continue ;
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
        sse += parms->l_external ;
    }
  }

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
#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
static int nsurfaces = 3 ;  // inflated.h, sulc, curv

static int   mrisLogStatus(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                           FILE *fp, float dt) ;
static int   mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                               int t) ;
double MRIScongealComputeSurfaceError(MRI_SURFACE *mris_ico, MRI_SURFACE *mris, 
                                      MHT *mht, INTEGRATION_PARMS *parms) ;
int MRIScongealRigidBodyAlignGlobal(MRI_SURFACE *mris_ico, MRI_SURFACE *mris, 
                                    MHT *mht, INTEGRATION_PARMS *parms,
                                    float min_degrees, float max_degrees, int nangles) ;
double
MRIScongealComputeSSE(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, MHT **mht_array, int nsubjects, 
                      INTEGRATION_PARMS *parms)
{
  int     vno, i ;
  VERTEX  *v, *vsurf ;
  double  mean, var, sse ;

  // for each vertex in the ico, find matching vertices and compute mean of them
  for (sse = 0.0, vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    for  (mean = var = 0.0, i = 0 ; i < nsubjects ; i++)
    {
      vsurf = MHTfindClosestVertexInTable(mht_array[i], mris_array[i], v->x, v->y, v->z, 1) ;
      if (vsurf == NULL)
        continue ;  // should never happen 
      mean += vsurf->curv ;
      var += (vsurf->curv * vsurf->curv) ;
    }
    mean /= nsubjects ;
    v->curv = mean ;
    v->val = var / nsubjects - mean*mean;
    if (vno == Gdiag_no)
      printf("ico v %d: %2.2f +- %2.2f\n", vno, v->curv, v->val) ;
    sse += v->val ;

  }

  return(sse) ;
}
double
MRIScongealEstimateTemplate(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, MHT **mht_array, 
                            int nsubjects, INTEGRATION_PARMS *parms, int exclude_i)
{
  int     vno, i ;
  VERTEX  *v, *vsurf ;
  double  mean, var, sse ;

  // for each vertex in the ico, find matching vertices and compute mean of them
  for (sse = 0.0, vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    for  (mean = var = 0.0, i = 0 ; i < nsubjects ; i++)
    {
      if (i == exclude_i)
        continue ;
      vsurf = MHTfindClosestVertexInTable(mht_array[i], mris_array[i], v->x, v->y, v->z, 1) ;
      if (vsurf == NULL)
        continue ;  // should never happen 
      mean += vsurf->curv ;
      var += (vsurf->curv * vsurf->curv) ;
    }
    if (exclude_i >= 0 && exclude_i < nsubjects)
    {
      mean /= (nsubjects-1) ;
      v->val = var / (nsubjects-1) - mean*mean;
    }
    else
    {
      mean /= nsubjects ;
      v->val = var / nsubjects - mean*mean;
    }
    v->curv = mean ;
    sse += v->val ;
    if (vno == Gdiag_no)
      printf("ico v %d: %2.2f +- %2.2f\n", vno, v->curv, v->val) ;

  }

  return(sse) ;
}

double
MRIScongealUpdateRegistration(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, 
                              MHT **mht_array, int nsubjects, 
                              INTEGRATION_PARMS *parms)
{
  double sse ;

  sse = MRIScongealComputeSSE(mris_ico, mris_array, mht_array, nsubjects, parms) ;
  return(sse) ;
}

double
MRIScongealUpdateRigidRegistration(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, 
                                   MHT **mht_array, MHT *mht_ico, int nsubjects, double degrees,
                                   INTEGRATION_PARMS *parms)
{
  double sse ;
  int    i ;
  static int start_t = 0 ;

  MRIScopyCurvatureToImagValues(mris_ico) ;
  for (i = 0 ; i < nsubjects ; i++)
  {
    MRIScongealEstimateTemplate(mris_ico, mris_array, mht_array, nsubjects, parms, i) ;
    parms->start_t = start_t ;
    MRIScongealRigidBodyAlignGlobal(mris_ico, mris_array[i], 
                                    mht_ico, parms,
                                    degrees, degrees, 6) ;
    MHTfree(&mht_array[i]) ;
    mht_array[i] = MHTcreateVertexTable(mris_array[i], CURRENT_VERTICES) ;
  }
  MRIScopyImagValuesToCurvature(mris_ico) ;

  start_t++ ;

  sse = MRIScongealComputeSSE(mris_ico, mris_array, mht_array, nsubjects, parms) ;
  return(sse) ;
}


MRI_SP *
MRIScongeal(MRI_SURFACE *mris_ico, MRI_SURFACE **mris_array, int nsubjects, 
            MRI_SP *mrisp, INTEGRATION_PARMS *parms)
{
  int done = 0, i, sno, iter ;
  double sse, last_sse = -1, pct_change, angle ;
  MHT  *mht_array[MAX_SUBJECTS], *mht_ico ;



  if (Gdiag & DIAG_WRITE)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.log", parms->base_name) ;
    INTEGRATION_PARMS_openFp(parms, fname, "a") ;
  }

  MRISscaleBrain(mris_ico, mris_ico, mris_array[0]->radius / mris_ico->radius) ;
  mris_ico->hemisphere = mris_array[0]->hemisphere ;

  mht_ico = MHTcreateVertexTable(mris_ico, CURRENT_VERTICES) ;
  if (mrisp == NULL)
  {
    mrisp = MRISPalloc(scale, nsurfaces * IMAGES_PER_SURFACE );
  }
  
  printf("creating hash tables....\n") ;
  for (i = 0 ; i < nsubjects ; i++)
    mht_array[i] = MHTcreateVertexTable(mris_array[i], CURRENT_VERTICES) ;

  for (sno = parms->flags & IP_USE_INFLATED ?  0 : 1 ; sno < nsurfaces ; sno++)
  {
    Gsno = sno ;  // diagnostics
    Giter = iter = 0 ;
    if (curvature_names[sno])
    {
      fprintf(stderr, "reading source curvature from %s\n",curvature_names[sno]) ;
      for (i = 0 ; i < nsubjects ; i++)
        MRISreadCurvatureFile(mris_array[i], curvature_names[sno]) ;
    }
    else
    {
      for (i = 0 ; i < nsubjects ; i++)
      {
        MRISsaveVertexPositions(mris_array[i], TMP_VERTICES) ;
        if (MRISreadVertexPositions(mris_array[i], surface_names[sno]) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read surface %s for %s",
                    Progname, surface_names[sno], mris_array[i]->fname) ;
        MRIScomputeSecondFundamentalForm(mris_array[i]) ;
        MRISuseMeanCurvature(mris_array[i]) ;
        MRISrestoreVertexPositions(mris_array[i], TMP_VERTICES) ;
      }
    }
    for (i = 0 ; i < nsubjects ; i++)
    {
      MRISaverageCurvatures(mris_array[i], navgs) ;
      MRISnormalizeCurvature(mris_array[i], which_norm) ;
    }
      
    last_sse = MRIScongealEstimateTemplate(mris_ico, mris_array, mht_array, nsubjects, parms,-1) ;
    for (angle = max_angle ; angle >= 0.5 ; angle /= 2)
    {
      if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s.template.sno%d.iter%d", parms->base_name, sno, iter) ;
        MRISwriteCurvature(mris_ico, fname) ;
      }
      MRIScongealUpdateRigidRegistration(mris_ico, mris_array, mht_array, mht_ico, nsubjects, angle, parms) ;
      sse = MRIScongealEstimateTemplate(mris_ico, mris_array, mht_array, nsubjects, parms,-1) ;
      pct_change = 100 * (last_sse - sse) / last_sse ;
      done = (last_sse >= 0) && (pct_change < parms->tol) ;
      iter++ ; 
      printf("iter %d: sse = %2.2f, last_sse = %2.2f, (%2.2f%%)\n",iter, sse, last_sse, pct_change) ;
      last_sse = sse ;
    }
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s.template.sno%d.iter%d", parms->base_name, sno, iter) ;
      MRISwriteCurvature(mris_ico, fname) ;
    }
    do
    {
      MRIScongealUpdateRegistration(mris_ico, mris_array, mht_array, nsubjects, parms) ;
      sse = MRIScongealEstimateTemplate(mris_ico, mris_array, mht_array, nsubjects, parms,-1) ;
      pct_change = 100 * (last_sse - sse) / last_sse ;
      done = (last_sse >= 0) && (pct_change < parms->tol) ;
      for (i = 0 ; i < nsubjects ; i++) {
        MHTfree(& mht_array[i]);
        MHTcreateVertexTable(mris_array[i], CURRENT_VERTICES) ;
      }
      iter++ ;
    } while (!done) ;
  } 
  for (i = 0 ; i < nsubjects ; i++)
    MHTfree(&mht_array[i]) ;
  MHTfree(&mht_ico) ;
  if (Gdiag & DIAG_WRITE)
  {
    INTEGRATION_PARMS_closeFp(parms);
  }
  return(mrisp) ;
}


#define STARTING_ANGLE   RADIANS(16.0f)
#define ENDING_ANGLE     RADIANS(4.0f)
#define NANGLES          8

double
MRIScongealComputeSurfaceError(MRI_SURFACE *mris_ico, MRI_SURFACE *mris, MHT *mht_ico, INTEGRATION_PARMS *parms) 
{
  double sse, error  ;
  int    vno ;
  VERTEX *v, *vsurf ;

  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    vsurf = MHTfindClosestVertexInTable(mht_ico, mris_ico, v->x, v->y, v->z, 1) ;
    if (vsurf == NULL)
      continue ;  // should never happen 
    if (vsurf-mris_ico->vertices == Gdiag_no)
      DiagBreak() ;
    error = (v->curv-vsurf->curv) ;
    sse += error*error ;
  }
  return(sse) ;
}

int
MRIScongealRigidBodyAlignGlobal(MRI_SURFACE *mris_ico, MRI_SURFACE *mris, 
                                MHT *mht, INTEGRATION_PARMS *parms,
                                float min_degrees, float max_degrees, int nangles)
{
  double   alpha, beta, gamma, degrees, delta, mina, minb, ming,
  sse, min_sse, ext_sse ;
  auto const old_status = mris->status;
  auto const old_norm   = parms->abs_norm ;
  
  parms->abs_norm = 1 ;
  min_degrees = RADIANS(min_degrees) ;
  max_degrees = RADIANS(max_degrees) ;
#if 0
  mrisOrientSurface(mris) ;
  mrisOrientSurface(mris_ico) ;
#endif
  mris->status = MRIS_RIGID_BODY ;
  if (!parms->start_t)
  {
    mrisLogStatus(mris, parms, stdout, 0.0f) ;
    if (Gdiag & DIAG_WRITE)
    {
      mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
      if (parms->write_iterations > 0)
        mrisWriteSnapshot(mris, parms, 0) ;
    }
  }
  for (degrees = max_degrees ; degrees >= min_degrees ; degrees /= 2.0f)
  {
    mina = minb = ming = 0.0 ;
    min_sse = MRIScongealComputeSurfaceError(mris_ico, mris, mht, parms) ; 
    if (gMRISexternalSSE)
    {
      ext_sse = (*gMRISexternalSSE)(mris, parms) ;
      min_sse += ext_sse ;
    }
    delta = 2*degrees / (float)nangles ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "scanning %2.2f degree nbhd, min sse = %2.2f\n",
              (float)DEGREES(degrees), (float)min_sse) ;
    for (alpha = -degrees ; alpha <= degrees ; alpha += delta)
    {
      for (beta = -degrees ; beta <= degrees ; beta += delta)
      {
        if (Gdiag & DIAG_SHOW)
        {
          fprintf(stdout, "\r(%+2.2f, %+2.2f, %+2.2f), "
                  "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                  (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                  DEGREES(-degrees), (float)DEGREES(mina),
                  (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);
          fflush(stdout) ;
        }

        for (gamma = -degrees ; gamma <= degrees ; gamma += delta)
        {
          MRISsaveVertexPositions(mris, TMP_VERTICES) ;
          MRISrotate(mris, mris, alpha, beta, gamma) ;
          sse = MRIScongealComputeSurfaceError(mris_ico, mris, mht, parms) ;
          if (gMRISexternalSSE)
          {
            ext_sse = (*gMRISexternalSSE)(mris, parms) ;
            sse += ext_sse ;
          }
          MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
          if (sse < min_sse)
          {
            mina = alpha ;
            minb = beta ;
            ming = gamma ;
            min_sse = sse ;
          }
#if 0
          if (Gdiag & DIAG_SHOW)
            fprintf(stdout, "\r(%+2.2f, %+2.2f, %+2.2f), "
                    "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                    (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                    DEGREES(gamma), (float)DEGREES(mina),
                    (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);
#endif
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "\n") ;
    if (!FZERO(mina) || !FZERO(minb) || !FZERO(ming))
    {
      MRISrotate(mris, mris, mina, minb, ming) ;
      sse = MRIScongealComputeSurfaceError(mris_ico, mris, mht, parms) ; 
      if (gMRISexternalSSE)
        sse += (*gMRISexternalSSE)(mris, parms) ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "min sse = %2.2f at (%2.2f, %2.2f, %2.2f)\n",
                sse, (float)DEGREES(mina), (float)DEGREES(minb),
                (float)DEGREES(ming)) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "rotating brain by (%2.2f, %2.2f, %2.2f), sse: %2.2f\n",
                (float)DEGREES(mina), (float)DEGREES(minb),
                (float)DEGREES(ming), (float)sse) ;
      if (Gdiag & DIAG_WRITE)
        mrisLogStatus(mris, parms, parms->fp, 0.0f) ;
      if (Gdiag & DIAG_SHOW)
        mrisLogStatus(mris, parms, stdout, 0.0f) ;
    }
    parms->start_t += 1.0f ;
    parms->t += 1.0f ;
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0)
      mrisWriteSnapshot(mris, parms, parms->start_t) ;
  }

  mris->status = old_status ;
  parms->abs_norm = old_norm ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisLogStatus(MRI_SURFACE *mris,INTEGRATION_PARMS *parms,FILE *fp, float dt)
{
#if 0
  float  area_rms, angle_rms, curv_rms, sse, dist_rms, corr_rms ;
  int    n,negative ;
  float nv;
  int fyi;

  if (!(Gdiag & DIAG_SHOW))
    return(NO_ERROR) ;

  negative = MRIScountNegativeTriangles(mris) ;

  fyi=0;
  if (parms->flags & IP_USE_MULTIFRAMES)
  {
    if (FZERO(parms->l_corr))
    { /* just for your information */
      /* check if we can load curvature information */
      for (n=0;n<parms->nfields;n++)
        if (parms->fields[n].field==CURVATURE_CORR_FRAME)
        {
          parms->frame_no=parms->fields[n].frame*IMAGES_PER_SURFACE;
          fyi=1;
          parms->l_corr=1.0f;
          break;
        }
    }
  }

  sse  = mrisComputeError(mris, parms,&area_rms,&angle_rms,&curv_rms,&dist_rms,
                          &corr_rms);

  if (fyi)
    parms->l_corr=0.0f;

#if  0
  sse = MRIScomputeSSE(mris, parms) ;
#endif
#if 0
  sse /= (float)MRISvalidVertices(mris) ;
  sse = sqrt(sse) ;
#endif

  if (mris->status==MRIS_SPHERICAL_PATCH) return NO_ERROR;

  if (FZERO(parms->l_corr) &&
      FZERO(parms->l_pcorr) &&
      ((parms->flags & IP_USE_MULTIFRAMES) == 0 ))
  {
    fprintf(fp, "%3.3d: dt: %2.2f, sse: %2.1f (%2.3f, %2.1f, %2.3f), "
            "neg: %d (%%%2.3f:%%%2.2f), avgs: %d\n",
            parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
            negative, 100.0*mris->neg_area/(mris->neg_area+mris->total_area),
            100.0*mris->neg_orig_area/(mris->orig_area),
            parms->n_averages);
    if (dist_rms > 20)
      DiagBreak() ;
  }
  else
  {
    if (parms->flags & IP_USE_MULTIFRAMES)
    {
      nv = (float) MRISvalidVertices(mris) ;

      fprintf(fp, "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
              "neg: %d (%%%2.2f:%%%2.2f), avgs: %d\n",
              parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
              corr_rms, negative,
              100.0*mris->neg_area/(mris->neg_area+mris->total_area),
              100.0*mris->neg_orig_area/(mris->orig_area),
              parms->n_averages);
      for ( n = 0 ; n < parms->nfields ; n++ )
      {
        if (FZERO(parms->fields[n].l_corr+parms->fields[n].l_pcorr)) continue;
        fprintf(stdout,"  (%d: %2.3f : %2.3f)",
                n,parms->fields[n].sse,sqrt(parms->fields[n].sse/nv));
      }
      fprintf(stdout,"\n");
    }
    else
      fprintf(fp, "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
              "neg: %d (%%%2.2f:%%%2.2f), avgs: %d\n",
              parms->t, dt, sse, area_rms, (float)DEGREES(angle_rms), dist_rms,
              corr_rms, negative,
              100.0*mris->neg_area/(mris->neg_area+mris->total_area),
              100.0*mris->neg_orig_area/(mris->orig_area),
              parms->n_averages);
  }
  fflush(fp) ;
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char fname[STRLEN], path[STRLEN], base_name[STRLEN], *cp ;

  FileNamePath(mris->fname, path) ;
  sprintf(base_name, "%s/%s.%s.sno%d.iter%d", path,
          mris->hemisphere == LEFT_HEMISPHERE ? "lh":"rh", parms->base_name, Gsno, Giter);
  if ((cp = strstr(base_name, ".geo")) != NULL)
  {
    *cp = 0;
    sprintf(fname, "%s%4.4d.geo", base_name, t) ;
    *cp = '.' ;
  }
  else
    sprintf(fname, "%s.%4.4d", base_name, t) ;
#if 1
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "writing %s...", fname) ;
  if (mris->patch)
    MRISwritePatch(mris, fname) ;
  else
    MRISwrite(mris, fname) ;

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stdout, "done.\n") ;
    fflush(stderr) ;
  }
#endif

  if (mris->status == MRIS_PARAMETERIZED_SPHERE && DIAG_VERBOSE_ON)
  {
    MRI_SP *mrisp = (MRI_SP *)mris->vp ;

    sprintf(fname, "%s%4.4d.hipl", parms->base_name, t) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "writing %s\n", fname) ;
    MRIStoParameterization(mris, mrisp, 1, 0) ;
    MRISPwrite(mrisp, fname) ;
  }

  return(NO_ERROR) ;
}

int
MRISblurCurvature(MRI_SURFACE *mris, double sigma)
{
  MRI_SP *mrisp ;
  
  mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
  MRISPblur(mrisp, mrisp, sigma, 0) ;
  MRISfromParameterization(mrisp, mris, 0);
  MRISPfree(&mrisp) ;
  return(NO_ERROR) ;
}
