/**
 * @brief cmd line utility for registering a subject's surface with an atlas,
 *
 * Command line utility for registering a subject's surface with an atlas.
 * Can also be used with -1 to register a surface to another individual surfce.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)

static char *starting_reg_fname = NULL ;
static int multi_scale = 0 ;
static int which_norm = NORM_MEAN ;
static int single_surf = 0 ;
static double l_ocorr = 1.0 ;
static char *annot_name = NULL ;
static int max_passes = 4 ;
static float min_degrees = 0.5 ;
static float max_degrees = 64.0 ;
static int nangles = 8 ;
static int nbrs = 1 ;
static float scale = 1.0f ;

static int target_hemi = LEFT_HEMISPHERE ;
static int reverse_flag = 0 ;

static float dalpha = 0.0f ;
static float dbeta = 0.0f ;
static float dgamma = 0.0f ;

static int navgs = 0 ;

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

static int use_defaults = 1 ;

static INTEGRATION_PARMS  parms ;
static int remove_negative = 1 ;

int
main(int argc, char *argv[])
{
  char         **av, *lh_surf_fname, *rh_surf_fname, *out_fname, fname[STRLEN],*cp ;
  int          ac, nargs, msec, sno, h, i ;
  MRI_SURFACE  *mris_lh, *mris_rh, *mris_template, *mris_mov ;
  MRI_SP       *mrisp_template ;
  float        *lh_coords[3], *rh_coords[3], **coords ;

  nargs = handleVersionOption(argc, argv, "mris_left_right_register");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Timer start;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

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
  parms.niterations = 100 ;
  parms.n_averages = 1024 ;   // used to be 256
  parms.first_pass_averages = 1024*16 ;  // only used in first pass
  parms.write_iterations = 100 ;
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

  lh_surf_fname = argv[1] ;
  rh_surf_fname = argv[2] ;
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

  fprintf(stderr, "reading lh surface from %s...\n", lh_surf_fname) ;
  mris_lh = MRISread(lh_surf_fname) ;
  if (!mris_lh)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh surface file %s",
              Progname, lh_surf_fname) ;

  fprintf(stderr, "reading rh surface from %s...\n", rh_surf_fname) ;
  mris_rh = MRISread(rh_surf_fname) ;
  if (!mris_rh)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh surface file %s",
              Progname, rh_surf_fname) ;


  MRISresetNeighborhoodSize(mris_lh, 1) ;
  MRISresetNeighborhoodSize(mris_rh, 1) ;
    
  if (MRISreadOriginalProperties(mris_lh, orig_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh original properties from %s",Progname, orig_name) ;
  if (MRISreadCanonicalCoordinates(mris_lh, canon_name) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not read lh canon surface %s", Progname, canon_name) ;
  if (MRISreadOriginalProperties(mris_rh, orig_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh original properties from %s",Progname, orig_name) ;
  if (MRISreadCanonicalCoordinates(mris_rh, canon_name) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not read rh canon surface %s", Progname, canon_name) ;


  if (target_hemi == LEFT_HEMISPHERE)
  {
    MRISreverse(mris_rh, REVERSE_X, 1) ;
    MRISreverseCoords(mris_rh, REVERSE_X, 0, CANONICAL_VERTICES) ;  // only reverse faces once
    MRISreverseCoords(mris_rh, REVERSE_X, 0, ORIGINAL_VERTICES) ;  // only reverse faces once
  }
  else
  {
    MRISreverse(mris_lh, REVERSE_X, 1) ;
    MRISreverseCoords(mris_lh, REVERSE_X, 0, CANONICAL_VERTICES) ;
    MRISreverseCoords(mris_lh, REVERSE_X, 0, ORIGINAL_VERTICES) ;
  }
  for (i = 0 ; i < 3 ; i++)
  {
    lh_coords[i] = (float *)calloc(mris_lh->nvertices, sizeof(float)) ;
    if (lh_coords[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d lh coords", Progname, mris_lh->nvertices) ;
    rh_coords[i] = (float *)calloc(mris_rh->nvertices, sizeof(float)) ;
    if (rh_coords[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d rh coords", Progname, mris_rh->nvertices) ;
  }

  for (h = LEFT_HEMISPHERE ; h <= RIGHT_HEMISPHERE ; h++)   // register left to right and right to left
  {
    char surf_dir[STRLEN], *template_fname;
    const char *template_hemi ;

    mrisp_template = MRISPalloc(scale, PARAM_IMAGES);
    if (h == LEFT_HEMISPHERE)   // we are moving the lh
    {
      mris_template = mris_rh ; template_hemi = "rh" ; template_fname = rh_surf_fname ; mris_mov = mris_lh ; coords = lh_coords ;
    }
    else
    {
      mris_template = mris_lh ; template_hemi = "lh" ; template_fname = lh_surf_fname ; mris_mov = mris_rh ; coords = rh_coords ;
    }
    for (sno = 0; sno < SURFACES ; sno++)
    {
      FileNamePath(template_fname, surf_dir) ;
      if (curvature_names[sno])  /* read in precomputed curvature file */
      {
	int req = snprintf(fname, STRLEN, "%s/%s.%s", surf_dir, template_hemi, curvature_names[sno]) ;
  if (req >= STRLEN) {
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
	int req = snprintf(fname, STRLEN, "%s/%s.%s", surf_dir, template_hemi, surface_names[sno]) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
	if (MRISreadVertexPositions(mris_template, fname) != NO_ERROR)
	  ErrorExit(ERROR_NOFILE,
		    "%s: could not read surface file %s",
		    Progname, fname) ;
	
	MRIScomputeMetricProperties(mris_template) ;
	MRIScomputeSecondFundamentalForm(mris_template) ;
	MRISuseMeanCurvature(mris_template) ;
	MRISaverageCurvatures(mris_template, navgs) ;
	MRISrestoreVertexPositions(mris_template, CANONICAL_VERTICES) ;
	MRISnormalizeCurvature(mris_template, which_norm) ;
      }
      MRIStoParameterization(mris_template, mrisp_template, scale, sno*IMAGES_PER_SURFACE) ;
      MRISPsetFrameVal(mrisp_template, sno*IMAGES_PER_SURFACE+1, 1.0) ;
    }
    //niter = parms.niterations ; parms.niterations = 0 ;  // only rigid
    MRISregister(mris_mov, mrisp_template,
		 &parms, max_passes,
		 min_degrees, max_degrees, nangles) ;
    MRISextractVertexCoords(mris_mov, coords, CURRENT_VERTICES) ;
    MRISrestoreVertexPositions(mris_mov, CANONICAL_VERTICES) ;
    //parms.niterations = niter ;
    parms.start_t = 0 ;
  }

  MRISfreeDistsButNotOrig(mris_lh);
  MRISfreeDistsButNotOrig(mris_rh);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  // average lh and rh warps and projects back onto sphere
  for (h = LEFT_HEMISPHERE ; h <= RIGHT_HEMISPHERE ; h++)   // register left to right and right to left
  {
    int    vno ;
    if (h == LEFT_HEMISPHERE)
    {
      mris_mov = mris_lh ; coords = lh_coords ;
    }
    else
    {
      mris_mov = mris_rh ; coords = rh_coords ;
    }

    float const radius = MRISaverageRadius(mris_mov);
    for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
    {
      VERTEX *v = &mris_mov->vertices[vno] ;
      float 
        x = (v->cx + coords[0][vno])/2,
        y = (v->cy + coords[1][vno])/2,
        z = (v->cz + coords[2][vno])/2;
        
      float norm = sqrt(x*x + y*y + z*z) ;
      
      MRISsetXYZ(mris_mov, vno,
        (x/norm)*radius,
        (y/norm)*radius,
        (z/norm)*radius);
    }
    
    printf("writing left/right registered hemi to %s\n", argv[3+h]) ;
    MRISwrite(mris_mov, argv[3+h]) ;
  }

  for (i = 0 ; i < 3 ; i++)
  {
    free(lh_coords[i]) ;
    free(rh_coords[i]) ;
  }
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
  int    nargs = 0 ;
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
  else if (!stricmp(option, "annot"))
  {
    annot_name = argv[2] ;
    fprintf(stderr,"zeroing medial wall in %s\n", annot_name) ;
    nargs=1;
  }
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
    strcpy(const_cast<char*>(curvature_names[0]), fname) ; // strcpy???? And the const_cast is.... ugly
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
    fprintf(stderr, "using %s as curvature function for surface 0.\n",
            curvature_names[0]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "curv1"))
  {
    curvature_names[1]  = argv[2];
    fprintf(stderr, "using %s as curvature function for surface 1.\n",
            curvature_names[1]) ;
    nargs = 1 ;
  }

  else if (!stricmp(option, "curv2"))
  {
    curvature_names[2]  = argv[2];
    fprintf(stderr, "using %s as curvature function for surface 2.\n",
            curvature_names[2]) ;
    nargs = 1 ;
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
//      gMRISexternalSSE = gcsaSSE ;
      break ;
    case 'E':
      parms.l_external = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting l_external = %2.1f\n", parms.l_external) ;
      break ;
    case 'C':
      strncpy(curvature_fname, argv[2], STRLEN) ; // Convert to strncpy at least
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
  printf("usage: mris_left_right_register lh.sphere rh.sphere lh.sphere.left_right rh.sphere.left_right\n") ;
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



