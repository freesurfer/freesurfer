#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"

static char vcid[] = "$Id: mris_register.c,v 1.3 1998/02/27 20:57:42 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int max_passes = 4 ;
static int nbrs = 2 ;
static float scale = 1.0f ;

char *Progname ;
static char curvature_fname[100] = "" ;

static INTEGRATION_PARMS  parms ;

int
main(int argc, char *argv[])
{
  char         **av, *surf_fname, *template_fname, *out_fname, fname[100], *cp;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp_template ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = PROJECT_SPHERE ;
  parms.tol = 1e-2 ;
  parms.min_averages = 0 ;
  parms.l_dist = 0.0 ;
  parms.l_area = 0.1 ;
  parms.l_parea = 0.0f ;
  parms.l_corr = 1.0f ;
  parms.l_pcorr = 0.0f ;
  parms.niterations = 25 ;
  parms.n_averages = 128 ;
  parms.write_iterations = 100 ;
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.99 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.dt_increase = 1.0 ;
  parms.dt_decrease = 1.0 ;
  parms.error_ratio = 1.1 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_ADAPTIVE ;
  parms.integration_type = INTEGRATE_MOMENTUM /*INTEGRATE_LINE_MINIMIZE*/ ;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.dt = 0.9 ;
  parms.momentum = 0.95 ;
  parms.desired_rms_height = -1.0 ;
  parms.nbhd_size = 7 ;
  parms.max_nbrs = 8 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  surf_fname = argv[1] ;
  template_fname = argv[2] ;
  out_fname = argv[3] ;
  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "sphere") ;
  }

  fprintf(stderr, "reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if (curvature_fname[0])
  {
    fprintf(stderr, "reading source curvature from %s\n",curvature_fname) ;
    MRISreadCurvatureFile(mris, curvature_fname) ;
  }
  fprintf(stderr, "reading template parameterization from %s...\n",
          template_fname) ;
#if 1
  mrisp_template = MRISPread(template_fname) ;
  if (!mrisp_template)
    ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
                Progname, template_fname) ;
#else
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISclearCurvature(mris) ;
  mris->vertices[0].curv = 1.0f ;
  MRISnormalizeCurvature(mris) ;
  mrisp_template = MRISPalloc(1, 6) ;
  MRISrotate(mris, mris, RADIANS(30.0f), RADIANS(30.0f), 0.0f) ;
  MRIStoParameterization(mris, mrisp_template, 1, 0) ;
  *IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2) = 1 ;   /* dof */
  /*  MRISPwrite(mrisp_template, "temp1.hipl") ;*/
  MRIStoParameterization(mris, mrisp_template, 1, 3) ;
  *IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 5) = 1 ;   /* dof */
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  /*  MRISPwrite(mrisp_template, "temp.hipl") ;*/
#endif
  
  
  MRISsetNeighborhoodSize(mris, nbrs) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  mris->status = MRIS_PARAMETERIZED_SPHERE ;
  MRIScomputeMetricProperties(mris) ;
  if (!FZERO(parms.l_dist))
    MRISscaleDistances(mris, scale) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISzeroNegativeAreas(mris) ;
  MRISstoreMetricProperties(mris) ;
  MRISstoreMeanCurvature(mris) ;  /* use curvature from file */

  MRISregister(mris, mrisp_template, &parms, max_passes) ;
  fprintf(stderr, "writing registered surface to %s...\n", out_fname) ;
  MRISwrite(mris, out_fname) ;

  MRISPfree(&mrisp_template) ;
  MRISfree(&mris) ;
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
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
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "parea"))
  {
    sscanf(argv[2], "%f", &parms.l_parea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_parea = %2.3f\n", parms.l_parea) ;
  }
  else if (!stricmp(option, "spring"))
  {
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "corr"))
  {
    sscanf(argv[2], "%f", &parms.l_corr) ;
    nargs = 1 ;
    fprintf(stderr, "using l_corr = %2.3f\n", parms.l_corr) ;
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
  else switch (toupper(*option))
  {
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
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
    fprintf(stderr, "using write iterations = %d\n", parms.write_iterations) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'P':
    max_passes = atoi(argv[2]) ;
    fprintf(stderr, "limitting unfolding to %d passes\n", max_passes) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
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
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
       "usage: %s [options] <input surface> <average surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program register a surface with  an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

