
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

static char vcid[] = "$Id: mris_inflate.c,v 1.10 1998/02/01 00:56:28 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static INTEGRATION_PARMS  parms ;
static int desired_curvature_set = 0 ;
static int talairach_flag = 0 ;

#define DESIRED_RADIUS   50.0   /* 10 cm radius */

#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;
int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, fname[100], *cp, path[100] ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  double       radius, scale ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.base_name[0] = 0 ;
  parms.projection = NO_PROJECTION ;
  parms.fi_desired = 25 ;
  parms.ici_desired = 3.5 ;
  parms.tol = 1e-4 ;
  parms.epsilon = EPSILON ;
  parms.dt = 0.35 ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;
  parms.n_averages = 32 ; /*N_AVERAGES*/ ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_dist = .1 ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_spring = 1.0 ;
  parms.l_curv = 0.0 ;
  parms.niterations = 80 ;   
  parms.write_iterations = 50 /*WRITE_ITERATIONS */;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.9 ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 1.0 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "inflated") ;
  }

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  MRISsetNeighborhoodSize(mris, 2) ;
  MRISreadOriginalProperties(mris, NULL) ;

  if (talairach_flag)
    MRIStalairachTransform(mris, mris) ;

  radius = MRISaverageRadius(mris) ;
  scale = DESIRED_RADIUS/radius ;
  MRISscaleBrain(mris, mris, scale) ;

  MRIScomputeSecondFundamentalForm(mris) ;
  MRISinflateBrain(mris, &parms) ;
  fprintf(stderr, "writing inflated surface to %s\n", out_fname) ;
  radius = MRISaverageRadius(mris) ;
  scale = DESIRED_RADIUS/radius ;
  MRISscaleBrain(mris, mris, scale) ;
  MRISwrite(mris, out_fname) ;
  FileNamePath(out_fname, path) ;
  sprintf(fname, "%s/%s.sulc", path, 
          mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;
  fprintf(stderr, "writing sulcal depths to %s\n", fname) ;
  MRISwriteCurvature(mris, fname) ;

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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  }
  else if (!stricmp(option, "angle"))
  {
    sscanf(argv[2], "%f", &parms.l_angle) ;
    nargs = 1 ;
    fprintf(stderr, "l_angle = %2.3f\n", parms.l_angle) ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "spring"))
  {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "no_proj"))
  {
    parms.projection = NO_PROJECTION ;
    fprintf(stderr, "disabling ellipsoid projection\n") ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  }
  else if (!stricmp(option, "hvariable"))
  {
    parms.flags |= IPFLAG_HVARIABLE ;
    fprintf(stderr, "variable Hdesired to drive integration\n") ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
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
  else if (!stricmp(option, "dt_dec"))
  {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else switch (toupper(*option))
  {
  case 'T':
    talairach_flag = 1 ;
    fprintf(stderr, "applying talairach transform to brain before inflation\n") ;
    break ;
  case 'S':
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
    break ;
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
    break ;
  case 'F':
    parms.fi_desired = atof(argv[2]) ;
    fprintf(stderr, "desired fi=%2.2f\n", parms.fi_desired) ;
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
  case 'E':
    parms.epsilon = atof(argv[2]) ;
    fprintf(stderr, "using epsilon=%2.4f\n", parms.epsilon) ;
    nargs = 1 ;
    break ;
  case 'H':
    desired_curvature_set = 1 ;
    parms.Hdesired = atof(argv[2]) ;
    fprintf(stderr, "inflating to desired curvature %2.4f\n", parms.Hdesired);
    nargs = 1 ;
    break ;
  case 'W':
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'A':
    sscanf(argv[2], "%d", &parms.n_averages) ;
    nargs = 1 ;
    fprintf(stderr, "n_averages = %d\n", parms.n_averages) ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "niterations = %d\n", parms.niterations) ;
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
          "usage: %s [options] <input surface file> <output surface file>"
          "\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will inflate a cortical surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

