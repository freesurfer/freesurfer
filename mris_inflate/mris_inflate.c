
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

static char vcid[] = "$Id: mris_inflate.c,v 1.4 1997/12/15 20:34:32 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static INTEGRATION_PARMS  parms ;
static int patch_flag = 0 ;
static int desired_curvature_set = 0 ;

#define DESIRED_RADIUS   50.0   /* 10 cm radius */

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, path[100], fname[100],*cp,hemi[10];
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  double       radius, scale ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  strcpy(parms.base_name, "inflated") ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = parms.base_dt = 0.01 ;
  parms.n_averages = 0 ; /*N_AVERAGES*/ ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_curv = 1.0 ;
  parms.niterations = 1000 ;
  parms.write_iterations = 25 /*WRITE_ITERATIONS */;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = MOMENTUM ;
  parms.dt_increase = DT_INCREASE ;
  parms.dt_decrease = DT_DECREASE ;
  parms.error_ratio = ERROR_RATIO ;
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

  if (patch_flag)
  {
    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
    {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    }
    else
      strcpy(hemi, "lh") ;
    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;

    FileNameOnly(in_fname, fname) ;
    if (MRISreadPatch(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;
  }
  else
  {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }

  MRIScomputeMetricProperties(mris) ;
  MRIStalairachTransform(mris, mris) ;
  radius = MRISaverageRadius(mris) ;
  scale = DESIRED_RADIUS/radius ;
  MRISscaleBrain(mris, mris, scale) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "scaling brain by %2.2f\n", scale) ;

  if (parms.niterations > 0)
  {
    MRISunfold(mris, &parms) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "writing inflated surface to %s\n", out_fname) ;
    MRISwrite(mris, out_fname) ;
  }


#if 0
  sprintf(fname, "%s.area_error", out_fname) ;
  printf("writing area errors to %s\n", fname) ;
  MRISwriteAreaError(mris, fname) ;
  sprintf(fname, "%s.angle_error", out_fname) ;
  printf("writing angle errors to %s\n", fname) ;
  MRISwriteAngleError(mris, fname) ;
  MRISfree(&mris) ;
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
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
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
    parms.tol = TOL ;
    fprintf(stderr, "integrating with line minimization (tol = %2.2e)\n", 
            parms.tol) ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.dt = parms.base_dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.base_dt) ;
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
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
    break ;
  case 'P':
    patch_flag = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
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

