
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
#include "utils.h"

static char vcid[] = "$Id: mris_flatten.c,v 1.13 1998/03/25 16:24:31 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISscaleUp(MRI_SURFACE *mris) ;

char *Progname ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE     1.0
static float base_dt_scale = BASE_DT_SCALE ;
static int nbrs = 2 ;
static int inflate = 0 ;
static double disturb = 0 ;
static int mrisDisturbVertices(MRI_SURFACE *mris, double amount) ;
static int randomly_flatten = 0 ;
static int   nospring = 0 ;
static float scale = 3 ;
static int   max_passes = 1 ;

#define SMOOTHWM_FNAME  "smoothwm" 

int
main(int argc, char *argv[])
{
  char         **av, in_surf_fname[100], *in_patch_fname, *out_patch_fname, 
               fname[100], path[100], *cp, hemi[10] ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;
  parms.dt = .1 ;
  parms.projection = PROJECT_PLANE ;
  parms.tol = 1e-1 ;
  parms.n_averages = 1024 /* 1024 */ ;
  parms.min_averages = 0 ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_neg = 0.0 ;
  parms.l_dist = 1.0 ;
  parms.l_spring = 0.0 ;
  parms.l_area = 1.0 ;
  parms.l_boundary = 0.0 ;
  parms.l_curv = 0.0 ;
  parms.niterations = 25 ;
  parms.write_iterations = 1000 ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.98 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.momentum = 0.9 ;
  parms.desired_rms_height = -1.0 ;
  parms.base_name[0] = 0 ;
  parms.Hdesired = 0.0 ;   /* a flat surface */
  parms.nbhd_size = 4 ;    /* out to 4-connected neighbors */
  parms.max_nbrs = 20 ;    /* 20 at each distance */

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

  parms.base_dt = base_dt_scale * parms.dt ;
  in_patch_fname = argv[1] ;
  out_patch_fname = argv[2] ;
  FileNamePath(in_patch_fname, path) ;
  cp = strrchr(in_patch_fname, '/') ;
  cp = strchr(cp, '.') ;
  if (cp)
  {
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
  }
  else
    strcpy(hemi, "lh") ;
  sprintf(in_surf_fname, "%s/%s.%s", path, hemi, SMOOTHWM_FNAME) ;

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_patch_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "flattened") ;
  }

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;


  if (MRISreadPatch(mris, in_patch_fname) != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "%s: could not read patch file %s",
              Progname, in_patch_fname) ;
  fprintf(stderr, "reading original vertex positions...\n") ;
  if (!FZERO(disturb))
    mrisDisturbVertices(mris, disturb) ;
  if (parms.niterations > 0)
  {
    if (randomly_flatten)
      MRISflattenPatchRandomly(mris) ;
    else
      MRISflattenPatch(mris) ;
    
    
    MRISreadOriginalProperties(mris, "smoothwm") ;
    MRISsetNeighborhoodSize(mris, nbrs) ;

    /* optimize metric properties of flat map */
    fprintf(stderr,"minimizing metric distortion induced by projection...\n");
    MRISscaleBrain(mris, mris, scale) ;
    MRIScomputeMetricProperties(mris) ;
    MRISunfold(mris, &parms, max_passes) ;
    fprintf(stderr, "writing flattened patch to %s\n", out_patch_fname) ;
    MRISwritePatch(mris, out_patch_fname) ;
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
  float f ;
  
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
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
  }
  else if (!stricmp(option, "curv"))
  {
    sscanf(argv[2], "%f", &parms.l_curv) ;
    nargs = 1 ;
    fprintf(stderr, "using l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "nospring"))
    nospring = 1 ;
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "boundary"))
  {
    sscanf(argv[2], "%f", &parms.l_boundary) ;
    nargs = 1 ;
    fprintf(stderr, "using l_boundary = %2.3f\n", parms.l_boundary) ;
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
  else if (!stricmp(option, "spring"))
  {
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using base name = %s\n", parms.base_name) ;
  }

  else if (!stricmp(option, "angle"))
  {
    sscanf(argv[2], "%f", &parms.l_angle) ;
    nargs = 1 ;
    fprintf(stderr, "using l_angle = %2.3f\n", parms.l_angle) ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
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
  case 'P':
    max_passes = atoi(argv[2]) ;
    fprintf(stderr, "limitting unfolding to %d passes\n", max_passes) ;
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
  case 'D':
    disturb = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "disturbing vertex positions by %2.3f\n",(float)disturb) ;
    break ;
  case 'R':
    randomly_flatten = !randomly_flatten ;
    fprintf(stderr, "using random placement for flattening.\n") ;
    break ;
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
    break ;
  case 'S':
    scale = atof(argv[2]) ;
    fprintf(stderr, "scaling brain by = %2.3f\n", (float)scale) ;
    nargs = 1 ;
    break ;
  case 'W':
    Gdiag |= DIAG_WRITE ;
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "using write iterations = %d\n", parms.write_iterations) ;
    break ;
  case 'I':
    inflate = 1 ;
    fprintf(stderr, "inflating brain...\n") ;
    break ;
  case 'A':
    sscanf(argv[2], "%d", &parms.n_averages) ;
    nargs = 1 ;
    fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
    break ;
  case '?':
  case 'U':
    print_help() ;
    exit(1) ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
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
          "usage: %s [options] <input patch> <output patch>\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will flatten a surface patch\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, " -w <# iterations>  "
          "write out the surface every # of iterations.\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
mrisDisturbVertices(MRI_SURFACE *mris, double amount)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->x += randomNumber(-amount, amount) ;
    v->y += randomNumber(-amount, amount) ;
  }

  MRIScomputeMetricProperties(mris) ;
  return(NO_ERROR) ;
}
int
MRISscaleUp(MRI_SURFACE *mris)
{
  int     vno, n, max_v, max_n ;
  VERTEX  *v ;
  float   ratio, max_ratio ;

  max_ratio = 0.0f ; max_v = max_n = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      if (FZERO(v->dist[n]))   /* would require infinite scaling */
        continue ;
      ratio = v->dist_orig[n] / v->dist[n] ;
      if (ratio > max_ratio)
      {
        max_v = vno ; max_n = n ;
        max_ratio = ratio ;
      }
    }
  }

  fprintf(stderr, "max @ (%d, %d), scaling brain by %2.3f\n", 
          max_v, max_n, max_ratio) ;
#if 0
  MRISscaleBrain(mris, mris, max_ratio) ;
#else
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
      v->dist_orig[n] /= max_ratio ;
  }
#endif
  return(NO_ERROR) ;
}

