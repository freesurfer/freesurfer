
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

static char vcid[] = "$Id: mris_flatten.c,v 1.3 1997/12/15 19:50:20 fischl Exp $";

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
static float min_neg_pct = 0.05f/100.0f ;  /* less than 0.05% negative */
static int   min_neg = 20 ;
static int   nospring = 0 ;

int
main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *in_patch_fname, *out_patch_fname, 
               fname[100], path[100], *cp ;
  int          ac, nargs, niterations, naverages, int_type, write_iterations ;
  MRI_SURFACE  *mris ;
  float        l_spring, l_dist, momentum, dt, l_area ;
  double       tol ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.dt = .1 ;
  parms.projection = PROJECT_PLANE ;
  parms.tol = 1e-2 ;
  parms.n_averages = 1024 ;
  parms.min_averages = 0 ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_neg = 0.0 ;
  parms.l_dist = 1.0 ;
  parms.l_spring = 0.0 ;
  parms.l_area = 1.0 ;
  parms.shrink = 1.0 ;
  parms.l_boundary = 0.0 ;
  parms.l_curv = 0.0 ;
  parms.niterations = 1 ;
  parms.write_iterations = 25 ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.98 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.momentum = 0.9 ;
  parms.fi_desired = -1.0 ;
  parms.ici_desired = -1.0 ;
  parms.base_name[0] = 0 ;
  parms.Hdesired = 0.0 ;   /* a flat surface */
  parms.nbhd_size = 8 ;
  parms.max_nbrs = 4 ;

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

  parms.base_dt = base_dt_scale * parms.dt ;
  in_surf_fname = argv[1] ;
  in_patch_fname = argv[2] ;
  out_patch_fname = argv[3] ;
  FileNamePath(in_surf_fname, path) ;

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
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading original vertex positions...") ;
  if (!FZERO(disturb))
    mrisDisturbVertices(mris, disturb) ;
#if 0
  MRISstoreCurrentPositions(mris) ;
  MRISreadVertexPositions(mris, "smoothwm") ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n"
            "Computing metric properties of original surface...") ;
  MRISsetNeighborhoodSize(mris, nbrs) ;
  MRIScomputeMetricProperties(mris) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "sampling distances on original surface...") ;
  MRISsampleDistances(mris, nbr_size, max_nbrs) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done, with avg # of neighbors = %2.1f\n", mris->avg_nbrs);
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeTriangleProperties(mris, 0) ;  /* hack */
  MRISstoreMetricProperties(mris) ;
  MRISrestoreOldPositions(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISupdateSurface(mris) ;
#endif
  if (parms.niterations > 0)
  {
    if (inflate)
      MRISinflateBrain(mris, &parms) ;
    else
    {
      if (randomly_flatten)
        MRISflattenPatchRandomly(mris) ;
      else
        MRISflattenPatch(mris) ;

      naverages = parms.n_averages ; parms.n_averages = 0 ;
      l_spring = parms.l_spring ; l_dist = parms.l_dist ;
      l_area = parms.l_area ;
      dt = parms.dt ; momentum = parms.momentum ; 
      int_type = parms.integration_type ;
      niterations = parms.niterations ; tol = parms.tol ;
      write_iterations = parms.write_iterations ;

      if (write_iterations > 0)
        parms.write_iterations = 50 ;
      /*      parms.momentum = 0.95 ; parms.dt = .95 ; */
      parms.integration_type = INTEGRATE_MOMENTUM /*INTEGRATE_ADAPTIVE */;
      parms.l_spring = 0.0 ; parms.l_dist = 0.01 ; parms.l_boundary = 0.0 ;
      parms.l_area = 1.0 ;
      parms.niterations = 1 ;
      parms.tol = 1e-5 ;
      MRISstoreCurrentPositions(mris) ;
      mris->status = MRIS_PATCH ;
      MRISreadVertexPositions(mris, "smoothwm") ;
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeTriangleProperties(mris, 0) ;  /* hack */
      MRISstoreMetricProperties(mris) ;
      MRISrestoreOldPositions(mris) ;
      mris->status = MRIS_PLANE ;
      MRISupdateSurface(mris) ;
#if 0
#if 1
      if (!nospring)
        MRISremoveNegativeVertices(mris, &parms, min_neg, min_neg_pct) ;
#else
      if (!nospring)
        MRISunfold(mris, &parms) ;  /* use spring force to remove neg. vert. */
#endif
#endif
      /* restore user-specified parameters */
      parms.l_spring = l_spring ; parms.l_dist = l_dist ;
      parms.l_area = l_area ;
      parms.niterations = niterations ; parms.n_averages = naverages ;
      parms.momentum = momentum ; parms.dt = dt ; 
      parms.integration_type = int_type ;
      parms.tol = tol ;
      parms.write_iterations = write_iterations ;

      /* read in original positions and calculate distances */
      if (nbrs > 1)
      {
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        mris->status = MRIS_PATCH ;  /* so no orientating will be done */
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
        MRISsetNeighborhoodSize(mris, nbrs) ;
        MRIScomputeMetricProperties(mris) ;
        MRIScomputeTriangleProperties(mris, 0) ;  /* hack */
        MRISstoreMetricProperties(mris) ;
        
        /* restore the current positions and properties */
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        mris->status = MRIS_PLANE ;
        MRIScomputeMetricProperties(mris) ;
        MRIScomputeTriangleProperties(mris, 0) ;  /* hack */
        MRISupdateSurface(mris) ;
      }

      if (Gdiag & DIAG_SHOW)
        fprintf(stderr,"surface unfolded - minimizing metric distortion...\n");
      /*      MRISscaleUp(mris) ;*/
      MRISunfold(mris, &parms) ;  /* optimize metric properties of flat map */
    }
    if (Gdiag & DIAG_SHOW)
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
#if 0
  else if (!stricmp(option, "neg"))
  {
    sscanf(argv[2], "%f", &parms.l_neg) ;
    nargs = 1 ;
    fprintf(stderr, "using l_neg = %2.3f\n", parms.l_neg) ;
  }
#endif
  else if (!stricmp(option, "nospring"))
    nospring = 1 ;
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "neg"))
  {
    min_neg = atoi(argv[2]) ;
    min_neg_pct = (float)atof(argv[3])/100.0f ;
    nargs = 2 ;
    fprintf(stderr,"negative vertex thresholds: count: %d or area: %2.2f%%\n",
            min_neg, 100.0f*min_neg_pct) ;
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
    parms.shrink = atof(argv[2]) ;
    fprintf(stderr, "using shrink = %2.3f\n", (float)parms.shrink) ;
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
    print_usage() ;
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
          "usage: %s [options] <surface file> <patch file name> <output patch>"
          "\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
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

