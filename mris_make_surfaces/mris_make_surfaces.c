
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "mrishash.h"
#include "macros.h"

static char vcid[] = "$Id: mris_make_surfaces.c,v 1.10 1999/01/21 00:18:56 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  mrisFindMiddleOfGray(MRI_SURFACE *mris) ;

char *Progname ;

static int navgs = 10 ;
static int create = 0 ;
static float sigma = 0.5f ;  /* should be around 1 */

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;


static int smooth = 2 ;
static int nwhite = 15 /*25*/ ;
static int ngray = 45 /*100*/ ;

static int nowhite = 0 ;
static float gray_surface = 55.0f ;
static int nbrs = 2 ;
static int write_vals = 0 ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm, *mri_kernel = NULL, *mri_smooth, 
                *mri_filled, *mri_T1 ;
  int           label_val, replace_val ;
  int           msec ;
  struct timeb  then ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = 0.75f ;
  parms.base_dt = BASE_DT_SCALE*parms.dt ;
  parms.n_averages = 2 /*8*/ /* 32*/ ; /*N_AVERAGES*/ ;
#if 1
  parms.l_tspring = 1.0f ; parms.l_curv = 0.5 ; parms.l_intensity = 0.1 ;

#else
#if 0
  parms.l_spring = 1.0f ;
#else
  parms.l_nspring = 0.25f ;  /* will be changed for gray matter surface */
  parms.l_tspring = 1.0f ;
#endif
  parms.l_intensity = 0.075 /* 0.025f*/ ;
  parms.l_grad = 1.0f ;
#endif


  parms.niterations = 0 ;
  parms.write_iterations = 5 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = -0.75 ;
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

  if (argc < 2)
    usage_exit() ;

  /* set default parameters for white and gray matter surfaces */
  parms.niterations = nwhite ;
  if (parms.momentum < 0.0)
    parms.momentum = 0.0 /*0.75*/ ;

  TimerStart(&then) ;
  sname = argv[1] ;
  hemi = argv[2] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, 
              "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
  strcpy(sdir, cp) ;
  
  if (sigma > 0.0)
    mri_kernel = MRIgaussian1d(sigma, 100) ;
  sprintf(fname, "%s/%s/mri/filled", sdir, sname) ;
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_filled = MRIread(fname) ;
  if (!mri_filled)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;
  if (!stricmp(hemi, "lh"))
  { label_val = LH_LABEL ; replace_val = RH_LABEL ; }
  else
  { label_val = RH_LABEL ; replace_val = LH_LABEL ; }

  sprintf(fname, "%s/%s/mri/T1", sdir, sname) ;
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1 = MRIread(fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  MRImask(mri_T1, mri_filled, mri_T1, replace_val) ; /* remove other hemi */
  MRIfree(&mri_filled) ;

  sprintf(fname, "%s/%s/mri/wm", sdir, sname) ;
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  if (mri_kernel)
  {
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n", sigma) ;
    mri_smooth = MRIclone(mri_T1, NULL) ;
    MRIconvolveGaussian(mri_T1, mri_smooth, mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) 
      MRIwrite(mri_smooth, "/tmp/brain_smooth.mnc") ;
  }
  else
    mri_smooth = mri_T1 ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, ORIG_NAME) ;
  fprintf(stderr, "reading original surface position from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (nbrs > 1)
    MRISsetNeighborhoodSize(mris, nbrs) ;

  if (!nowhite)
  {
    fprintf(stderr, "computing target intensity values...\n") ;
    if (smooth)
      MRISsmoothSurfaceNormals(mris, smooth);
    MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
    MRIScomputeWhiteSurfaceValues(mris, mri_T1, mri_smooth, mri_wm) ;
    if (write_vals)
    {
      sprintf(fname, "./%s-white.w", hemi) ;
      MRISwriteValues(mris, fname) ;
    }
    fprintf(stderr,"repositioning cortical surface to gray/white boundary\n");
    strcpy(parms.base_name, WHITE_MATTER_NAME) ;
    
    MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);

    sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, WHITE_MATTER_NAME) ;
    fprintf(stderr, "writing white matter surface to %s...\n", fname) ;
    MRISwrite(mris, fname) ;
    if (create)
    {
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s.curv",mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh");
      fprintf(stderr, "writing smoothed curvature to %s\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      sprintf(fname, "%s.area",mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh");
      fprintf(stderr, "writing smoothed area to %s\n", fname) ;
      MRISwriteArea(mris, fname) ;
    }

    if (parms.flags & IPFLAG_NO_SELF_INT_TEST)
    {
      msec = TimerStop(&then) ;
      fprintf(stderr,
              "refinement took %2.1f minutes\n", (float)msec/(60*1000.0f));
      exit(0) ;
    }
  }
  else
  {
    if (MRISreadVertexPositions(mris, "white") != NO_ERROR)
      ErrorExit(Gerror, "%s: could not read white matter surfaces.",
                Progname) ;
    MRIScomputeMetricProperties(mris) ;
  }

  if (mri_kernel)
  {
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n", sigma) ;
    MRIconvolveGaussian(mri_smooth, mri_smooth, mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) 
      MRIwrite(mri_smooth, "/tmp/brain_smooth2.mnc") ;
    MRIfree(&mri_kernel) ;
  }
  else
    mri_smooth = mri_T1 ;

  fprintf(stderr, "repositioning cortical surface to gray/csf boundary.\n") ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
  parms.niterations = ngray ;
  parms.t = parms.start_t = 0.0 ;
  strcpy(parms.base_name, GRAY_MATTER_NAME) ;
  MRIScomputeGraySurfaceValues(mris, mri_T1, mri_smooth, gray_surface) ;
  if (write_vals)
  {
    sprintf(fname, "./%s-gray.w", hemi) ;
    MRISwriteValues(mris, fname) ;
  }
  MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAY_MATTER_NAME) ;
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;
  MRISmeasureCorticalThickness(mris) ;
  fprintf(stderr, 
          "writing cortical thickness estimate to 'thickness' file.\n") ;
  MRISwriteCurvature(mris, "thickness") ;

  /* at this point, the v->curv slots contain the cortical surface. Now
     move the white matter surface out by 1/2 the thickness as an estimate
     of layer IV.
     */
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  mrisFindMiddleOfGray(mris) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAYMID_NAME) ;
  fprintf(stderr, "writing layer IV surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  msec = TimerStop(&then) ;
  fprintf(stderr,"positioning took %2.1f minutes\n", (float)msec/(60*1000.0f));
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
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "write_vals"))
  {
    write_vals = 1 ;
    fprintf(stderr,  "writing gray and white surface targets to w files\n") ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using dt = %2.1e\n", parms.dt) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "spring"))
  {
    parms.l_spring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "grad"))
  {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_grad = %2.3f\n", parms.l_grad) ;
  }
  else if (!stricmp(option, "tspring"))
  {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  }
  else if (!stricmp(option, "nspring"))
  {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "smooth"))
  {
    smooth = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing for %d iterations\n", smooth) ;
  }
  else if (!stricmp(option, "intensity"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "nwhite"))
  {
    nwhite = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, 
           "integrating gray/white surface positioning for %d time steps\n",
           nwhite) ;
  }
  else if (!stricmp(option, "nowhite"))
  {
    nowhite = 1 ;
    fprintf(stderr, "reading previously compute gray/white surface\n") ;
  }
  else if (!stricmp(option, "pial"))
  {
    gray_surface = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, 
           "settting pial surface target value to %2.1f\n", gray_surface) ;
  }
  else if (!stricmp(option, "ngray"))
  {
    ngray = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"integrating pial surface positioning for %d time steps\n",
            ngray) ;
  }
  else if (!stricmp(option, "sigma"))
  {
    sigma = atof(argv[2]) ;
    fprintf(stderr,  "smoothing volume with Gaussian sigma = %2.1f\n", sigma) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'Q':
    parms.flags |= IPFLAG_NO_SELF_INT_TEST ;
    fprintf(stderr, 
            "doing quick (no self-intersection) white matter refinement.\n") ;
    break ;
  case 'A':
    parms.n_averages = atoi(argv[2]) ;
    fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
    nargs = 1 ;
    break ;
  case 'M':
    parms.integration_type = INTEGRATE_MOMENTUM ;
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
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
  case 'C':
    create = 1 ;
    fprintf(stderr, "creating area and curvature files for wm surface...\n") ;
    break ;
  case 'W':
    sscanf(argv[2], "%d", &parms.write_iterations) ;
    nargs = 1 ;
    fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'N':
    sscanf(argv[2], "%d", &parms.niterations) ;
    nargs = 1 ;
    fprintf(stderr, "niterations = %d\n", parms.niterations) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    print_help() ;
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
  fprintf(stderr, "usage: %s [options] <subject name> <hemisphere>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program positions the tessellation of the cortical surface\n"
          "at the white matter surface, then the gray matter surface\n"
          "and generate surface files for these surfaces as well as a\n"
          "'curvature' file for the cortical thickness, and a surface file\n"
          "which approximates layer IV of the cortical sheet.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, 
          "-q    omit self-intersection and only generate "
          "gray/white surface.\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
mrisFindMiddleOfGray(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;
  float   nx, ny, nz, thickness ;

  MRISaverageCurvatures(mris, 3) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ; ny = v->ny ; nz = v->nz ;
    thickness = 0.5 * v->curv ;
    v->x = v->origx + thickness * nx ;
    v->y = v->origy + thickness * ny ;
    v->z = v->origz + thickness * nz ;
  }
  return(NO_ERROR) ;
}

