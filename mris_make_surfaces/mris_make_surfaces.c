
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

static char vcid[] = "$Id: mris_make_surfaces.c,v 1.3 1998/10/28 15:35:31 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  mrisFindMiddleOfGray(MRI_SURFACE *mris) ;

char *Progname ;

static float nsigma = 0.0f ;  /* not used for gray matter currently */
static float sigma = 0.75f ;  /* should be around 1 */

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE    1.0
static float base_dt_scale = BASE_DT_SCALE ;

static int sigma_set = 0 ;
static int nsigma_set = 0 ;

static int nwhite = 15 /*25*/ ;
static int ngray = 65 /*100*/ ;

static float gray_surface = 65.0f ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm, *mri_brain, *mri_kernel = NULL, *mri_smooth, 
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
#if 0
  parms.l_spring = 1.0f ;
#else
  parms.l_nspring = 0.25f ;  /* will be changed for gray matter surface */
  parms.l_tspring = 1.0f ;
#endif
  parms.l_intensity = 0.075 /* 0.025f*/ ;
  parms.l_grad = 1.0f ;
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
  if (!sigma_set)
    sigma = 0.25f ;
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
  sprintf(fname, "%s/%s/mri/brain", sdir, sname) ;
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_brain = MRIread(fname) ;
  if (!mri_brain)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

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
  MRImask(mri_brain, mri_filled, mri_brain, replace_val) ;
  /*  MRIwrite(mri_brain, "/tmp/brain.mnc") ;*/
  MRIfree(&mri_filled) ;


  sprintf(fname, "%s/%s/mri/T1", sdir, sname) ;
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_T1 = MRIread(fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

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
    MRIfree(&mri_kernel);
  }
  else
    mri_smooth = mri_T1 ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, ORIG_NAME) ;
  fprintf(stderr, "reading spherical surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  MRIScomputeWhiteSurfaceValues(mris, mri_brain, mri_wm, nsigma) ;
  fprintf(stderr, "repositioning cortical surface to gray/white boundary.\n");
  if (!nsigma_set)
    nsigma = 0.1f ;
  MRISaverageVertexPositions(mris, 4) ;  /* start with a smooth surface */
  strcpy(parms.base_name, WHITE_MATTER_NAME) ;


  MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, WHITE_MATTER_NAME) ;
  fprintf(stderr, "writing white matter surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;

  if (!nsigma_set)
    nsigma = 1.5f ;
  fprintf(stderr, "repositioning cortical surface to gray/csf boundary.\n") ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */
  parms.niterations = ngray ;
  parms.t = parms.start_t = 0.0 ;
  parms.l_nspring = 0.5 ;
  strcpy(parms.base_name, GRAY_MATTER_NAME) ;
  MRIScomputeGraySurfaceValues(mris, mri_brain, mri_wm, nsigma, gray_surface) ;
  MRISpositionSurface(mris, mri_T1, mri_smooth,&parms);

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAY_MATTER_NAME) ;
  fprintf(stderr, "writing pial surface to %s...\n", fname) ;
  MRISwrite(mris, fname) ;
  MRISmeasureCorticalThickness(mris, mri_brain) ;
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
  else if (!stricmp(option, "nsigma"))
  {
    nsigma_set = 1 ;
    nsigma = atof(argv[2]) ;
    fprintf(stderr,  "searching for boundary %2.1f standard deviations "
            "from mean\n", nsigma) ;
    nargs = 1 ;
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
            nwhite) ;
  }
  else if (!stricmp(option, "sigma"))
  {
    sigma_set = 1 ;
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
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
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

