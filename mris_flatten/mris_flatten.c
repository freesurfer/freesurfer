/**
 * @file  mris_flatten.c
 * @brief flatten a surface patch
 *
 * "Cortical Surface-Based Analysis II: Inflation, Flattening, and a
 * Surface-Based Coordinate System", Fischl, B., Sereno, M.I., Dale, A.M.
 * (1999) NeuroImage, 9(2):195-207.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:32 $
 *    $Revision: 1.37 $
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

#include "macros.h"
#include "cma.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "utils.h"
#include "version.h"
#include "fastmarching.h"

static char vcid[] =
  "$Id: mris_flatten.c,v 1.37 2011/03/02 00:04:32 nicks Exp $";

int main(int argc, char *argv[]) ;


static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISscaleUp(MRI_SURFACE *mris) ;

char *Progname ;

static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE     1.0
static float base_dt_scale = BASE_DT_SCALE ;
static char *label_fname = NULL ;
static int nbrs = 2 ;
static int inflate = 0 ;
static double disturb = 0 ;
static int mrisDisturbVertices(MRI_SURFACE *mris, double amount) ;
static int randomly_flatten = 0 ;
static int   nospring = 0 ;
static float scale = 3 ;
static int   max_passes = 1 ;

static int sphere_flag = 0 ;
static int plane_flag = 0 ;
static int dilate = 0 ;
static int dilate_label = 0 ; // how many times to dilate label after reading

static int one_surf_flag = 0 ;
static char *original_surf_name = SMOOTH_NAME ;
static char *original_unfold_surf_name = ORIG_NAME ;
static float rescale = 1.0f ;


int
main(int argc, char *argv[])
{
  char         **av, in_surf_fname[STRLEN], *in_patch_fname, *out_patch_fname,
  fname[STRLEN], path[STRLEN], *cp, hemi[10] ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_flatten.c,v 1.37 2011/03/02 00:04:32 nicks Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= (DIAG_SHOW | DIAG_WRITE) ;
  memset(&parms, 0, sizeof(parms)) ;
  parms.dt = .1 ;
  parms.projection = PROJECT_PLANE ;
  parms.tol = 0.2 ;
  parms.n_averages = 1024 ;
  parms.l_dist = 1.0 ;
  parms.l_nlarea = 1.0 ;
  parms.niterations = 40 ;
  parms.area_coef_scale = 1.0 ;
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.98 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.momentum = 0.9 ;
  parms.desired_rms_height = -1.0 ;
  parms.base_name[0] = 0 ;
  parms.nbhd_size = 7 ;    /* out to 7-connected neighbors */
  parms.max_nbrs = 12 ;    /* 12 at each distance */
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    print_help() ;

  parms.base_dt = base_dt_scale * parms.dt ;
  in_patch_fname = argv[1] ;
  out_patch_fname = argv[2] ;
  FileNamePath(in_patch_fname, path) ;
  cp = strrchr(in_patch_fname, '/') ;
  if (!cp)
    cp = in_patch_fname ;
  cp = strchr(cp, '.') ;
  if (cp)
  {
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
  }
  else
    strcpy(hemi, "lh") ;
  if (one_surf_flag)
    sprintf(in_surf_fname, "%s", in_patch_fname) ;
  else
    sprintf(in_surf_fname, "%s/%s.%s", path, hemi, original_surf_name) ;

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

  if (sphere_flag)
  {
    MRIScenter(mris, mris) ;
    mris->radius = MRISaverageRadius(mris) ;
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  }

  if (Gdiag_no >= 0)
  {
    int n ;
    printf("vertex %d has %d nbrs before patch:\n",
           Gdiag_no, mris->vertices[Gdiag_no].vnum) ;
    for (n = 0 ; n < mris->vertices[Gdiag_no].vnum ; n++)
      printf("\t%d\n", mris->vertices[Gdiag_no].v[n]) ;
  }
  if (one_surf_flag)  /* only have the 1 surface - no patch file */
  {
    mris->patch = 1 ;
    mris->status = MRIS_PATCH ;
    if (!FEQUAL(rescale,1))
    {
      MRISscaleBrain(mris, mris, rescale) ;
      MRIScomputeMetricProperties(mris) ;
    }
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

  } 
  else
  {
    MRISresetNeighborhoodSize(mris, mris->vertices[0].nsize) ; // set back to max
    if (label_fname) // read in a label instead of a patch
    {
      LABEL *area ;
      area = LabelRead(NULL, label_fname) ;
      if (area == NULL)
        ErrorExit(ERROR_BADPARM, "%s: could not read label file %s",
                  Progname, label_fname) ;

      LabelDilate(area, mris, dilate_label) ;
      MRISclearMarks(mris) ;
      LabelMark(area, mris) ;
      MRISripUnmarked(mris) ;
      MRISripFaces(mris);
      mris->patch = 1 ;
      mris->status = MRIS_CUT ;
      LabelFree(&area) ;
      printf("%d valid vertices (%2.1f %% of total)\n",
             MRISvalidVertices(mris), 
             100.0*MRISvalidVertices(mris)/mris->nvertices) ;
    }
    else
    {
      if (MRISreadPatch(mris, in_patch_fname) != NO_ERROR)
        ErrorExit(ERROR_BADPARM, "%s: could not read patch file %s",
                  Progname, in_patch_fname) ;
      if (dilate)
      {
        printf("dilating patch %d times\n", dilate) ;
        MRISdilateRipped(mris, dilate) ;
        printf("%d valid vertices (%2.1f %% of total)\n",
               MRISvalidVertices(mris), 100.0*MRISvalidVertices(mris)/mris->nvertices) ;
      }
    }
    MRISremoveRipped(mris) ;
    MRISupdateSurface(mris) ;
#if 0
    mris->nsize = 1 ; // before recalculation of 2 and 3-nbrs
    {
      int vno ;
      VERTEX *v ;
      for (vno= 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        v->vtotal = v->vnum ;
        v->nsize = 1 ;
      }
    }
    MRISsetNeighborhoodSize(mris, nbrs) ;
#endif
  }

  if (Gdiag_no >= 0)
    printf("vno %d is %sin patch\n", Gdiag_no,
           mris->vertices[Gdiag_no].ripflag ? "NOT " : "") ;

  if (Gdiag_no >= 0 && mris->vertices[Gdiag_no].ripflag == 0)
  {
    int n ;
    printf("vertex %d has %d nbrs after patch:\n",
           Gdiag_no, mris->vertices[Gdiag_no].vnum) ;
    for (n = 0 ; n < mris->vertices[Gdiag_no].vnum ; n++)
      printf("\t%d\n", mris->vertices[Gdiag_no].v[n]) ;
  }
  fprintf(stderr, "reading original vertex positions...\n") ;
  if (!FZERO(disturb))
    mrisDisturbVertices(mris, disturb) ;
  if (parms.niterations > 0)
  {
    MRISresetNeighborhoodSize(mris, nbrs) ;

    if (!FZERO(parms.l_unfold) || !FZERO(parms.l_expand))
    {
      static INTEGRATION_PARMS p2 ;
      sprintf(in_surf_fname, "%s/%s.%s", path, hemi, original_surf_name) ;
      if (stricmp(original_unfold_surf_name,"none") == 0)
      {
        printf("using current position of patch as initial position\n") ;
        MRISstoreMetricProperties(mris) ;  /* use current positions */
      }
      else if (!sphere_flag && !one_surf_flag)
        MRISreadOriginalProperties(mris, original_unfold_surf_name) ;
      *(&p2) = *(&parms) ;
      p2.l_dist = 0 ;
      p2.niterations = 100 ;
      p2.nbhd_size = p2.max_nbrs = 1 ;
      p2.n_averages = 0 ;
      p2.write_iterations = parms.write_iterations > 0 ? 25 : 0 ;
      p2.tol = -1 ;
      p2.dt = 0.5 ;
      p2.l_area = 0.0 ;
      p2.l_spring = 0.9 ;
      p2.l_convex = 0.9 ;
      p2.momentum = 0 ;
      p2.integration_type = INTEGRATE_MOMENTUM ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
#if 0
      p2.flags |= IPFLAG_NO_SELF_INT_TEST ;
      printf("expanding surface....\n") ;
      MRISexpandSurface(mris, 4.0, &p2) ;  // push it away from fissure
#endif
      p2.niterations = 100 ;
      MRISunfold(mris, &p2, 1) ;
      p2.niterations = 300 ;
      p2.l_unfold *= 0.25 ;
      MRISunfold(mris, &p2, 1) ;
      p2.l_unfold *= 0.25 ;
      MRISunfold(mris, &p2, 1) ;
#if 0
      printf("smoothing unfolded surface..\n");
      p2.niterations = 200 ;
      p2.l_unfold = 0 ;  // just smooth it
      MRISunfold(mris, &p2, max_passes) ;
#endif
      parms.start_t = p2.start_t ;
      parms.l_unfold = parms.l_convex = parms.l_boundary = parms.l_expand=0 ;
      MRIfree(&parms.mri_dist) ;
    }

    sprintf(in_surf_fname, "%s/%s.%s", path, hemi, original_surf_name) ;
    if (!sphere_flag && !one_surf_flag)
      MRISreadOriginalProperties(mris, original_surf_name) ;
    if (randomly_flatten)
      MRISflattenPatchRandomly(mris) ;
    else
      MRISflattenPatch(mris) ;

    /* optimize metric properties of flat map */
    fprintf(stderr,"minimizing metric distortion induced by projection...\n");
    MRISscaleBrain(mris, mris, scale) ;
    MRIScomputeMetricProperties(mris) ;
    MRISunfold(mris, &parms, max_passes) ;
    MRIScenter(mris, mris) ;
    fprintf(stderr, "writing flattened patch to %s\n", out_patch_fname) ;
    MRISwritePatch(mris, out_patch_fname) ;
  }

  if (plane_flag || sphere_flag)
  {
    char fname[STRLEN] ;
    FILE *fp ;

#if 0
    sprintf(fname, "%s.%s.out",
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms.base_name);
#else
    sprintf(fname, "flatten.log") ;
#endif
    fp = fopen(fname, "a") ;

    if (plane_flag)
      MRIScomputeAnalyticDistanceError(mris, MRIS_PLANE, fp) ;
    else if (sphere_flag)
      MRIScomputeAnalyticDistanceError(mris, MRIS_SPHERE, fp) ;
    fclose(fp) ;
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
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "norand"))
  {
    setRandomSeed(0L) ;
  }
  else if (!stricmp(option, "sphere"))
  {
    sphere_flag = 1 ;
  }
  else if (!stricmp(option, "plane"))
  {
    plane_flag = 1 ;
  }
  else if (!stricmp(option, "rescale"))
  {
    rescale = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "rescaling brain by %2.3f\n", rescale) ;
  }
  else if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dilating cuts %d times\n", dilate) ;
  }
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  }
  else if (!stricmp(option, "expand"))
  {
    sscanf(argv[2], "%f", &parms.l_expand) ;
    nargs = 1 ;
    printf("setting l_expand = %2.3f\n", parms.l_expand) ;
  }
  else if (!stricmp(option, "unfold"))
  {
    MRI *mri_kernel, *mri_tmp ;

    sscanf(argv[2], "%f", &parms.l_unfold) ;
    mri_tmp = MRIread(argv[3]) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read distance map %s...\n",
                argv[3]) ;

    mri_kernel = MRIgaussian1d(1.0, -1) ;
    parms.mri_dist = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
    MRIfree(&mri_kernel) ;
    MRIfree(&mri_tmp) ;
    nargs = 2 ;
    fprintf(stderr, "l_unfold = %2.3f\n", parms.l_unfold) ;
  }
  else if (!stricmp(option, "boundary"))
  {
    sscanf(argv[2], "%f", &parms.l_boundary) ;
    nargs = 1 ;
    fprintf(stderr, "l_boundary = %2.3f\n", parms.l_boundary) ;
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
  {
    nospring = 1 ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "nlarea"))
  {
    sscanf(argv[2], "%f", &parms.l_nlarea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_nlarea = %2.3f\n", parms.l_nlarea) ;
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
  else if (!stricmp(option, "complete"))
  {
    parms.complete_dist_mat = 1 ;
    fprintf(stderr, "using complete distance matrix\n") ;
  }
  else if (!stricmp(option, "vnum") || (!stricmp(option, "distances")))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "sampling %d neighbors out to a distance of %d mm\n",
            parms.max_nbrs, parms.nbhd_size) ;
  }
  else if (!stricmp(option, "dt_dec"))
  {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else if (!stricmp(option, "ou"))
  {
    original_unfold_surf_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"reading original unfolding surface from %s...\n",
            original_unfold_surf_name);
  }
  else if (!stricmp(option, "as"))
  {
    parms.area_coef_scale = atof(argv[2]) ;
    printf("setting area coef scale to %2.3f\n",
           parms.area_coef_scale) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'P':
      max_passes = atoi(argv[2]) ;
      fprintf(stderr, "limitting unfolding to %d passes\n", max_passes) ;
      nargs = 1 ;
      break ;
    case '1':
      one_surf_flag = 1 ;  /* patch is only surface file */
      break ;
    case 'O':
      original_surf_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr,"reading original surface from %s...\n",
              original_surf_name);
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
    case 'L':
      label_fname = argv[2] ;
      dilate_label = atof(argv[3]) ;
      nargs = 2;
      printf("loading label %s and dilating it %d times before flattening\n",
             label_fname, dilate_label) ;
      break ;
    case 'D':
      disturb = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "disturbing vertex positions by %2.3f\n",
              (float)disturb) ;
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
      fprintf(stderr, "using write iterations = %d\n",
              parms.write_iterations) ;
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
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
      break ;
    default:
    case 'H':
    case '?':
    case 'U':
      print_help() ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void)
{
  fprintf(stderr,
          "Usage: %s [options] <input patch> <output patch>\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\nThis program will flatten a surface patch\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, " -w <# iterations>\n\t"
          "write out the surface every # of iterations.\n") ;
  fprintf(stderr,
          " -distances <nbhd size> <# of vertices at each distance>\n\t"
          "specify size of neighborhood and number of vertices at each\n\t"
          "distance to be used in the optimization.\n") ;
  fprintf(stderr,
          " -dilate <# of dilations>\n\t"
          "specify the number of times to dilate the ripped edges to ensure a clean cut\n") ;
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

  max_ratio = 0.0f ;
  max_v = max_n = 0 ;
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
        max_v = vno ;
        max_n = n ;
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
