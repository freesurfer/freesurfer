/*
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

#include "version.h"
#include "macros.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"



static float sigmas[] = {
                          4.0f, 2.0f, 1.0f, 0.5f
                        } ;

const char *surface_names[] = {
                          "hippocampus"
                        } ;

static const char *curvature_names[] = {
                                   "hippocampus.curv",
                                 } ;



int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  compute_area_ratios(MRI_SURFACE *mris) ;

static int  mrisRegister(MRI_SURFACE *mris, MRI_SP *mrisp_template, INTEGRATION_PARMS *parms, int max_passes) ;
static int  mrisLogIntegrationParms2(FILE *fp, MRI_SURFACE *mris,INTEGRATION_PARMS *parms);
static int  mrisClearMomentum(MRI_SURFACE *mris);
static int  mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,int base_averages);

static int which_norm = NORM_MEAN;
static int max_passes = 4 ;
static int nbrs = 1 ;
static float scale = 1.0f ;

static int reverse_flag = 0 ;

static float dalpha = 0.0f ;
static float dbeta = 0.0f ;
static float dgamma = 0.0f ;


const char *Progname ;
static const char curvature_fname[STRLEN] = "" ;
static const char *orig_name = "hippocampus" ;
static char *jacobian_fname = NULL ;

static int use_defaults = 1 ;


#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])
#define NSIGMAS  (sizeof(sigmas)  / sizeof(sigmas[0]))
#define MAX_NBHD_SIZE  200


static INTEGRATION_PARMS  parms ;

int
main(int argc, char *argv[]) {
  char         **av, *surf_fname, *template_fname, *out_fname, fname[STRLEN],*cp;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp_template ;

  Gdiag = DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.projection = PROJECT_SPHERE ;
  parms.flags |= IP_USE_CURVATURE ;
  parms.tol = 1e-0*10 ;
  parms.min_averages = 0 ;
  parms.l_area = 0.0 ;
  parms.l_parea = 0.2f ;
  parms.l_dist = 0.1 ;
  parms.l_corr = 1.0f ;
  parms.l_nlarea = 1 ;
  parms.l_pcorr = 0.0f ;
  parms.niterations = 25 ;
  parms.n_averages = 256 ;
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
  parms.nbhd_size = 0 ;
  parms.max_nbrs = 0 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  surf_fname = argv[1] ;
  template_fname = argv[2] ;
  out_fname = argv[3] ;
  if (parms.base_name[0] == 0) {
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

  if (!FZERO(dalpha) || !FZERO(dbeta) || !FZERO(dgamma))
    MRISrotate(mris, mris, RADIANS(dalpha), RADIANS(dbeta),
               RADIANS(dgamma)) ;

  if (curvature_fname[0]) {
    fprintf(stderr, "reading source curvature from %s\n",curvature_fname) ;
    MRISreadCurvatureFile(mris, curvature_fname) ;
  }
  fprintf(stderr, "reading template parameterization from %s...\n",
          template_fname) ;
  mrisp_template = MRISPread(template_fname) ;
  if (!mrisp_template)
    ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
              Progname, template_fname) ;
  if (use_defaults) {
    if (*IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2) <= 1.0)  /* 1st time */
    {
      parms.l_dist = 0.5 ;
      parms.l_corr = 1.0 ;
      parms.l_parea = 0.1 ;
    } else   /* subsequent alignments */
    {
      parms.l_dist = 0.1 ;
      parms.l_corr = 1.0 ;
      parms.l_parea = 0.2 ;
    }
  }

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  if (reverse_flag)
    MRISreverse(mris, REVERSE_X, 1) ;
  mris->status = MRIS_PARAMETERIZED_SPHERE ;
  MRIScomputeMetricProperties(mris) ;
  if (!FZERO(parms.l_dist))
    MRISscaleDistances(mris, scale) ;
#if 0
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISzeroNegativeAreas(mris) ;
  MRISstoreMetricProperties(mris) ;
#endif
  MRISstoreMeanCurvature(mris) ;  /* use curvature from file */
  /*  MRISsetOriginalFileName(mris, orig_name) ;*/
  MRISreadOriginalProperties(mris, orig_name) ;
  mrisRegister(mris, mrisp_template, &parms, max_passes) ;
  fprintf(stderr, "writing registered surface to %s...\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  if (jacobian_fname) {
    MRIScomputeMetricProperties(mris) ;
    compute_area_ratios(mris) ;  /* will put results in v->curv */
#if 0
    MRISwriteArea(mris, jacobian_fname) ;
#else
    MRISwriteCurvature(mris, jacobian_fname) ;
#endif
  }

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
get_option(int argc, char *argv[]) {
  int    nargs = 0 ;
  char   *option ;
  float  f ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "vnum") || !stricmp(option, "distances")) {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  } else if (!stricmp(option, "rotate")) {
    dalpha = atof(argv[2]) ;
    dbeta = atof(argv[3]) ;
    dgamma = atof(argv[4]) ;
    fprintf(stderr, "rotating brain by (%2.2f, %2.2f, %2.2f)\n",
            dalpha, dbeta, dgamma) ;
    nargs = 3 ;
  } else if (!stricmp(option, "reverse")) {
    reverse_flag = 1 ;
    fprintf(stderr, "mirror image reversing brain before morphing...\n") ;
  } else if (!stricmp(option, "jacobian")) {
    jacobian_fname = argv[2] ;
    nargs = 1 ;
    printf("writing out jacobian of mapping to %s\n", jacobian_fname) ;
  } else if (!stricmp(option, "dist")) {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    use_defaults = 0 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  } else if (!stricmp(option, "norot")) {
    fprintf(stderr, "disabling initial rigid alignment...\n") ;
    parms.flags |= IP_NO_RIGID_ALIGN ;
  } else if (!stricmp(option, "lm")) {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  } else if (!stricmp(option, "search")) {
    parms.integration_type = INTEGRATE_LM_SEARCH ;
    fprintf(stderr, "integrating with binary search line minimization\n") ;
  } else if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = .2*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
  } else if (!stricmp(option, "area")) {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  } else if (!stricmp(option, "parea")) {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_parea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_parea = %2.3f\n", parms.l_parea) ;
  } else if (!stricmp(option, "nlarea")) {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_nlarea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_nlarea = %2.3f\n", parms.l_nlarea) ;
  } else if (!stricmp(option, "spring")) {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  } else if (!stricmp(option, "corr")) {
    use_defaults = 0 ;
    sscanf(argv[2], "%f", &parms.l_corr) ;
    nargs = 1 ;
    fprintf(stderr, "using l_corr = %2.3f\n", parms.l_corr) ;
  } else if (!stricmp(option, "curv")) {
    parms.flags |= IP_USE_CURVATURE ;
    fprintf(stderr, "using smoothwm curvature for final alignment\n") ;
  } else if (!stricmp(option, "nocurv")) {
    parms.flags &= ~IP_USE_CURVATURE ;
    fprintf(stderr, "using smoothwm curvature for final alignment\n") ;
  } else if (!stricmp(option, "adaptive")) {
    parms.integration_type = INTEGRATE_ADAPTIVE ;
    fprintf(stderr, "using adaptive time step integration\n") ;
  } else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  } else if (!stricmp(option, "tol")) {
    if (sscanf(argv[2], "%e", &f) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan tol from %s",
                Progname, argv[2]) ;
    parms.tol = (double)f ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.2e\n", (float)parms.tol) ;
  } else if (!stricmp(option, "error_ratio")) {
    parms.error_ratio = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "error_ratio=%2.3f\n", parms.error_ratio) ;
  } else if (!stricmp(option, "dt_inc")) {
    parms.dt_increase = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_increase=%2.3f\n", parms.dt_increase) ;
  } else if (!stricmp(option, "vnum")) {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "nbr size = %d, max neighbors = %d\n",
            parms.nbhd_size, parms.max_nbrs) ;
  } else if (!stricmp(option, "dt_dec")) {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  } else switch (toupper(*option)) {
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
      break ;
    case 'C':
      strncpy(const_cast<char*>(curvature_fname), argv[2], STRLEN) ; // Well... at least it's strncpy
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
    case 'O':
      orig_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "using %s for original properties...\n", orig_name) ;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <average surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program register a surface with  an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
compute_area_ratios(MRI_SURFACE *mris) {
  VERTEX  *v ;
  int     vno ;
  float   area_scale ;

  area_scale = mris->total_area / mris->orig_area  ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    v->curv = v->area / (v->origarea*area_scale) ;
  }

  return(NO_ERROR) ;
}

static int
mrisRegister(MRI_SURFACE *mris, MRI_SP *mrisp_template,
             INTEGRATION_PARMS *parms, int max_passes) {
  float   sigma ;
  int     i, /*steps,*/ done, sno, ino, msec ;
  MRI_SP  *mrisp ;
  char    fname[STRLEN], base_name[STRLEN], path[STRLEN] ;
  double  base_dt ;
  static  int first = 1 ;

  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;
  Timer start;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  FileNamePath(mris->fname, path) ;
  sprintf(base_name, "%s/%s.%s", path,
          mris->hemisphere == LEFT_HEMISPHERE ? "lh":"rh", parms->base_name);

  base_dt = parms->dt ;
  if (Gdiag & DIAG_WRITE) {
    sprintf(fname, "%s.%s.out",
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",parms->base_name);
    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w") ;
      if (!parms->fp)
        ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                  Progname, fname) ;
    }
    mrisLogIntegrationParms2(parms->fp, mris,parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms2(stderr, mris,parms) ;

  MRISuseMeanCurvature(mris) ;
  MRISnormalizeCurvature(mris, which_norm) ;
  MRISstoreMeanCurvature(mris) ;

  if (parms->nbhd_size > 0)  /* compute long-range distances */
  {
    int i, nbrs[MAX_NBHD_SIZE] ;
    for (i = mris->nsize+1 ; i <= parms->nbhd_size ; i++)
      nbrs[i] = parms->max_nbrs ;
  }

  for (sno = 0 ; sno < SURFACES ; sno++) {
    if (!first && ((parms->flags & IP_USE_CURVATURE) == 0))
      break ;

    ino = parms->frame_no = sno*IMAGES_PER_SURFACE ;
    if (curvature_names[sno])  /* read in precomputed curvature file */
    {
      sprintf(fname, "%s.%s",
              mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
              curvature_names[sno]) ;
      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read curvature file '%s'\n",
                  "mrisRegister", fname) ;
      MRISnormalizeCurvature(mris, which_norm) ;
    } else                       /* compute curvature of surface */
    {
      sprintf(fname, "%s", surface_names[sno]) ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  "mrisRegister", fname) ;

      MRISsetNeighborhoodSizeAndDist(mris, -1) ;  /* back to max */
      MRIScomputeMetricProperties(mris) ;
      MRIScomputeSecondFundamentalForm(mris) ;
      MRISuseMeanCurvature(mris) ;
      MRISnormalizeCurvature(mris, which_norm) ;
      MRISresetNeighborhoodSize(mris,1);/*only use nearest neighbor distances*/
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    }
    MRISstoreMeanCurvature(mris) ;

    if (Gdiag & DIAG_SHOW) {
      if (curvature_names[sno])
        fprintf(stdout, "reading precomputed curvature from %s\n",fname) ;
      else
        fprintf(stdout, "calculating curvature of %s surface\n",fname) ;
    }

    if (Gdiag & DIAG_WRITE)
      fprintf(parms->fp,"calculating curvature of %s surface\n",fname);

    if (!first && parms->flags & IP_USE_CURVATURE) {
      /* only small adjustments needed after 1st time around */
      parms->tol *= 2.0f ;
      parms->l_corr /= 20.0f ;  /* should be more adaptive */
      if (Gdiag & DIAG_WRITE)
        mrisLogIntegrationParms2(parms->fp, mris, parms) ;
      if (Gdiag & DIAG_SHOW)
        mrisLogIntegrationParms2(stderr, mris, parms) ;
    } else
      if (!first) /* don't do curvature alignment */
        break ;   /* finished */

    for (i = 0 ; i < NSIGMAS ; i++)  /* for each spatial scale (blurring) */
    {
      parms->sigma = sigma = sigmas[i] ;
      parms->dt = base_dt ;
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "\nblurring surfaces with sigma=%2.2f...", sigma) ;
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,"\ncorrelating surfaces with with sigma=%2.2f\n",
                sigma) ;
      if (Gdiag & DIAG_WRITE && !i && !parms->start_t) {
        MRISfromParameterization(mrisp_template, mris, ino);
        sprintf(fname, "%s/%s.target", path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh") ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }
      MRISuseMeanCurvature(mris) ;
      mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
      parms->mrisp = MRISPblur(mrisp, NULL, sigma, 0) ;
      parms->mrisp_template = MRISPblur(mrisp_template, NULL, sigma, ino) ;
      MRISPblur(parms->mrisp_template, NULL, sigma, ino+1) ; /* variances */
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "done.\n") ;
      /* normalize curvature intensities for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino);
      MRISnormalizeCurvature(mris,which_norm) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino) ;

#if 0
      /* normalize variances for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino+1);
      MRISnormalizeCurvature(mris, which_norm) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino+1) ;
#endif

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        sprintf(fname, "%s/%s.%4.4dtarget%2.2f",
                path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->start_t, sigma) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }

      MRISfromParameterization(parms->mrisp, mris, 0);
      MRISnormalizeCurvature(mris, which_norm) ;
      MRIStoParameterization(mris, parms->mrisp, 1, 0) ;
      MRISPfree(&mrisp) ;

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        MRISPwrite(parms->mrisp, "mrisp_blur.hipl") ;
        MRISPwrite(parms->mrisp_template, "mrisp_template_blur.hipl") ;
      }
      mris->vp = (void *)parms->mrisp ;  /* hack to get it to projectSurface */

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        sprintf(fname, "%s/%s.%4.4dblur%2.2f",
                path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->start_t, sigma) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing curvature file %s...", fname) ;
        MRISwriteCurvature(mris, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
        sprintf(fname, "target.%s.%4.4d.hipl",parms->base_name,parms->start_t);
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "writing parameterization file %s...", fname) ;
        MRISPwrite(parms->mrisp_template, fname) ;
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout, "done.\n") ;
      }
      if (first)  /* only do rigid alignment first time through */
      {
        first = 0 ;
        if ((parms->flags & IP_NO_RIGID_ALIGN) == 0) {
          if (Gdiag & DIAG_SHOW)
            fprintf(stdout, "finding optimal rigid alignment\n") ;
          if (Gdiag & DIAG_WRITE)
            fprintf(parms->fp, "finding optimal rigid alignment\n") ;
          MRISrigidBodyAlignGlobal(mris, parms, 0.5f, 32.0f, 8) ;
          if (Gdiag & DIAG_WRITE && parms->write_iterations != 0)
            MRISwrite(mris, "rotated") ;
        }
      }

      mrisClearMomentum(mris) ;
      done = 0 ;
      mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
    }
  }

  parms->tol /= 10 ;  /* remove everything possible pretty much */
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "removing remaining folds...\n") ;
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "removing remaining folds...\n") ;
#if 1
  parms->l_nlarea *= 5 ;
  mrisIntegrationEpoch(mris, parms, parms->n_averages) ;
#else
  parms->l_nlarea = 1 ;
  parms->l_corr /= 10.0 ;
  parms->l_area = parms->l_parea = 0 ;
  mrisRemoveNegativeArea(mris,parms,parms->n_averages,MAX_NEG_AREA_PCT,3);
#endif
  MRISPfree(&parms->mrisp) ;
  MRISPfree(&parms->mrisp_template) ;
  msec = start.milliseconds() ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "registration took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "registration took %2.2f hours\n",
            (float)msec/(1000.0f*60.0f*60.0f));
  return(NO_ERROR) ;
}
static int
mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,int base_averages) {
  int   total_steps, done, steps, n_averages, old_averages ;
  const char* snum;
  const char* sdenom ;
  float ratio, *pdenom, *pnum ;

  if (!FZERO(parms->l_corr)) {
    sdenom = "corr" ;
    pdenom = &parms->l_corr  ;
  } else {
    sdenom = "dist" ;
    pdenom = &parms->l_dist  ;
  }

  if (!FZERO(parms->l_area)) {
    snum = "area" ;
    pnum = &parms->l_area ;
  } else if (!FZERO(parms->l_parea)) {
    snum = "parea" ;
    pnum = &parms->l_parea  ;
  } else if (!FZERO(parms->l_nlarea)) {
    snum = "nlarea" ;
    pnum = &parms->l_nlarea  ;
  } else {
    snum = "spring" ;
    pnum = &parms->l_spring  ;
  }

  if (Gdiag & DIAG_SHOW)
    mrisLogIntegrationParms2(stderr, mris, parms) ;
  if (Gdiag & DIAG_WRITE)
    mrisLogIntegrationParms2(parms->fp, mris, parms) ;
  if (!FZERO(*pdenom)) {
    ratio = *pnum / *pdenom ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN] ;
      if (!parms->fp) {
        sprintf(fname, "%s.%s.out",
                mris->hemisphere == RIGHT_HEMISPHERE ? "rh":"lh",
                parms->base_name);
        if (!parms->start_t)
          INTEGRATION_PARMS_openFp(parms, fname, "w") ;
        else
          INTEGRATION_PARMS_openFp(parms, fname, "a") ;
        if (!parms->fp)
          ErrorExit(ERROR_NOFILE, "%s: could not open log file %s",
                    Progname, fname) ;
      }
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio) ;
    }
  }

  old_averages = parms->n_averages ;
  for (done = total_steps = 0, n_averages = base_averages ; !done ;
       n_averages /= 4) {
    parms->n_averages = n_averages ;
    steps = MRISintegrate(mris, parms, n_averages) ;
    if (n_averages > 0 && parms->flags & IP_RETRY_INTEGRATION &&
        ((parms->integration_type == INTEGRATE_LINE_MINIMIZE) ||
         (parms->integration_type == INTEGRATE_LM_SEARCH))) {
      int niter = parms->niterations ;
      int integration_type = parms->integration_type ;

      fprintf(stdout, "taking momentum steps...\n") ;
      parms->integration_type = INTEGRATE_MOMENTUM ;
      parms->niterations = 10 ;
      parms->start_t += steps ;
      total_steps += steps ;
      steps = MRISintegrate(mris, parms, n_averages) ;
      parms->integration_type = integration_type ;
      parms->niterations = niter ;
      parms->start_t += steps ;
      total_steps += steps ;
      steps = MRISintegrate(mris, parms, n_averages) ;
    }
    parms->start_t += steps ;
    total_steps += steps ;
    done = n_averages == parms->min_averages ;
    if (mris->status == MRIS_SPHERE) {
      if (Gdiag & DIAG_SHOW)
        MRISprintTessellationStats(mris, stderr) ;
      parms->scale *= parms->dt_decrease ;
      if (parms->scale < 1.0f)
        parms->scale = 1.0f ;
    }
  }
#if 0
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
#endif
  parms->n_averages = old_averages  ;  /* hack, but no time to clean up now */
  return(total_steps) ;
}


static int
mrisClearMomentum(MRI_SURFACE *mris) {
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->odx = 0 ;
    v->ody = 0 ;
    v->odz = 0 ;
  }
  return(NO_ERROR) ;
}

static int
mrisLogIntegrationParms2(FILE *fp, MRI_SURFACE *mris,INTEGRATION_PARMS *parms) {
  char  *cp, host_name[STRLEN] ;

  if (!fp)
    return(NO_ERROR) ;

  cp = getenv("HOST") ;
  if (cp)
    strcpy(host_name, cp) ;
  else
    strcpy(host_name, "unknown") ;

  fprintf(fp, "tol=%2.1e, host=%5.5s, nav=%d, nbrs=%d",
          (float)parms->tol, host_name, parms->n_averages, mris->nsize) ;
  if (!FZERO(parms->l_area))
    fprintf(fp, ", l_area=%2.3f", parms->l_area) ;
  if (!FZERO(parms->l_external))
    fprintf(fp, ", l_extern=%2.3f", parms->l_external) ;
  if (!FZERO(parms->l_parea))
    fprintf(fp, ", l_parea=%2.3f", parms->l_parea) ;
  if (!FZERO(parms->l_nlarea))
    fprintf(fp, ", l_nlarea=%2.3f", parms->l_nlarea) ;
  if (!FZERO(parms->l_nldist))
    fprintf(fp, ", l_nldist=%2.3f", parms->l_nldist) ;
  if (!FZERO(parms->l_angle))
    fprintf(fp, ", l_angle=%2.3f", parms->l_angle) ;
  if (!FZERO(parms->l_repulse))
    fprintf(fp, ", l_repulse=%2.3f", parms->l_repulse) ;
  if (!FZERO(parms->l_repulse_ratio))
    fprintf(fp, ", l_repulse_ratio=%2.3f", parms->l_repulse_ratio) ;
  if (!FZERO(parms->l_surf_repulse))
    fprintf(fp, ", l_surf_repulse=%2.3f", parms->l_surf_repulse) ;
  if (!FZERO(parms->l_corr))
    fprintf(fp, ", l_corr=%2.3f", parms->l_corr) ;
  if (!FZERO(parms->l_spring))
    fprintf(fp, ", l_spring=%2.3f", parms->l_spring) ;
  if (!FZERO(parms->l_spring_norm))
    fprintf(fp, ", l_spring_norm=%2.3f", parms->l_spring_norm) ;
  if (!FZERO(parms->l_tspring))
    fprintf(fp, ", l_tspring=%2.3f", parms->l_tspring) ;
  if (!FZERO(parms->l_nspring))
    fprintf(fp, ", l_nspring=%2.3f", parms->l_nspring) ;
  if (!FZERO(parms->l_dist))
    fprintf(fp, ", l_dist=%2.3f", parms->l_dist) ;
  if (!FZERO(parms->l_intensity))
    fprintf(fp, ", l_intensity=%2.3f", parms->l_intensity) ;
  if (!FZERO(parms->l_grad))
    fprintf(fp, ", l_grad=%2.3f", parms->l_grad) ;
  if (!FZERO(parms->l_sphere))
    fprintf(fp, ", l_sphere=%2.3f", parms->l_sphere) ;
  if (!FZERO(parms->l_expand))
    fprintf(fp, ", l_expand=%2.3f", parms->l_expand) ;
  if (!FZERO(parms->l_curv))
    fprintf(fp, ", l_curv=%2.3f", parms->l_curv) ;
  if (!FZERO(parms->l_convex))
    fprintf(fp, ", l_convex=%2.3f", parms->l_convex) ;
  if (!FZERO(parms->l_boundary))
    fprintf(fp, ", l_boundary=%2.3f", parms->l_boundary) ;
  if (!FZERO(parms->l_neg))
    fprintf(fp, ", l_neg=%2.3f", parms->l_neg) ;
  if (!FZERO(parms->l_tsmooth))
    fprintf(fp, ", l_tsmooth=%2.3f", parms->l_tsmooth) ;
  fprintf(fp, "\n") ;
  switch (parms->integration_type) {
  case INTEGRATE_LM_SEARCH:
    fprintf(fp, "using binary search line minimization\n") ;
    break ;
  case INTEGRATE_LINE_MINIMIZE:
    fprintf(fp, "using quadratic fit line minimization\n") ;
    break ;
  case INTEGRATE_ADAPTIVE:
    fprintf(fp,
            "mom=%2.2f, dt=%2.2f, base_dt=%2.3f, dt_inc=%2.2f, "
            "dt_dec=%2.2f, err_rat=%2.2f\n",
            (float)parms->momentum, (float)parms->dt,
            (float)parms->base_dt, (float)parms->dt_increase,
            (float)parms->dt_decrease, (float)parms->error_ratio) ;
    break ;
  default:
  case INTEGRATE_MOMENTUM:
    fprintf(fp,
            "mom=%2.2f, dt=%2.2f\n",(float)parms->momentum, (float)parms->dt);
    break ;
  }
#if 0
  fprintf(fp, "nbhd_size=%d, max_nbrs=%d ", parms->nbhd_size,parms->max_nbrs);
#endif
  if (parms->desired_rms_height > 0.0)
    fprintf(fp, "desired rms height=%2.3f", parms->desired_rms_height) ;
  fprintf(fp, "\n") ;
  fflush(fp) ;
  return(NO_ERROR) ;
}
