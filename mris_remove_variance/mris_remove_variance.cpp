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

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISremoveValueVarianceFromCurvature(MRI_SURFACE *mris) ;
double MRIScomputeCurvatureValueCorrelationCoefficient(MRI_SURFACE *mris) ;

const char *Progname ;

static int navgs = 0 ;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *in_curv, *var_curv, *out_curv, fname[200],
  hemi[10], path[200], name[200],*cp ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  double       r ;

  nargs = handleVersionOption(argc, argv, "mris_remove_variance");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  in_fname = argv[1] ;
  in_curv = argv[2] ;
  var_curv = argv[3] ;
  out_curv = argv[4] ;

  FileNamePath(in_fname, path) ;
  FileNameOnly(in_fname, name) ;
  cp = strchr(name, '.') ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from '%s'",
              Progname, fname) ;
  strncpy(hemi, cp-2, 2) ;
  hemi[2] = 0 ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading surface file %s...\n", in_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  MRISreadCurvatureFile(mris, var_curv) ;
  MRIScopyCurvatureToValues(mris) ;
  MRISreadCurvatureFile(mris, in_curv) ;
  r = MRIScomputeCurvatureValueCorrelationCoefficient(mris) ;
  fprintf(stderr, "correlation coefficient = %2.3f\n", (float)r) ;
  MRISremoveValueVarianceFromCurvature(mris) ;
  if (Gdiag & DIAG_SHOW) {
    r = MRIScomputeCurvatureValueCorrelationCoefficient(mris) ;
    fprintf(stderr, "after decorrelation coefficient = %2.3f\n", (float)r) ;
  }
  MRISaverageCurvatures(mris, navgs) ;
  MRISwriteCurvature(mris, out_curv) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
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
  fprintf(stderr, "usage: %s [options] <input surface file> <curv. file>\n"
          "\t<curv. file to remove> <output curvature file>.\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will remove the linear component of the variance \n"
          "\taccounted for by one curvature vector from another.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

int
MRISremoveValueVarianceFromCurvature(MRI_SURFACE *mris) {
  int      vno ;
  VERTEX   *v ;
  double   len, dot, mean_curv, mean_val, n, val ;
  float    min_curv, max_curv ;

  /* compute means of two vectors */
  n = (double)mris->nvertices ;
  for (mean_val = mean_curv = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean_val += v->val ;
    mean_curv += v->curv ;
  }

  mean_val /= n ;
  mean_curv /= n ;

  /* first build normalize vector in vertex->val field */

  /* measure length of val vector */
  for (len = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    val = (double)v->val - mean_val ;
    len += val*val ;
  }

  if (FZERO(len))
    return(NO_ERROR) ;

  len = sqrt(len) ;

  /* normalize it and compute dot product with curvature */
  for (dot = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->d = (v->val-mean_val) / len ;
    dot += (double)(v->curv-mean_curv) * (double)v->d ;
  }

  /* now subtract the the scaled val vector from the curvature vector */
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv -= dot * (double)v->d ;
  }

  /* dot should now be 0 */
  for (dot = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dot += (double)(v->curv-mean_curv) * (double)(v->val-mean_val) ;
  }

  /* recompute min, max and mean curvature */
  min_curv = 10000.0f ;
  max_curv = -min_curv ;
  mean_curv = 0.0f ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean_curv += v->curv ;
    if (v->curv < min_curv)
      min_curv = v->curv ;
    if (v->curv > max_curv)
      max_curv = v->curv ;
  }

  mris->min_curv = min_curv ;
  mris->max_curv = max_curv ;
  mean_curv /= (float)mris->nvertices ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,
            "after variance removal curvature %2.1f --> %2.1f, mean %2.1f.\n",
            min_curv, max_curv, mean_curv) ;
  return(NO_ERROR) ;
}

double
MRIScomputeCurvatureValueCorrelationCoefficient(MRI_SURFACE *mris) {
  double   r, mean_val, mean_curv, s_val, s_curv, n ;
  int      vno ;
  VERTEX   *v ;

  /* first build normalized vector in vertex->val field */
  n = (double)mris->nvertices ;
  for (mean_val = mean_curv = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mean_val += v->val ;
    mean_curv += v->curv ;
  }

  mean_val /= n ;
  mean_curv /= n ;

  /* compute standard deviations */
  for (s_val = s_curv = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    s_val += SQR(v->val-mean_val) ;
    s_curv += SQR(v->curv-mean_curv) ;
  }

  s_val /= (n-1.0) ;
  s_curv /= (n-1.0) ;
  s_val = sqrt(s_val) ;
  s_curv = sqrt(s_curv) ;

  /* measure length of val vector */
  for (r = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    r += (v->val - mean_val) * (v->curv - mean_curv) ;
  }
  r /= ((n-1.0) * s_val * s_curv) ;

  return(r) ;
}
