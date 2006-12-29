/**
 * @file  mris_curvature.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2006/12/29 20:32:17 $
 *    $Revision: 1.21 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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

static char vcid[] = "$Id: mris_curvature.c,v 1.21 2006/12/29 20:32:17 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int write_flag = 0 ;
static int nbrs = 2 ;
static int navgs = 0 ;
static char *param_file = NULL ;
static int normalize = 0 ;
static int diff_flag = 0 ;
static int max_flag = 0 ;
static int min_flag = 0 ;
static int stretch_flag = 0 ;
static int patch_flag = 0 ;
static int neg_flag = 0 ;
static int param_no = 0 ;
static int normalize_param = 0 ;
static int ratio_flag = 0 ;
static int contrast_flag = 0 ;

#define MAX_NBHD_SIZE 100
static int nbhd_size = 0 ;
static int nbrs_per_distance = 0 ;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname,fname[STRLEN],hemi[10], path[STRLEN],
  name[STRLEN],*cp ;
  int          ac, nargs, nhandles ;
  MRI_SURFACE  *mris ;
  double       ici, fi, var ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_curvature.c,v 1.21 2006/12/29 20:32:17 fischl Exp $", "$Name:  $");
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

  if (argc < 2)
    usage_exit() ;

  in_fname = argv[1] ;

  FileNamePath(in_fname, path) ;
  FileNameOnly(in_fname, name) ;
  cp = strchr(name, '.') ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from '%s'",
              Progname, fname) ;
  strncpy(hemi, cp-2, 2) ;
  hemi[2] = 0 ;

  if (patch_flag)  /* read the orig surface, then the patch file */
  {
    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISfastRead(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "reading patch file %s...\n", in_fname) ;
    if (MRISreadPatch(mris, in_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;

  } else   /* just read the surface normally */
  {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }

  MRISsetNeighborhoodSize(mris, nbrs) ;

  if (nbhd_size > 0)
    MRISsampleAtEachDistance(mris, nbhd_size, nbrs_per_distance) ;

  if (param_file) {
    MRI_SP *mrisp ;
    mrisp = MRISPread(param_file) ;
    if (normalize_param)
      MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
    else
      MRISfromParameterization(mrisp, mris, param_no) ;
    MRISPfree(&mrisp) ;
    if (normalize)
      MRISnormalizeCurvature(mris) ;
    sprintf(fname, "%s/%s.param", path,name) ;
    fprintf(stderr, "writing parameterized curvature to %s...", fname) ;
    MRISwriteCurvature(mris, fname) ;
    fprintf(stderr, "done.\n") ;
  } else {
    MRIScomputeSecondFundamentalForm(mris) ;
    nhandles = nint(1.0 - mris->Ktotal / (4.0*M_PI)) ;
    fprintf(stderr, "total integrated curvature = %2.3f*4pi (%2.3f) --> "
            "%d handles\n", (float)(mris->Ktotal/(4.0f*M_PI)),
            (float)mris->Ktotal, nhandles) ;

#if 0
    fprintf(stderr, "0: k1 = %2.3f, k2 = %2.3f, H = %2.3f, K = %2.3f\n",
            mris->vertices[0].k1, mris->vertices[0].k2,
            mris->vertices[0].H, mris->vertices[0].K) ;
    fprintf(stderr, "0: vnum = %d, v2num = %d, total=%d, area=%2.3f\n",
            mris->vertices[0].vnum, mris->vertices[0].v2num,
            mris->vertices[0].vtotal,mris->vertices[0].area) ;
#endif
    MRIScomputeCurvatureIndices(mris, &ici, &fi);
    var = MRIStotalVariation(mris) ;
    fprintf(stderr,"ICI = %2.1f, FI = %2.1f, variation=%2.3f\n", ici, fi, var);

    if (diff_flag) {
      MRISuseCurvatureDifference(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s/%s.diff", path,name) ;
      fprintf(stderr, "writing curvature difference to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (ratio_flag) {
      MRISuseCurvatureRatio(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.ratio", path,name) ;
      fprintf(stderr, "writing curvature ratio to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (contrast_flag) {
      MRISuseCurvatureContrast(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.contrast", path,name) ;
      fprintf(stderr, "writing curvature contrast to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (neg_flag) {
      int neg ;
      if (mris->patch)
        mris->status = MRIS_PLANE ;
      MRIScomputeMetricProperties(mris) ;
      neg = MRIScountNegativeTriangles(mris) ;
      MRISuseNegCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s/%s.neg", path,name) ;
      fprintf(stderr, "writing negative vertex curvature to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "%d negative triangles\n", neg) ;
      fprintf(stderr, "done.\n") ;
      {
        int    vno, fno ;
        VERTEX *v ;
        FACE   *f ;
        for (vno = 0 ; vno < mris->nvertices ; vno++) {
          v = &mris->vertices[vno] ;
          if (v->ripflag)
            continue ;
          neg = 0 ;
          for (fno = 0 ; fno < v->num ; fno++) {
            f = &mris->faces[v->f[fno]] ;
            if (f->area < 0.0f)
              neg = 1 ;
          }
          if (neg)
            fprintf(stdout, "%d\n", vno) ;
        }
      }
    }

    if (max_flag) {
      MRISuseCurvatureMax(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.max", path,name) ;
      fprintf(stderr, "writing curvature maxima to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (min_flag) {
      MRISuseCurvatureMin(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.min", path,name) ;
      fprintf(stderr, "writing curvature minima to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (stretch_flag) {
      MRISreadOriginalProperties(mris, NULL) ;
      MRISuseCurvatureStretch(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.stretch", path,name) ;
      fprintf(stderr, "writing curvature stretch to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (write_flag) {
      MRISuseGaussianCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      sprintf(fname, "%s/%s.K", path,name) ;
      fprintf(stderr, "writing Gaussian curvature to %s...", fname) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      MRISwriteCurvature(mris, fname) ;
      MRISuseMeanCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
        MRISnormalizeCurvature(mris) ;
      sprintf(fname, "%s/%s.H", path,name) ;
      fprintf(stderr, "done.\nwriting mean curvature to %s...", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
  }
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
  else if (!stricmp(option, "distances") || !stricmp(option, "vnum")) {
    nbhd_size = atoi(argv[2]) ;
    nbrs_per_distance = atoi(argv[3]) ;
    fprintf(stderr, "sampling %d neighbors out to a distance of %d mm\n",
            nbrs_per_distance, nbhd_size) ;
    nargs = 2 ;
  } else if (!stricmp(option, "diff"))
    diff_flag = 1 ;
  else if (!stricmp(option, "ratio") || !stricmp(option, "defect"))
    ratio_flag = 1 ;
  else if (!stricmp(option, "contrast"))
    contrast_flag = 1 ;
  else if (!stricmp(option, "neg"))
    neg_flag = 1 ;
  else if (!stricmp(option, "max"))
    max_flag = 1 ;
  else if (!stricmp(option, "min"))
    min_flag = 1 ;
  else if (!stricmp(option, "stretch"))
    stretch_flag = 1 ;
  else if (!stricmp(option, "param")) {
    param_file = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using parameterization file %s\n", param_file) ;
  } else if (!stricmp(option, "nparam")) {
    char *cp ;
    cp = strchr(argv[2], '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    } else
      param_no = 0 ;
    normalize_param = 1 ;
  } else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  } else switch (toupper(*option)) {
    case 'N':
      fprintf(stderr, "normalizing curvature values.\n") ;
      normalize = 1 ;
      break ;
    case 'P':
      patch_flag = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      write_flag = 1 ;
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
  fprintf(stderr, "usage: %s [options] <input surface file>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will the second fundamental form of a cortical surface."
          "\nIt will create two new files <hemi>.H and <hemi>.K with the\n"
          "mean and Gaussian curvature respectively.");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

