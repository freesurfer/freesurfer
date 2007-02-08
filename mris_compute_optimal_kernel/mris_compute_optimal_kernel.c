/**
 * @file  mris_compute_optimal_kernel.c
 * @brief program for computing the optimal blurring kernel between an individual label and a group.
 *
 * computes the isotropic gaussian blurring kernel that is optimal in the lms sense between a
 * group average label and an individual. Outputs the standard deviation of the kernel, which
 * can be used as a measure of the accuracy of a coordinate system.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/02/08 20:59:50 $
 *    $Revision: 1.1 $
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
#include "label.h"

static int MRIScomputeOptimalGaussianKernel(MRI_SURFACE *mris, LABEL *subject_label, 
                                            LABEL *group_label, int step_size, int max_avgs, double *psigma) ;

static char vcid[] = "$Id: mris_compute_optimal_kernel.c,v 1.1 2007/02/08 20:59:50 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char orig_name[STRLEN] = ORIG_NAME ;

static char sdir[STRLEN] = "" ;
static int step_size = 10 ;
static int max_avgs = 1000 ;

int
main(int argc, char *argv[]) {
  char          **av, *subject, *cp, fname[STRLEN], *hemi, *subject_label_name, *group_label_name ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  LABEL         *subject_label, *group_label ;
  double        sigma ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_compute_optimal_kernel.c,v 1.1 2007/02/08 20:59:50 fischl Exp $", "$Name:  $");
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

  if (argc < 5)
    usage_exit() ;

  subject = argv[1] ;
  hemi = argv[2] ;
  subject_label_name = argv[3] ;
  group_label_name = argv[4] ;
  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }


  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, orig_name) ;
  fprintf(stderr, "reading surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  subject_label = LabelRead(subject, subject_label_name) ;
  if (!subject_label)
    ErrorExit(ERROR_NOFILE, "%s: could not read label file %s", Progname, subject_label) ;
  group_label = LabelRead(subject, group_label_name) ;
  if (!group_label)
    ErrorExit(ERROR_NOFILE, "%s: could not read label file %s", Progname, group_label) ;
  

  MRIScomputeOptimalGaussianKernel(mris, subject_label, group_label, step_size, max_avgs, &sigma) ;
  printf("%f\n", sigma) ;
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
  else if (!stricmp(option, "step")) 
  {
    step_size = atoi(argv[2]) ;
    fprintf(stderr,  "using step size %d\n", step_size) ;
    nargs = 1 ;
  } 
  else if (!stricmp(option, "max")) 
  {
    max_avgs = atoi(argv[2]) ;
    fprintf(stderr,  "setting max avgs to %d\n", max_avgs) ;
    nargs = 1 ;
  } 
  else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "orig")) {
    strcpy(orig_name, argv[2]) ;
    fprintf(stderr,  "reading surface from file named %s\n", orig_name) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
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
          "usage: %s [options] <subject name> <hemi> <individual label> <group label>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program measures the thickness of the cortical surface\n"
          "and writes the resulting measurement into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
MRIScomputeOptimalGaussianKernel(MRI_SURFACE *mris, LABEL *subject_label, 
                                 LABEL *group_label, int step_size, int max_avgs, double *psigma)
{
  int          n, vno ;
  double       rms, min_rms, min_sigma ;
  VERTEX       *v ;
  LABEL_VERTEX *lv ;
  
  // write individual subject label into v->val field
  for (n = 0 ; n < subject_label->n_points ; n++)
  {
    lv = &subject_label->lv[n] ;
    if (lv->vno < 0 || lv->vno >= mris->nvertices)
      ErrorExit(ERROR_BADFILE, "%s: group label has invalid vertex %d at index %d",
                Progname, lv->vno, n) ;
    v = &mris->vertices[lv->vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->val = 1 ;
  }

  // write group stats into v->val2 field
  for (n = 0 ; n < group_label->n_points ; n++)
  {
    lv = &group_label->lv[n] ;
    if (lv->vno < 0 || lv->vno >= mris->nvertices)
      ErrorExit(ERROR_BADFILE, "%s: group label has invalid vertex %d at index %d",
                Progname, lv->vno, n) ;
    v = &mris->vertices[lv->vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->val2 = lv->stat ;
  }

  min_sigma = 0 ; min_rms = mris->nvertices*1e6;

  for (n = 0 ; n < max_avgs ; n += step_size)
  {
    for (rms = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      rms += SQR(v->val - v->val2);
    }
    if (rms < min_rms)
    {
      min_rms = rms ;
      min_sigma = n ;
      fprintf(stderr, "new optimum found at avgs = %d\n", n) ;
    }
    MRISaverageVals(mris, step_size) ;
  }

  *psigma = min_sigma ;
  return(NO_ERROR) ;
}


