/**
 * @file  label_area.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:03 $
 *    $Revision: 1.7 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void print_usage(void) ;

char *Progname ;

static int verbose = 0 ;


#define NAME_LEN      100

static float compute_label_area(MRI_SURFACE *mris, char *subject_name,
                                char *area_name) ;

static char subjects_dir[NAME_LEN] = "" ;



int
main(int argc, char *argv[]) {
  int          ac, nargs, i ;
  char         **av, *cp, surf_name[100], *hemi, *subject_name, *area_name ;
  MRI_SURFACE  *mris ;
  float        area ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: label_area.c,v 1.7 2006/12/29 02:09:03 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* read in command-line options */
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    print_usage() ;

  hemi = argv[1] ;
  subject_name = argv[2] ;

  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "no subjects directory in environment.\n") ;
  strcpy(subjects_dir, cp) ;
  sprintf(surf_name,"%s/%s/surf/%s.smoothwm",subjects_dir,subject_name,hemi);
  fprintf(stderr, "reading %s...\n", surf_name) ;
  mris = MRISread(surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s\n",
              surf_name) ;
  MRIScomputeMetricProperties(mris) ;

  /*  read in the area names */
  for (i = 3 ; i < argc ; i++) {
    area_name = argv[i] ;
    if (verbose > 3)  /* never */
      fprintf(stderr, "reading area %s\n", area_name) ;
    area = compute_label_area(mris, subject_name, area_name) ;
    if (mris->group_avg_surface_area > 0) {
      MRIScomputeMetricProperties(mris) ;
      printf("adjusting for group surface area %2.2f (%2.0f / %2.0f)\n",
             mris->group_avg_surface_area/mris->total_area, mris->group_avg_surface_area,mris->total_area);
      area *= mris->group_avg_surface_area/mris->total_area ;
    }
    fprintf(stderr, "%s:  %s - %2.3f square mm\n",
            subject_name, area_name, area) ;
  }


  if (verbose)
    fprintf(stderr, "done.\n") ;
  exit(0) ;
  return(0) ;  /* ansi */
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
  switch (toupper(*option)) {
  case 'V':
    verbose = !verbose ;
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
#if 0
static Transform *
load_transform(char *subject_name, General_transform *transform) {
  char xform_fname[100] ;

  sprintf(xform_fname, "%s/%s/mri/transforms/talairach.xfm",
          subjects_dir, subject_name) ;
  if (input_transform_file(xform_fname, transform) != OK)
    ErrorExit(ERROR_NOFILE, "%s: could not load transform file '%s'",
              Progname, xform_fname) ;

  if (verbose == 2)
    fprintf(stderr, "transform read successfully from %s\n", xform_fname) ;
  return(get_linear_transform_ptr(transform)) ;
}
#endif

static float
compute_label_area(MRI_SURFACE *mris, char *subject_name, char *area_name) {
  int      nlines, vno ;
  char     *cp, line[200], fname[200] ;
  FILE     *fp ;
  float    total_area ;
  VERTEX  *v ;

  sprintf(fname, "%s/%s/label/%s.label",subjects_dir,subject_name,area_name);
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
              Progname, fname) ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorExit(ERROR_BADFILE, "%s: empty label file %s", Progname, fname) ;
  if (!sscanf(cp, "%d", &nlines))
    ErrorExit(ERROR_BADFILE, "%s: could not scan # of lines from %s",
              Progname, fname) ;

  nlines = 0 ;
  total_area = 0.0f ;
  while ((cp = fgetl(line, 199, fp)) != NULL) {
    if (sscanf(cp, "%d %*f %*f %*f", &vno) != 1)
      ErrorExit(ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, nlines, fname) ;

    v = &mris->vertices[vno] ;
    total_area += v->area ;
    nlines++ ;
  }

  if (!nlines)
    ErrorExit(ERROR_BADFILE, "%s: no data in label file %s", Progname,fname);

  fclose(fp) ;
  return(total_area) ;
}

static void
print_usage(void) {
  printf("usage: %s <hemi> <subject name> <label file name>...\n",Progname);
  exit(1) ;
}
