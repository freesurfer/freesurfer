/**
 * @file  label_border.c
 * @brief 
 *
 * compute the boundary of a label
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:11 $
 *    $Revision: 1.2 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "const.h"
#include "utils.h"
#include "mrisurf.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static int  unmark_interior(MRI_SURFACE *mris) ;

char *Progname ;

static int verbose = 0 ;
static int print_radius = 0 ;

static char subjects_dir[STRLEN] = "" ;



int
main(int argc, char *argv[]) {
  int          ac, nargs ;
  char         **av, *cp, surf_name[STRLEN], *hemi, *subject_name, *area_name,
               *out_name ;
  MRI_SURFACE  *mris ;
  LABEL        *area ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: label_border.c,v 1.2 2011/03/02 00:04:11 nicks Exp $", "$Name:  $");
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

  if (argc < 5)
    print_usage() ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  area_name = argv[3] ;
  out_name = argv[4] ;

  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "no subjects directory in environment.\n") ;
  strcpy(subjects_dir, cp) ;
  sprintf(surf_name,"%s/%s/surf/%s.white",subjects_dir,subject_name,hemi);
  fprintf(stderr, "reading %s...\n", surf_name) ;
  mris = MRISread(surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s\n",
              surf_name) ;

  /*  read in the area names */
  area = LabelRead(subject_name, area_name) ;
  fprintf(stderr, "reading area %s\n", area_name) ;
  MRISclearMarks(mris) ;
  LabelMarkSurface(area, mris) ;
  unmark_interior(mris) ;
  LabelFree(&area) ;
  area = LabelFromMarkedSurface(mris) ;
  LabelWrite(area, out_name) ;

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
  case 'R':
    print_radius = 1 ;
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

static void
print_usage(void) {
  printf("usage: %s <subject name> <hemi> <input label file> < output label file>\n",Progname);
  exit(1) ;
}
static int
unmark_interior(MRI_SURFACE *mris)
{
  int    vno, n, unmarked ;
  VERTEX *v, *vn ;

  MRIScopyMarkedToMarked2(mris) ;
  for (vno = 0 ; vno <mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == 0)
      continue ;
    for (unmarked = n = 0 ; unmarked == 0 && n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked == 0)
        unmarked = 1 ;
    }
    if (unmarked == 0)
      v->marked2 = 0 ;
  }
  MRIScopyMarked2ToMarked(mris) ;

  if (print_radius)
  {
    double xc, yc ;
    int    num ;

    for (xc = yc = 0.0, num = vno = 0 ; vno <mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked == 0)
        continue ;
      num++ ;
    }
  }

  return(NO_ERROR) ;
}

