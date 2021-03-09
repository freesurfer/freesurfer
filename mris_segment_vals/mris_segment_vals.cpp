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
#include "label.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int thresh = 0 ;
static int area_thresh = 0 ;

int
main(int argc, char *argv[]) {
  char               **av, *in_fname, *out_fname, *surf_fname ;
  int                ac, nargs, vno, nlabels, lno ;
  MRI_SURFACE        *mris ;
  VERTEX             *v ;
  LABEL              **label_array ;

  nargs = handleVersionOption(argc, argv, "mris_segment_vals");
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
    print_help() ;

  surf_fname = argv[1] ;
  in_fname = argv[2] ;
  out_fname = argv[3] ;

  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if (MRISreadValues(mris, in_fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read val file %s",
              Progname, in_fname) ;

  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || fabs(v->val)<=thresh)
      continue ;
    v->marked = 1 ;
  }

  MRISsegmentMarked(mris, &label_array, &nlabels, area_thresh) ;
  printf("%d labels found...\n", nlabels) ;
  MRISsetVals(mris, 0) ;
  for (lno = 0 ; lno < nlabels ; lno++) {
    for (vno = 0 ; vno < label_array[lno]->n_points ; vno++) {
      v = &mris->vertices[label_array[lno]->lv[vno].vno] ;
      v->curv = lno+1 ;
    }
  }

  fprintf(stderr, "writing segmented vals to %s\n", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;
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
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      thresh = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'A':
      area_thresh = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <input w/curv file> <output w/curv file>\n",
          Progname) ;
  fprintf(stderr, "where options are:\n\t-T <threshold> (default is 0)\n"
          "\t-A <area thresh> ignore segments smaller than <area thresh> mm(default 0)\n");
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program segments an input val file into connected components\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

