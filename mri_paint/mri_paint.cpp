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
#include "transform.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int imageoffset = 0 ;
static char *thickness_fname = NULL ;

static int coords = TALAIRACH_COORDS ;

/*
usage: stat_paint [options] <input> <subject> <hemisphere> <surface> <output>
*/
int
main(int argc, char *argv[]) {
  char        **av, *in_vol, *in_surf, *xform_fname, *out_fname ;
  int         ac, nargs ;
  MRI_SURFACE *mris ;
  MRI         *mri ;
  LTA         *lta ;
  MATRIX      *m ;

  nargs = handleVersionOption(argc, argv, "mri_paint");
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

  in_vol = argv[1] ;
  in_surf = argv[2] ;
  xform_fname = argv[3] ;
  out_fname = argv[4] ;

  fprintf(stderr, "reading volume from %s and painting it onto %s\n", in_vol, in_surf) ;
  mri = MRIread(in_vol) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, in_vol) ;

  mris = MRISread(in_surf) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, in_surf) ;

  if (thickness_fname) {
    if (MRISreadCurvature(mris, thickness_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read thickness file %s", Progname, thickness_fname) ;
  }

  if ((stricmp(xform_fname, "identity") == 0) ||
      (stricmp(xform_fname, "I") == 0))
    lta = LTAalloc(1, NULL) ;
  else
    lta = LTAread(xform_fname) ;
  if (!lta)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_fname) ;

  m = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  MatrixFree(&lta->xforms[0].m_L) ;
  lta->xforms[0].m_L = m ;
  MRISpaintVolume(mris, lta, mri) ;
  MRIScopyValuesToCurvature(mris) ;    /* write them out in new curv format */
  fprintf(stderr, "writing output to %s.\n", out_fname) ;
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
  else if (!stricmp(option, "imageoffset")) {
    imageoffset = atoi(argv[2]) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'T':
      thickness_fname = argv[2] ;
      nargs = 1 ;
      printf("using thickness file %s to sample with...\n", thickness_fname) ;
      break ;
    case 'S':
    case 'C':
      coords = SPHERICAL_COORDS ;
      fprintf(stderr, "using spherical coordinates\n") ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input volume> <input surface> <registration file> <output .float file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will paint a average Talairach stats onto a surface\n");
  fprintf(stderr, "-imageoffset <image offset> - set offset to use\n") ;
  fprintf(stderr, "-S                          - paint using surface "
          "coordinates\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

