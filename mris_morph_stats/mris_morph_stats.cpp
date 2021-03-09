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
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "mrishash.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int write_vals = 0 ;

static int nbrs = 1 ;

int
main(int argc, char *argv[]) {
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500],
  *morph_name, *out_name ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;

  nargs = handleVersionOption(argc, argv, "mris_morph_stats");
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

  sname = argv[1] ;
  hemi = argv[2] ;
  morph_name = argv[3] ;
  out_name = argv[4] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
  strcpy(sdir, cp) ;

  /* read in white matter surface and use it to compute original properties */
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, WHITE_MATTER_NAME) ;
  fprintf(stderr, "reading original surface position from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */
  MRISstoreMetricProperties(mris) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ; /* save white-matter */

  fprintf(stderr, "reading morphed surface position from %s...\n",morph_name);
  MRISreadVertexPositions(mris, morph_name) ;
  MRIScomputeMetricProperties(mris) ;    /* recompute surface normals */

  fprintf(stderr, "writing morph areal errors to %s\n",out_name);
  if (write_vals) {
#if 0
    char *cp ;
    cp = strrchr(out_name, '.') ;
    if (!cp || *(cp+1) != 'w') {
      sprintf(fname, "%s.w", out_name) ;
      out_name = fname ;
    }
#endif
    MRISwriteAreaErrorToValFile(mris, out_name) ;
  } else
    MRISwriteAreaError(mris, out_name) ;

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
  else if (!stricmp(option, "nbrs")) {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  } else if (!stricmp(option, "write_vals")) {
    write_vals = 1 ;
    fprintf(stderr,  "writing morph errors to w file\n") ;
  } else if (!stricmp(option, "name")) {
#if 0
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "base name = %s\n", parms.base_name) ;
#endif
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr, "usage: %s [options] <subject name> <hemisphere> <morphed surface> <output name>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program generates statistics which characterize\n"
          "a surface-based deformation field\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


