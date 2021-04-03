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

const char *Progname ;

static int patch_flag = 0 ;
static int which = REVERSE_X ;

int
main(int argc, char *argv[]) {
  char               **av, *in_fname, *out_fname, path[STRLEN], fname[STRLEN],
  hemi[STRLEN], *cp ;
  int                ac, nargs ;
  MRI_SURFACE        *mris ;

  nargs = handleVersionOption(argc, argv, "mris_reverse");
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

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (patch_flag) {
    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, hemi) ;
    cp = strchr(hemi, '.') ;
    if (cp)
      *cp = 0 ;
    else
      ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from %s\n",
                in_fname) ;
    int req = snprintf(fname, STRLEN, "%s/%s.%s", path, hemi, ORIG_NAME) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    if (MRISreadPatch(mris, in_fname) != NO_ERROR)
      ErrorExit(Gerror, "%s: could not read patch\n", Progname) ;
  } else {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }

  FileNamePath(out_fname, path) ;
  MRISreverse(mris, which, 1) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing reversed surface to %s\n", out_fname) ;
  mris->type = MRIS_TRIANGULAR_SURFACE ;
  if (patch_flag)
    MRISwritePatch(mris, out_fname) ;
  else
    MRISwrite(mris, out_fname) ;
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
    case 'P':
      patch_flag = 1 ;
      break ;
    case 'Y':
      which = REVERSE_Y ;
      break ;
    case 'Z':
      which = REVERSE_Z ;
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
          "usage: %s [options] <input surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr, "\nThis reverses a cortical surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

