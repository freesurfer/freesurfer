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
static void translate_annotation(MRI_SURFACE *mris, char *trans_name) ;


const char *Progname ;

static char subjects_dir[STRLEN] = "" ;

int
main(int argc, char *argv[]) {
  char         *cp, fname[STRLEN], **av, *subject, *hemi, *in_annot, *trans_name, *out_annot ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;

  nargs = handleVersionOption(argc, argv, "mris_translate_annotation");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag = DIAG_SHOW ;
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

  if (argc < 6)
    usage_exit() ;

  subject = argv[1] ;
  hemi = argv[2] ;
  in_annot = argv[3] ;
  trans_name = argv[4] ;
  out_annot = argv[5] ;

  if (strlen(subjects_dir) == 0) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment or cmd line", Progname) ;
    strcpy(subjects_dir, cp) ;
  }
  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.orig", subjects_dir, subject, hemi) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading surface from %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
  if (MRISreadAnnotation(mris, in_annot) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annot file %s", Progname, in_annot) ;

  translate_annotation(mris, trans_name) ;
  printf("writing translated annotation to %s...\n", out_annot) ;
  MRISwriteAnnotation(mris, out_annot) ;

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

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
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
  fprintf(stderr,
          "usage: %s [options] <subject> <hemi> <in annot> <translation file> <out annot>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;


  fprintf(stderr,
          "\nThis program applies a translation table to an annotation file") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static void
translate_annotation(MRI_SURFACE *mris, char *trans_name) {
  FILE   *fp ;
  int    vno, in_annot, out_annot, rin, gin, bin, rout, gout, bout ;
  char   *cp, line[STRLEN] ;
  VERTEX *v ;

  fp = fopen(trans_name, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not read translation file %s...\n", trans_name) ;

  while ((cp = fgetl(line, STRLEN-1, fp)) != NULL) {
    sscanf(cp, "%d %d %d %d %d %d", &rin, &gin, &bin, &rout, &gout, &bout) ;
    in_annot = rin +   (gin  << 8) + (bin  << 16) ;
    out_annot = rout + (gout << 8) + (bout << 16) ;

    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->annotation == in_annot)
        v->annotation = out_annot ;
    }
  }
  fclose(fp) ;
}

