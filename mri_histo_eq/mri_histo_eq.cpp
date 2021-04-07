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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "mrinorm.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "transform.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static int   adaptive_normalize = 0 ;

static char *xform_fname = NULL ;

int
main(int argc, char *argv[]) {
  char   **av, *out_fname ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI    *mri_src, *mri_template ;

  nargs = handleVersionOption(argc, argv, "mri_histo_eq");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;


  mri_src = MRIread(argv[1]) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[1]) ;
  mri_template = MRIread(argv[2]) ;
  if (!mri_template)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[2]) ;
  out_fname = argv[3] ;

  if (xform_fname) {
    char   path[STRLEN], fname[STRLEN] ;
    LTA    *lta_src, *lta_template ;
    MATRIX *m_L, *m_inv ;
    MRI    *mri_tmp ;

    FileNameOnly(xform_fname, xform_fname) ;

    FileNamePath(argv[1], path) ;
    int req = snprintf(fname, STRLEN, "%s/%s", path, xform_fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    lta_src = LTAread(fname) ;
    if (!lta_src)
      ErrorExit(ERROR_BADPARM, "%s: could not read xform from %s",
                Progname, fname) ;

    /* now put template in source space */
    m_inv = MatrixInverse(lta_src->xforms[0].m_L, NULL) ;
    if (!m_inv)
      ErrorExit(ERROR_BADPARM, "%s: non-invertible transform in %s!",
                Progname, fname) ;

    FileNamePath(argv[2], path) ;
    req = snprintf(fname, STRLEN, "%s/%s", path, xform_fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    lta_template = LTAread(fname) ;
    if (!lta_template)
      ErrorExit(ERROR_BADPARM, "%s: could not read xform from %s",
                Progname, fname) ;

    m_L = MatrixMultiply(m_inv, lta_template->xforms[0].m_L, NULL) ;
    LTAfree(&lta_src) ;
    LTAfree(&lta_template) ;
    MatrixFree(&m_inv) ;
    mri_tmp = MRIlinearTransform(mri_template, NULL, m_L) ;
    MRIfree(&mri_template) ;
    mri_template = mri_tmp ;
    MatrixFree(&m_L);
  }

  MRI *mri_eq=NULL;
  if (adaptive_normalize)
    mri_eq = MRIadaptiveHistoNormalize(mri_src, NULL, mri_template,
                                       8, 32, 30) ;
  else
    mri_eq = MRIhistoNormalize(mri_src, NULL, mri_template, 30, 170) ;
  MRIwrite(mri_eq, out_fname) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOGRAM *histo,*hsmooth ;

    histo = MRIhistogram(mri_src, 0) ;
    HISTOplot(histo, "src.plt") ;
    histo = MRIhistogram(mri_template, 0) ;
    HISTOplot(histo, "template.plt") ;
    histo = MRIhistogram(mri_eq, 0) ;
    HISTOplot(histo, "eq.plt") ;
    hsmooth = HISTOsmooth(histo, NULL, 2) ;
    HISTOplot(hsmooth, "eqs.plt") ;
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  fprintf(stderr, "histogram normalization took %d minutes and %d seconds.\n",
          minutes, seconds) ;

  exit(0) ;
  return(0) ;
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
  case 'A':
    adaptive_normalize = 1 ;
    break ;
  case 'T':
    xform_fname = argv[2] ;
    nargs = 1 ;
    break;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <volume 1> <volume 2>",
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  exit(code) ;
}
