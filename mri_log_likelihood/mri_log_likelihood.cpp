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


//
// mri_log_likelihood
// written by Bruce Fischl
//
// Warning: Do not edit the following four lines.  CVS maintains them.
//
////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "tags.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"

const char *Progname ;

static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;

static char *orig_fname = NULL ;

/*
   command line consists of three inputs:

   argv[1]  - input volume
   argv[2]  - atlas (gca)
   argv[3]  - transform (lta/xfm/m3d)
*/


int
main(int argc, char *argv[]) {
  char         *gca_fname, *in_fname, **av, *xform_fname ;
  MRI          *mri_in, *mri_tmp, *mri_orig = NULL ;
  GCA          *gca ;
  int          ac, nargs, input, ninputs ;
  TRANSFORM    *transform = NULL ;
  double       ll ;

  std::string cmdline = getAllInfo(argc, argv, "mri_log_likelihood");

  nargs = handleVersionOption(argc, argv, "mri_log_likelihood");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  setRandomSeed(-1L) ;
  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit
    (ERROR_BADPARM,
     "usage: %s [<options>] <inbrain1> <inbrain2> ... "
     "<atlas> <transform file> ...\n",
     Progname) ;

  ninputs = (argc - 1) / 2 ;
  if (DIAG_VERBOSE_ON)
    printf("reading %d input volume%ss\n", ninputs, ninputs > 1 ? "s" : "") ;
  in_fname = argv[1] ;
  gca_fname = argv[1+ninputs] ;
  xform_fname = argv[2+ninputs] ;
  transform = TransformRead(xform_fname) ;
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read input transform from %s",
              Progname, xform_fname) ;

  if (DIAG_VERBOSE_ON)
    printf("reading atlas from '%s'...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read input atlas from %s",
              Progname, gca_fname) ;

  fflush(stdout) ;
  for (input = 0 ; input < ninputs ; input++) {
    in_fname = argv[1+input] ;
    if (DIAG_VERBOSE_ON)
      printf("reading input volume from %s...\n", in_fname) ;
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, in_fname) ;
    MRImakePositive(mri_tmp, mri_tmp) ;
    if (input == 0) {
      mri_in =
        MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                         mri_tmp->type, ninputs) ;
      if (!mri_in)
        ErrorExit(ERROR_NOMEMORY,
                  "%s: could not allocate input volume %dx%dx%dx%d",
                  mri_tmp->width,mri_tmp->height,mri_tmp->depth,ninputs) ;
      MRIcopyHeader(mri_tmp, mri_in) ;
    }

    MRIcopyFrame(mri_tmp, mri_in, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  MRIaddCommandLine(mri_in, cmdline) ;

  TransformInvert(transform, mri_in) ;

  if (orig_fname) {
    mri_orig = MRIread(orig_fname) ;
    if (mri_orig == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read orig volume from %s", Progname, orig_fname) ;
  }
  ll = GCAimageLogLikelihood(gca, mri_in, transform, 1, mri_orig) ;
  printf("%2.0f\n", 10000*ll) ;

  MRIfree(&mri_in) ;

  if (gca)
    GCAfree(&gca) ;
  if (mri_in)
    MRIfree(&mri_in) ;
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
  StrUpper(option) ;
  if (!strcmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "orig")) {
    orig_fname = argv[2] ;
    nargs = 1 ;
  } else switch (*option) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(int code) {
  printf("usage: %s <in volume> <atlas> <transform>\n\n",
         Progname) ;
  exit(code) ;
}
