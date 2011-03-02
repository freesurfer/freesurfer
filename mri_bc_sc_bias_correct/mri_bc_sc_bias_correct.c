/**
 * @file  mri_bc_sc_bias_correct.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
 *    $Revision: 1.5 $
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
#include <string.h>
#include <ctype.h>

#include "mri.h"
#include "diag.h"
#include "version.h"
#include "error.h"
#include "mrisegment.h"
#include "mrinorm.h"
#include "proto.h"
#include "macros.h"

static char vcid[] =
  "$Id: mri_bc_sc_bias_correct.c,v 1.5 2011/03/02 00:04:13 nicks Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static float sigma = 10 ;

char *Progname ;

int
main(int argc, char *argv[]) {
  char             **av, *out_fname ;
  int              ac, nargs, b, segno ;
  MRI              *mri_bc, *mri_sc, *mri_bias, *mri_in, *mri_mask, *mri_ctrl, *mri_kernel,
  *mri_smooth_bias, *mri_tmp, *mri_out  ;
  float            thresh ;
  HISTOGRAM        *histo ;
  MRI_SEGMENTATION *mriseg ;

  /* rkt: check for and handle version tag */
  Progname = argv[0] ;
  nargs = handle_version_option (argc, argv, "$Id: mri_bc_sc_bias_correct.c,v 1.5 2011/03/02 00:04:13 nicks Exp $", "$Name:  $");
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
  mri_in = MRIread(argv[1]) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read input image from %s", argv[1]) ;
  mri_tmp = MRIchangeType(mri_in, MRI_FLOAT, 0, 1, 1);
  MRIfree(&mri_in);
  mri_in = mri_tmp;

  mri_bc = MRIread(argv[2]) ;
  if (mri_bc == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read body coil PD map from %s", argv[2]) ;
  mri_sc = MRIread(argv[3]) ;
  if (mri_sc == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface coil PD map from %s", argv[3]) ;

  histo = MRIhistogram(mri_sc, 256) ;

  out_fname = argv[4] ;
  for (b = 255 ; b >= 0 ; b--)
    if (histo->counts[b] > 0)
      break ;
  thresh = histo->bins[b]*.25 ;
  printf("largest bin found at %d = %2.1f, setting thresh = %2.1f\n",
         b, histo->bins[b], thresh) ;

  mriseg = MRIsegment(mri_sc, thresh, histo->bins[b]) ;
  segno = MRIfindMaxSegmentNumber(mriseg) ;
  printf("max seg %d found with %d voxels\n", segno, mriseg->segments[segno].nvoxels) ;
  mri_mask = MRIsegmentToImage(mri_sc, NULL, mriseg, segno) ;
  MRIsegmentFree(&mriseg) ;
  MRIerode(mri_mask, mri_mask) ;
  mri_ctrl = MRIbinarize(mri_mask, NULL, 1, 0, CONTROL_MARKED) ;
  mri_tmp = MRIchangeType(mri_ctrl, MRI_UCHAR, 0, 1, 1);
  MRIfree(&mri_ctrl);
  mri_ctrl = mri_tmp;
  MRImask(mri_bc, mri_mask, mri_bc, 0, 0) ;
  MRImask(mri_sc, mri_mask, mri_sc, 0, 0) ;
  mri_tmp = MRIchangeType(mri_bc, MRI_FLOAT, 0, 1, 1) ;
  MRIfree(&mri_bc) ;
  mri_bc = mri_tmp ;
  mri_tmp = MRIchangeType(mri_sc, MRI_FLOAT, 0, 1, 1) ;
  MRIfree(&mri_sc) ;
  mri_sc = mri_tmp ;
  mri_bias = MRIdivide(mri_bc, mri_sc, NULL) ;
  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_mask, "mask.mgz") ;
    MRIwrite(mri_ctrl, "ctrl.mgz") ;
    MRIwrite(mri_bc, "bc.mgz") ;
    MRIwrite(mri_sc, "sc.mgz") ;
    MRIwrite(mri_bias, "bias1.mgz") ;
    MRIwrite(mri_bias, "bias_v.mgz") ;
  }

  MRIfree(&mri_ctrl) ;
  MRIfree(&mri_mask) ;
  MRIfree(&mri_bc) ;
  MRIfree(&mri_sc) ;

  mri_kernel = MRIgaussian1d(sigma, 100) ;
  mri_smooth_bias = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_smooth_bias, "bias.mgz") ;

  mri_out = MRIapplyBiasCorrection(mri_in, mri_smooth_bias, NULL) ;

  MRIfree(&mri_kernel) ;
  MRIfree(&mri_bias) ;
  MRIfree(&mri_smooth_bias) ;
  printf("writing output to %s.\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;

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
    case 'S':
      printf("using sigma = %2.2f (default=%2.2f)\n", atof(argv[2]), sigma) ;
      sigma = atof(argv[2]) ;
      nargs = 1 ;
      break ;
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
      print_usage();
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
  printf("usage: %s [options] <input volume> <BC PD map> <SC PD map> <output volume>\n",
         Progname) ;
}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program will intensity correct a volume based on surface coil and body coil images\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
