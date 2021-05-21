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
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "image.h"
#include "density.h"
#include "mrisegment.h"
#include "proto.h"
#include "macros.h"

const char *Progname ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int main(int argc, char *argv[]) ;

static int nbins = 256 ;
static int valid1_min = -1 ;
static int valid1_max = -1 ;
static int valid2_min = -1 ;
static int valid2_max = -1 ;
static float sigma = 2.0 ;

static int nlevels = 1 ;

#define MAX_LEVELS 10

int
main(int argc, char *argv[]) {
  char        **av ;
  int         ac, nargs, i, level ;
  MRI         *mri1, *mri2, *mri_tmp ;
  DENSITY     *pdf ;
  float       min_val1, min_val2, max_val1, max_val2 ;
  int *valid1 = NULL ;
  int *valid2 = NULL ;
  MRI_SEGMENTATION *mriseg ;
  char fname[STRLEN] ;
  MRI         *mri1_pyramid[MAX_LEVELS], *mri2_pyramid[MAX_LEVELS];

  nargs = handleVersionOption(argc, argv, "histo_compute_joint_density");
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

  mri1 = MRIread(argv[1]) ;
  if (!mri1)
    ErrorExit(ERROR_BADPARM, "%s: could not read volume 1 from %s\n", argv[1]) ;
  mri2 = MRIread(argv[2]) ;
  if (!mri2)
    ErrorExit(ERROR_BADPARM, "%s: could not read volume 2 from %s\n", argv[2]) ;

  MRIreplaceValues(mri1, mri1, 0, 1) ;  /* reserve 0 for background */
  MRIreplaceValues(mri2, mri2, 0, 1) ;  /* reserve 0 for background */
  if (valid1_min >= 0) {
    MRIvalRange(mri1, &min_val1, &max_val1) ;
    valid1_min = MAX(valid1_min, min_val1) ;
    valid1_max = MIN(valid1_max, max_val1) ;
    valid1 = (int *)calloc(256, sizeof(int)) ;
    for (i = valid1_min ; i <= valid1_max ; i++)
      valid1[i] = 1 ;

    mriseg = MRImaxsegment(mri1, valid1_min, valid1_max) ;
    mri_tmp = MRIsegmentToImage(mri1, NULL, mriseg, 0) ;
    MRIfree(&mri1) ;
    mri1 = mri_tmp ;
    mri_tmp = MRIsegmentToImage(mri2, NULL, mriseg, 0) ;
    MRIfree(&mri2) ;
    mri2 = mri_tmp ;
    MRIsegmentFree(&mriseg) ;
  } else {
    valid1 = (int *)calloc(256, sizeof(int)) ;
    for (i = 0 ; i <= 255 ; i++)
      valid1[i] = 1 ;
  }

  if (valid2_min >= 0) {
    MRIvalRange(mri2, &min_val2, &max_val2) ;
    valid2_min = MAX(valid2_min, min_val2) ;
    valid2_max = MIN(valid2_max, max_val2) ;
    valid2 = (int *)calloc(256, sizeof(int)) ;
    for (i = valid2_min ; i <= valid2_max ; i++)
      valid2[i] = 1 ;
  } else {
    valid2 = (int *)calloc(256, sizeof(int)) ;
    for (i = 0 ; i <= 255 ; i++)
      valid2[i] = 1 ;
  }

  mri1_pyramid[0] = mri1 ;
  mri2_pyramid[0] = mri2 ;
  for (level = 1 ; level < MAX_LEVELS ; level++) {
    mri1_pyramid[level] = MRIreduce2D(mri1_pyramid[level-1], NULL) ;
    mri2_pyramid[level] = MRIreduce2D(mri1_pyramid[level-1], NULL) ;
    sprintf(fname, "h%d.mgz", level) ;
    MRIwrite(mri1_pyramid[level], fname) ;
    sprintf(fname, "b%d.mgz", level) ;
    MRIwrite(mri2_pyramid[level], fname) ;
  }

  for (level = nlevels-1 ; level >= 0 ; level--) {
    pdf = DensityHistogramEstimate(mri1_pyramid[level], mri2_pyramid[level], nbins, sigma, valid1, valid2) ;
    sprintf(fname, "%s_level%d.den", argv[3], level) ;
    printf("writing estimated density to %s\n", fname);
    DensityWrite(pdf, fname) ;
    /*  DensityFree(&pdf) ;*/
  }

  exit(0) ;
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
  else if (!stricmp(option, "valid1")) {
    valid1_min = atoi(argv[2]) ;
    valid1_max = atoi(argv[3]) ;
    nargs = 2 ;
    printf("limiting histogram to intensity in range [%d, %d]\n", valid1_min, valid1_max) ;
  } else if (!stricmp(option, "valid2")) {
    valid2_min = atoi(argv[2]) ;
    valid2_max = atoi(argv[3]) ;
    nargs = 2 ;
    printf("limiting histogram to intensity in range [%d, %d]\n", valid2_min, valid2_max) ;
  } else if (!stricmp(option, "nlevels")) {
    nlevels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("computing %d levels of gaussian pyramid\n", nlevels) ;
  } else switch (toupper(*option)) {
    case 'S':
      sigma = atof(argv[2]) ;
      nargs = 1 ;
      printf("using Gaussian smoothing of histogram with sigma=%2.2f\n", sigma) ;
      break ;
    case 'N':
      nbins = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting # of bins in joint histogram to %d\n", nbins) ;
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
  //  print_usage() ; // print_help _calls print_usage
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <volume1> <volume2> <joint density file>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;

  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

