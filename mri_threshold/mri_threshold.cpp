/**
 * @brief Function to lower (or upper) threshold the input volume.
 *
 * This is a function that allows for intensity thresholding the input volume. By default the threshold is a lower threshold, but with the -u flag this can be changed to be an upper one.
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
#include <ctype.h>

#include "mri.h"
#include "diag.h"
#include "tags.h"
#include "version.h"
#include "error.h"
#include "proto.h"
#include "macros.h"


int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int binarize = 0 ;
const char *Progname ;
int specificframe = 0;
int frame = 0;
int upperthreshold = 0;

int
main(int argc, char *argv[]) {
  char        **av, *in_fname, *out_fname ;
  int         ac, nargs ;
  MRI         *mri_in, *mri_out ;
  float       thresh ;


  std::string cmdline = getAllInfo(argc, argv, "mri_threshold");
  Progname = argv[0] ;
  nargs = handleVersionOption(argc, argv, "mri_threshold");
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
  in_fname = argv[1] ;
  thresh = atof(argv[2]) ;
  out_fname = argv[3] ;

  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname,
              in_fname) ;

  if (! upperthreshold) // Keep everything above the threshold value
    {
      printf("lower thresholding volume at %2.4f....\n", thresh) ;
      //mri_out = MRIthreshold(mri_in, NULL, thresh) ;
      if (specificframe)
	mri_out = MRIthresholdFrame(mri_in, NULL, thresh, frame) ;
      else // if value in any frame < thresh ==> reject in all frames
	mri_out = MRIthresholdAllFrames(mri_in, NULL, thresh) ;
    }
  else // Keep everything below the given threshold value
    {
      printf("upper thresholding volume at %2.4f....\n", thresh) ;
      if (specificframe)
	mri_out = MRIupperthresholdFrame(mri_in, NULL, thresh, frame) ;
      else // if value in any frame < thresh ==> reject in all frames
	mri_out = MRIupperthresholdAllFrames(mri_in, NULL, thresh) ;
    }

  if (binarize > 0)
    {
      if (upperthreshold == 0)
	MRIbinarize(mri_out, mri_out, thresh, 0, binarize) ;
      else
	MRIbinarizeNoThreshold(mri_out, mri_out) ;
    }
  
  printf("writing output to %s.\n", out_fname) ;
  MRIaddCommandLine(mri_out, cmdline) ;
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
    case 'B':
      binarize = atoi(argv[2]) ;
      printf("binarizing output volume to be [0, %d]\n", binarize) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'F':
      frame = atoi(argv[2]) ;
      specificframe = 1;
      nargs = 1 ;
      break ;
    case 'U':
      upperthreshold = 1;
      break ;
    case '?':
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
  printf("usage: %s [options] <in vol> <thresh> <out_vol>\n",
         Progname) ;
}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program will lower threshold an input volume\n") ;
  printf("\n") ;
  printf("Options:\n\n") ;
  printf("  -B bval: It will binarize the output volume to be bval.\n");
  printf("  -U Instead of lower thresholding the volume it will upper threshold it.\n");
  printf("  -F fnum: Apply thresholding to a specific frame indexed by fnum.\n");
  printf("  -? : print usage\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
