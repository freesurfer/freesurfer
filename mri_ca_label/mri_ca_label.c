#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gca.h"
#include "transform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;
MRI *MRImodeFilter(MRI *mri_src, MRI *mri_dst, int niter) ;

static int filter = 0 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  char   *in_fname, *out_fname,  *gca_fname, *xform_fname ;
  MRI    *mri_in, *mri_labeled ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  GCA     *gca ;
  LTA     *lta ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;
  in_fname = argv[1] ;
  xform_fname = argv[2];
  gca_fname = argv[3] ;
  out_fname = argv[4] ;

  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
              Progname, in_fname) ;

  /*  fprintf(stderr, "mri_in read: xform %s\n", mri_in->transform_fname) ;*/
  printf("reading classifier array from %s...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier array from %s",
              Progname, gca_fname) ;

  lta = LTAread(xform_fname) ;
  if (!lta)
    ErrorExit(ERROR_NOFILE, "%s: could not open transform", xform_fname) ;

  printf("labeling volume...\n") ;
  mri_labeled = GCAlabel(mri_in, gca, NULL, lta) ;
  GCAfree(&gca) ; MRIfree(&mri_in) ;
  if (filter)
  {
    MRI *mri_tmp ;

    mri_tmp = MRImodeFilter(mri_labeled, NULL, filter) ;
    MRIfree(&mri_labeled) ;
    mri_labeled = mri_tmp ;
  }

  MRIwrite(mri_labeled, out_fname) ;

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "auto-labeling took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  switch (toupper(*option))
  {
  case 'F':
    filter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("applying mode filter %d times to output of labelling\n",filter);
    break ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <input volume> <xform> <gca file>"
         " <output volume>\n", Progname) ;
  exit(code) ;
}
