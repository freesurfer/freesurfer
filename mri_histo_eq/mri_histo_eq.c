#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[])
{
  char   **av, *out_fname ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI    *mri_src, *mri_template, *mri_eq ;

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

  mri_eq = MRIhistoEqualize(mri_src, mri_template, NULL, 30, 170) ;
  MRIwrite(mri_eq, out_fname) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
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

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "overlap calculation took %d minutes and %d seconds.\n", 
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
  printf("usage: %s [options] <volume 1> <volume 2>",
         Progname) ;
  printf(
         "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
         );
  exit(code) ;
}
