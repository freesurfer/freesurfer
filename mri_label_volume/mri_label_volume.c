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
static char *log_fname = NULL ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, msec, minutes, label, volume, seconds ;
  struct timeb start ;
  MRI    *mri ;
  FILE   *log_fp ;

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


  mri = MRIread(argv[1]) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[1]) ;
  label = atoi(argv[2]) ;

  volume = MRIvoxelsInLabel(mri, label) ;
  if (log_fname)
  {
    log_fp = fopen(log_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
  }
  else
    log_fp = NULL ;

  printf("%d voxels in label\n", volume) ;
  if (log_fp)
  {
    fprintf(log_fp,"%d\n", volume) ;
    fclose(log_fp) ;
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
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "logging results to %s\n", log_fname) ;
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
  printf("usage: %s [options] <volume 1> <volume 2>",
         Progname) ;
  printf(
         "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
         );
  exit(code) ;
}
