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
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static char *log_fname = NULL ;
static char *brain_fname = NULL ;
static void usage_exit(int code) ;

static int in_label = -1 ;
static int out_label = -1 ;

static int compute_pct = 0 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, msec, minutes, label, volume, seconds, i, brain_volume ;
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

  if (in_label >= 0)
    MRIreplaceValues(mri, mri, in_label, out_label) ;

  if (compute_pct)
  {
    for (brain_volume = label = 0 ; label <= MAX_CMA_LABEL ; label++)
    {
      if (!IS_BRAIN(label))
        continue ;
      brain_volume += MRIvoxelsInLabel(mri, label) ;
    }
  }
  else if (brain_fname)
  {
    MRI *mri_brain ;

    mri_brain = MRIread(brain_fname) ;
    if (!mri_brain)
      ErrorExit(ERROR_NOFILE, "%s: could not read brain volume from %s",
                Progname, brain_fname) ;
    MRIbinarize(mri_brain, mri_brain, 1, 0, 255) ;
    brain_volume = MRIvoxelsInLabel(mri_brain, 255) ;
    MRIfree(&mri_brain) ;
    compute_pct = 1 ; /* use this brain volume to normalize */
  }
  else
    brain_volume = 1 ;

  for (i = 2 ; i < argc ; i++)
  {
    label = atoi(argv[i]) ;
    printf("processing label %d...\n", label) ;


    volume = MRIvoxelsInLabel(mri, label) ;
    if (log_fname)
    {
      char fname[STRLEN] ;

      sprintf(fname, log_fname, label) ;
      printf("logging to %s...\n", fname) ;
      log_fp = fopen(fname, "a+") ;
      if (!log_fp)
        ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                  Progname, fname) ;
    }
    else
      log_fp = NULL ;

    if (compute_pct)
    {
      printf("%d voxels in label %d, %%%2.6f of brain volume (%d)\n", 
             volume, label, 100.0*(float)volume/(float)brain_volume,
             brain_volume) ;
      if (log_fp)
      {
        fprintf(log_fp,"%2.6f\n", 100.0*(float)volume/(float)brain_volume) ;
        fclose(log_fp) ;
      }
    }
    else
    {
      printf("%d voxels in label %d\n", volume, label) ;
      if (log_fp)
      {
        fprintf(log_fp,"%d\n", volume) ;
        fclose(log_fp) ;
      }
    }
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
  case 'T':
    in_label = atoi(argv[2]) ;
    out_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("translating label %d to label %d\n", in_label, out_label) ;
    break ;
  case 'B':
    brain_fname = argv[2] ;
    nargs = 1 ;
    break ;
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    /*    fprintf(stderr, "logging results to %s\n", log_fname) ;*/
    break ;
  case 'P':
    compute_pct = 1 ;
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
  printf("usage: %s [options] <volume> <label 1> <label 2> ...\n", Progname) ;
  printf(
       "\tp             - compute brain volume as a pct of all brain labels\n"
       "\tl <log fname> - log results to file (note %%d will include label #)\n"
       "\tb <brain vol> - load brain vol and use it to normalize volumes\n"
         );
  exit(code) ;
}
