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
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static char *log_fname = NULL ;
static void usage_exit(int code) ;

static int in_label = -1 ;
static int out_label = -1 ;

static int all_flag = 0 ;
static int compute_pct = 0 ;
static char *brain_fname = NULL ;
static char *icv_fname = NULL ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, msec, minutes, label, volume, seconds, i ;
  struct timeb start ;
  MRI    *mri ;
  FILE   *log_fp ;
	double  vox_volume, brain_volume ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_label_volume.c,v 1.12 2003/04/15 21:07:30 kteich Exp $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  if ((all_flag && argc < 2) || (all_flag == 0 && argc < 3))
    usage_exit(1) ;


  mri = MRIread(argv[1]) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[1]) ;

	vox_volume = mri->xsize * mri->ysize * mri->zsize ;

  if (in_label >= 0)
    MRIreplaceValues(mri, mri, in_label, out_label) ;

	if (all_flag)
	{
		int   nvox ;
		float volume ;
		
		nvox = MRItotalVoxelsOn(mri, WM_MIN_VAL) ;
		volume = nvox * mri->xsize * mri->ysize * mri->zsize ;
		printf("total volume = %d voxels, %2.1f mm^3\n", nvox, volume) ;
		exit(0) ;
	}

  if (brain_fname)
	{
		MRI *mri_brain = MRIread(brain_fname) ;
		if (mri_brain == NULL)
			ErrorExit(ERROR_BADPARM, "%s: could not read brain volume from %s\n", Progname,brain_fname) ;

		brain_volume = (double)MRItotalVoxelsOn(mri_brain, WM_MIN_VAL) ;
		MRIfree(&mri_brain) ;
		brain_volume *= (mri->xsize * mri->ysize * mri->zsize) ;
	}
  else if (compute_pct)
  {
    for (brain_volume = 0.0, label = 0 ; label <= MAX_CMA_LABEL ; label++)
    {
      if (!IS_BRAIN(label))
        continue ;
      brain_volume += (double)MRIvoxelsInLabel(mri, label) ;
    }
		brain_volume *= (mri->xsize * mri->ysize * mri->zsize) ;
  }
	else if (icv_fname)
	{
		FILE *fp ;
		fp = fopen(icv_fname, "r") ;
		if (fp == NULL)
			ErrorExit(ERROR_NOFILE, "%s: could not open ICV file %s\n", Progname, icv_fname) ;
		if (fscanf(fp, "%lf", &brain_volume) != 1)
			ErrorExit(ERROR_NOFILE, "%s: could not read ICV from %s\n", Progname, icv_fname) ;
		fclose(fp) ;
		printf("using intra-cranial volume = %2.1f\n", brain_volume) ;
	}
	else
    brain_volume = 1.0 ;

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

    if (compute_pct || icv_fname)
    {
      printf("%d voxels (%2.1f mm^3) in label %d, %%%2.6f of brain volume (%2.0f)\n", 
             volume, volume*vox_volume,label, 100.0*(float)volume/(float)brain_volume,
             brain_volume) ;
      if (log_fp)
      {
        fprintf(log_fp,"%2.6f\n", 100.0*(float)volume*vox_volume/(float)brain_volume) ;
        fclose(log_fp) ;
      }
    }
    else
    {
      printf("%d (%2.1f mm^3) voxels in label %d\n", volume, 
						 volume*vox_volume, label) ;
      if (log_fp)
      {
        fprintf(log_fp,"%2.1f\n", vox_volume*(float)volume) ;
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
	if (!stricmp(option, "ICV"))
	{
		icv_fname = argv[2] ;
		printf("reading ICV from %s\n", icv_fname) ;
		nargs = 1 ;
	}
	else switch (toupper(*option))
  {
	case 'A':
		all_flag = 1 ;
		printf("computing volume of all non-zero voxels\n") ;
		break ;
  case 'T':
    in_label = atoi(argv[2]) ;
    out_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("translating label %d to label %d\n", in_label, out_label) ;
    break ;
	case 'B':
		brain_fname = argv[2] ;
    compute_pct = 1 ;
		nargs = 1 ;
		printf("reading brain volume from %s...\n", brain_fname) ;
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
  printf("usage: %s [options] <volume 1> <volume 2>",
         Progname) ;
  printf(
         "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
         );
  exit(code) ;
}
