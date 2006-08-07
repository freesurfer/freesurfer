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
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
                

char *Progname ;

static int use_log = 0 ;
static float sigma = 0 ;
static int write_areas = 0 ;

static void usage_exit(int code) ;


int
main(int argc, char *argv[])
{
  char      **av, *out_fname ;
  int       ac, nargs ;
	GCA_MORPH *gcam ;
  int       msec, minutes, seconds ;
  struct timeb start ;
	MRI       *mri, *mri_jacobian, *mri_area, *mri_orig_area ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_jacobian.c,v 1.1 2006/08/07 21:03:33 fischl Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit(1) ;

	out_fname = argv[argc-1] ;
	gcam = GCAMread(argv[1]) ;
	if (gcam == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read input morph %s\n", Progname,argv[1]);
	mri = MRIread(argv[2]) ;
	if (gcam == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read template volume %s\n", Progname,argv[2]);

	mri_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_AREA, 0, 0);
	mri_orig_area = GCAMmorphFieldFromAtlas(gcam, mri, GCAM_ORIG_AREA, 0, 0);
	if (FZERO(sigma) == 0)
	{
		MRI *mri_kernel, *mri_smooth ;
		mri_kernel = MRIgaussian1d(sigma, 100) ;
		mri_smooth = MRIconvolveGaussian(mri_area, NULL, mri_kernel) ;
		MRIfree(&mri_area) ; mri_area = mri_smooth ;
		mri_smooth = MRIconvolveGaussian(mri_orig_area, NULL, mri_kernel) ;
		MRIfree(&mri_orig_area) ; mri_orig_area = mri_smooth ;

		MRIfree(&mri_kernel) ; 
	}
	mri_jacobian = MRIdivide(mri_area, mri_orig_area, mri_jacobian) ;
	if (use_log)
		MRIlog10(mri_jacobian, mri_jacobian, 0) ;
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_jacobian, out_fname) ;
  MRIfree(&mri_jacobian) ;
	if (write_areas)
	{
		MRIwrite(mri_area, "area.mgz") ;
		MRIwrite(mri_orig_area, "orig_area.mgz") ;
	}
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ; minutes = seconds / 60 ;seconds = seconds % 60 ;
  fprintf(stderr, "jacobian calculation took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "dt"))
  {
  }
  else switch (toupper(*option))
  {
	case 'W':
		write_areas = 1 ;
		printf("writing area volumes\n") ;
		break ;
	case 'L':
		use_log = 1 ;
		printf("taking log of jacobian values before saving\n") ;
		break ;
	case 'S':
		sigma = atof(argv[2]) ;
		printf("smoothing jacobian volume with sigma=%2.2f\n", sigma) ;
		nargs = 1 ;
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
  printf("usage: %s [options] <3d morph> <template volume> <output volume>\n", Progname) ;
  exit(code) ;
}

