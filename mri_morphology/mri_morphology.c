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
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
                

char *Progname ;

#define ERODE          1
#define DILATE         2
#define CLOSE          3
#define OPEN           4
#define DILATE_LABEL   5

static void usage_exit(int code) ;



int
main(int argc, char *argv[])
{
  char   *out_fname, **av ;
  int    ac, nargs, niter, operation ;
  MRI    *mri_src, *mri_dst ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_morphology.c,v 1.1 2005/05/04 19:28:42 fischl Exp $", "$Name:  $");
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
	mri_src = MRIread(argv[1]) ;
	if (mri_src == NULL)
		ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n", 
							Progname,argv[1]);

	if (!stricmp(argv[2], "dilate"))
		operation = DILATE ;
	else 	if (!stricmp(argv[2], "dilatelabel"))
		operation = DILATE_LABEL ;
	else 	if (!stricmp(argv[2], "open"))
		operation = OPEN ;
	else 	if (!stricmp(argv[2], "close"))
		operation = CLOSE ;
	else 	if (!stricmp(argv[2], "open"))
		operation = OPEN ;
	else
	{
		operation = 0 ;
		ErrorExit(ERROR_UNSUPPORTED, "morphological operation '%s'  is not supported", argv[2]) ;
	}

	niter = atoi(argv[3]) ;
	switch (operation)
	{
	case DILATE_LABEL:
		{
			int label ;

			label = atoi(argv[4]) ;
			printf("dilating label %s (%d) %d times\n", 
						 cma_label_to_name(label), label, niter) ;
			mri_dst = MRIdilateLabel(mri_src, NULL, label, niter) ;
			break ;
		}
	default:
		break ;
	}
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  MRIfree(&mri_dst) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphological processing took %d minutes and %d seconds.\n", 
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
  if (!stricmp(option, "dt"))
  {
  }
  else switch (toupper(*option))
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
  printf("usage: %s [options] <volume> <operation> <# iter> <out volume>\n", Progname) ;
  printf("\t-u              where which can be [dilate,erode]\n");
  exit(code) ;
}

