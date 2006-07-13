#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mri.h"
#include "diag.h"
#include "tags.h"
#include "version.h"
#include "error.h"

static char vcid[] = "$Id: mri_threshold.c,v 1.4 2006/07/13 16:07:44 fischl Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int binarize = 0 ;
char *Progname ;

int
main(int argc, char *argv[])
{
  char        **av, *in_fname, *out_fname ;
  int         ac, nargs ;
  MRI         *mri_in, *mri_out ;
  float       thresh ;

	char cmdline[CMD_LINE_LEN] ;
	
  make_cmd_version_string (argc, argv, "$Id: mri_threshold.c,v 1.4 2006/07/13 16:07:44 fischl Exp $", "$Name:  $", cmdline);
  /* rkt: check for and handle version tag */
	Progname = argv[0] ;
  nargs = handle_version_option (argc, argv, "$Id: mri_threshold.c,v 1.4 2006/07/13 16:07:44 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
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

	printf("thresholding volume at %2.4f....\n", thresh) ;
	mri_out = MRIthreshold(mri_in, NULL, thresh) ;

	if (binarize > 0)
		MRIbinarize(mri_out, mri_out, thresh, 0, binarize) ;

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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
	case 'B':
		binarize = atoi(argv[2]) ;
		printf("binarizing output volume to be [0, %d]\n", binarize) ;
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
usage_exit(void)
{
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] <in vol> <thresh> <out_vol>\n",
				 Progname) ;
}

static void
print_help(void)
{
  fprintf(stderr, 
          "\nThis program will threshold an input volume\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
