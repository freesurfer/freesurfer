#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "mriclass.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "classify.h"
#include "gcarray.h"

#define SCALE     16

static int scale = SCALE ;
static int features = FEATURE_INTENSITY | FEATURE_ZSCORE3 | FEATURE_MEAN3 | FEATURE_DIRECTION ;

static int extract = 0 ;
static int classifier = CLASSIFIER_GAUSSIAN ;
static char priors_fname[100] = "priors.mnc" ;
static int  verbose = 0 ;

char *Progname ;

void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

void 
main(int argc, char *argv[])
{
  MRIC    *mric ;
  char    *training_file_name, *output_file_name ;
  int     nargs ;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,"usage: %s <training file name>  <mric file>");

  training_file_name = argv[1] ;
  output_file_name = argv[2] ;

  mric = MRICalloc(classifier, features, NULL) ;
  MRICtrain(mric, training_file_name, priors_fname) ;
  MRICwrite(mric, output_file_name) ;
  MRICfree(&mric) ;
  exit(0) ;
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
  case 'V':
    verbose = !verbose ;
    break ;
  case 'F':
    if (sscanf(argv[2], "0x%x", &features) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using features 0x%x\n", features) ;
    break ;
  case 'P':
    strcpy(priors_fname, argv[2]) ;
    nargs = 1 ;
    if (verbose)
      fprintf(stderr, "using priors file %s\n", priors_fname) ;
    break ;
  case 'S':
    if (sscanf(argv[2], "%d", &scale) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'X':
    if (sscanf(argv[2], "%d", &extract) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'N':
#if 0
    if (sscanf(argv[2], "%d", &ninputs) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
#endif
    break ;
  case '?':
  case 'U':
    printf("usage: %s <training file> <output file>\n", Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
