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

static int features = FEATURE_INTENSITY | FEATURE_MEAN3 | FEATURE_DIRECTION |
                        FEATURE_CPOLV_MEDIAN5 ;

static int extract = 0 ;
static int classifier = CLASSIFIER_RBF ;
static char priors_fname[100] = "none" ;
static int  verbose = 0 ;

char *Progname ;

void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

#define NCLUSTERS  6

static int nclusters = 0 ;
static int train_cpolv = 0 ;

static RBF_PARMS rbf_parms =
{
{ 2, NCLUSTERS, NCLUSTERS, NCLUSTERS, 2}
} ;
void 
main(int argc, char *argv[])
{
  MRIC    *mric ;
  char    *training_file_name, *output_file_name ;
  int     nargs, error ;

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
    ErrorExit(ERROR_BADPARM,"usage: %s <training file name>  <mric file>",
              Progname);

  training_file_name = argv[1] ;
  output_file_name = argv[2] ;

  if (nclusters > 0)
  {
    rbf_parms.max_clusters[GRAY_MATTER] = nclusters ;
    rbf_parms.max_clusters[WHITE_MATTER] = nclusters ;
    rbf_parms.max_clusters[BORDER_MATTER] = nclusters ;
  }
  else
    nclusters = NCLUSTERS ;

  if (train_cpolv)
  {
    rbf_parms.max_clusters[CSF] = nclusters ;
    rbf_parms.max_clusters[GRAY_MATTER] = 0 ;
    rbf_parms.max_clusters[WHITE_MATTER] = nclusters/2 ;
    rbf_parms.max_clusters[BORDER_MATTER] = nclusters/2 ;
    rbf_parms.max_clusters[BRIGHT_MATTER] = 0 ;
  }

  mric = MRICalloc(1, &classifier, &features, (void *)&rbf_parms) ;
  if ((strlen(priors_fname) > 1) && stricmp(priors_fname, "none"))
    error = MRICtrain(mric, training_file_name, priors_fname) ;
  else
    error = MRICtrain(mric, training_file_name, NULL) ;

  if (error != NO_ERROR)
    ErrorExit(error, "training failed.\n") ;

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
  if (!stricmp(option, "cpolv"))
    train_cpolv = 1 ;
  else switch (toupper(*option))
  {
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
    if (sscanf(argv[2], "%d", &nclusters) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
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
  case 'X':
    if (sscanf(argv[2], "%d", &extract) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
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
