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

static int extract = 0 ;
static int  verbose = 0 ;

char *Progname ;

void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

void 
main(int argc, char *argv[])
{
  MRIC    *mric ;
  MRI     *mri_input, *mri_output ;
  char    *classifier_file_name, *input_file_name, *output_file_name ;
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

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <classifier file> <input volume> <output volume>",
              Progname);

  classifier_file_name = argv[1] ;
  input_file_name = argv[2] ;
  output_file_name = argv[3] ;

  mric = MRICread(classifier_file_name) ;
  if (!mric)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier from %s",
              Progname, classifier_file_name) ;
  mri_input = MRIread(input_file_name) ;
  if (!mri_input)
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
              Progname, input_file_name) ;

  mri_output = MRICclassify(mric, mri_input, NULL, 0.0f, NULL, NULL) ;
  MRIwrite(mri_output, output_file_name) ;
  MRICfree(&mric) ;
  MRIfree(&mri_input) ;
  MRIfree(&mri_output) ;
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
  case 'X':
    if (sscanf(argv[2], "%d", &extract) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    fprintf(stderr,
            "usage: %s <classifier file> <input volume> <output volume>\n",
            Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
