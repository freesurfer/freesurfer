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

#define SCALE     16
#define NINPUTS   5

static int scale = SCALE ;
static int ninputs = NINPUTS ;
static int extract = 0 ;

char *Progname ;

void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

void 
main(int argc, char *argv[])
{
  FILE    *fp ;
  MRIC    *mric ;
  char    *training_file_name, *output_file_name, line[300], *cp, fname[100] ;
  int     nargs ;
  MRI     *mri ;

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

  /* figure out proper dimensions */
  fp = fopen(training_file_name, "r") ;
  if (!fp)
    ErrorExit(ERROR_NO_FILE,
                 "%s:  could not open file %s", Progname, training_file_name);
  
  cp = fgetl(line, 299, fp) ;
  if (!cp)
    ErrorExit(ERROR_BADFILE, "%s: could not scan any file names from %s",
              Progname, training_file_name) ;
  sscanf(cp, "%s", fname) ;
  mri = MRIreadInfo(fname) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not open info file in %s",
              Progname, fname) ;
  fclose(fp) ;
  mric =  MRIclassAlloc(mri->width, mri->height, mri->depth, scale, ninputs) ;
  MRIfree(&mri) ;
  MRIclassTrainAll(mric, training_file_name) ;
  MRIclassWrite(mric, output_file_name) ;
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
    if (sscanf(argv[2], "%d", &ninputs) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    if (ninputs != 2 && ninputs != 5)
      ErrorExit(ERROR_BADPARM, "%s: # of inputs must be 2 (default) or 5",
                Progname) ;
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
