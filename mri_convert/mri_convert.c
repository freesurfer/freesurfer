#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"


void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static int verbose = 0 ;

void
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri ;
  char   *in_fname, *out_fname ;

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

  if (argc < 1)
    argc = 1 ;

  if (argc < 1)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[1] ;

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  out_fname = argv[2] ;

  if (verbose)
    fprintf(stderr, "reading from %s...", in_fname) ;
  mri = MRIread(in_fname) ;
  if (!mri)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, in_fname) ;

  if (!mri->imnr0)
  {
    mri->imnr0++ ;
    mri->imnr1++ ;
  }
  if (verbose)
    fprintf(stderr, "done.\nwriting to %s...", out_fname) ;
  MRIwrite(mri, out_fname) ;
  if (verbose)
    fprintf(stderr, "done.\n") ;
  MRIfree(&mri) ;
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
  case 'N':
#if 0
    sscanf(argv[2], "%d", &reductions) 
    fprintf(stderr, "reducing %d times\n", reductions) ;
#endif
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
