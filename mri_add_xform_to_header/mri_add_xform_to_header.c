#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static int verbose = 0 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri ;
  char   *xform_fname, *in_fname, *out_fname ;

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

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no transform name specified", Progname) ;
  xform_fname = argv[1] ;

  if (argc < 3)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[2] ;

  if (argc < 4)
    out_fname = in_fname ;
  else
    out_fname = argv[3] ;

  if (verbose)
    fprintf(stderr, "reading from %s...", in_fname) ;
  mri = MRIreadInfo(in_fname) ;
  if (!mri)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, xform_fname) ;

  if (input_transform_file(xform_fname, &mri->transform) != OK)
    ErrorPrintf(ERROR_NO_MEMORY, 
                "%s: could not read xform file '%s'\n", Progname, xform_fname);

  mri->linear_transform = get_linear_transform_ptr(&mri->transform) ;
  mri->inverse_linear_transform = 
    get_inverse_linear_transform_ptr(&mri->transform) ;
  mri->free_transform = 1 ;
  strcpy(mri->transform_fname, xform_fname) ;
  if (verbose)
    fprintf(stderr, "done.\nwriting to %s...", out_fname) ;
  MRIwriteInfo(mri, out_fname) ;
  if (verbose)
    fprintf(stderr, "done.\n") ;
  MRIfree(&mri) ;
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
   printf("usage: %s [xform file name] [input directory] [output directory]\n",
          argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
