#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "classify.h"
#include "mri.h"

static int extract = 0 ;
static int  verbose = 0 ;
static int wsize = 5 ;
static float pct = 0.8 ;

char *Progname ;

void main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

void 
main(int argc, char *argv[])
{
  MRI     *mri_src, *mri_dst, *mri_cpolv, *mri_tmp ;
  char    *input_file_name, *output_file_name ;
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
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input volume> <output volume>", Progname);

  input_file_name = argv[1] ;
  output_file_name = argv[2] ;

  mri_src = MRIread(input_file_name) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
              Progname, input_file_name) ;

  fprintf(stderr, "computing orientation of gray/white interface...\n") ;
  mri_cpolv = MRIcentralPlaneOfLeastVarianceNormal(mri_src, NULL, wsize) ;
  fprintf(stderr, "labelling white matter...\n") ;
  mri_tmp = MRIwmfilter(mri_src, mri_cpolv, NULL) ;
  mri_dst = MRIremoveHoles(mri_tmp, NULL, 3, pct) ;

  fprintf(stderr, "writing output to %s...", output_file_name) ;
  MRIwrite(mri_dst, output_file_name) ;
  fprintf(stderr, "done.\n") ;

  MRIfree(&mri_cpolv) ;
  MRIfree(&mri_src) ;
  MRIfree(&mri_dst) ;
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
  case 'P':
    pct = atof(argv[2]) ;
    fprintf(stderr, "using %2.0f%% threshold\n", pct*100.0f) ;
    break ;
  case 'X':
    if (sscanf(argv[2], "%d", &extract) != 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                Progname, argv[2]) ;
    nargs = 1 ;
    break ;
  case 'W':
    wsize = atoi(argv[2]) ;
    fprintf(stderr, "using wsize = %d\n", wsize) ;
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
