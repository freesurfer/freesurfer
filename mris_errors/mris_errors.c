#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"

static char vcid[]="$Id: mris_errors.c,v 1.1 1997/09/23 21:55:27 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int    MRISareaErrors(MRI_SURFACE *mris) ;

char *Progname ;
static MRI_SURFACE  *mris ;


int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, fname[100] ;
  int          ac, nargs ;

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
    usage_exit() ;

  in_fname = argv[1] ;
#if 0
  out_fname = argv[2] ;
  cp = strrchr(out_fname, '.') ;
#endif

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  MRISareaErrors(mris) ;

  if (0)
  {
    sprintf(fname, "%s.area_error", out_fname) ;
    printf("writing area errors to %s\n", fname) ;
    MRISwriteAreaError(mris, fname) ;
    sprintf(fname, "%s.angle_error", out_fname) ;
    printf("writing angle errors to %s\n", fname) ;
    MRISwriteAngleError(mris, fname) ;
  }

  MRISfree(&mris) ;

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
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will unfold an MRI on the surface of an ellipsoid.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

int
MRISareaErrors(MRI_SURFACE *mris)
{
  int      vno, fno, max_v = -1, n ;
  VERTEX   *v ;
  float    verror, max_verror, total_error, total_sq_error,
           error, mean_error, std_error, pct_error ;

  MRISupdateEllipsoidSurface(mris) ;
  total_error = total_sq_error = max_verror = 0.0f ;
  n = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    verror = 0.0f ;
    for (fno = 0 ; fno < v->num ; n++, fno++)  /* for each neighboring face */
    {
      error = v->tri_area[fno] - v->orig_tri_area[fno] ;
      pct_error = error / v->orig_tri_area[fno] * 100.0f ;
      printf("%d %d %d %2.3f %2.3f %2.3f %2.1f\n", n, vno, fno, 
             v->orig_tri_area[fno], v->tri_area[fno], error, pct_error) ;
      total_sq_error += (error * error) ;
      verror += fabs(error) ;
      total_error += error ;
    }
    if (verror >= max_verror)
    {
      max_verror = verror ;
      max_v = vno ;
    }
  }

  mean_error = total_error / (float)n ;
  std_error = sqrt(total_sq_error / (float)n - mean_error*mean_error) ;
  fprintf(stderr, "max error occurs at %d, error = %2.3f\n",max_v, max_verror);
  fprintf(stderr, "mean error = %2.3f, std = %2.3f\n", mean_error, std_error);
  return(NO_ERROR) ;
}

