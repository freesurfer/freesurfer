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

static char vcid[]="$Id: mris_errors.c,v 1.2 1997/09/25 22:32:15 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISareaErrors(MRI_SURFACE *mris) ;

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

  MRISreadTriangleProperties(mris, in_fname) ;
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
  int      fno, tno, max_f = -1 ;
  FACE     *face ;
  float    ferror, max_ferror, total_error, total_sq_error,
           error, mean_error, std_error, pct_error, n ;

  MRISupdateEllipsoidSurface(mris) ;
  total_error = total_sq_error = max_ferror = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    ferror = 0.0f ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      error = face->area[tno] - face->orig_area[tno] ;
      pct_error = error / face->orig_area[tno] * 100.0f ;
      printf("%d %2.3f %2.3f %2.3f %2.1f\n", fno, face->orig_area[tno], 
             face->area[tno], error, pct_error) ;
      total_sq_error += (error * error) ;
      ferror += fabs(error) ;
      total_error += error ;
    }
    if (ferror >= max_ferror)
    {
      max_ferror = ferror ;
      max_f = fno ;
    }
  }

  n = (float)(2*mris->nfaces) ;
  mean_error = total_error / n ;
  std_error = sqrt(total_sq_error / (float)n - mean_error*mean_error) ;
  fprintf(stderr, "max error occurs at %d, error = %2.3f\n",max_f, max_ferror);
  fprintf(stderr, "mean error = %2.3f, std = %2.3f\n", mean_error, std_error);
  return(NO_ERROR) ;
}

