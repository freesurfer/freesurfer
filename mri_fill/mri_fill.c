 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <memory.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mri.h"

static char vcid[] = "$Id: mri_fill.c,v 1.2 1997/09/15 19:34:57 fischl Exp $";

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

#define CSF_THRESHOLD        50 
#define VENTRICLE_FILL_VAL   200

static Real tal_seed_x = -2 ;
static Real tal_seed_y = -2 ;
static Real tal_seed_z = 19 ;

static int threshold = CSF_THRESHOLD ;
static int fill_val = VENTRICLE_FILL_VAL ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  Real   xv, yv, zv ;
  char   *in_fname, *wm_fname, *out_fname ;
  MRI    *mri_src, *mri_wm ;

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
  wm_fname = argv[2] ;
  out_fname = argv[3] ;

  mri_src = MRIread(in_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read '%s'", Progname, in_fname) ;
  mri_wm = MRIread(wm_fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read '%s'", Progname, wm_fname) ;

  MRItalairachToVoxel(mri_src, tal_seed_x, tal_seed_y, tal_seed_z,&xv,&yv,&zv);
  MRIfill(mri_src, mri_wm, nint(xv), nint(yv), nint(zv), threshold, fill_val) ;
  MRIwrite(mri_wm, out_fname) ;
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
  case 'S':
    if (sscanf(argv[2], "%lf", &tal_seed_x) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan x seed from %s\n",
                Progname, argv[2]) ;
    if (sscanf(argv[3], "%lf", &tal_seed_y) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan y seed from %s\n",
                Progname, argv[3]) ;
    if (sscanf(argv[4], "%lf", &tal_seed_z) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan z seed from %s\n",
                Progname, argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "using seed = (%2.1f, %2.1f, %2.1f)\n", 
            tal_seed_x, tal_seed_y, tal_seed_z) ;
    break ;
  case 'T':
    if (sscanf(argv[2], "%d", &threshold) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan threshold from %s\n",
                Progname, argv[2]) ;
    fprintf(stderr, "using threshold = %d\n", threshold) ;
    nargs = 1 ;
    break ;
  case 'F':
    if (sscanf(argv[2], "%d", &fill_val) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan fill val from %s\n",
                Progname, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using fill val = %d\n", fill_val) ;
    break ;
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
          "usage: %s [options] <input image file> <output image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, "\nThis program will fill a region of an MR image\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, 
          "  --version         - display version #.\n") ;
  fprintf(stderr, 
          "  --help            - display this help message.\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

