#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "transform.h"

static char vcid[] = "$Id: mri_synthesize.c,v 1.1 2002/04/09 17:42:08 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;

static MRI *MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE) ;

char *Progname ;

int
main(int argc, char *argv[])
{
  char        **av, *out_fname, *T1_fname, *PD_fname ;
  int         ac, nargs ;
  MRI         *mri_T1, *mri_PD, *mri_out ;
  float       TR, TE, alpha ;

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

  if (argc < 6)
    usage_exit() ;

  TR = atof(argv[1]) ;
  alpha = atof(argv[2]) ;
  TE = atof(argv[3]) ;
  T1_fname = argv[4] ;
  PD_fname = argv[5] ;
  out_fname = argv[6] ;

  printf("reading T1 volume from %s...\n", T1_fname) ;
  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume %s", Progname, 
              T1_fname) ;

  printf("reading PD volume from %s...\n", PD_fname) ;
  mri_PD = MRIread(PD_fname) ;
  if (!mri_PD)
    ErrorExit(ERROR_NOFILE, "%s: could not read PD volume %s", Progname, 
              PD_fname) ;

  printf("synthesizing volume with TR=%2.1f msec, TE=%2.1f msec, and alpha=%2.2f degrees...\n",
         TR, TE, alpha) ;
  mri_out = MRIsynthesize(mri_T1, mri_PD, NULL, TR, RADIANS(alpha), TE) ;
  printf("writing output to %s.\n", out_fname) ;
  MRIwrite(mri_out, out_fname) ;

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
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <TR> <alpha (deg)> <TE> <T1 volume> <PD volume> <output volume>\n",
          Progname) ;
}

static void
print_help(void)
{
  fprintf(stderr, 
          "\nThis program will synthesize a flash acquisition based on previously computed T1/PD maps\n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double  CFA = 1, SFA = 0 ;
  double  E1, FLASH ;

  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;

  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}

static MRI *
MRIsynthesize(MRI *mri_T1, MRI *mri_PD, MRI *mri_dst, double TR, double alpha, double TE)
{
  int   x, y, z, width, height, depth ;
  double flash, T1, PD ;
  

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  width = mri_T1->width ; height = mri_T1->height ; depth = mri_T1->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        T1 = MRISvox(mri_T1, x, y, z) ;
        PD = MRISvox(mri_PD, x, y, z) ;
        flash = FLASHforwardModel(alpha, TR, PD, T1) ;
        MRISvox(mri_dst, x, y, z) = (short)nint(flash) ;
      }
    }
  }

  return(mri_dst) ;
}

