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
#include "transform.h"

static char vcid[] = "$Id: mris_transform.c,v 1.1 2002/07/02 16:41:58 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static MRI_SURFACE  *mris ;
static char *inverse_xform_fname = NULL ;
static int invert = 0 ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, *xform_fname ;
  int          ac, nargs ;
  LTA          *lta ;

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

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  xform_fname = argv[2] ;
  out_fname = argv[3] ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  lta = LTAread(xform_fname) ;
  if (!lta)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
              Progname, xform_fname) ;

  if (lta->type == LINEAR_VOX_TO_VOX)
    LTAvoxelTransformToCoronalRasTransform(lta) ;

  if (invert)
  {
    MATRIX *m_tmp = lta->xforms[0].m_L ;
    lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
    MatrixFree(&m_tmp) ;
  }

  if (inverse_xform_fname)
  {
    LTA          *inv_lta ;
    MATRIX       *m_tmp, *m_inv ;

    inv_lta = LTAread(inverse_xform_fname) ;
    if (!inv_lta)
    ErrorExit(ERROR_NOFILE, "%s: could not read inverse transform file %s",
              Progname, inverse_xform_fname) ;

    if (inv_lta->type == LINEAR_VOX_TO_VOX)
      LTAvoxelTransformToCoronalRasTransform(inv_lta) ;
    m_inv = MatrixInverse(inv_lta->xforms[0].m_L, NULL) ;
    m_tmp = MatrixMultiply(m_inv, lta->xforms[0].m_L, NULL) ;
    MatrixCopy(m_tmp, lta->xforms[0].m_L) ;
    MatrixFree(&m_tmp) ; MatrixFree(&m_inv) ;
  }
  MRIStransform(mris, NULL, lta) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;

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
  else if (!stricmp(option, "invert"))
    invert = 1 ;
  else switch (toupper(*option))
  {
  case 'I':
    inverse_xform_fname = argv[2] ;
    nargs = 1 ;
    printf("left multiplying by inverse transform %s\n",
           inverse_xform_fname) ;
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
          "usage: %s [options] <input surf> <transform file> <output surf>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will transform an MRI surface into Talairach space.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

