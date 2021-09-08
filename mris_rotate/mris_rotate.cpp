/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */



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
#include "mri.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
char *regfile = NULL;
int VolGeomValid = 1;

int main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  float        alpha, beta, gamma ;

  nargs = handleVersionOption(argc, argv, "mris_rotate");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  in_fname = argv[1] ;
  if (sscanf(argv[2], "%f", &alpha) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan alpha from %s",
              Progname, argv[2]) ;
  if (sscanf(argv[3], "%f", &beta) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan beta from %s",
              Progname, argv[3]) ;
  if (sscanf(argv[4], "%f", &gamma) != 1)
    ErrorExit(ERROR_BADPARM, "%s: could not scan gamma from %s",
              Progname, argv[4]) ;
  out_fname = argv[5] ;

  mris = MRISfastRead(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  if(regfile){
    printf("Applying rotational components from regfile, ignoring alpha, beta, gamma\n");
    MRIScenter(mris, mris) ; // not sure this is needed if this is a sphere
    printf("Reading in reg file %s\n",regfile);
    LTA *lta = LTAread(regfile);
    if(lta==NULL) exit(1);
    printf("Extracting rotational components\n");
    LTAmat2RotMat(lta);
    printf("Applying rotation matrix to surface\n");
    MatrixPrint(stdout,lta->xforms[0].m_L);
    int err = MRISltaMultiply(mris, lta);
    if(err) exit(1);
    LTAfree(&lta);
    // This may be needed because lta is changed to tkreg and it gets an offset
    MRIScenter(mris, mris) ;
  }
  else {
    printf("alpha = %g  beta = %g gamma = %g\n",alpha,beta,gamma);
    alpha = RADIANS(alpha) ;
    beta = RADIANS(beta) ;
    gamma = RADIANS(gamma) ;
    MRIScenter(mris, mris) ;
    MRISrotate(mris, mris, alpha, beta, gamma) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not rotate surface", Progname) ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing rotated surface to %s\n", out_fname) ;

  if(VolGeomValid == 0){
    // This can be useful with spherical surfs when viewing in freeview
    printf("Making vg invalid\n");
    mris->vg.valid = 0;
  }

  MRISwrite(mris, out_fname) ;

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'R':
      regfile = argv[2];
      nargs++;
      break ;
    case 'N':
      VolGeomValid = 0;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  printf("%s [options] input_surface alphaDeg betaDeg gammaDeg output_surface\n", Progname) ;
  printf("-r regfile : extract angles from regfile (ignores alpha, beta, gamma)\n");
  printf("-n : invalidate volume geometry in output\n");
}

static void
print_help(void) {
  print_usage() ;
  printf("This program will rotate a surface given the three angles\n");
  exit(1) ;
}

static void
print_version(void) {
  printf("%s\n", getVersion().c_str()) ;
  exit(1) ;
}

