/**
 * @brief program for deforming a surface to lie at the gray/white or pial boundary from 
 *  ultra-high res data
 *
 * Fit a generative piecewise constant model to the data to determine
 * target locations and deform the surface to match them.
 */
/*
 * Original Author: Bruce Fischl
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "transform.h"
#include "mrisurf.h"
#include "label.h"
#include "tritri.h"
#include "filter.h"


const char *Progname = NULL ;

static void usage_exit(int code) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

#define MAX_VERTICES 10000
static int nvertices = 0;
static int target_vnos[MAX_VERTICES] ;
static float target_vals[MAX_VERTICES] ;
static double sigma = 2.0 ;

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs, nsize ;
  MRI_SURFACE  *mris ;
  MRI          *mri ;

  nargs = handleVersionOption(argc, argv, "mris_nudge");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gx = Gy = Gz = -1 ;
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

  if (argc < 7)
    usage_exit(1) ;

  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, argv[1]) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;

  mri = MRIread(argv[2]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s", Progname, argv[2]) ;

  target_vnos[nvertices] = atoi(argv[3]) ;
  target_vals[nvertices] = atof(argv[4]) ;
  nsize = atoi(argv[5]) ;

  printf("nudging %d vertex region around vertex %d to target val %2.1f\n", 
         nsize, target_vnos[nvertices], target_vals[nvertices]) ;
  nvertices++ ;

  MRISerodeRipped(mris, nsize) ;
  MRISrepositionSurface(mris, mri, target_vnos, target_vals, nvertices, nsize, sigma, 0) ;

  MRISunrip(mris) ;
  printf("writing repositioned surface to %s\n", argv[6]) ;
  MRISwrite(mris, argv[6]) ;
  return(0) ;
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
  if (!stricmp(option, "vavgs")) {
  }
  else switch (toupper(*option)) {
  case 'S':
    sigma = atof(argv[2]) ;
    printf("using sigma = %2.1f\n", sigma) ;
    nargs = 1 ;
    break ;
  case 'V':
    target_vnos[nvertices] = atoi(argv[2]) ;
    target_vals[nvertices] = atof(argv[2]) ;
    printf("moving vertex %d to value %2.0f\n", target_vnos[nvertices], target_vals[nvertices]) ;
    nvertices++ ;
    nargs = 2 ;
    break ;
  default:
    usage_exit(-1) ;
    break ;
  }
  return(nargs) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <input surface> <input volume> <vertex> <target val> <nbhd> <output surf>\n",
         Progname) ;
  printf(
         "\t\n") ;
  exit(code) ;
}


