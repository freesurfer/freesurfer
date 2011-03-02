/**
 * @file  mris_density.c
 * @brief program for generating a surface map of the density of the interior voxels.
 *
 * Compute the # of voxels that are interior to the surface in a user-specified radius
 * at each point on the surface.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:31 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

static char vcid[] = "$Id: mris_density.c,v 1.7 2011/03/02 00:04:31 nicks Exp $";


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static double resolution = 1.0/8.0 ;
static double radius = 20 ;
static char *density_fname = NULL ;
static char *translate_fname = NULL ;

int
main(int argc, char *argv[]) {
  char          **av, *out_fname, *in_fname ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  MRI           *mri_density ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_density.c,v 1.7 2011/03/02 00:04:31 nicks Exp $", "$Name:  $");
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

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  fprintf(stderr, "reading surface from %s...\n", in_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",Progname, in_fname) ;
  if (translate_fname && Gdiag_no >= 0)
  {
    MRI_SURFACE *mris2 ;
    float       x0, y0, z0, dist, min_dist, x, y, z ;
    int         vno, min_vno=0;
    char        surf_name[STRLEN] ;
    VERTEX      *v ;

    mris2 = MRISread(translate_fname) ;
    if (!mris2)
      ErrorExit(ERROR_NOFILE, "%s: could not load translation surface %s", Progname, translate_fname) ;
    FileNameOnly(translate_fname, surf_name) ;
    if (MRISreadCanonicalCoordinates(mris, surf_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read canonical coords from %s", Progname, surf_name) ;
    x0 = mris2->vertices[Gdiag_no].x ;
    y0 = mris2->vertices[Gdiag_no].y ;
    z0 = mris2->vertices[Gdiag_no].z ;
    min_dist = 1e10 ; min_vno = -1 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      x = v->cx ; y = v->cy ; z = v->cz ;
      dist = sqrt(SQR(x-x0) + SQR(y-y0) + SQR(z-z0)) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_vno = vno ;
      }
    }
    printf("translating to vertex %d (min dist = %2.1f)\n", min_vno, min_dist) ;
    Gdiag_no = min_vno ;
  }
    
  MRISmakeDensityMap(mris, resolution, radius, Gdiag_no, &mri_density) ;
  if (density_fname != NULL)
  {
    printf("writing density volume for %d to %s\n", Gdiag_no, density_fname) ;
    MRIwrite(mri_density, density_fname) ;
    MRIfree(&mri_density) ;
  }

  fprintf(stderr, "writing density map to curvature file %s...\n", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;
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
  else if (!stricmp(option, "radius")) {
    radius = atof(argv[2]) ;
    nargs =  1 ;
    printf("using radius = %2.3f\n", radius) ;
  } else if (!stricmp(option, "debug")) {
    Gdiag_no = atoi(argv[2]) ;
    density_fname = argv[3] ;
    nargs =  2 ;
    printf("debugging vertex %d, and writing density map for it to %s\n", Gdiag_no, density_fname) ;
  } else switch (toupper(*option)) {
  case 'R':
    resolution = (double)atof(argv[2]) ;
    nargs = 1 ;
    printf("setting resolution for intermediate calculations to %2.4f\n", resolution) ;
    break ;
  case 'T':
    translate_fname = argv[2] ;
    printf("translating vertex %d on surface %s to current surface\n", Gdiag_no, translate_fname) ;
    nargs = 1 ;
    break ;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <output map>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program computes a density map and attaches it to each vertex on the surface.\n"
          "The resulting measurement are written into a 'curvature' file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


