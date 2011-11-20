/**
 * @file  mris_curvature_vectors.c
 * @brief returns surface normals and directions of principal curvature directions
 *
 */
/*
 * Original Author:
 * CVS Revision Info:
 *    $Author: jonp $
 *    $Date: 2011/11/20 01:17:34 $
 *    $Revision: 1.1 $
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

#include <assert.h>
#include <errno.h>

#include <getopt.h>
#include <stdarg.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "xDebug.h"
#include "label.h"

#define  STRBUF         65536


static char vcid[] =
  "$Id: mris_curvature_vectors.c,v 1.1 2011/11/20 01:17:34 jonp Exp $";

int main(int argc, char *argv[]) ;

static void  usage_exit(int code) ;
static void  print_usage(void) ;
static void  print_help(void) ;
static void  print_version(void) ;

static int get_option(int argc, char *argv[]) ;

int MRISsortPrincipalDirectionsByCurvatures(MRI_SURFACE *mris) ;
int MRISorderPrincipalDirectionsConsistentWithNormal(MRI_SURFACE *mris) ;


// Global variables

//char*  Progname ;
char  Progname[64] ;
char*  hemi ;

short b_sortBySignedPrincipalCurv = 0 ;  // boolean flag, user option for sorting by signed or unsigned k1 and k2

static int      G_nbrs                          = 2 ;

int
main(int argc, char *argv[]) {
  char          output_filename[STRBUF] ;
  char          *surf_name, *output_filename_stem ;
  int           nargs ;
  MRI_SURFACE   *mris ;

  InitDebugging( "mris_curvature_vectors" ) ;
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv,
                                 "$Id: mris_curvature_vectors.c,v 1.1 2011/11/20 01:17:34 jonp Exp $", "$Name:  $") ;

  //Progname = argv[0] ;
  sprintf(Progname, "mris_curvature_vectors") ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc < 2 || argc > 4)
    usage_exit(1) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  surf_name  = argv[1] ;

  if ( argc <= 2 )
    {
      // if no output stem is specified, use the input file name as a stem
      output_filename_stem  = surf_name ;
    }
  else
    {
      output_filename_stem  = argv[2] ;
    }

  fprintf(stdout, "reading surface file:  %s\n", surf_name) ;


  mris = MRISread(surf_name) ;

  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_name) ;

  // borrowed the 'cprints' scheme from mris_curvature_stats  :)
  cprints("calculating surface normals...", "") ;
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeNormals(mris) ;
  cprints("", "ok") ;

  sprintf(output_filename, "%s_normals.mgz", output_filename_stem) ;
  fprintf(stdout, "writing normals file:  %s\n", output_filename) ;
  MRISwriteNormals(mris, output_filename) ;
  MRISrestoreRipFlags(mris) ;

  // bruce says: these routines will "use a quadratic patch fit using
  // the 2-nbrs (probably 3 is the most you would want. 1 might not be
  // enough to do the fit and would be noisy)".

  MRISsetNeighborhoodSize(mris, G_nbrs) ;

  cprints("calculating second fundamental form...", "") ;
  MRIScomputeSecondFundamentalForm(mris) ;
  cprints("", "ok") ;

  // TODO: sort PCDs so that max and min are saved properly
  // TODO: smooth surface, or vectors themselves, to reduce noise?
  // DONE: flip PCDs to be consistent with outward surface normal, cross(k1,k2) == n

  MRISsortPrincipalDirectionsByCurvatures(mris) ;

  MRISorderPrincipalDirectionsConsistentWithNormal(mris) ;

  sprintf(output_filename, "%s_principal1.mgz", output_filename_stem) ;
  fprintf(stdout, "writing principal curvature directions file:  %s\n", output_filename) ;
  MRISwritePrincipalDirection(mris, 1, output_filename) ;

  sprintf(output_filename, "%s_principal2.mgz", output_filename_stem) ;
  fprintf(stdout, "writing principal curvature directions file:  %s\n", output_filename) ;
  MRISwritePrincipalDirection(mris, 2, output_filename) ;


  MRISfree(&mris) ;
  fprintf(0, "\n\n") ;
  exit(0) ;
  return(0) ;  /* for ansi */
}


int
MRISsortPrincipalDirectionsByCurvatures(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  double holder[3] ;
  float k1, k2 ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;

      k1 = v->k1 ;
      k2 = v->k2 ;

      // it is unlikely that the user will want this, but just in case...
      if ( b_sortBySignedPrincipalCurv )
        {
        if ( k1 < k2 )
          {
            // swap the two directions if the curvature corresponding to
            // the first direction is smaller than that of the second
            // direction

            holder[0] = v->e1x ;
            holder[1] = v->e1y ;
            holder[2] = v->e1z ;

            v->e1x = v->e2x ;
            v->e1y = v->e2y ;
            v->e1z = v->e2z ;

            v->e2x = holder[0] ;
            v->e2y = holder[1] ;
            v->e2z = holder[2] ;
          }
        }
      else
        {
        /*
         * in general, for directions to be consistently oriented over
         * large extents of folded surface, the appropriate comparison
         * will be the absolute or unsigned curvature. the sign
         * denotes the orientation of the curving relative to the
         * normal, but here we only care about the radius of
         * curvature. so, e.g., if we consider a chunk of surface from
         * the crest of a gyrus down the bank to the fundus of a
         * sulcus that is flat in the orthogonal direction, the first
         * PCD should always point from the gyrus to the sulcus (maybe
         * with a singularity somewhere in the bank) and the second
         * PCD should always point in the orthogonal direction---which
         * would not be the case if we used the SIGNED curvatures
         * since the sign would flip at the transition from gyrus to
         * sulcus.
         */
        if ( abs(k1) < abs(k2) )
          {
            // swap the two directions if the curvature corresponding to
            // the first direction is smaller than that of the second
            // direction

            holder[0] = v->e1x ;
            holder[1] = v->e1y ;
            holder[2] = v->e1z ;

            v->e1x = v->e2x ;
            v->e1y = v->e2y ;
            v->e1z = v->e2z ;

            v->e2x = holder[0] ;
            v->e2y = holder[1] ;
            v->e2z = holder[2] ;
          }
        }
    }

  return(NO_ERROR) ;
}



#define F_DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define F_CROSS(a,b,d) (d[0]=a[1]*b[2]-b[1]*a[2],       \
                        d[1]=a[2]*b[0]-b[2]*a[0],       \
                        d[2]=a[0]*b[1]-b[0]*a[1])

int
MRISorderPrincipalDirectionsConsistentWithNormal(MRI_SURFACE *mris)
{
  int     vno ;
  VERTEX  *v ;

  double v1[3], v2[3], v3[3], n0[3], dot ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;

      v1[0] = v->e1x ;
      v1[1] = v->e1y ;
      v1[2] = v->e1z ;

      v2[0] = v->e2x ;
      v2[1] = v->e2y ;
      v2[2] = v->e2z ;

      n0[0] = v->nx ;
      n0[1] = v->ny ;
      n0[2] = v->nz ;

      F_CROSS(v1,v2,v3) ;

      // check to see if normal is pointed in same direction as v1 x v2
      dot = F_DOT(n0, v3) ;

      // flip orientation if dot product is negative
      if ( dot < 0 )
        {
          v->e2x = -v->e2x ;
          v->e2y = -v->e2y ;
          v->e2z = -v->e2z ;
        }
    }

  return(NO_ERROR) ;
}


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-version")) {
    print_version() ;
  } else if (!stricmp(option, "-all-info")) {
    // do nothing
    exit(0) ;
  }
  else switch (toupper(*option)) {
  case 'A':
    b_sortBySignedPrincipalCurv = 0 ;
    break ;
  case 'S':
    b_sortBySignedPrincipalCurv = 1 ;
    break ;
  case '?':
    print_help() ;
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(int code) {
  print_usage() ;
  exit(code) ;
}

static void
print_usage(void) {
  fprintf(stdout, "usage:  %s surface [outputstem]\n", Progname) ;
}

static void
print_help(void) {
  // borrowed from mris_curvature_stats
  char  pch_synopsis[65536] ;

  sprintf(pch_synopsis, "\n\
 \n\
 \n\
    NAME \n\
 \n\
          mris_curvature_vectors \n\
 \n\
    SYNOPSIS \n\
 \n\
          mris_curvature_vectors [OPTIONS] surface [outputstem] \n\
 \n\
    DESCRIPTION \n\
 \n\
         wrapper around MRIScomputeSecondFundamentalForm, outputs \n\
         surface normals and the directions of the principal curvatures. \n\
         output files are three-frame MGZ files containing the vectors \n\
         for each vertex in the input surface file. \n\
 \n\
    OPTIONS \n\
 \n\
       -A   order Principal Curvatures by UNSIGNED (abs) curvature values [default] \n\
       -S   order Principal Curvatures by SIGNED curvature values \n\
       -U   print usage \n\
       -?   print this help message \n\
       --version  print version control info \n\
 \n\
    EXAMPLES \n\
 \n\
        the command \n\
 \n\
            mris_curvature_vectors rh.white \n\
 \n\
        will output files rh.white_normals.mgz, rh.white_principal1.mgz, \n\
        and rh.white_principal2.mgz. \n\
 \n\
 \n\
    contact: jonathan polimeni <jonp@nmr.mgh.harvard.edu> \n\
\n") ;

  fprintf(stdout, "%s", pch_synopsis) ;
  exit(0) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(0) ;
}

