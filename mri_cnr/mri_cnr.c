/**
 * @file  mri_cnr.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:06 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
#include "mri.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri_conform.h"
#include "version.h"

static char vcid[] = "$Id: mri_cnr.c,v 1.6 2006/12/29 02:09:06 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static double compute_volume_cnr(MRI_SURFACE *mris, MRI *mri) ;
char *Progname ;

int
main(int argc, char *argv[]) {
  char        **av, *mri_name,  fname[STRLEN], *hemi, *path ;
  int         ac, nargs, i, j ;
  MRI         *mri, *mri_template = NULL, *mri_tmp ;
  MRI_SURFACE *mris ;
  double      cnr_total, cnr = 0.0 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_cnr.c,v 1.6 2006/12/29 02:09:06 nicks Exp $", "$Name:  $");
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

  path = argv[1] ;
  for (cnr_total = 0.0, j = 0 ; j <= 1 ; j++) {
    if (j == 0)
      hemi = "lh" ;
    else
      hemi = "rh" ;
    sprintf(fname, "%s/%s.white", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    sprintf(fname, "%s/%s.pial", path, hemi) ;
    if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;

    for (i = 2 ; i < argc ; i++) {
      mri_name = argv[i] ;
      if (!j)
        printf("processing MRI volume %s...\n", mri_name) ;
      mri = MRIread(mri_name) ;
      if (!mri)
        ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, mri_name) ;
      if (!j && !FZERO(mri->tr))
        printf("TR = %2.1f msec, flip angle = %2.0f degrees, TE = %2.1f msec\n",
               mri->tr, DEGREES(mri->flip_angle), mri->te) ;

      if (!mri_template)
        mri_template = MRIconform(mri) ;
      mri_tmp = MRIresample(mri, mri_template, SAMPLE_NEAREST) ;
      MRIfree(&mri) ;
      mri = mri_tmp ;

      if (!j)
        cnr = compute_volume_cnr(mris, mri) ;
      else {
        cnr += compute_volume_cnr(mris, mri) ;
        cnr /= 2.0 ;
        printf("CNR = %2.3f\n", cnr) ;
      }

      MRIfree(&mri) ;
    }
    cnr_total += cnr ;
    MRISfree(&mris) ;
  }
  printf("total CNR = %2.3f\n", cnr/(double)((argc-2))) ;

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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <surf directory> <vol 1> <vol 2> ...\n",
          Progname) ;
}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program will compute the gray/white/csf contrast-to-noise ratio\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
/*
  orig vertex positions = white
  current vertex positions = pial
*/
static double
compute_volume_cnr(MRI_SURFACE *mris, MRI *mri) {
  double  gray_white_cnr, gray_csf_cnr, gray_var, white_var, csf_var, gray_mean, white_mean, csf_mean ;
  float   thickness ;
  Real    x, y, z, gray, white, csf ;
  int     vno ;
  VERTEX  *v ;

  MRIScomputeMetricProperties(mris) ;
  gray_mean = gray_var = white_mean = white_var = csf_mean = csf_var = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    //MRIworldToVoxel(mri, v->x+v->nx, v->y+v->ny, v->z+v->nz, &x, &y, &z) ;
    MRIsurfaceRASToVoxel(mri, v->x+v->nx, v->y+v->ny, v->z+v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &csf) ;
    thickness = v->curv*.5 ;
    // MRIworldToVoxel(mri, v->x-thickness*v->nx, v->y-thickness*v->ny, v->z-thickness*v->nz, &x, &y, &z) ;
    MRIsurfaceRASToVoxel(mri, v->x-thickness*v->nx, v->y-thickness*v->ny, v->z-thickness*v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &gray) ;

    gray_var += (gray*gray) ;
    gray_mean += gray ;
    csf_var += (csf*csf) ;
    csf_mean += csf ;
  }

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  for (white_var = gray_white_cnr = 0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    // MRIworldToVoxel(mri, v->x-v->nx, v->y-v->ny, v->z-v->nz, &x, &y, &z) ;
    MRIsurfaceRASToVoxel(mri, v->x-v->nx, v->y-v->ny, v->z-v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &white) ;
    thickness = v->curv*.5 ;
    // MRIworldToVoxel(mri, v->x+thickness*v->nx, v->y+thickness*v->ny, v->z+thickness*v->nz, &x, &y, &z) ;
    MRIsurfaceRASToVoxel(mri, v->x+thickness*v->nx, v->y+thickness*v->ny, v->z+thickness*v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &gray) ;

    gray_var += (gray*gray) ;
    gray_mean += gray ;
    white_var += (white*white) ;
    white_mean += white ;
  }

  white_mean /= (double)mris->nvertices ;
  csf_mean /= (double)mris->nvertices ;
  gray_mean /= (double)mris->nvertices*2.0 ;

  white_var = white_var / (double)mris->nvertices - white_mean*white_mean ;
  gray_var = gray_var / ((double)mris->nvertices*2.0) - gray_mean*gray_mean ;
  csf_var = csf_var / (double)mris->nvertices - csf_mean*csf_mean ;

  printf("\twhite = %2.1f+-%2.1f, gray = %2.1f+-%2.1f, csf = %2.1f+-%2.1f\n",
         white_mean, sqrt(white_var), gray_mean, sqrt(gray_var),
         csf_mean, sqrt(csf_var)) ;

  gray_white_cnr = SQR(gray_mean - white_mean) / (gray_var+white_var) ;
  gray_csf_cnr = SQR(gray_mean - csf_mean) / (gray_var+csf_var) ;

  printf("\tgray/white CNR = %2.3f, gray/csf CNR = %2.3f\n",
         gray_white_cnr, gray_csf_cnr) ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  return((gray_white_cnr + gray_csf_cnr)/2.0) ;
}

