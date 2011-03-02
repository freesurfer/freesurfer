/**
 * @file  mri_relabel_hypointensities.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
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
#include "fio.h"
#include "mrishash.h"
#include "cma.h"

static char vcid[] = "$Id: mri_relabel_hypointensities.c,v 1.7 2011/03/02 00:04:24 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int relabel_hypointensities(MRI *mri, MRI_SURFACE *mris, int right) ;
static int relabel_hypointensities_neighboring_gray(MRI *mri) ;


char *Progname ;

static char *surf_name = "white" ;

int
main(int argc, char *argv[]) {
  char          **av, *hemi, fname[STRLEN],
  *in_aseg_name, *out_aseg_name, *surf_dir ;
  int           ac, nargs, h ;
  MRI_SURFACE   *mris ;
  MRI           *mri_aseg ;

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

  in_aseg_name = argv[1] ;
  surf_dir = argv[2] ;
  out_aseg_name = argv[3] ;

  mri_aseg = MRIread(in_aseg_name) ;
  if (!mri_aseg)
    ErrorExit(ERROR_NOFILE, "%s: could not read input segmentation %s", Progname, in_aseg_name) ;


  for (h = 0 ; h <= 1 ; h++) {
    if (h == 0)
      hemi = "lh" ;
    else
      hemi = "rh" ;
    sprintf(fname, "%s/%s.%s", surf_dir, hemi, surf_name)  ;
    printf("reading input surface %s...\n", fname) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;

    printf("relabeling %s hypointensities...\n", hemi) ;
    relabel_hypointensities(mri_aseg, mris, h) ;
    MRISfree(&mris) ;
  }
  relabel_hypointensities_neighboring_gray(mri_aseg) ;

  MRIwrite(mri_aseg, out_aseg_name) ;
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
  else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (toupper(*option)) {
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input aseg> <surface  directory> <output aseg>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program measures a variety of anatomical properties\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "-l <label file>              - limit calculations to specified "
          "label\n") ;
  fprintf(stderr,
          "-a <annotation file>         - compute properties for each label\n"
          "                               in the annotation file separately"
          "\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
relabel_hypointensities(MRI *mri, MRI_SURFACE *mris, int right) {
  int   x, y, z, label, changed ;
  MRIS_HASH_TABLE *mht ;
  VERTEX           *v ;
  float            dx, dy, dz, dot, dist ;
  Real             xw, yw, zw ;

  mht = MHTfillVertexTableRes(mris,NULL, CURRENT_VERTICES, 8.0f) ;
  for (changed = x = 0 ; x < mri->width ; x++) {
    for (y = 0 ; y < mri->height ; y++) {
      for (z = 0 ; z < mri->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri, x, y, z) ;
        if (label == Left_WM_hypointensities)
          MRIvox(mri, x, y, z) = WM_hypointensities ;
        else if (label == Right_WM_hypointensities)
          MRIvox(mri, x, y, z) = WM_hypointensities ;
        if ((!right && (label != Left_Cerebral_Cortex)) ||
            (right && (label != Right_Cerebral_Cortex)))
          continue ;

        // MRIvoxelToWorld(mri, x, y, z, &xw, &yw, &zw) ;
        MRIvoxelToSurfaceRAS(mri, x, y, z, &xw, &yw, &zw);
        v = MHTfindClosestVertexInTable(mht, mris, xw, yw, zw, 0) ;
        if (v == NULL)  /* no vertices within range - assume it is hypointensity */
        {
          dot = -1 ;
          dist = 1000 ;
        } else {
          dx = xw - v->x ;
          dy = yw - v->y ;
          dz = zw - v->z ;
          dot = v->nx*dx + v->ny*dy + v->nz*dz ;
          dist = sqrt(dx*dx+dy*dy+dz*dz) ;
        }
        if (dot < 0 && dist > 1) {
          changed++ ;
          MRIvox(mri, x, y, z) = WM_hypointensities ;
        }
      }
    }
  }

  printf("%d voxels changed to hypointensity...\n", changed) ;
  MHTfree(&mht) ;
  return(NO_ERROR) ;
}
int
relabel_hypointensities_neighboring_gray(MRI *mri) {
  int    x, y, z, label, changed, i ;
  MRI    *mri_tmp = NULL ;

  for (changed = i = 0 ; i < 2 ; i++) {
    mri_tmp = MRIcopy(mri, mri_tmp) ;
    for (x = 0 ; x < mri->width ; x++) {
      for (y = 0 ; y < mri->height ; y++) {
        for (z = 0 ; z < mri->depth ; z++) {
          label = MRIvox(mri_tmp, x, y, z) ;
          if (label != WM_hypointensities)
            continue ;
          if (MRIneighbors(mri_tmp, x, y, z, Left_Cerebral_Cortex) > 0) {
            MRIvox(mri, x, y, z) = Left_Cerebral_Cortex ;
            changed++ ;
          } else  if (MRIneighbors(mri_tmp, x, y, z, Right_Cerebral_Cortex) > 0) {
            MRIvox(mri, x, y, z) = Right_Cerebral_Cortex ;
            changed++ ;
          }
        }
      }
    }
  }
  printf("%d hypointense voxels neighboring cortex changed\n", changed) ;
  return(NO_ERROR) ;
}
