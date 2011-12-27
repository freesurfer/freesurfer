/**
 * @file  mris_compute_layer_intensities.c
 * @brief compute the % of gm layers 1-6, wm and CSF in each voxel in a volume
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/12/27 18:43:06 $
 *    $Revision: 1.2 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "mrisurf.h"
#include "registerio.h"
#include "cma.h"

#define LAYER1         1
#define LAYER2         2
#define LAYER3         3
#define LAYER4         4
#define LAYER5         5
#define LAYER6         6
#define WM_VAL         7
#define CSF_VAL        8
#define SUBCORT_GM_VAL 9

#define NLAYERS        6
#define NLABELS        9
#define MAX_LAYERS     50

static int nlayers = NLAYERS ;
static char *LAMINAR_NAME = "gwdist";

static char *subject_name = NULL ;
static char *hemi = "lh" ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static int Gwhalf = 3 ;
static MRI *compute_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, MRI_SURFACE **mris, int nlayers, int whalf0, MRI *mri_layer_intensities) ;

int
main(int argc, char *argv[]) 
{
  char   **av ;
  int    ac, nargs ;
  int    msec, minutes, seconds, i ;
  struct timeb start ;
  MRI_SURFACE *mris[MAX_LAYERS];
  char        fname[STRLEN] ;
  MRI         *mri_intensities, *mri_volume_fractions, *mri_layer_intensities ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_compute_layer_intensities.c,v 1.2 2011/12/27 18:43:06 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  mri_intensities = MRIread(argv[1]) ;
  if (mri_intensities == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load intensity volume from %s", Progname, argv[1]) ;
  mri_volume_fractions = MRIread(argv[2]) ;
  if (mri_volume_fractions == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load volume fractions from %s", Progname, argv[2]) ;
  for (i = 0 ; i <= nlayers ; i++)
  {
    sprintf(fname, "%s%d", argv[3], i) ;
    printf("reading laminar surface %s\n", fname) ;
    mris[i] = MRISread(fname) ;
    if (mris[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load surface from %s", Progname, fname) ;
  }
  
  mri_layer_intensities = compute_layer_intensities(mri_intensities, mri_volume_fractions, mris, nlayers, Gwhalf, NULL) ;
  printf("writing layer intensities to %s\n", argv[4]) ;
  MRIwrite(mri_layer_intensities, argv[4]) ;

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("layer intensities calculation took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
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
  if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  } else if (!stricmp(option, "nlayers")) {
    nlayers = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d input layers for laminar analysis\n", nlayers) ;
  } else if (!stricmp(option, "rh") || !stricmp(option, "lh")) {
    hemi = option ;
    printf("processing %s hemisphere\n", hemi) ;
  } else switch (toupper(*option)) {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    break ;
  case 'W':
    Gwhalf = atoi(argv[2]) ;
    printf("using half window size = %d\n", Gwhalf) ;
    nargs = 1 ;
    break ;
    break ;
  case 'N':
    LAMINAR_NAME = argv[2] ;
    printf("using %s as layer surface name\n", LAMINAR_NAME) ;
    nargs = 1 ;
    break ;
  case 'S':
    subject_name = argv[2] ;
    nargs = 1 ;
    printf("overriding subject name in .dat file with %s\n", subject_name) ;
    break ;
  case '?':
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
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <input intensity volume> <layer volume fractions file> <input surface> <output overlay>\n",
         Progname) ;
  printf(
         "\t\n");
  exit(code) ;
}

static MRI *
compute_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, MRI_SURFACE **mris, int nlayers, int whalf0, MRI *mri_layer_intensities)
{
  MATRIX  *mF, *mL, *mI, *mFinv ;
  int     whalf, vno, n, t1, t2, nvals, xvi, yvi, zvi, n1, max_vals ;
  VERTEX  *v ;
  double  step, xs, ys, zs, xv, yv, zv ;
  MRI     *mri_visited ;

  step = mri_intensities->xsize/2 ;
  if (mri_layer_intensities == NULL)
    mri_layer_intensities = MRIallocSequence(mris[0]->nvertices, 1, 1, MRI_FLOAT, nlayers) ;

  mri_visited = MRIcloneDifferentType(mri_intensities, MRI_UCHAR) ;
  for (n = 0 ; n <= nlayers ; n++)
  {
    MRISsetNeighborhoodSize(mris[n],3);
    MRIScomputeSecondFundamentalForm(mris[n]) ;
  }

  mL = MatrixAlloc(nlayers+2, 1, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris[0]->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;

    whalf = whalf0 ;
    do
    {
      // find how many unique vals are within whalf
      for (max_vals = 0, n = 0 ; n < nlayers ; n++)  // in normal direction
      {
	v = &mris[n]->vertices[vno] ;
	
	for (t1 = -whalf ; t1 <= whalf ; t1++)    // in tangent plane
	  for (t2 = -whalf ; t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y;
	    zs = v->z + step*t1*v->e1z + step*t2*v->e2z;
	    MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	      continue ;
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) == 0)
	    {
	      max_vals++ ;
	      MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;
	    }
	  }
      }

      if (max_vals == 0)
	continue ;
      mF = MatrixAlloc(max_vals, nlayers+2, MATRIX_REAL) ;
      mI = MatrixAlloc(max_vals, 1, MATRIX_REAL) ;

      // unmark the vertices
      for (n = 0 ; n <= nlayers ; n++)  // in normal direction
      {
	v = &mris[n]->vertices[vno] ;
	
	for (t1 = -whalf ; t1 <= whalf ; t1++)    // in tangent plane
	  for (t2 = -whalf ; t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y;
	    zs = v->z + step*t1*v->e1z + step*t2*v->e2z;
	    MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	      continue ;
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 0) ;
	  }
      }
      // build the matrices
      for (nvals = 0, n = 0 ; n <= nlayers ; n++)  // in normal direction
      {
	v = &mris[n]->vertices[vno] ;
	
	for (t1 = -whalf ; t1 <= whalf ; t1++)    // in tangent plane
	  for (t2 = -whalf ; t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y;
	    zs = v->z + step*t1*v->e1z + step*t2*v->e2z;
	    MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	      continue ;
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    if (nvals >= max_vals)
	      continue ;
	    if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) == 0)
	    {
	      MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;
	      *MATRIX_RELT(mI, nvals+1, 1) = MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ;
	      for (n1 = 0 ; n1 < nlayers ; n1++)
		*MATRIX_RELT(mF, nvals+1, n1+1) = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, nlayers-n1) ;
	      // add wm and csf estimates
	      *MATRIX_RELT(mF, nvals+1, nlayers+1) = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, nlayers+1) ;
	      *MATRIX_RELT(mF, nvals+1, nlayers+2) = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, nlayers+2) ;
	      nvals++ ;
	    }
	  }
      }
      // unmark the vertices
      for (n = 0 ; n <= nlayers ; n++)  // in normal direction
      {
	v = &mris[n]->vertices[vno] ;
	
	for (t1 = -whalf ; t1 <= whalf ; t1++)    // in tangent plane
	  for (t2 = -whalf ; t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y;
	    zs = v->z + step*t1*v->e1z + step*t1*v->e2z;
	    MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	      continue ;
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 0) ;
	  }
      }
      mFinv = MatrixPseudoInverse(mF, NULL) ;
      MatrixMultiply(mFinv, mI, mL) ;
      for (n = 0 ; n < nlayers ; n++) // fill in outputs
	MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n, *MATRIX_RELT(mL, n+1, 1)) ;

      MatrixFree(&mF) ; MatrixFree(&mI) ; MatrixFree(&mFinv) ;
    } while (whalf == 0) ;
  }

  MatrixFree(&mL) ;
  MRIfree(&mri_intensities) ;
  return(mri_layer_intensities) ;
}

