/**
 * @brief compute the % of gm layers 1-6, wm and CSF in each voxel in a volume
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021
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
#include "mrinorm.h"


#define NLAYERS        6
#define MAX_LAYERS     50

static int nlayers = NLAYERS ;

static double vfrac_thresh = -1 ;

static char *subject_name = NULL ;
static char *hemi = NULL ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int FS_names = 0 ;
static int Gwhalf = 3 ;
static MRI *compute_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, 
				      MRI_SURFACE **mris, int nlayers, int whalf0, MRI *mri_layer_intensities,
				      int curv_bins) ;
static MRI *compute_thresholded_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, 
						  MRI_SURFACE **mris, int nlayers, int whalf0, 
						  MRI *mri_layer_intensities,
						  int curv_bins, double vfrac_thresh) ;

static int curv_bins = 0 ;
int
main(int argc, char *argv[]) 
{
  char   **av ;
  int    ac, nargs ;
  int    msec, minutes, seconds, i ;
  Timer start ;
  MRI_SURFACE *mris[MAX_LAYERS];
  char        fname[STRLEN] ;
  MRI         *mri_intensities, *mri_volume_fractions, *mri_layer_intensities ;

  nargs = handleVersionOption(argc, argv, "mris_compute_layer_intensities");
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

  start.reset() ;

  if (hemi == NULL)
    ErrorExit(ERROR_BADPARM, "%s: must specify -rh or -lh", Progname) ;

  mri_intensities = MRIread(argv[1]) ;
  if (mri_intensities == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load intensity volume from %s", Progname, argv[1]) ;
  mri_volume_fractions = MRIread(argv[2]) ;
  if (mri_volume_fractions == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load volume fractions from %s", Progname, argv[2]) ;
  if (mri_volume_fractions->nframes < nlayers+1)
    ErrorExit(ERROR_BADFILE, "%s: volume fraction input has fewer frames (%d) than needed (%d)\n", mri_volume_fractions->nframes, nlayers+1);
  if (FS_names && subject_name == NULL)
    ErrorExit(ERROR_UNSUPPORTED, "%s: if specifying FS_names must use -s <subject>", Progname) ;
  for (i = 0 ; i <= nlayers ; i++)
  {
    if (FS_names && (i == 0 || i == nlayers))
    {
      char *sdir = getenv("SUBJECTS_DIR") ;
      if (i == 0)
	sprintf(fname, "%s/%s/surf/%s.white", sdir, subject_name,hemi) ;
      else
	sprintf(fname, "%s/%s/surf/%s.pial", sdir, subject_name,hemi) ;
    }
    else
    {
      if (i == 10)
	argv[3][strlen(argv[3])-1] = 0 ;  // make it layer010 not layer0010
      sprintf(fname, "%s%3.3d", argv[3], i) ;
    }
    printf("reading laminar surface %s\n", fname) ;
    mris[i] = MRISread(fname) ;
    if (mris[i] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load surface from %s", Progname, fname) ;
  }

  if (vfrac_thresh > 0)
  mri_layer_intensities = 
    compute_thresholded_layer_intensities(mri_intensities, mri_volume_fractions, mris, nlayers, 
					  Gwhalf, NULL, curv_bins, vfrac_thresh) ;
  else
    mri_layer_intensities = compute_layer_intensities(mri_intensities, mri_volume_fractions, mris, nlayers, 
						      Gwhalf, NULL, curv_bins) ;
  {
    sprintf(fname, "wsize.mgz") ;
    printf("writing half window sizes to %s\n", fname) ;
    MRISwriteValues(mris[0], fname) ;
  }

  printf("writing layer intensities to %s\n", argv[4]) ;
  MRIwrite(mri_layer_intensities, argv[4]) ;

  msec = start.milliseconds() ;
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
  } else if (!stricmp(option, "FS_names")) {
    printf("using standard FS names white and pial\n") ;
    FS_names = 1 ;
  } else if (!stricmp(option, "thresh")) {
    vfrac_thresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("only using voxels with at least %2.2f volume fraction to estimate intensities\n",
	   vfrac_thresh) ;
  } else if (!stricmp(option, "curv")) {
    curv_bins = atoi(argv[2]) ;
    nargs = 1 ;
    printf("binning curvature into %d bins\n", curv_bins) ;
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

#define CURV_THRESH .02

static MRI *
compute_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, MRI_SURFACE **mris, int nlayers, 
			  int whalf0, MRI *mri_layer_intensities, int curv_bins)
{
  MATRIX  *mF, *mL, *mI, *mFinv ;
  int     whalf, vno, n, nd, t1, t2, nvals, xvi, yvi, zvi, n1, out_of_fov ;
  VERTEX  *v ;
  double  step, xs, ys, zs, xv, yv, zv, vfrac, nx, ny, nz ;
  MRI     *mri_visited, *mri_curv_bins ;
  double  bin_size, Hmin, whalf_total = 0, val, vfrac_thresh ;
  int     bin0, bin, found_layer[MAX_LAYERS], estimable, nwindows = 0 ;

  printf("computing partial volume corrected layer intensities\n") ;
/*
  try to find voxels with a lot of each layer, but the more layers there are the less likely we are to be able
  to do so. The higher res the data the higher the threshold should be (as the more likely we are to be able
  to find voxels with a lot of one layer.
  xsize = 1,  nlayers = 2 --> vfrac = .75
  xsize = .5,  nlayers = 2 --> vfrac = 1
  xsize = .5, nlayers = 6 --> vfrac = .667
  xsize = .5, nlayers = 12 --> vfrac = .583
*/

  vfrac_thresh = .4 + .5*1/(mri_intensities->xsize*nlayers) ;
  if (vfrac_thresh > .995)
    vfrac_thresh = .995 ;

  step = mri_intensities->xsize/4 ;
  if (mri_layer_intensities == NULL)
    mri_layer_intensities = MRIallocSequence(mris[0]->nvertices, 1, 1, MRI_FLOAT, nlayers+2) ;

  mri_visited = MRIcloneDifferentType(mri_intensities, MRI_UCHAR) ;
  for (n = 0 ; n <= nlayers ; n++)
  {
    MRISsaveVertexPositions(mris[n], TMP_VERTICES) ;
    MRISaverageVertexPositions(mris[n], 5) ;
    MRISsetNeighborhoodSizeAndDist(mris[n],3);
    MRIScomputeSecondFundamentalForm(mris[n]) ;
    MRISsmoothCurvatures(mris[n], 5) ;
    MRISrestoreVertexPositions(mris[n], TMP_VERTICES) ;
  }
  mri_curv_bins = MRIcloneDifferentType(mri_intensities, MRI_UCHAR) ;
  if (curv_bins > 1)
  {
    MRI     *mri_ctrl ;
    int     x, y, z ;

    mri_ctrl = MRIclone(mri_curv_bins, NULL) ;
    printf("constructing %d curvature bins\n", curv_bins) ;
    Hmin = mris[0]->Hmin ;
    bin_size = (mris[0]->Hmax - mris[0]->Hmin) / (float)(curv_bins-1) ;
    for (vno = 0 ; vno < mris[0]->nvertices ; vno++)
    {
      v = &mris[0]->vertices[vno] ;
      bin = nint((v->H - Hmin) / bin_size) ;
      bin = v->H < -CURV_THRESH ? 1 : (v->H > CURV_THRESH ? 3 : 2) ;
      for (n = 0 ; n <= nlayers ; n++)
      {
	v = &mris[n]->vertices[vno] ;
	MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, v->x, v->y, v->z, &xv,&yv,&zv) ;
	if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	  continue ;
	xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	MRIsetVoxVal(mri_ctrl, xvi, yvi, zvi, 0, CONTROL_MARKED) ;
	MRIsetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0, bin) ;
      }
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_ctrl, "ctrl.mgz") ;
      MRIwrite(mri_curv_bins, "cbins.mgz") ;
    }
    MRIbuildVoronoiDiagram(mri_curv_bins, mri_ctrl, mri_curv_bins);
    /*
      set bin to unused value for each spot that has no non-zero fractions so that that voxel
      it isn't used below
    */
    for (x = 0 ; x < mri_curv_bins->width ; x++)
      for (y = 0 ; y < mri_curv_bins->height ; y++)
	for (z = 0 ; z < mri_curv_bins->depth ; z++)
	{
	  for (nvals = 0, n = 0 ; n <= nlayers ; n++)  // in normal direction
	  {
	    if (MRIgetVoxVal(mri_volume_fractions, x, y, z, n) > 0)
	    {
	      nvals++ ;
	      break ;
	    }
	  }
	  if (nvals == 0)
	  {
	    DiagBreak() ;
	    MRIsetVoxVal(mri_curv_bins, x, y, z, 0, -10000) ;
	  }
	}

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_curv_bins, "cbins.mgz") ;
    }
    MRIfree(&mri_ctrl) ;
  }
  else
    bin_size = 0 ;

#define MAX_VALS 2000
  mF = MatrixAlloc(MAX_VALS, nlayers+2, MATRIX_REAL) ;
  mI = MatrixAlloc(MAX_VALS, 1, MATRIX_REAL) ;
  mL = MatrixAlloc(nlayers+2, 1, MATRIX_REAL) ; // wm+nlayers+csf+ (subcortical gray maybe)
  for (vno = 0 ; vno < mris[0]->nvertices ; vno++)
  {
    if ((vno % 500) == 0)
      printf("processing vno %d of %d (%2.1f%%)\n", vno, mris[0]->nvertices, 100.0f*vno/mris[0]->nvertices) ;
    v = &mris[0]->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (FEQUAL(mris[0]->vertices[vno].x, mris[nlayers]->vertices[vno].x) &&
	FEQUAL(mris[0]->vertices[vno].y, mris[nlayers]->vertices[vno].y) &&
	FEQUAL(mris[0]->vertices[vno].z, mris[nlayers]->vertices[vno].z))
    {
      DiagBreak() ;
#if 0
      continue ;   // pial and white in same place - vertex is not cortical and can't be estimated
#endif
    }
    if (bin_size == 0)
      bin0 = 0 ;
    else
    {
      bin0 = (v->H - Hmin) / (bin_size-1) ;
      bin0 = v->H < -CURV_THRESH ? 1 : (v->H > CURV_THRESH ? 3 : 2) ;
    }

    whalf = whalf0/step ;
    out_of_fov = 0 ;
    do
    {
      // build the matrices
      memset(found_layer, 0, sizeof(found_layer)) ;
      for (nvals = 0, n = 0 ; nvals < MAX_VALS && n <= nlayers ; n++)  // in normal direction
      {
	for (t1 = -whalf ; nvals < MAX_VALS && t1 <= whalf ; t1++)    // in tangent plane
	{
	  for (t2 = -whalf ; nvals < MAX_VALS && t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y;
	    zs = v->z + step*t1*v->e1z + step*t2*v->e2z;
	    MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	    {
	      out_of_fov = 1 ;
	      break ;
	    }
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    bin = MRIgetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0) ;
	    if (bin != bin0)  
	      continue ;  // not the same curvature
	    if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) == 0)  // new voxel found
	    {
	      MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;
	      val = MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ; 
	      if (val < 0)
		DiagBreak() ;
	      *MATRIX_RELT(mI, nvals+1, 1) = val ;
	      for (n1 = 0 ; n1 <= nlayers+1 ; n1++)
	      {
		vfrac = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, n1) ;
		if (vfrac >= vfrac_thresh)
		  found_layer[n1] = 1 ;
		*MATRIX_RELT(mF, nvals+1, n1+1) = vfrac ;
		if (vno == Gdiag_no)
		  printf("layer %d v %d: val %d with fraction %f found at (%d, %d, %d)\n",
			 n1, vno, (int)val, vfrac, xvi, yvi, zvi) ;
	      }
	      // add wm and csf estimates
	      *MATRIX_RELT(mF, nvals+1, 1) = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, 0) ; // wm
	      *MATRIX_RELT(mF, nvals+1, nlayers+2) = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, nlayers+1) ;//csf
	      nvals++ ;
	    }
	  }
	}
      }
/*
  now go whalf into the wm and whalf outside the last surface to get good estimates of WM and CSF
*/
      for (n = 0 ; nvals < MAX_VALS && n <= nlayers ; n += nlayers)
      {
	v = &mris[n]->vertices[vno] ;
	if (n == 0)  // look inwards from inside surface
	{
	  nx = -v->nx ; ny = -v->ny ; nz = -v->nz ;
	}
	else   // look outwards from outside surface
	{
	  nx = v->nx ; ny = v->ny ; nz = v->nz ;
	}

	for (nd = -1 ; nvals < MAX_VALS && nd <= whalf+1 ; nd++)
	{
	  for (t1 = -whalf ; nvals < MAX_VALS && t1 <= whalf ; t1++)    // in tangent plane
	  {
	    for (t2 = -whalf ; nvals < MAX_VALS && t2 <= whalf ; t2++)
	    {
	      xs = v->x + step*t1*v->e1x + step*t2*v->e2x + nd*step*nx;
	      ys = v->y + step*t1*v->e1y + step*t2*v->e2y + nd*step*ny;
	      zs = v->z + step*t1*v->e1z + step*t2*v->e2z + nd*step*nz;
	      MRISsurfaceRASToVoxelCached(mris[nlayers], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	      if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	      {
		out_of_fov = 1 ;
		break ;
	      }
	      xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	      bin = MRIgetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0) ;
	      if (bin != bin0)  
		continue ;  // not the same curvature
	      if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) == 0)  // new voxel found
	      {
		MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;
		val = MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ;
		if (val < 0)
		  DiagBreak() ;
		*MATRIX_RELT(mI, nvals+1, 1) = val ;
		for (n1 = 0 ; n1 <= nlayers+1 ; n1++)
		{
		  vfrac = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, n1) ;
		  if (vfrac >= vfrac_thresh)
		    found_layer[n1] = 1 ;
		  *MATRIX_RELT(mF, nvals+1, n1+1) = vfrac ;
		}
		nvals++ ;
	      }
	    }
	  }
	}
      }

      if (vno == Gdiag_no)
	DiagBreak() ;
      mF->rows = nvals ; mI->rows = nvals ;  // temporary - will be reset below (avoids allocating frequently)
      MRIclear(mri_visited) ;
      for (estimable = 1, n1 = 0 ; n1 <= nlayers+1 ; n1++)
	if (found_layer[n1] == 0)
	{
	  estimable = 0 ;
	  break ;
	}
      if (estimable)
      {
	if (MatrixNSConditionNumber(mF) > 100000) // can't reliably estimate parameters
	  estimable = 0 ;
	else  // matrix is well-conditioned and we have found voxels with majority of each layer
	{
	  mFinv = MatrixPseudoInverse(mF, NULL) ;
	  if (mFinv == NULL)  // shouldn't happen since we check condition # above
	    estimable = 0 ;
	  else   // inversion ok - compute parameter estimates
	  {
	    MatrixMultiply(mFinv, mI, mL) ;
	    for (n = 0 ; n <= nlayers+1 ; n++) // fill in outputs
	    {
	      val = *MATRIX_RELT(mL, n+1, 1) ;
	      if (val < 0 || val > 150)
		DiagBreak() ;
	      if (vno == Gdiag_no)
		printf("v %d layer %d: %2.1f\n", vno, n, val) ;
	      MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n, val) ;
	    }
	    MatrixFree(&mFinv) ;
	  }
	}
      }
      whalf++ ;  // if couldn't estimate inverse use a bigger neighborhood and more data
      mF->rows = MAX_VALS ; mI->rows = MAX_VALS ;
      if (whalf > 10*(whalf0+1)/step || (!estimable && nvals >= MAX_VALS-1))
	break ;  // failed to find enough vals
    } while (!estimable) ;      // most values are 0 - make sure we have enough to estimate
    if (estimable)
    {
      whalf_total += (step*whalf) ; nwindows++ ;
      mris[0]->vertices[vno].val = whalf/step ;
    }

    mF->rows = MAX_VALS ; mI->rows = MAX_VALS ;  // in case loop was broken out of
  }

  printf("%d vertices estimated, average window size = %2.1fmm\n", nwindows, whalf_total/nwindows) ;
  MatrixFree(&mF) ; MatrixFree(&mI) ; 
  MatrixFree(&mL) ;
  return(mri_layer_intensities) ;
}

static MRI *
compute_thresholded_layer_intensities(MRI *mri_intensities, MRI *mri_volume_fractions, 
				      MRI_SURFACE **mris, int nlayers, 
				      int whalf0, MRI *mri_layer_intensities, int curv_bins,
				      double vfrac_thresh)
{
  int     whalf, vno, n, t1, t2, xvi, yvi, zvi, n1, out_of_fov, num_found[MAX_LAYERS], xv0, yv0, zv0 ;
  VERTEX  *v ;
  double  step, xs, ys, zs, xv, yv, zv, vfrac, nx, ny, nz ;
  MRI     *mri_curv_bins ;
  double  bin_size, Hmin, val ;
  int     bin0, bin, nd ;
  MRI     *mri_visited ;

  printf("computing thresholded layer intensities\n") ;
  step = mri_intensities->xsize/4 ;
  if (mri_layer_intensities == NULL)
    mri_layer_intensities = MRIallocSequence(mris[0]->nvertices, 1, 1, MRI_FLOAT, nlayers+2) ;

  mri_visited = MRIcloneDifferentType(mri_volume_fractions, MRI_UCHAR) ;
  for (n = 0 ; n <= nlayers ; n++)
  {
    MRISsaveVertexPositions(mris[n], TMP_VERTICES) ;
    MRISaverageVertexPositions(mris[n], 5) ;
    MRISsetNeighborhoodSizeAndDist(mris[n],3);
    MRIScomputeSecondFundamentalForm(mris[n]) ;
    MRISsmoothCurvatures(mris[n], 5) ;
    MRISrestoreVertexPositions(mris[n], TMP_VERTICES) ;
  }
  mri_curv_bins = MRIcloneDifferentType(mri_intensities, MRI_UCHAR) ;
  if (curv_bins > 1)
  {
    MRI     *mri_ctrl ;

    mri_ctrl = MRIclone(mri_curv_bins, NULL) ;
    printf("constructing %d curvature bins\n", curv_bins) ;
    Hmin = mris[0]->Hmin ;
    bin_size = (mris[0]->Hmax - mris[0]->Hmin) / (float)(curv_bins-1) ;
    for (vno = 0 ; vno < mris[0]->nvertices ; vno++)
    {
      v = &mris[0]->vertices[vno] ;
      bin = nint((v->H - Hmin) / bin_size) ;
      bin = v->H < -CURV_THRESH ? 1 : (v->H > CURV_THRESH ? 3 : 2) ;
      for (n = 0 ; n <= nlayers ; n++)
      {
	v = &mris[n]->vertices[vno] ;
	MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, v->x, v->y, v->z, &xv,&yv,&zv) ;
	if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	  continue ;
	xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	MRIsetVoxVal(mri_ctrl, xvi, yvi, zvi, 0, CONTROL_MARKED) ;
	MRIsetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0, bin) ;
      }
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_ctrl, "ctrl.mgz") ;
      MRIwrite(mri_curv_bins, "cbins.mgz") ;
    }
    MRIbuildVoronoiDiagram(mri_curv_bins, mri_ctrl, mri_curv_bins);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_curv_bins, "cbins.mgz") ;
    }
    MRIfree(&mri_ctrl) ;
  }
  else
    bin_size = 0 ;

  for (vno = 0 ; vno < mris[0]->nvertices ; vno++)
  {
    MRIclear(mri_visited) ;
    if ((vno % 500) == 0)
      printf("processing vno %d of %d (%2.1f%%)\n", vno, mris[0]->nvertices, 100.0f*vno/mris[0]->nvertices) ;
    v = &mris[0]->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (FEQUAL(mris[0]->vertices[vno].x, mris[nlayers]->vertices[vno].x) &&
	FEQUAL(mris[0]->vertices[vno].y, mris[nlayers]->vertices[vno].y) &&
	FEQUAL(mris[0]->vertices[vno].z, mris[nlayers]->vertices[vno].z))
    {
      DiagBreak() ;
      continue ;   // pial and white in same place - vertex is not cortical and can't be estimated
    }
    if (bin_size == 0)
      bin0 = 0 ;
    else
    {
      bin0 = (v->H - Hmin) / (bin_size-1) ;
      bin0 = v->H < -CURV_THRESH ? 1 : (v->H > CURV_THRESH ? 3 : 2) ;
    }

    whalf = whalf0/step ;
    out_of_fov = 0 ;

    // look for supra-threshold volume fraction voxels
    memset(num_found, 0, sizeof(num_found)) ;
    for ( n = 0 ; n <= nlayers ; n++)  // in normal direction
    {
      v = &mris[n]->vertices[vno] ;
      MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, v->x, v->y, v->z, &xv,&yv,&zv) ;
      xv0 = nint(xv) ; yv0 = nint(yv) ; zv0 = nint(zv) ;
      if ((MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0) ||
	  MRIindexNotInVolume(mri_intensities, xv0, yv0, zv0) != 0)
      {
	out_of_fov = 1 ;
	break ;
      }
      
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
	  if (xvi == Gx && yvi == Gy && zvi == Gz)
	    DiagBreak() ;
	  bin = MRIgetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0) ;
	  if (bin != bin0)
	    continue ;
	  if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) > 0)
	    continue ;
	  MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;

	  for (n1 = 0 ; n1 <= nlayers ; n1++) // see if any are above threshold
	  {
	    vfrac = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, n1) ;
	    if (vfrac > vfrac_thresh)
	    {
	      float vnew ;

	      num_found[n1]++ ;
	      val = MRIgetVoxVal(mri_layer_intensities, vno, 0, 0, n1) ;
	      vnew = MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ;
	      if (vno == Gdiag_no)
		printf("layer %d, v %d: val %d with fraction %f found at (%d, %d, %d)\n",
		       n1, vno, (int)vnew, vfrac, xvi, yvi, zvi) ;
	      
	      MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n1, val+vnew) ;
	    }
	  }
	}
    }
/*
  now go whalf into the wm and whalf outside the last surface to get good estimates of WM and CSF
*/
    if (num_found[0] == 0)
      DiagBreak() ;
    for (n = 0 ; n <= nlayers ; n += nlayers)
    {
      v = &mris[n]->vertices[vno] ;
      MRISsurfaceRASToVoxelCached(mris[n], mri_intensities, v->x, v->y, v->z, &xv,&yv,&zv) ;
      xv0 = nint(xv) ; yv0 = nint(yv) ; zv0 = nint(zv) ;
      if (n == 0)  // look inwards from inside surface
      {
	nx = -v->nx ; ny = -v->ny ; nz = -v->nz ;
      }
      else   // look outwards from outside surface
      {
	nx = v->nx ; ny = v->ny ; nz = v->nz ;
      }
      
      for (nd = -1 ; nd <= whalf+1 ; nd++)
      {
	for (t1 = -whalf ; t1 <= whalf ; t1++)    // in tangent plane
	{
	  for (t2 = -whalf ; t2 <= whalf ; t2++)
	  {
	    xs = v->x + step*t1*v->e1x + step*t2*v->e2x + nd*step*nx;
	    ys = v->y + step*t1*v->e1y + step*t2*v->e2y + nd*step*ny;
	    zs = v->z + step*t1*v->e1z + step*t2*v->e2z + nd*step*nz;
	    MRISsurfaceRASToVoxelCached(mris[nlayers], mri_intensities, xs, ys, zs, &xv,&yv,&zv) ;
	    if (MRIindexNotInVolume(mri_intensities, xv, yv, zv) != 0)
	    {
	      out_of_fov = 1 ;
	      break ;
	    }
	    xvi = nint(xv) ; yvi = nint(yv) ; zvi = nint(zv) ;
	    bin = MRIgetVoxVal(mri_curv_bins, xvi, yvi, zvi, 0) ;
	    if (bin != bin0)  
	      continue ;  // not the same curvature
	    if (MRIgetVoxVal(mri_visited, xvi, yvi, zvi, 0) > 0)
	      continue ;
	    MRIsetVoxVal(mri_visited, xvi, yvi, zvi, 0, 1) ;
	    val = MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ;
	    if (val < 0)
	      DiagBreak() ;
	    for (n1 = 0 ; n1 <= nlayers+1 ; n1++)
	    {
	      vfrac = MRIgetVoxVal(mri_volume_fractions, xvi, yvi, zvi, n1) ;
	      if (vfrac >= vfrac_thresh)
	      {
		num_found[n1]++ ;
		val = MRIgetVoxVal(mri_layer_intensities, vno, 0, 0, n1) ;
		val += MRIgetVoxVal(mri_intensities, xvi, yvi, zvi, 0) ;
		MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n1, val) ;
	      }
	    }
	  }
	}
      }
    }

    for (n = 0 ; n <= nlayers+1 ; n++) // fill in outputs
    {
      if (num_found[n] == 0)
	MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n, -1) ;
      else
      {
	val = MRIgetVoxVal(mri_layer_intensities, vno, 0, 0, n) ;
	MRIsetVoxVal(mri_layer_intensities, vno, 0, 0, n, val/num_found[n]) ;
	if (vno == Gdiag_no)
	  printf("v %d, layer %d: %2.1f (%d found)\n", vno, n, val/num_found[n], 
		 num_found[n]) ;
      }
    }
  }

  MRIfree(&mri_visited) ;
  return(mri_layer_intensities) ;
}

