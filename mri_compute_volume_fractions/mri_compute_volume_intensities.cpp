/**
 * @brief compute the unpartial-volumed intensities given an input volume and volume fracion maps
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


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int whalf = 4 ;
static double sigma = 1 ;
static int separate_frames = 0;

static MRI *compute_unpartial_volumed_intensities(MRI *mri_src, MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort, 
						 MRI *mri_vfrac_csf, 
						  int whalf0,  double sigma, MRI *mri_dst, int separate_frames) ;

static void
patch_csf_vol(MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort,  MRI *mri_vfrac_csf)
{
  int x, y, z ;
  double v ;

  for (x = 0 ; x < mri_vfrac_wm->width ; x++)
    for (y = 0 ; y < mri_vfrac_wm->height ; y++)
      for (z = 0 ; z < mri_vfrac_wm->depth ; z++)
      {
	v = MRIgetVoxVal(mri_vfrac_wm, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_cortex, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_subcort, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_csf, x, y, z, 0) ;
	if (FZERO(v))
	  MRIsetVoxVal(mri_vfrac_csf, x, y, z, 0, 1.0) ;
      }
}
int
main(int argc, char *argv[]) 
{
  char   **av ;
  int    ac, nargs ;
  int    msec, minutes, seconds ;
  Timer start ;
  char        fname[STRLEN], *stem ;
  MRI         *mri_src, *mri_vfrac_wm, *mri_vfrac_cortex, *mri_vfrac_subcort, *mri_vfrac_csf, *mri_unpv_intensities ;

  nargs = handleVersionOption(argc, argv, "mri_compute_volume_intensities");
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

  if (argc < 4)
    usage_exit(1) ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;
  
  mri_src = MRIread(argv[1]) ;
  if (mri_src == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load intensity volume from %s", Progname, argv[1]) ;

  stem = argv[2] ;
  sprintf(fname, "%s.cortex.mgz", stem) ;
  mri_vfrac_cortex = MRIread(fname) ;
  if (mri_vfrac_cortex == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read vfrac volume from %s", Progname, fname) ;
  sprintf(fname, "%s.subcort_gm.mgz", stem) ;
  mri_vfrac_subcort = MRIread(fname) ;
  if (mri_vfrac_subcort == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read vfrac volume from %s", Progname, fname) ;
  sprintf(fname, "%s.csf.mgz", stem) ;
  mri_vfrac_csf = MRIread(fname) ;
  if (mri_vfrac_csf == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read vfrac volume from %s", Progname, fname) ;
  sprintf(fname, "%s.wm.mgz", stem) ;
  mri_vfrac_wm = MRIread(fname) ;
  if (mri_vfrac_wm == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read vfrac volume from %s", Progname, fname) ;
  patch_csf_vol(mri_vfrac_wm, mri_vfrac_cortex, mri_vfrac_subcort,  mri_vfrac_csf) ;

  mri_unpv_intensities =   
    compute_unpartial_volumed_intensities(mri_src, mri_vfrac_wm,  mri_vfrac_cortex, mri_vfrac_subcort, 
					  mri_vfrac_csf, whalf,  sigma, NULL, separate_frames) ;
  printf("writing unpartial-volumed intensities to %s\n", argv[3]) ;
  MRIwrite(mri_unpv_intensities, argv[3]) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("unpartial-volume intensities calculation took %d minutes"
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
  } else if (!stricmp(option, "SEPARATE_FRAMES")) {
    separate_frames = 1 ;
    printf("outputting multi-frame volume instead of merging\n") ;
  } else switch (toupper(*option)) {
  case 'W':
    whalf = atoi(argv[2]) ;
    printf("using half window size = %d\n", whalf) ;
    if (whalf < 1)
      ErrorExit(ERROR_BADPARM, "half window size must be >= 1 to estimate unpartial-volumed intensities") ;
    nargs = 1 ;
    break ;
  case 'S':
    sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting Gaussian smoothing sigma to %2.3fmm\n", sigma) ;
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
  printf("usage: %s [options] <input intensity volume> <volume fraction stem> <output volume>\n",
         Progname) ;
  printf(
         "\t\n");
  exit(code) ;
}


static MRI *
compute_unpartial_volumed_intensities(MRI *mri_src, MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort, MRI *mri_vfrac_csf, 
				      int whalf0,  double sigma, MRI *mri_dst, int separate_frames) 
{
  int     xi, yi, zi, xk, yk, zk, x, y, z, nvox, whalfx, whalfy, whalfz, num, row ;
  double  w, distsq, dx, dy, dz, norm, total_gm, total_csf, total_wm, vwm, vgm, vcsf ;
  MATRIX   *m_A_pinv, *m_A3, *m_A2, *m_A1, *m_A ;
  VECTOR  *v_I, *v_s3, *v_s2, *v_s1, *v_s ;
  float   wm, gm, csf ;

  whalfx = (int)ceil(whalf0 / mri_src->xsize) ;
  whalfy = (int)ceil(whalf0 / mri_src->ysize) ;
  whalfz = (int)ceil(whalf0 / mri_src->zsize) ;
  nvox = (2*whalfx + 1) * (2*whalfy + 1) * (2*whalfz + 1) ;
  m_A3 = MatrixAlloc(nvox, 3, MATRIX_REAL) ;
  m_A2 = MatrixAlloc(nvox, 2, MATRIX_REAL) ;
  m_A1 = MatrixAlloc(nvox, 1, MATRIX_REAL) ;
  v_s3 = VectorAlloc(3, MATRIX_REAL) ;
  v_s2 = VectorAlloc(2, MATRIX_REAL) ;
  v_s1 = VectorAlloc(1, MATRIX_REAL) ;
  v_I = VectorAlloc(nvox, MATRIX_REAL) ;


  if (mri_dst == NULL)
  {
    if (separate_frames)
    {
      mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT, 3) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else
      mri_dst = MRIcloneDifferentType(mri_src, MRI_FLOAT) ;
  }

  for (x = 0 ; x < mri_src->width ; x++)
  {
    if (!(x % 10))
      printf("%d of %d\n", x, mri_src->width) ;
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;

	total_gm = total_csf = total_wm = 0.0 ;
	for (norm = 0.0, num = 0, xk = -whalfx ; xk <= whalfx ; xk++)
	{
	  dx = xk*mri_src->xsize ;
	  xi = mri_src->xi[x+xk] ;
	  for (yk = -whalfy ; yk <= whalfy ; yk++)
	  {
	    dy = yk*mri_src->ysize ;
	    yi = mri_src->yi[y+yk] ;
	    for (zk = -whalfz ; zk <= whalfz ; zk++)
	    {
	      dz = zk*mri_src->zsize ;
	      zi = mri_src->zi[z+zk] ;
	      distsq = dx*dx + dy*dy + dz*dz ;
	      w = exp(-.5*distsq/(sigma*sigma)) ;
	      norm += w ;
	      VECTOR_ELT(v_I, num+1) = w*MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
	      vwm = MRIgetVoxVal(mri_vfrac_wm, xi, yi, zi, 0) ;
	      vgm = MRIgetVoxVal(mri_vfrac_cortex, xi, yi, zi, 0) + MRIgetVoxVal(mri_vfrac_subcort, xi, yi, zi, 0) ;
	      vcsf = MRIgetVoxVal(mri_vfrac_csf, xi, yi, zi, 0)  ; ;
	      *MATRIX_RELT(m_A3, num+1, 1) = w*vwm ; *MATRIX_RELT(m_A3, num+1, 2) = w*vgm ; *MATRIX_RELT(m_A3, num+1, 3) = w*vcsf ;
	      total_gm += w*vgm ; total_csf += w*vcsf ; total_wm += w*vwm ;
	      num++ ;
	      if (num  == Gdiag_no)
		DiagBreak() ;
	    }
	  }
	}

	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	if (total_csf < total_gm/100 || total_csf < total_wm/100)
	  total_csf = 0 ;
	if (total_gm < total_csf/100 || total_gm < total_wm/100)
	  total_gm = 0 ;
	if (total_wm < total_gm/100 || total_wm < total_csf/100)
	  total_wm = 0 ;

	if (!FZERO(total_gm))   // some gm in this voxel
	{
	  if (!FZERO(total_wm) && !FZERO(total_csf))
	  {
	    v_s = v_s3 ;
	    m_A = m_A3 ;
	  }

	  else if (!FZERO(total_wm))   // estimate wm and gm, but not csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 2) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else if (!FZERO(total_csf))  // estimate gm and csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 2) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 3) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else   // only gm in this voxel
	  {
	    for (row = 1 ; row < m_A3->rows ; row++)
	      *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 2) ;

	    v_s = v_s1 ;
	    m_A = m_A1 ;
	  }
	}
	else if (!FZERO(total_wm))  // some wm in voxel, but no gm
	{
	  if (!FZERO(total_csf))  // estimate wm and csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 3) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else   // only wm in this voxel
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++) {
	      *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;
	    }

	    v_s = v_s1 ;
	    m_A = m_A1 ;
	  }
	}
	else  // only csf in this region
	{
	  for (row = 1 ; row <= m_A3->rows ; row++) {
	    *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 3) ;
	  }

	  v_s = v_s1 ;
	  m_A = m_A1 ;
	}

	m_A_pinv = MatrixPseudoInverse(m_A, NULL) ;
	if (m_A_pinv == NULL)
	  continue ;

	MatrixMultiply(m_A_pinv, v_I, v_s) ;
	vwm = MRIgetVoxVal(mri_vfrac_wm, x, y, z, 0) ;
	vcsf = MRIgetVoxVal(mri_vfrac_csf, x, y, z, 0) ;
	vgm = MRIgetVoxVal(mri_vfrac_cortex, x, y, z, 0) + MRIgetVoxVal(mri_vfrac_subcort, x, y, z, 0) ;
	wm = gm = csf = 0.0 ;
	if (!FZERO(total_wm))  // wm in 1st col
	{
	  wm = (float)VECTOR_ELT(v_s, 1) ; 
	  if (!FZERO(total_gm))
	  {
	    gm = (float)VECTOR_ELT(v_s, 2) ; 
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 3) ; 
	  }
	  else   // wm but no gm, csf in col 2
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 2) ; 
	}
	else  // no wm in this region
	{
	  if (!FZERO(total_gm))  // some gm in this voxel
	  {
	    gm = (float)VECTOR_ELT(v_s, 1) ; 
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 2) ; 
	  }
	  else  // only csf in this voxel
	    csf = (float)VECTOR_ELT(v_s, 1) ; 
	}

	if (separate_frames)
	{
	  MRIsetVoxVal(mri_dst, x, y, z, 0, wm) ;
	  MRIsetVoxVal(mri_dst, x, y, z, 1, gm) ;
	  MRIsetVoxVal(mri_dst, x, y, z, 2, csf) ;
	}
	else  // store the value for the class with the largest volume fraction
	{
	  if (vwm > vgm && vwm > vcsf) // wm always in first col
	    MRIsetVoxVal(mri_dst, x, y, z, 0, wm) ;  // mostly wm
	  else if (vgm > vcsf)   // mostly gm
	    MRIsetVoxVal(mri_dst, x, y, z, 0, gm) ;  
	  else   // wm not estimated, in col 1
	    MRIsetVoxVal(mri_dst, x, y, z, 0, csf) ;  

	  if (!devFinite(MRIgetVoxVal(mri_dst, x, y, z, 0)))
	  {
	    DiagBreak() ;
	    MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
	  }
	}
	MatrixFree(&m_A_pinv) ;
      }
  }


  MatrixFree(&m_A3) ; VectorFree(&v_s3) ; VectorFree(&v_I) ;
  MatrixFree(&m_A2) ; VectorFree(&v_s2) ; 
  MatrixFree(&m_A1) ; VectorFree(&v_s1) ; 
  return(mri_dst) ;
}

