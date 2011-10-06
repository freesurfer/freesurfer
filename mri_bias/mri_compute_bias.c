/**
 * @file  mri_compute_bias.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/10/06 01:19:35 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "mrinorm.h"
#include "version.h"
#include "transform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int update_bias(MRI *mri_orig, MRI *mri_T1, MRI *mri_brain, MRI *mri_bias, MRI *mri_counts, int normalize) ;
static int normalize_bias(MRI *mri_bias, MRI *mri_counts, int normalize);

char *Progname ;

static void usage_exit(int code) ;

static char sdir[STRLEN] = "" ;
static float pad=20;
static double sigma = 4.0 ;
static int normalize = 0 ;
static char *xform_name ;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN], *subject, *cp ;
  int          ac, nargs, i ;
  MRI          *mri_orig, *mri_T1, *mri_bias, *mri_counts, *mri_kernel, *mri_smooth, *mri_brain ;
  char         *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;
	double       c_r, c_a, c_s ;
  MATRIX       *m_vox2vox ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mri_compute_bias.c,v 1.7 2011/10/06 01:19:35 fischl Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
	{
		nargs = get_option(argc, argv) ;
		argc -= nargs ;
		argv += nargs ;
	}

  if (argc < 3)
    usage_exit(1) ;

	if (strlen(sdir) == 0)
	{
		cp = getenv("SUBJECTS_DIR");
		if (cp)
			strcpy(sdir, cp) ;
		else
			ErrorExit(ERROR_UNSUPPORTED, "%s: must specify SUBJECTS_DIR in env or with -sdir on cmd line", Progname);
	}
	
  out_fname = argv[argc-1] ;

	mri_bias = MRIalloc(256+2*pad, 256+2*pad, 256+2*pad, MRI_FLOAT) ;
	mri_counts = MRIalloc(256+2*pad, 256+2*pad, 256+2*pad, MRI_FLOAT) ;
#if 0
	mri_bias->c_r = (double)mri_bias->width/2.0 ;
	mri_bias->c_a = (double)mri_bias->height/2.0 ;
	mri_bias->c_s = (double)mri_bias->depth/2.0 ;
	mri_counts->c_r = (double)mri_counts->width/2.0 ;
	mri_counts->c_a = (double)mri_counts->height/2.0 ;
	mri_counts->c_s = (double)mri_counts->depth/2.0 ;
#endif

  for (i = 1 ; i < argc-1 ; i++)
	{
		subject = argv[i] ;
		fprintf(stderr, "subject %s, %d of %d...\n", subject, i, argc-2) ;

		sprintf(fname, "%s/%s/mri/orig.mgz", sdir, subject) ;
		mri_orig = MRIread(fname) ;
		if (!mri_orig)
			ErrorExit(Gerror, "%s: could not read orig volume %s", Progname, fname) ;
		sprintf(fname, "%s/%s/mri/T1.mgz", sdir, subject) ;
		mri_T1 = MRIread(fname) ;
		if (!mri_T1)
			ErrorExit(Gerror, "%s: could not read T1 volume %s", Progname, fname) ;
		sprintf(fname, "%s/%s/mri/brainmask.mgz", sdir, subject) ;
		mri_brain = MRIread(fname) ;
		if (!mri_T1)
			ErrorExit(Gerror, "%s: could not read T1 volume %s", Progname, fname) ;
		if (i == 1)
		{
			mri_bias->xstart = mri_orig->xstart - pad*mri_orig->xsize ;
			mri_bias->ystart = mri_orig->ystart - pad*mri_orig->ysize ;
			mri_bias->zstart = mri_orig->zstart - pad*mri_orig->zsize ;
			mri_bias->xend   = mri_orig->xend   + pad*mri_orig->xsize ;
			mri_bias->yend   = mri_orig->yend   + pad*mri_orig->ysize ;
			mri_bias->zend   = mri_orig->zend   + pad*mri_orig->zsize ;

			mri_counts->xstart = mri_orig->xstart - pad*mri_orig->xsize ;
			mri_counts->ystart = mri_orig->ystart - pad*mri_orig->ysize ;
			mri_counts->zstart = mri_orig->zstart - pad*mri_orig->zsize ;
			mri_counts->xend   = mri_orig->xend   + pad*mri_orig->xsize ;
			mri_counts->yend   = mri_orig->yend   + pad*mri_orig->ysize ;
			mri_counts->zend   = mri_orig->zend   + pad*mri_orig->zsize ;
			MRIreInitCache(mri_bias) ; MRIreInitCache(mri_counts) ;
			MRIcalcCRASforExtractedVolume
				(mri_orig, mri_bias, pad, pad, pad, 
				 pad+mri_orig->width, pad+mri_orig->height, pad+mri_orig->depth, 
				 &c_r, &c_a, &c_s);
			mri_bias->c_r = mri_bias->c_a = mri_bias->c_s = 0 ;
			mri_counts->c_r = mri_counts->c_a = mri_counts->c_s = 0;
			MRIreInitCache(mri_bias) ; MRIreInitCache(mri_counts) ;
		}
    if (xform_name)
    {
      TRANSFORM    *transform ;
      MRI          *mri ;
      sprintf(fname, "%s/%s/mri/transforms/%s", sdir, subject, xform_name) ;
      transform = TransformRead(fname) ;
      if (transform == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not load transform from %s", Progname, fname) ;
      m_vox2vox = MRIgetVoxelToVoxelXform(mri_bias, mri_orig) ;
      //      TransformCompose(transform, m_vox2vox, NULL, transform) ;
      mri = MRIclone(mri_bias, NULL) ;
      //      mri->c_r = mri->c_a = mri->c_s = 0 ;
      TransformApplyType(transform, mri_orig, mri, SAMPLE_TRILINEAR) ;
      MRIfree(&mri_orig) ; mri_orig = mri ;
      mri = MRIclone(mri_bias, NULL) ;
      //      mri->c_r = mri->c_a = mri->c_s = 0 ;
      TransformApplyType(transform, mri_T1, mri, SAMPLE_TRILINEAR) ;
      MRIfree(&mri_T1) ; mri_T1 = mri ;

      mri = MRIclone(mri_bias, NULL) ;
      //      mri->c_r = mri->c_a = mri->c_s = 0 ;
      TransformApplyType(transform, mri_brain, mri, SAMPLE_TRILINEAR) ;
      MRIfree(&mri_brain) ; mri_brain = mri ;
      TransformFree(&transform) ;
    }
		update_bias(mri_orig, mri_T1, mri_brain, mri_bias, mri_counts, 0) ;
	}

	normalize_bias(mri_bias, mri_counts, normalize);
  mri_kernel = MRIgaussian1d(sigma, 0) ;
	mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
	MRIfree(&mri_kernel) ; 

  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_smooth, out_fname) ;
  MRIfree(&mri_smooth) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "bias calculation took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "sdir"))
	{
		strcpy(sdir, argv[2]) ;
		printf("using SUBJECTS_DIR %s\n", sdir);
		nargs = 1 ;
	}
  else if (!stricmp(option, "debug_voxel")) 
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else switch (toupper(*option))
	{
  case 'T':
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("applying xform %s to input datasets\n", xform_name) ;
    break ;
	case 'S':
		sigma = atof(argv[2]) ;
		nargs = 1 ;
		printf("using sigma = %2.2f\n", sigma) ;
		break ;
	case 'N':
		normalize = 1 ;
		printf("normalizing bias maps before combining.\n") ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <subject 1> <subject 2> ... <output volume>\n", Progname) ;
  exit(code) ;
}
static int
update_bias(MRI *mri_orig, MRI *mri_T1, MRI *mri_brain, MRI *mri_bias, MRI *mri_counts, int normalize)
{
	MATRIX   *m_vox2vox;
	VECTOR   *v1, *v2;
	int      x, y, z, num ;
	double   xd, yd, zd, bias, val_orig, val_T1, mn, val_brain ;
	MRI      *mri_tmp ;
	
	m_vox2vox = MRIgetVoxelToVoxelXform(mri_orig, mri_bias) ;
	v1 = VectorAlloc(4, MATRIX_REAL);
	v2 = VectorAlloc(4, MATRIX_REAL);
	VECTOR_ELT(v1, 4) = 1.0 ;
	VECTOR_ELT(v2, 4) = 1.0 ;

	mri_tmp = MRIalloc(mri_orig->width, mri_orig->height, mri_orig->depth, MRI_FLOAT) ;
	MRIcopyHeader(mri_orig, mri_tmp) ;

	for (mn = 0.0, num = x = 0 ; x < mri_orig->width ; x++)
	{
		for (y = 0 ; y < mri_orig->height ; y++)
		{
			for (z = 0 ; z < mri_orig->depth ; z++)
			{
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
				val_orig = MRIgetVoxVal(mri_orig, x, y, z, 0) ;
				val_T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
				if (FZERO(val_orig))
					continue ;
				val_brain = MRIgetVoxVal(mri_brain, x, y, z, 0) ;
				if (FZERO(val_brain))  // don't use it for mean calculation
					continue ;
				bias = val_T1/val_orig ;
				MRIsetVoxVal(mri_tmp, x, y, z, 0, bias) ;
				mn += bias ;
				num++ ;
			}
		}
	}
	if (num > 0)
		mn /= num ;
	if (normalize)
	{
		printf("normalizing bias by %2.3f\n", mn) ;
		for (x = 0 ; x < mri_orig->width ; x++)
		{
			for (y = 0 ; y < mri_orig->height ; y++)
			{
				for (z = 0 ; z < mri_orig->depth ; z++)
				{
          val_brain = MRIgetVoxVal(mri_brain, x, y, z, 0) ;
          if (FZERO(val_brain))  // don't use it for mean calculation
            continue ;
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
					bias = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;
					bias /= mn ;
					MRIsetVoxVal(mri_tmp, x, y, z, 0, bias) ;
				}
			}
		}
	}


	for (x = 0 ; x < mri_orig->width ; x++)
	{
		V3_X(v1) = x ;
		for (y = 0 ; y < mri_orig->height ; y++)
		{
			V3_Y(v1) = y ;
			for (z = 0 ; z < mri_orig->depth ; z++)
			{
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
				val_orig = MRIgetVoxVal(mri_orig, x, y, z, 0) ;
				if (FZERO(val_orig))
					continue ;  // not set
				val_brain = MRIgetVoxVal(mri_brain, x, y, z, 0) ;
				if (FZERO(val_brain))  // don't use it for mean calculation
					continue ;
				V3_Z(v1) = z ;
				bias = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;
				MatrixMultiply(m_vox2vox, v1, v2) ;
				xd = V3_X(v2) ; yd = V3_Y(v2) ;zd = V3_Z(v2);
        if (nint(xd) == Gx && nint(yd) == Gy && nint(zd) == Gz)
          DiagBreak() ;
				MRIinterpolateIntoVolume(mri_bias, xd, yd, zd, bias);
				MRIinterpolateIntoVolume(mri_counts, xd, yd, zd, 1);
			}
		}
	}

	MRIfree(&mri_tmp) ;
	MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
	return(NO_ERROR) ;
}

static int
normalize_bias(MRI *mri_bias, MRI *mri_counts, int normalize)
{
	int    x, y, z, num ;
	double count, bias, mn ;
	MRI    *mri_ctrl ;

  mri_ctrl = MRIalloc(mri_bias->width, mri_bias->height, mri_bias->depth, MRI_UCHAR) ;
	for (mn = 0.0, num = x = 0 ; x < mri_bias->width ; x++)
	{
		for (y = 0 ; y < mri_bias->height ; y++)
		{
			for (z = 0 ; z < mri_bias->depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				count = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
				if (FZERO(count))
					continue ;
				bias = MRIgetVoxVal(mri_bias, x, y, z, 0) ;
				MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_MARKED) ;
				bias /= count ;
        mn += bias ;
        num++ ;
				MRIsetVoxVal(mri_bias, x, y, z, 0, bias) ;
			}
		}
	}

  if (normalize && num > 0)
  {
    mn /= num ;
    MRIscalarMul(mri_bias, mri_bias, 1.0/mn) ;
  }
  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
  MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, 20, 1) ;

  MRIfree(&mri_ctrl) ;
  return(NO_ERROR) ;
}
