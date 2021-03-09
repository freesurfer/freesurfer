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
#include "label.h"
#include "transform.h"
#include "mri_conform.h"
#include "ctrpoints.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int update_bias(MRI *mri_orig, MRI *mri_T1, MRI *mri_brain, MRI *mri_bias, MRI *mri_counts, int normalize) ;
static int normalize_bias(MRI *mri_bias, MRI *mri_counts, int normalize);

const char *Progname ;

static void usage_exit(int code) ;

static LABEL *label = NULL ;
static char *control_fname = NULL ;
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
  Timer start ;
	double       c_r, c_a, c_s ;
  MATRIX       *m_vox2vox ;

  nargs = handleVersionOption(argc, argv, "mri_compute_bias");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  out_fname = argv[argc-1] ;
  
  if (label)
  {
    MRI *mri, *mri_control, *mri_bias, *mri_kernel, *mri_smooth ;
    double mean, val ;
    int   n, xv, yv, zv ;
    
    printf("reading label %s and volume %s to compute bias field\n", label->name, argv[1]) ;
    mri = MRIread(argv[1]) ;
    if (mri == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", Progname, argv[1]) ;
    mean = LabelMeanIntensity(label, mri) ;
    printf("mean in WM label is %2.1f\n", mean) ;

    mri_control = MRIcloneDifferentType(mri, MRI_UCHAR) ;
    mri_bias = MRIcloneDifferentType(mri, MRI_FLOAT) ;

    if (control_fname)
    {
      int   n, total ;
      MRI *mri_conf ;
      MATRIX *m_vox2vox ;
      LABEL *ltmp, *lconf ;
      MPoint *pointArray ;

      mri_conf = MRIconform(mri) ;
      m_vox2vox = MRIgetVoxelToVoxelXform(mri, mri_conf) ;
      printf("writing output control points to %s\n", control_fname) ;
      MatrixPrint(stdout, m_vox2vox) ;

      ltmp = LabelToVoxel(label, mri, NULL) ;

      lconf = LabelApplyMatrix(ltmp, m_vox2vox, NULL) ;
      LabelVoxelToSurfaceRAS(lconf, mri_conf, lconf) ;

      pointArray=(MPoint*) malloc(lconf->n_points * sizeof(MPoint));
      for (n = 0 ; n < lconf->n_points ; n++)
      {
	pointArray[n].x = (double)lconf->lv[n].x ;
	pointArray[n].y = (double)lconf->lv[n].y ;
	pointArray[n].z = (double)lconf->lv[n].z ;
      }

      if (FileExists(control_fname))
      {
	MPoint *existingPointArray, *newArray ;
	int    count, useRealRAS, i1, i2, duplicate ;

	existingPointArray = MRIreadControlPoints(control_fname, &count, &useRealRAS);
	printf("examining %d existing control points to remove duplicated\n",count) ;
	if (useRealRAS)
	  ErrorExit(ERROR_UNSUPPORTED, "%s: cannot combine existing RAS and current tkregRAS control points", Progname) ;
	duplicate = 0 ;
	for (i1 = 0 ; i1 < lconf->n_points ; i1++)
	  for (i2 = 0 ; i2 < count ; i2++)
	  {
	    if ((nint(pointArray[i1].x) == nint(existingPointArray[i2].x)) &&
		(nint(pointArray[i1].y) == nint(existingPointArray[i2].y)) &&
		(nint(pointArray[i1].z) == nint(existingPointArray[i2].z)))
	    {
	      duplicate++ ;
	      break ;
	    }
	  }
	total = count + lconf->n_points - duplicate ;
	printf("%d duplicate points detected,  saving %d\n", duplicate, total) ;
	if (total > lconf->n_points)  // if they are are not duplicate
	{
	  newArray = (MPoint *)calloc(total, sizeof(MPoint)) ;
	  if (newArray == NULL)
	    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d-len Point array", Progname, total) ;
	  memmove(newArray, pointArray, lconf->n_points*sizeof(MPoint)) ;
	  //  add only those existing points that are unique
	  for (total = lconf->n_points, i2 = 0 ; i2 < count ; i2++)
	  {
	    duplicate = 0 ;
	    for (i1 = 0 ; i1 < lconf->n_points ; i1++)
	    {
	      if ((nint(pointArray[i1].x) == nint(existingPointArray[i2].x)) &&
		  (nint(pointArray[i1].y) == nint(existingPointArray[i2].y)) &&
		  (nint(pointArray[i1].z) == nint(existingPointArray[i2].z)))
	      {
		duplicate = 1 ;
		break ;
	      }
	    }
	    if  (duplicate)
	      continue ;
	    newArray[total].x = existingPointArray[i2].x ;
	    newArray[total].y = existingPointArray[i2].y ;
	    newArray[total].z = existingPointArray[i2].z ;
	    total++ ;
	  }
	  free(pointArray) ;
	  pointArray = newArray ;
	}
	free(existingPointArray) ;
      }
      else
	total = lconf->n_points ;
      MRIwriteControlPoints(pointArray, total, 0, control_fname) ;
      MatrixFree(&m_vox2vox) ;
    }
    LabelToVoxel(label, mri, label) ;
    for (n = 0 ; n < label->n_points ; n++)
    {
      xv = nint(label->lv[n].x) ; yv = nint(label->lv[n].y) ; zv = nint(label->lv[n].z) ;
      MRIsetVoxVal(mri_control, xv, yv, zv, 0, CONTROL_MARKED) ;
      val = MRIgetVoxVal(mri, xv, yv, zv, 0) ;
      MRIsetVoxVal(mri_bias, xv, yv, zv, 0, mean/val) ;
    }
    MRIbuildVoronoiDiagram(mri_bias, mri_control, mri_bias) ;
    mri_kernel = MRIgaussian1d(sigma, 0) ;
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
    MRIfree(&mri_kernel) ; 

    printf("writing bias field to %s\n", out_fname) ;
    MRIwrite(mri_smooth, out_fname) ;
    exit(0) ;
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
  msec = start.milliseconds() ;
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
  else if (!stricmp(option, "label")) 
  {
    label = LabelRead(NULL, argv[2]) ;
    nargs = 1 ;
    printf("reading label %s and using it to compute bias field\n", argv[2]) ;
  }
  else 
    switch (toupper(*option))
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
    case 'C':
      control_fname = argv[2] ;
      nargs = 1;  
      printf("writing label to control point file %s\n", control_fname) ;
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
