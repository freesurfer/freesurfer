/**
 * @brief synthesize a histological volume from an MRI
 *
 * synthesize a histological volume from an MRI
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


#include "romp_support.h"

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

typedef struct
{
  int    len ;
  int    whalf ;
  double *vals ;
  int    flags ;
} FEATURE ;


static int dump_window(MRI *mri, const char *fname, int x0, int y0, int z0, int wsize)  ;
static int prune_indices(MRI *mri_mask, MRI *histo, short *xind, short *yind, short *zind, int nind, int use_val, float hthresh) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags, char *fname, MRI *mri_train_src, MRI *mri_train_dst) ;
static int extract_feature_vector(MRI *mri, int x, int y, int z, FEATURE *feature) ;
static int find_most_similar_location(MRI *mri, FEATURE *fsrc, int z, int *pxd, int *pyd,int *pzd,
                                      MRI *mri_mask, MRI_REGION *box, int flags,
                                      short *xind, short *yind, short *zind, int nind,
				      double tol,
                                      int num_notfound,
				      int x0, int y0, int z0, double min_dist) ;
static double feature_distance(FEATURE *f1, FEATURE *f2, int which) ;

const char *Progname ;
static void usage_exit(int code) ;

static char base_name[STRLEN] = "" ;
static int crop_width = 0 ;
static int test_slice = 20 ;  // in histo coords
static int train_slice = 30 ; // in MRI coords
static int wsize = 3 ;
static int downsample = 0 ;
static double tol = 0 ;
static int num_notfound = 1000 ;  // # of search voxels to terminate after
static double min_training_dist = 100 ;


static char *training_src_fname  = NULL ;
static char *training_dst_fname = NULL ;

#define SUBTRACT_CENTER  0x00001
#define L1_NORM          0x00002
#define MRI_SPACE        0x00004
#define TRAINING_PAIR    0x00008

#define L2_NORM_DIST     0
#define L1_NORM_DIST     1


static int flags = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI          *mri, *histo, *hsynth, *mri_train_src = NULL, *mri_train_dst = NULL ;

  setRandomSeed(-1L) ;

  nargs = handleVersionOption(argc, argv, "histo_synthesize");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  FileNameRemoveExtension(argv[3], base_name) ;

  if (argc < 3)
    usage_exit(1) ;

  mri = MRIread(argv[1]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, argv[1]) ;
  if (downsample > 0)
  {
    MRI *mri_tmp ;
    int n ;

    for (n = 0 ; n < downsample ; n++)
    {
      mri_tmp = MRIreduce(mri, NULL) ;
      MRIfree(&mri) ; mri = mri_tmp ;
      train_slice /= 2 ;
      if (flags & MRI_SPACE)
        test_slice /= 2 ;

    }
  }
  histo = MRIread(argv[2]) ;
  if (histo == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read histological volume %s", Progname, argv[2]) ;
  if (training_src_fname)
  {
    mri_train_src = MRIread(training_src_fname) ;
    if (!mri_train_src)
      ErrorExit(ERROR_NOFILE, "%s: could not read training src file %s", Progname, training_src_fname) ;
    mri_train_dst = MRIread(training_dst_fname) ;
    if (!mri_train_dst)
      ErrorExit(ERROR_NOFILE, "%s: could not read training dst file %s", Progname, training_dst_fname) ;
    if (train_slice >= mri_train_dst->depth)
      train_slice = MAX(0, (mri_train_dst->depth-1)/2) ;
  }
  else
  {
    if (train_slice >= histo->depth)
      train_slice = MAX(0, (histo->depth-1)/2) ;
  }
  if (test_slice >= mri->depth)
    test_slice = (mri->depth-1)/2 ;
  hsynth = HISTOsynthesize(mri, histo, test_slice, train_slice, NULL, wsize, flags, argv[3], mri_train_src, mri_train_dst) ;
  printf("writing output to %s\n", argv[3]) ;
  MRIwrite(hsynth, argv[3]) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "image synthesis took %d minutes"
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
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atof(argv[2]) ;
    Gy = atof(argv[3]) ;
    Gz = atof(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "test_slice"))
  {
    test_slice = atoi(argv[2]) ;
    nargs = 1 ;
    printf("testing on slice %d\n", test_slice) ;
  }
  else if (!stricmp(option, "train"))
  {
    training_src_fname = argv[2] ;
    training_dst_fname = argv[3] ;
    flags |= TRAINING_PAIR ;
    nargs = 2 ;
    printf("reading training pair %s->%s and appending to images\n", training_src_fname, training_dst_fname) ;
    min_training_dist = 0 ;
  }
  else if (!stricmp(option, "crop_width"))
  {
    crop_width = atoi(argv[2]) ;
    nargs = 1 ;
    printf("cropping width to %d in synthesis\n", crop_width) ;
  }
  else if (!stricmp(option, "tol"))
  {
    tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting search termination tol = %f\n", tol) ;
  }
  else if (!stricmp(option, "num"))
  {
    num_notfound = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting search termination num = %d\n", num_notfound) ;
  }
  else if (!stricmp(option, "train_slice"))
  {
    train_slice = atoi(argv[2]) ;
    nargs = 1 ;
    printf("training on slice %d\n", train_slice) ;
  }
  else if (!stricmp(option, "L1"))
  {
    flags |= L1_NORM ;
    printf("using L1 norm for distance measure\n") ;
  }
  else if (!stricmp(option, "MRI"))
  {
    flags |= MRI_SPACE ;
    printf("synthesizing data in MRI coords\n") ;
  }
  else switch (toupper(*option)) {
  case 'D':
    downsample = atoi(argv[2]) ;
    printf("downsampling %d times\n", downsample) ;
    nargs = 1 ;
    break ;
  case 'W':
    wsize = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using window size = %d\n",wsize) ;
    break ;
  case 'C':
    flags |= SUBTRACT_CENTER ;
    printf("subtracting mean from feature vectors\n") ;
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
  printf("usage: %s [options] <MRI volume> <HISTO volume> <synthetic histo>",
         Progname) ;
  printf(
    "\t\n"
  );
  exit(code) ;
}


MRI *
HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags, char *fname, MRI *mri_train_src, MRI *mri_train_dst)
{
  int     x, y, z, xd=0, yd=0, zd=0, xm, ym, zm, x1, y1, z1, identity, tid, nind, width ;
  short   *xind, *yind, *zind ;
  double  val, xh, yh, zh, sval ;
  FEATURE farray[_MAX_FS_THREADS], *f = &farray[0] ;
  MATRIX  *m_histo2mri, *m_mri2histo ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_mask ;
  MRI_REGION box_train, box_test ;

  xh = yh = zh = y = 0 ;  // compiler warning - not sure why
  for (x = 0 ; x < _MAX_FS_THREADS ; x++)
  {
    farray[x].len = wsize*wsize*wsize ;
    farray[x].whalf = (wsize-1)/2 ;
    farray[x].flags = flags ;
    farray[x].vals = (double *)calloc(farray[x].len, sizeof(double)) ;
    if (farray[x].vals == NULL)
      ErrorExit(ERROR_NOMEMORY, "HISTOsynthesize: could not allocate %d-len feature vector", farray[x].len) ;
  }

  if (hsynth == NULL)
  {
    if (flags & MRI_SPACE)
      hsynth = MRIclone(mri, NULL) ;
    else
      hsynth = MRIclone(histo, NULL) ;
    MRIsetValues(hsynth, 255) ;
  }

  if (mri_train_src == NULL)
    mri_train_src = mri ;
  if (mri_train_dst == NULL)
    mri_train_dst = histo ;

  m_histo2mri = MRIgetVoxelToVoxelXform(mri_train_dst, mri_train_src) ;
  m_mri2histo = MRIgetVoxelToVoxelXform(mri_train_src, mri_train_dst) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;
  identity = MatrixIsIdentity(m_histo2mri) ;
  if (identity)
    printf("image geometries are identical - disabling transforms\n") ;

  // compute region of MRI that maps to valid histo
  box_train.x = mri_train_src->width; box_train.y = mri_train_src->height ; box_train.z = mri_train_src->depth ;
  box_train.dx = 0 ; box_train.dy = 0 ; box_train.dz = 0 ;

  mri_mask = MRIclone(mri_train_src, NULL) ;
  for (x = 0 ; x < mri_train_src->width ; x++)
  {
    for (y = 0 ; y < mri_train_src->height ; y++)
    {
      for (z = 0 ; z < mri_train_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
	if (identity)
	{
	  xh = x ; yh = y ; zh = z ;
	  val = MRIgetVoxVal(mri_train_dst, x, y, z, 0) ;
	  sval = MRIgetVoxVal(mri_train_src, x, y, z, 0) ;
	}
	else
	{
	  V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
	  MatrixMultiply(m_mri2histo, v1, v2) ;
	  xh = nint(V3_X(v2)) ; yh = nint(V3_Y(v2)) ; zh = nint(V3_Z(v2)) ;  // MRI coord
	  MRIsampleVolume(mri_train_dst, xh, yh, zh, &val) ;
	  MRIsampleVolume(mri_train_src, xh, yh, zh, &sval) ;
	}
	if (FZERO(sval))
	  continue ;
        if ((val < 240) && (val != 4080) && MRIindexNotInVolume(mri_train_dst, xh, yh, zh) == 0)
        {
          x1 = box_train.x + box_train.dx-1 ;
          y1 = box_train.y + box_train.dy-1 ;
          z1 = box_train.z + box_train.dz-1 ;
          if (x < box_train.x)
            box_train.x = x ;
          if (y < box_train.y)
            box_train.y = y ;
          if (z < box_train.z)
            box_train.z = z ;
          if (x > x1)
            x1 = x ;
          if (y > y1)
            y1 = y ;
          if (z > z1)
            z1 = z ;
          box_train.dx = x1 - box_train.x + 1 ;
          box_train.dy = y1 - box_train.y + 1 ;
          box_train.dz = z1 - box_train.z + 1 ;
          MRIsetVoxVal(mri_mask, x, y, z, 0, 1) ;
        }
      }
    }
  }

  nind = mri_train_src->width*mri_train_src->height*mri_train_src->depth ;
  xind = (short *)calloc(nind, sizeof(short)) ;
  yind = (short *)calloc(nind, sizeof(short)) ;
  zind = (short *)calloc(nind, sizeof(short)) ;
  MRIcomputeVoxelPermutation(mri_train_src, xind, yind,zind) ;
  nind = prune_indices(mri_mask, mri_train_dst, xind, yind, zind, nind, 1, 240) ;
  if (mri_train_src != mri)
    MRIhistogramNormalize(mri_train_src,  mri, mri_train_src) ;
//  MRImask(mri_train_src, mri_mask, mri_train_src, 0, 255) ;
  if (flags & MRI_SPACE) 
  {
    MRIsetValues(hsynth, 4080) ;
    box_test.x = mri->width; box_test.y = mri->height ; box_test.z = mri->depth ;
    box_test.dx = 0 ; box_test.dy = 0 ; box_test.dz = 0 ;
    z = test_slice ;
    for (x = 0 ; x < mri->width ; x++)
    {
      for (y = 0 ; y < mri->height ; y++)
      {
        if (x == Gx && y == Gy)
          DiagBreak() ;
	if (identity)
	{
	  xh = x ; yh = y ; zh = test_slice ;
	  val = MRIgetVoxVal(histo, xh, yh, zh, 0) ;
	}
	else
	{
	  V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
	  MatrixMultiply(m_mri2histo, v1, v2) ;
	  xh = nint(V3_X(v2)) ; yh = nint(V3_Y(v2)) ; zh = nint(V3_Z(v2)) ;  // MRI coord
	  MRIsampleVolume(histo, xh, yh, zh, &val) ;
	}
        if ((val != 4080) && MRIindexNotInVolume(histo, xh, yh, zh) == 0)
        {
          x1 = box_test.x + box_test.dx-1 ;
          y1 = box_test.y + box_test.dy-1 ;
          z1 = box_test.z + box_test.dz-1 ;
          if (x < box_test.x)
            box_test.x = x ;
          if (y < box_test.y)
            box_test.y = y ;
          if (z < box_test.z)
            box_test.z = z ;
          if (x > x1)
            x1 = x ;
          if (y > y1)
            y1 = y ;
          if (z > z1)
            z1 = z ;
          box_test.dx = x1 - box_test.x + 1 ;
          box_test.dy = y1 - box_test.y + 1 ;
          box_test.dz = z1 - box_test.z + 1 ;
        }
      }
    }
  
    for (x = box_test.x ; x < box_test.x+box_test.dx ; x++)
    {
      if ((x % 10 == 0))
      {
        if (fname && (x-box_test.x) > 0)
        {
          printf("%03d of %03d - writing file %s\n", x-box_test.x, box_test.dx, fname) ;
          MRIwrite(hsynth, fname) ;
        }
        else
          printf("%03d of %03d\n", x-box_test.x, box_test.dx) ;
        fflush(stdout) ;
      }
      for (y = box_test.y ; y < box_test.y+box_test.dy ; y++)
      {
        if ( x == Gx && y == Gy && test_slice == Gz)
          DiagBreak() ;
        if (MRIgetVoxVal(mri_mask, x, y, test_slice,0) == 0)
          continue ;
        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
        extract_feature_vector(mri, x, y, test_slice, f) ;
        find_most_similar_location(mri, f, train_slice, &xd, &yd, &zd, mri_mask, &box_train, flags,
                                   xind, yind, zind, nind, tol, num_notfound,
				   x, y, test_slice, min_training_dist) ;
	if (identity)
	{
	  xh = xd ; yh = yd ; zh = zd ;
	  val =	MRIgetVoxVal(histo, xh, yh, zh, 0) ;
	}
	else
	{
	  V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ;
	  MatrixMultiply(m_mri2histo, v1, v2) ;
	  xh = V3_X(v2) ;  yh = V3_Y(v2) ;  zh = V3_Z(v2) ; 
	  MRIsampleVolume(histo, xh, yh, zh, &val) ;
	}
#if 0
        if (val == 4080)
          val = 0 ;  // make background black
#endif
        MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
      }
    }  
  }
  else    // go through each pixel in the histo space
  {
    width = crop_width > 0 ? crop_width : hsynth->width ;
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    xm = ym = zm = 0 ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(m_histo2mri, m_mri2histo, xind, yind, zind, hsynth, y, val, histo, test_slice, Gx, Gy, xm, ym, zm, identity, mri_mask, mri, train_slice, min_training_dist, tol, num_notfound, xh, yh, zh, xd, yd, zd, f, farray) shared(v1, v2) schedule(static,1)
#endif
  for (x = 0 ; x < width ; x++)   // synth is in histo coords
  {
    ROMP_PFLB_begin
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
    f = &farray[tid] ;
#else
    tid = 0;
#endif
    if ((x % 10 == 0))
    {
      printf("%03d of %03d\n", x, hsynth->width) ;
      fflush(stdout) ;
    }
    if (x % 50 == 0)
    {
      char fname[STRLEN] ;
      int req = snprintf(fname, STRLEN, "%s.%03d.mgz", base_name, x) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing snapshot to %s\n", fname) ;
      MRIwrite(hsynth, fname) ;
    }
    for (y = 0 ; y < hsynth->height ; y++)
    {
      val = MRIgetVoxVal(histo, x, y, test_slice, 0) ;
      if (val == 4080 || val > 240)
      {
        MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
        continue ;
      }
      if ( x == Gx && y == Gy && test_slice == Gz)
        DiagBreak() ;
      if (identity)
      {
	xm = x ; ym = y ; zm = test_slice ;
      }
      else
      {
	V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
	MatrixMultiply(m_histo2mri, v1, v2) ;
	xm = nint(V3_X(v2)) ; ym = nint(V3_Y(v2)) ; zm = nint(V3_Z(v2)) ;  // MRI coord
      }

      extract_feature_vector(mri, xm, ym, zm, f) ;
      find_most_similar_location(mri_train_src, f, train_slice, &xd, &yd, &zd, mri_mask, &box_train, flags,
				 xind, yind, zind, nind, tol, num_notfound,
				 xm, ym, zm, 0) ;
      if (x == Gx && y == Gy && test_slice == Gz)
      {
	dump_window(mri, "test.dat", x, y, Gz, wsize) ;
	dump_window(mri_train_src, "train.dat", xd, yd, zd, wsize) ;
      }
      if (identity)
      {
	xh = xd ; yh = yd ; zh = zd ;
	val = MRIgetVoxVal(mri_train_dst, xh, yh, zh, 0) ;
      }
      else
      {
	V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ;
	MatrixMultiply(m_mri2histo, v1, v2) ;
	xh = V3_X(v2) ;  yh = V3_Y(v2) ;  zh = V3_Z(v2) ; 
	MRIsampleVolume(mri_train_dst, xh, yh, zh, &val) ;
      }

      if (val > 240)
	DiagBreak() ;
      MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  }
  printf("\n") ;

  for (x = 0 ; x < _MAX_FS_THREADS ; x++)
    free(farray[x].vals) ;
  MatrixFree(&m_histo2mri) ; MatrixFree(&m_mri2histo) ;MatrixFree(&v1) ; MatrixFree(&v2);
  return(hsynth) ;
}

static int
extract_feature_vector(MRI *mri, int x, int y, int z, FEATURE *f)
{
  int xi, yi, zi, xk, yk, zk,i  ;
  double val ;

  if (f->flags & SUBTRACT_CENTER)
    val = MRIgetVoxVal(mri, x, y, z, 0) ;
  else
    val = 0 ;
  for (i = 0, xk = -f->whalf ; xk <= f->whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -f->whalf ; yk <= f->whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (zk = -f->whalf ; zk <= f->whalf ; zk++, i++)
      {
        zi = mri->zi[z+zk] ;
        f->vals[i] = MRIgetVoxVal(mri, xi, yi, zi, 0) ;
        if (f->flags & SUBTRACT_CENTER)
          f->vals[i] -= val ;
      }
    }
  }

  return(NO_ERROR) ;
}


static int
find_most_similar_location(MRI *mri, FEATURE *fsrc, int z, int *pxd, int *pyd, int *pzd, 
                           MRI *mri_mask, MRI_REGION *box, int flags,
                           short *xind, short *yind, short *zind, int nind, double tol, int num_notfound,
			   int x0, int y0, int z0, double min_training_dist)
{
  FEATURE f ;
  double  dist, min_dist ;
  int     x, y, x1, y1, num = 0, ind ;

  f.len = fsrc->len ;
  f.whalf = fsrc->whalf ; ;
  f.flags = fsrc->flags ;
  f.vals = (double *)calloc(f.len, sizeof(double)) ;
  if (f.vals == NULL)
    ErrorExit(ERROR_NOMEMORY, "find_most_similar_location: could not allocate %d-len feature vector", f.len) ;

  min_dist = 1e12 ;

  x1 = box->x + box->dx - 1 ; y1 = box->y + box->dy - 1 ;
#if 0
#ifdef HAVE_OPENMP
  dist = 0 ; x = y = 0 ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(nind, x, y, x1, y1, Gx, Gy, Gz, box, x0, y0, fsrc, f, flags, min_training_dist, dist, num_notfound) shared(min_dist, num) schedule(static,1)
#endif
#endif
  for (ind = 0 ; ind < nind ; ind++)
  {
    ind = randomNumber(0, nind-.1) ;
    x = xind[ind] ;
    y = yind[ind] ;
    if (x < box->x || x > x1 || y < box->y || y > y1)
      continue ;
    if ( x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
      continue ;
    if ((mri->depth == 1) && 
	((abs(x-x0)<min_training_dist) || (abs(y-y0)<min_training_dist)))
      continue ;  // too close to test voxel
    extract_feature_vector(mri, x, y, z, &f) ;
    dist = feature_distance(fsrc, &f, flags & L1_NORM ? L1_NORM_DIST : L2_NORM_DIST) ;
    if (dist < min_dist)
    {
      num = 0 ;
      min_dist = dist ;
      *pxd = x ; *pyd = y ; *pzd = z ;
    }
    else
      if (num++ > num_notfound)
#ifdef HAVE_OPENMP
	continue ;
#else
        break;
#endif
  }

  free(f.vals) ;
  return(NO_ERROR) ;
}


static double
feature_distance(FEATURE *f1, FEATURE *f2, int which)
{
  int  i,num=0 ;
  double sse = 0.0;
  float  v1, v2 ;

  switch (which)
  {
  default:
  case L2_NORM_DIST:
    for (i = 0 ; i < f1->len ; i++)
    {
      v1 = f1->vals[i] ; v2 = f2->vals[i] ;
      if (FZERO(v1) || FZERO(v2) || FEQUAL(v1,255) || FEQUAL(v2,255))
	continue ;
      num++ ;
      sse += SQR(v1 - v2) ;
    }
    sse = sqrt(sse/num) ;
    break ;
  case L1_NORM_DIST:
    for (i = 0 ; i < f1->len ; i++)
    {
      v1 = f1->vals[i] ; v2 = f2->vals[i] ;
      if (FZERO(v1) || FZERO(v2) || FEQUAL(v1,255) || FEQUAL(v2,255))
	continue ;
      num++ ;
      sse += fabs(v1 - v2) ;
    }
    sse /= num ;
    break ;
  }
  

  return(sse)  ;
}
static int
prune_indices(MRI *mri_mask, MRI *histo, short *xind, short *yind, short *zind, int nind, int use_val, float hthresh)
{
  int isrc, idst, x, y, z, val, deleted ;
  short *xi, *yi, *zi ;
  float hval ;

  xi = (short *)calloc(nind, sizeof(short)) ;
  yi = (short *)calloc(nind, sizeof(short)) ;
  zi = (short *)calloc(nind, sizeof(short)) ;
  
  for (isrc = idst = deleted = 0 ; isrc < nind ; isrc++)
  {
    x = xind[isrc] ;  y = yind[isrc] ;  z = zind[isrc] ; 
    val = (int)MRIgetVoxVal(mri_mask, x, y, z,0) ;
    hval = (int)MRIgetVoxVal(histo, x, y, z,0) ;
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    if (val == use_val && hval < hthresh)
    {
      xi[idst] = x ; yi[idst] = y ; zi[idst] = z ;
      idst++ ;
    }
    else
      deleted++ ;
  }

  nind -= deleted ;
  for (isrc = 0 ; isrc < nind ; isrc++)
  {
    xind[isrc] = xi[isrc] ;  
    yind[isrc] = yi[isrc] ;  
    zind[isrc] = zi[isrc] ;  
  }
  free(xi) ; free(yi) ; free(zi) ;
  return(nind) ;
}
static int
dump_window(MRI *mri, const char *fname, int x0, int y0, int z0, int wsize) 
{
  int   xk, yk, xi, yi, whalf ;
  float val ;
  FILE  *fp ;

  fp = fopen(fname, "w") ;
  if (fp == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "dump_window(%s): could not open file",fname));
  whalf = (wsize-1)/2 ;
  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      xi = mri->xi[xk+x0] ;
      yi = mri->yi[yk+y0] ;
      if (xi >= 0 && yi >= 0 && xi < mri->width && yi < mri->height)
	val = MRIgetVoxVal(mri, xi, yi, z0, 0) ;
      else
	val = 0 ;
      fprintf(fp, "%f ", val) ;
    }
    fprintf(fp, "\n") ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

