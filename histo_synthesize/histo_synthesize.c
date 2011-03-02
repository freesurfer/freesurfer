/**
 * @file  histo_synthesize.c
 * @brief synthesize a histological volume from an MRI
 *
 * synthesize a histological volume from an MRI
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.4 $
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


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags, char *fname) ;
static int extract_feature_vector(MRI *mri, int x, int y, int z, FEATURE *feature) ;
static int find_most_similar_location(MRI *mri, FEATURE *fsrc, int z, int *pxd, int *pyd,int *pzd,
                                      MRI *mri_mask, MRI_REGION *box, int flags,
                                      short *xind, short *yind, short *zind, double tol,
                                      int num_notfound) ;
static double feature_distance(FEATURE *f1, FEATURE *f2, int which) ;

char *Progname ;
static void usage_exit(int code) ;

static int test_slice = 20 ;  // in histo coords
static int train_slice = 30 ; // in MRI coords
static int wsize = 3 ;
static int downsample = 0 ;
static double tol = 0 ;
static int num_notfound = 1000 ;  // # of search voxels to terminate after


#define SUBTRACT_CENTER  0x00001
#define L1_NORM          0x00002
#define MRI_SPACE        0x00004

#define L2_NORM_DIST     0
#define L1_NORM_DIST     1


static int flags = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI          *mri, *histo, *hsynth ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: histo_synthesize.c,v 1.4 2011/03/02 00:04:09 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

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
      mri_tmp = MRIdownsample2(mri, NULL) ;
      MRIfree(&mri) ; mri = mri_tmp ;
      train_slice /= 2 ;
      if (flags & MRI_SPACE)
        test_slice /= 2 ;

    }
  }
  histo = MRIread(argv[2]) ;
  if (histo == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read histological volume %s", Progname, argv[2]) ;
  hsynth = HISTOsynthesize(mri, histo, test_slice, train_slice, NULL, wsize, flags, argv[3]) ;
  printf("writing output to %s\n", argv[3]) ;
  MRIwrite(hsynth, argv[3]) ;
  msec = TimerStop(&start) ;
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
HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags, char *fname)
{
  int     x, y, z, xd, yd, zd, xm, ym, zm, x1, y1, z1 ;
  short   *xind, *yind, *zind ;
  double  val, xh, yh, zh ;
  FEATURE f ;
  MATRIX  *m_histo2mri, *m_mri2histo ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_mask ;
  MRI_REGION box_train, box_test ;

  f.len = wsize*wsize*wsize ;
  f.whalf = (wsize-1)/2 ;
  f.flags = flags ;
  f.vals = (double *)calloc(f.len, sizeof(double)) ;
  if (f.vals == NULL)
    ErrorExit(ERROR_NOMEMORY, "HISTOsynthesize: could not allocate %d-len feature vector", f.len) ;

  if (hsynth == NULL)
  {
    if (flags & MRI_SPACE)
      hsynth = MRIclone(mri, NULL) ;
    else
      hsynth = MRIclone(histo, NULL) ;
  }

  m_histo2mri = MRIgetVoxelToVoxelXform(histo, mri) ;
  m_mri2histo = MRIgetVoxelToVoxelXform(mri, histo) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;

  // compute region of MRI that maps to valid histo
  box_train.x = mri->width; box_train.y = mri->height ; box_train.z = mri->depth ;
  box_train.dx = 0 ; box_train.dy = 0 ; box_train.dz = 0 ;

  mri_mask = MRIclone(mri, NULL) ;
  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
        MatrixMultiply(m_mri2histo, v1, v2) ;
        xh = nint(V3_X(v2)) ; yh = nint(V3_Y(v2)) ; zh = nint(V3_Z(v2)) ;  // MRI coord
        MRIsampleVolume(histo, xh, yh, zh, &val) ;
        if ((val != 4080) && MRIindexNotInVolume(histo, xh, yh, zh) == 0)
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
  xind = (short *)calloc(hsynth->width*hsynth->height*hsynth->depth, sizeof(short)) ;
  yind = (short *)calloc(hsynth->width*hsynth->height*hsynth->depth, sizeof(short)) ;
  zind = (short *)calloc(hsynth->width*hsynth->height*hsynth->depth, sizeof(short)) ;
  MRIcomputeVoxelPermutation(hsynth, xind, yind,zind) ;
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
        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
        MatrixMultiply(m_mri2histo, v1, v2) ;
        xh = nint(V3_X(v2)) ; yh = nint(V3_Y(v2)) ; zh = nint(V3_Z(v2)) ;  // MRI coord
        MRIsampleVolume(histo, xh, yh, zh, &val) ;
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
        extract_feature_vector(mri, x, y, test_slice, &f) ;
        find_most_similar_location(mri, &f, train_slice, &xd, &yd, &zd, mri_mask, &box_train, flags,
                                   xind, yind, zind, tol, num_notfound) ;
        V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ;
        MatrixMultiply(m_mri2histo, v1, v2) ;
        xh = V3_X(v2) ;  yh = V3_Y(v2) ;  zh = V3_Z(v2) ; 
        MRIsampleVolume(histo, xh, yh, zh, &val) ;
#if 0
        if (val == 4080)
          val = 0 ;  // make background black
#endif
        MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
      }
    }  
  }
  else for (x = 0 ; x < hsynth->width ; x++)   // synth is in histo coords
  {
    if ((x % 10 == 0))
    {
      printf("%03d of %03d\n", x, hsynth->width) ;
      fflush(stdout) ;
    }
    for (y = 0 ; y < hsynth->height ; y++)
    {
      val = MRIgetVoxVal(histo, x, y, test_slice, 0) ;
      if (val == 4080)
      {
        MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
        continue ;
      }
      if ( x == Gx && y == Gy && test_slice == Gz)
        DiagBreak() ;
      V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
      MatrixMultiply(m_histo2mri, v1, v2) ;
      xm = nint(V3_X(v2)) ; ym = nint(V3_Y(v2)) ; zm = nint(V3_Z(v2)) ;  // MRI coord
      
      extract_feature_vector(mri, xm, ym, zm, &f) ;
      find_most_similar_location(mri, &f, train_slice, &xd, &yd, &zd, mri_mask, &box_train, flags,
                                 xind, yind, zind, tol, num_notfound) ;
      V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ;
      MatrixMultiply(m_mri2histo, v1, v2) ;
      xh = V3_X(v2) ;  yh = V3_Y(v2) ;  zh = V3_Z(v2) ; 
      MRIsampleVolume(histo, xh, yh, zh, &val) ;
      MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
    }
  }
  printf("\n") ;

  free(f.vals) ;
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
                           short *xind, short *yind, short *zind, double tol, int num_notfound)
{
  FEATURE f ;
  double  dist, min_dist ;
  int     x, y, x1, y1, num, ind, nind ;

  f.len = fsrc->len ;
  f.whalf = fsrc->whalf ; ;
  f.flags = fsrc->flags ;
  f.vals = (double *)calloc(f.len, sizeof(double)) ;
  if (f.vals == NULL)
    ErrorExit(ERROR_NOMEMORY, "find_most_similar_location: could not allocate %d-len feature vector", f.len) ;

  min_dist = 1e12 ;

  nind = mri->width*mri->height*mri->depth ;
  x1 = box->x + box->dx - 1 ; y1 = box->y + box->dy - 1 ;
  for (num = ind = 0 ; ind < nind ; ind++)
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
        break;
  }

  free(f.vals) ;
  return(NO_ERROR) ;
}


static double
feature_distance(FEATURE *f1, FEATURE *f2, int which)
{
  int  i ;
  double sse ;

  switch (which)
  {
  default:
  case L2_NORM_DIST:
    for (sse = 0.0, i = 0 ; i < f1->len ; i++)
      sse += SQR(f1->vals[i] - f2->vals[i]) ;
    sse = sqrt(sse/f1->len) ;
    break ;
  case L1_NORM_DIST:
    for (sse = 0.0, i = 0 ; i < f1->len ; i++)
      sse += fabs(f1->vals[i] - f2->vals[i]) ;
    sse /= f1->len ;
    break ;
  }
  

  return(sse)  ;
}
