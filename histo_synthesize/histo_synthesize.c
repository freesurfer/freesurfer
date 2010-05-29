/**
 * @file  histo_synthesize.c
 * @brief synthesize a histological volume from an MRI
 *
 * synthesize a histological volume from an MRI
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2010/05/29 23:30:21 $
 *    $Revision: 1.1 $
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

typedef struct
{
  int    len ;
  int    whalf ;
  double *vals ;
  int    flags ;
} FEATURE ;


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags) ;
static int extract_feature_vector(MRI *mri, int x, int y, int z, FEATURE *feature) ;
static int find_most_similar_location(MRI *mri, FEATURE *fsrc, int z, int *pxd, int *pyd, int *pzd, MRI *mri_mask,
                                      MRI_REGION *box, int flags) ;
static double feature_distance(FEATURE *f1, FEATURE *f2, int which) ;

char *Progname ;
static void usage_exit(int code) ;

static int test_slice = 20 ;  // in histo coords
static int train_slice = 51 ; // in MRI coords
static int wsize = 5 ;


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
  nargs = handle_version_option (argc, argv, "$Id: histo_synthesize.c,v 1.1 2010/05/29 23:30:21 fischl Exp $", "$Name:  $");
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
  histo = MRIread(argv[2]) ;
  hsynth = HISTOsynthesize(mri, histo, test_slice, train_slice, NULL, wsize, flags) ;
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
HISTOsynthesize(MRI *mri, MRI *histo, int test_slice, int train_slice, MRI *hsynth, int wsize, int flags)
{
  int     x, y, z, xd, yd, zd, xm, ym, zm, x1, y1, z1 ;
  double  val, xh, yh, zh ;
  FEATURE f ;
  MATRIX  *m_histo2mri, *m_mri2histo ;
  VECTOR  *v1, *v2 ;
  MRI     *mri_mask ;
  MRI_REGION box ;

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
  box.x = mri->width; box.y = mri->height ; box.z = mri->depth ;
  box.dx = 0 ; box.dy = 0 ; box.dz = 0 ;

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
          x1 = box.x + box.dx-1 ;
          y1 = box.y + box.dy-1 ;
          z1 = box.z + box.dz-1 ;
          if (x < box.x)
            box.x = x ;
          if (y < box.y)
            box.y = y ;
          if (z < box.z)
            box.z = z ;
          if (x > x1)
            x1 = x ;
          if (y > y1)
            y1 = y ;
          if (z > z1)
            z1 = z ;
          box.dx = x1 - box.x + 1 ;
          box.dy = y1 - box.y + 1 ;
          box.dz = z1 - box.z + 1 ;
          MRIsetVoxVal(mri_mask, x, y, z, 0, 1) ;
        }
      }
    }
  }
  if (flags & MRI_SPACE) for (x = 0 ; x < hsynth->width ; x++)
  {
    if ((x % 10 == 0))
    {
      printf("%03d of %03d\n", x, hsynth->width) ;
      fflush(stdout) ;
    }
    for (y = 0 ; y < hsynth->height ; y++)
    {
      if ( x == Gx && y == Gy && test_slice == Gz)
        DiagBreak() ;
      if (MRIgetVoxVal(mri_mask, x, y, test_slice,0) == 0)
        continue ;
      V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = test_slice ;
      extract_feature_vector(mri, x, y, test_slice, &f) ;
      find_most_similar_location(mri, &f, train_slice, &xd, &yd, &zd, mri_mask, &box, flags) ;
      V3_X(v1) = xd ; V3_Y(v1) = yd ; V3_Z(v1) = zd ;
      MatrixMultiply(m_mri2histo, v1, v2) ;
      xh = V3_X(v2) ;  yh = V3_Y(v2) ;  zh = V3_Z(v2) ; 
      MRIsampleVolume(histo, xh, yh, zh, &val) ;
      if (val == 4080)
        val = 0 ;  // make background black
      MRIsetVoxVal(hsynth, x, y, test_slice, 0, val) ;
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
      find_most_similar_location(mri, &f, train_slice, &xd, &yd, &zd, mri_mask, &box, flags) ;
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
find_most_similar_location(MRI *mri, FEATURE *fsrc, int z, int *pxd, int *pyd, int *pzd, MRI *mri_mask, MRI_REGION *box, int flags)
{
  FEATURE f ;
  double  dist, min_dist ;
  int     x, y, x1, y1 ;

  f.len = fsrc->len ;
  f.whalf = fsrc->whalf ; ;
  f.flags = fsrc->flags ;
  f.vals = (double *)calloc(f.len, sizeof(double)) ;
  if (f.vals == NULL)
    ErrorExit(ERROR_NOMEMORY, "find_most_similar_location: could not allocate %d-len feature vector", f.len) ;

  min_dist = 1e12 ;

  x1 = box->x + box->dx - 1 ; y1 = box->y + box->dy - 1 ;
  for (x = box->x ; x <= x1 ; x++)
    for (y = box->y ; y <= y1 ; y++)
    {
      if ( x == Gx && y == Gy && z == Gz)
        DiagBreak() ;
      if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
        continue ;
      extract_feature_vector(mri, x, y, z, &f) ;
      dist = feature_distance(fsrc, &f, flags & L1_NORM ? L1_NORM_DIST : L2_NORM_DIST) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        *pxd = x ; *pyd = y ; *pzd = z ;
      }
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
