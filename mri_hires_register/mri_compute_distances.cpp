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


//
// mri_compute_distances.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
//
// Warning: Do not edit the following four lines.  CVS maintains them.
//
////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gcamorph.h"
#include "mri.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "utils.h"
#include "matrix.h"
#include "gca.h"
#include "mrisegment.h"
#include "cma.h"
#include "version.h"

const char *Progname ;
static void usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;
static float MRIcomputeAverageHausdorffDistance(MRI *mri1,
    MRI *mri2,
    int label,
    double *psigma) ;
static float MRIcomputeSegmentPairHausdorffDistance(MRI *mri1,
    MRI *mri2,
    MRI_SEGMENT *mseg1,
    MRI_SEGMENT *mseg2) ;
static float MRIcomputeSegmentPairCentroidDistance(MRI *mri1,
    MRI *mri2,
    MRI_SEGMENT *mseg1,
    MRI_SEGMENT *mseg2);

#define MAX_TRANS 1000
static int translate_in[MAX_TRANS] ;
static int translate_out[MAX_TRANS] ;
static int ntrans = 0 ;
static int centroid = 0 ;

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  MRI          *mri1, *mri2 ;
  Timer start ;
  int          msec, minutes, seconds, i, label ;
  double       dist, sigma = 0;

  start.reset() ;
  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

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

  fprintf(stderr, "reading %s...\n", argv[1]) ;
  mri1 = MRIread(argv[1]) ;
  if (!mri1)
    ErrorExit(ERROR_NOFILE, "%s: could not read first label volume %s",
              Progname, argv[1]) ;
  fprintf(stderr, "reading %s...\n", argv[2]) ;
  mri2 = MRIread(argv[2]) ;
  if (!mri2)
    ErrorExit(ERROR_NOFILE, "%s: could not read second label volume %s",
              Progname, argv[2]) ;


  for (i = 0 ; i < ntrans ; i++) {
    fprintf(stderr, "translating label %s (%d) to %s (%d)\n",
            cma_label_to_name(translate_in[i]), translate_in[i],
            cma_label_to_name(translate_out[i]), translate_out[i]) ;
    MRIreplaceValues(mri1, mri1, translate_in[i], translate_out[i]) ;
    MRIreplaceValues(mri2, mri2, translate_in[i], translate_out[i]) ;
  }
  for (i = 3 ; i < argc ; i++) {
    label = atoi(argv[i]) ;
    fprintf(stderr, "computing %s distance for label %s (%d)...\n",
            centroid ? "centroid" : "Hausdorff",
            cma_label_to_name(label), label) ;
    dist = MRIcomputeAverageHausdorffDistance(mri1, mri2, label, &sigma) ;
    printf("%d %s %f %f\n", label, cma_label_to_name(label), dist, sigma) ;
    fprintf(stderr,"%d %s %f +- %2.3f\n",
            label, cma_label_to_name(label), dist,sigma) ;
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "distance calculation took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
}


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (*option) {
    case 'L':
      printf("%s ", argv[2]) ;   // put it into stdout for the log file
      nargs = 1 ;
      break ;
    case 'T':
      translate_in[ntrans] = atoi(argv[2]) ;
      translate_out[ntrans] = atoi(argv[3]) ;
      nargs = 2 ;
      fprintf(stderr, "translating label %s (%d) to %s (%d)\n",
              cma_label_to_name(translate_in[ntrans]),
              translate_in[ntrans],
              cma_label_to_name(translate_out[ntrans]),
              translate_out[ntrans]) ;
      ntrans++ ;
      break ;
    case 'C':
      centroid = 1 ;
      fprintf(stderr, "computing centroid distances\n") ;
      break ;
    case '?':
    case 'U':
      usage_exit(1);
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      usage_exit(1) ;
      break ;
    }
  return(nargs) ;
}


static void
usage_exit(int ecode) {
  printf("usage: %s <source> <target> <output xform>\n",Progname) ;
  exit(ecode) ;
}


static float
MRIcomputeAverageHausdorffDistance(MRI *mri1,
                                   MRI *mri2,
                                   int label,
                                   double *psigma) {
  double           total_dist, min_dist, dist, dsq, sigma ;
  MRI_SEGMENTATION *mriseg1, *mriseg2 ;
  MRI_SEGMENT      *mseg1, *mseg2 ;
  int              s1, s2 ;

  mriseg1 = MRIsegment(mri1, label, label) ;
  mriseg2 = MRIsegment(mri2, label, label) ;
  MRIremoveSmallSegments(mriseg1, 100) ;
  MRIremoveSmallSegments(mriseg2, 100) ;
  if (mriseg1->nsegments == 0 || mriseg2->nsegments == 0)
    ErrorReturn(-1.0, (ERROR_BADPARM,
                       "MRIcomputeAverageHausdorffDistance: "
                       "invalid nseg %d or %d",
                       mriseg1->nsegments, mriseg2->nsegments)) ;
  fprintf(stderr, "vol 1 with %d segments, vol 2 with %d segments\n",
          mriseg1->nsegments, mriseg2->nsegments) ;


  for (dsq = 0.0, total_dist = 0.0, s1 = 1 ; s1 < mriseg1->nsegments ; s1++) {
    mseg1 = &mriseg1->segments[s1] ;
    for (min_dist = -1.0, s2 = 0 ; s2 < mriseg2->nsegments ; s2++) {
      mseg2 = &mriseg2->segments[s2] ;
      if (centroid)
        dist =
          MRIcomputeSegmentPairCentroidDistance(mri1, mri2, mseg1, mseg2);
      else
        dist =
          MRIcomputeSegmentPairHausdorffDistance(mri1, mri2, mseg1, mseg2);
      if (dist < min_dist || min_dist < 0)
        min_dist = dist ;
    }

    if (min_dist >= 0) {
      dsq += (min_dist * min_dist) ;
      total_dist += min_dist ;
    }
  }

  total_dist /= mriseg1->nsegments ;  // mean
  sigma =
    sqrt(dsq/mriseg1->nsegments - total_dist*total_dist) ; // standard dev
  if (psigma)
    *psigma = sigma / sqrt(mriseg1->nsegments) ;

  MRIsegmentFree(&mriseg1) ;
  MRIsegmentFree(&mriseg2) ;


  return(total_dist) ;
}


static float
MRIcomputeSegmentPairHausdorffDistance(MRI *mri1,
                                       MRI *mri2,
                                       MRI_SEGMENT *mseg1,
                                       MRI_SEGMENT *mseg2) {
  double hdist = 0.0, dist, min_dist, dx, dy, dz, x1, y1, z1, x2, y2, z2 ;
  int    i1, i2 ;
  MATRIX *m_vox2ras1, *m_vox2ras2, *vvox, *vras ;

  m_vox2ras1 = MRIgetVoxelToRasXform(mri1) ;
  m_vox2ras2 = MRIgetVoxelToRasXform(mri2) ;
  vvox = VectorAlloc(4, MATRIX_REAL) ;
  vras = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(vvox,4,1) = 1.0 ;
  *MATRIX_RELT(vras,4,1) = 1.0 ;

  for (i1 = 0 ; i1 < mseg1->nvoxels ; i1++) {
    min_dist = -1 ;
    V3_X(vvox) = mseg1->voxels[i1].x ;
    V3_Y(vvox) = mseg1->voxels[i1].y ;
    V3_Z(vvox) = mseg1->voxels[i1].z ;
    MatrixMultiply(m_vox2ras1, vvox, vras) ;
    x1 = V3_X(vras) ;
    y1 = V3_Y(vras) ;
    z1 = V3_Z(vras) ;
    for (i2 = 0 ; i2 < mseg2->nvoxels ; i2++) {
      V3_X(vvox) = mseg2->voxels[i2].x ;
      V3_Y(vvox) = mseg2->voxels[i2].y ;
      V3_Z(vvox) = mseg2->voxels[i2].z ;
      MatrixMultiply(m_vox2ras2, vvox, vras) ;
      x2 = V3_X(vras) ;
      y2 = V3_Y(vras) ;
      z2 = V3_Z(vras) ;
      dx = x2-x1 ;
      dy = y2-y1 ;
      dz = z2-z1 ;
      dist = dx*dx + dy*dy + dz*dz ;
      if (dist < min_dist || i2 == 0)
        min_dist = dist ;
    }
    if (min_dist > hdist)
      hdist = min_dist ;
  }
  MatrixFree(&m_vox2ras1) ;
  MatrixFree(&m_vox2ras2) ;
  MatrixFree(&vvox) ;
  MatrixFree(&vras) ;
  return(sqrt(hdist)) ;
}


static float
MRIcomputeSegmentPairCentroidDistance(MRI *mri1,
                                      MRI *mri2,
                                      MRI_SEGMENT *mseg1,
                                      MRI_SEGMENT *mseg2) {
  double dist, dx, dy, dz, x1, y1, z1, x2, y2, z2 ;
  MATRIX *m_vox2ras1, *m_vox2ras2, *vvox, *vras ;

  m_vox2ras1 = MRIgetVoxelToRasXform(mri1) ;
  m_vox2ras2 = MRIgetVoxelToRasXform(mri2) ;
  vvox = VectorAlloc(4, MATRIX_REAL) ;
  vras = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(vvox,4,1) = 1.0 ;
  *MATRIX_RELT(vras,4,1) = 1.0 ;

  V3_X(vvox) = mseg1->cx ;
  V3_Y(vvox) = mseg1->cy ;
  V3_Z(vvox) = mseg1->cz ;
  MatrixMultiply(m_vox2ras1, vvox, vras) ;
  x1 = V3_X(vras) ;
  y1 = V3_Y(vras) ;
  z1 = V3_Z(vras) ;

  V3_X(vvox) = mseg2->cx ;
  V3_Y(vvox) = mseg2->cy ;
  V3_Z(vvox) = mseg2->cz ;
  MatrixMultiply(m_vox2ras2, vvox, vras) ;
  x2 = V3_X(vras) ;
  y2 = V3_Y(vras) ;
  z2 = V3_Z(vras) ;
  dx = x2-x1 ;
  dy = y2-y1 ;
  dz = z2-z1 ;
  dist = dx*dx + dy*dy + dz*dz ;

  MatrixFree(&m_vox2ras1) ;
  MatrixFree(&m_vox2ras2) ;
  MatrixFree(&vvox) ;
  MatrixFree(&vras) ;
  return(sqrt(dist)) ;
}

