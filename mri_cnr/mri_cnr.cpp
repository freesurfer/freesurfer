/**
 * @brief use white and pial surfaces to compute the CNR of a volume
 *
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri_conform.h"
#include "version.h"
#include "label.h"


int main(int argc, char *argv[]) ;

static int MRIScomputeSlope(MRI_SURFACE *mris, MRI *mri, double dist_in, double dist_out,
                            double step_in, double step_out, MATRIX *m_linfit);
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_version(void) ;
static double compute_volume_cnr(MRI_SURFACE *mris, MRI *mri, char *log_fname) ;
const char *Progname ;
static char *log_fname = NULL ;

static char *slope_fname = NULL ;
static double dist_in, dist_out, step_in, step_out ;
static int interp = SAMPLE_TRILINEAR ;

static LABEL *lh_area, *rh_area ;

static int only_total = 0 ;

int
main(int argc, char *argv[]) {
  char        **av, *mri_name,  fname[STRLEN],*path ;
  const char* hemi;
  int         ac, nargs, i, j ;
  MRI         *mri, *mri_template = NULL, *mri_tmp ;
  MRI_SURFACE *mris ;
  double      cnr_total, cnr = 0.0 ;

  nargs = handleVersionOption(argc, argv, "mri_cnr");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  path = argv[1] ;
  for (cnr_total = 0.0, j = 0 ; j <= 1 ; j++) {
    if (j == LEFT_HEMISPHERE)
      hemi = "lh" ;
    else
      hemi = "rh" ;
    sprintf(fname, "%s/%s.white", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    sprintf(fname, "%s/%s.pial", path, hemi) ;
    if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;

    if (lh_area && j == LEFT_HEMISPHERE)
      LabelRipRestOfSurface(lh_area, mris) ;
    else if (rh_area && j == RIGHT_HEMISPHERE)
      LabelRipRestOfSurface(rh_area, mris) ;

    for (i = 2 ; i < argc ; i++) {
      mri_name = argv[i] ;
      if (!j)
        fprintf(stderr, "processing MRI volume %s...\n", mri_name) ;
      mri = MRIread(mri_name) ;
      if (!mri)
        ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, mri_name) ;
      if (!j && !FZERO(mri->tr))
        fprintf(stderr, "TR = %2.1f msec, flip angle = %2.0f degrees, TE = %2.1f msec\n",
               mri->tr, DEGREES(mri->flip_angle), mri->te) ;

      if (0)
      {
        if (!mri_template)
          mri_template = MRIconform(mri) ;
        mri_tmp = MRIresample(mri, mri_template, SAMPLE_NEAREST) ;
        MRIfree(&mri) ;
        mri = mri_tmp ;
      }

      if (j == LEFT_HEMISPHERE)
      {
        cnr = compute_volume_cnr(mris, mri, log_fname) ;
	if (only_total == 0)
	  printf("%s CNR = %2.3f\n", hemi, cnr) ;
      }
      else {
        double rh_cnr ;
        rh_cnr = compute_volume_cnr(mris, mri, log_fname) ;
        cnr = (cnr + rh_cnr) / 2.0 ;
	if (only_total == 0)
	  printf("%s CNR = %2.3f\n", hemi, rh_cnr) ;
      }

      if (slope_fname)
      {
        MATRIX *m_linfit ;

        m_linfit = MatrixAlloc(mris->nvertices, 2, MATRIX_REAL) ;
        if (m_linfit == NULL)
          ErrorExit(ERROR_NOMEMORY, "%s: could not allocate slope/offset matrix", Progname) ;
        MRIScomputeSlope(mris, mri, dist_in, dist_out, step_in, step_out, m_linfit);
        MRISimportValFromMatrixColumn(mris, m_linfit, 1) ;
        sprintf(fname, "%s/%s.%s.slope.mgz", path, hemi, slope_fname) ;
        MRISwriteValues(mris, fname) ;
        MRISimportValFromMatrixColumn(mris, m_linfit, 2) ;
        sprintf(fname, "%s/%s.%s.offset.mgz", path, hemi, slope_fname) ;
        MRISwriteValues(mris, fname) ;


        MatrixFree(&m_linfit) ;
      }
      MRIfree(&mri) ;
    }
    cnr_total += cnr ;
    MRISfree(&mris) ;
  }
  if (only_total == 0)
    printf("total CNR = %2.3f\n", cnr/(double)((argc-2))) ;
  else
    printf("%2.3f\n", cnr/(double)((argc-2))) ;

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
    usage_exit() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "label"))
  {
    fprintf(stderr, "reading lh and rh labels from %s and %s\n", argv[2], argv[3]) ;
    lh_area = LabelRead(NULL, argv[2]) ;
    if (lh_area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label %s", Progname, argv[2]) ;
    rh_area = LabelRead(NULL, argv[3]) ;
    if (rh_area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label %s", Progname, argv[3]) ;
    nargs = 2 ;
  }
  else switch (toupper(*option)) {
  case 'L':
    log_fname = argv[2] ;
    printf("logging g/w cnr to %s\n", log_fname) ;
    nargs = 1 ;
    break ;
  case 'S':
    slope_fname = argv[2];
    dist_in = atof(argv[3]) ;
    dist_out = atof(argv[4]) ;
    step_in = atof(argv[5]) ;
    step_out = atof(argv[6]) ;
    nargs = 5 ;
    break ;
    case 'T':
      only_total = 1 ;
      fprintf(stderr, "stdout will only have total CNR\n") ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      usage_exit() ;
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
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "%s -- compute the gray/white/csf contrast-to-noise ratio for volumes.\n",
          Progname) ;
  fprintf(stderr,
          "usage: %s [options] <surf directory> <vol 1> <vol 2> ...\n",
          Progname) ;
  fprintf(stderr,
          "usage example (assumes fs pipeline has finished for subject subj1): %s subj1/surf subj1/mri/orig.mgz\n",
          Progname) ;
  fprintf(stderr, "Available options:\n") ;
  fprintf(stderr,
          "\t-s <slope_fname> <dist in> <dist out> <step in> <step out>: compute slope based on given values, write it to slope and offset files labeled <slope_fname> (e.g., `lh.<slope_fname>.slope.mgz')\n") ;
  fprintf(stderr,
          "\t-t : print only the total CNR to stdout (stderr still contains more information)\n") ;
  fprintf(stderr,
          "\t-l <logfile>: log cnr to file <logfile>. Will contain 8 values in the following order: gray_white_cnr, gray_csf_cnr, white_mean, gray_mean, csf_mean, sqrt(white_var), sqrt(gray_var), sqrt(csf_var)\n") ;
  fprintf(stderr,
            "\tlabel <lh> <rh>: read hemisphere labels from <lh> and <rh>\n") ;
  fprintf(stderr,
          "\t-u, -?, -help : print usage information and quit\n") ;
  fprintf(stderr,
          "\t-version : print software version information and quit\n") ;
}



static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
/*
  orig vertex positions = white
  current vertex positions = pial
*/
static double
compute_volume_cnr(MRI_SURFACE *mris, MRI *mri, char *log_fname) {
  double  gray_white_cnr, gray_csf_cnr, gray_var, white_var, csf_var, gray_mean, white_mean, csf_mean ;
  float   thickness ;
  double  x, y, z, gray, white, csf ;
  int     vno ;
  VERTEX  *v ;

//  MRISsetVolumeForSurface(mris, mri);
  MRIScomputeMetricProperties(mris) ;
  gray_mean = gray_var = white_mean = white_var = csf_mean = csf_var = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    //MRIworldToVoxel(mri, v->x+v->nx,
    MRISsurfaceRASToVoxelCached(mris, mri, v->x+v->nx, v->y+v->ny, v->z+v->nz, &x, &y, &z) ;
//    MRIsurfaceRASToVoxel(mri, v->x+v->nx, v->y+v->ny, v->z+v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &csf) ;
    thickness = v->curv*.5 ;
    // MRIworldToVoxel(mri, v->x-thickness*v->nx, v->y-thickness*v->ny, v->z-thickness*v->nz, &x, &y, &z) ;
//    MRIsurfaceRASToVoxel(mri, v->x-thickness*v->nx, v->y-thickness*v->ny, v->z-thickness*v->nz, &x, &y, &z) ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x-thickness*v->nx, v->y-thickness*v->ny, v->z-thickness*v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &gray) ;

    gray_var += (gray*gray) ;
    gray_mean += gray ;
    csf_var += (csf*csf) ;
    csf_mean += csf ;
  }

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  for (white_var = gray_white_cnr = 0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    // MRIworldToVoxel(mri, v->x-v->nx, v->y-v->ny, v->z-v->nz, &x, &y, &z) ;
//    MRIsurfaceRASToVoxel(mri, v->x-v->nx, v->y-v->ny, v->z-v->nz, &x, &y, &z) ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x-v->nx, v->y-v->ny, v->z-v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &white) ;
    thickness = v->curv*.5 ;
    // MRIworldToVoxel(mri, v->x+thickness*v->nx, v->y+thickness*v->ny, v->z+thickness*v->nz, &x, &y, &z) ;
//    MRIsurfaceRASToVoxel(mri, v->x+thickness*v->nx, v->y+thickness*v->ny, v->z+thickness*v->nz, &x, &y, &z) ;
    MRISsurfaceRASToVoxelCached(mris, mri, v->x+thickness*v->nx, v->y+thickness*v->ny, v->z+thickness*v->nz, &x, &y, &z) ;
    MRIsampleVolume(mri, x, y, z, &gray) ;

    gray_var += (gray*gray) ;
    gray_mean += gray ;
    white_var += (white*white) ;
    white_mean += white ;
  }

  white_mean /= (double)mris->nvertices ;
  csf_mean /= (double)mris->nvertices ;
  gray_mean /= (double)mris->nvertices*2.0 ;

  white_var = white_var / (double)mris->nvertices - white_mean*white_mean ;
  gray_var = gray_var / ((double)mris->nvertices*2.0) - gray_mean*gray_mean ;
  csf_var = csf_var / (double)mris->nvertices - csf_mean*csf_mean ;

  if (only_total == 0)
    printf("\twhite = %2.1f+-%2.1f, gray = %2.1f+-%2.1f, csf = %2.1f+-%2.1f\n",
	   white_mean, sqrt(white_var), gray_mean, sqrt(gray_var),
	   csf_mean, sqrt(csf_var)) ;

  gray_white_cnr = SQR(gray_mean - white_mean) / (gray_var+white_var) ;
  gray_csf_cnr = SQR(gray_mean - csf_mean) / (gray_var+csf_var) ;

  if (only_total == 0)
    printf("\tgray/white CNR = %2.3f, gray/csf CNR = %2.3f\n",
	   gray_white_cnr, gray_csf_cnr) ;

  if (log_fname)
  {
    FILE    *fp ;
    fp = fopen(log_fname, "a") ;
    if (fp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open file %s", Progname, log_fname) ;
    fprintf(fp, "%f %f %f %f %f %f %f %f\n", gray_white_cnr, gray_csf_cnr, white_mean, gray_mean, csf_mean,
            sqrt(white_var), sqrt(gray_var), sqrt(csf_var)) ;

  }
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

  return((gray_white_cnr + gray_csf_cnr)/2.0) ;
}

static int
MRIScomputeSlope(MRI_SURFACE *mris, MRI *mri, double dist_in, double dist_out,
                 double step_in, double step_out, MATRIX *m_linfit)
{
  int      vno, nsamples, n ;
  VERTEX   *v ;
  MATRIX   *m_X, *m_Y, *m_P, *m_pinv = NULL ;
  double   d, xv, yv, zv, x, y, z, val ;

  nsamples = 1 + (int)(dist_in / step_in) + (int)(dist_out / step_out) ;
  m_X = MatrixAlloc(2, nsamples, MATRIX_REAL) ;
  m_Y = MatrixAlloc(1, nsamples, MATRIX_REAL) ;
  m_P = MatrixAlloc(1, 2, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0, d = -dist_in ; d < 0.0 ; d += step_in, n++)
    {
      x = v->whitex + d*v->nx ; y = v->whitey + d*v->ny ; z = v->whitez + d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
      MRIsampleVolumeType(mri, xv, yv, zv, &val, interp) ;
      *MATRIX_RELT(m_Y, 1, n+1) = val ;
      *MATRIX_RELT(m_X, 1, n+1) = d ;
      *MATRIX_RELT(m_X, 2, n+1) = 1 ;
    }
    n++ ;   // include val at 0
    x = v->whitex ; y = v->whitey ; z = v->whitez ;
    MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri, xv, yv, zv, &val, interp) ;
    *MATRIX_RELT(m_Y, 1, n+1) = val ;
    *MATRIX_RELT(m_X, 1, n+1) = 0.0 ;
    *MATRIX_RELT(m_X, 2, n+1) = 1 ;

    for (d = step_out ; d <= dist_out ; d += step_out, n++)
    {
      x = v->whitex + d*v->nx ; y = v->whitey + d*v->ny ; z = v->whitez + d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv) ;
      MRIsampleVolumeType(mri, xv, yv, zv, &val, interp) ;
      *MATRIX_RELT(m_Y, 1, n+1) = val ;
      *MATRIX_RELT(m_X, 1, n+1) = d ;
      *MATRIX_RELT(m_X, 2, n+1) = 1 ;
    }

    m_pinv = MatrixPseudoInverse(m_X, m_pinv) ;
    MatrixMultiply(m_Y, m_pinv, m_P) ;
    *MATRIX_RELT(m_linfit, vno+1, 1) = *MATRIX_RELT(m_P, 1, 1) ;
    *MATRIX_RELT(m_linfit, vno+1, 2) = *MATRIX_RELT(m_P, 1, 2) ;
    MatrixFree(&m_pinv) ;
  }

  return(NO_ERROR) ;
}
