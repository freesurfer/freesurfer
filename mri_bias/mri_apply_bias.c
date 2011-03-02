/**
 * @file  mri_apply_bias.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
 *    $Revision: 1.5 $
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
static MRI *apply_bias(MRI *mri_orig, MRI *mri_norm, MRI *mri_bias) ;

char *Progname ;

static void usage_exit(int code) ;
static char *xform_fname = NULL ;

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  MRI          *mri_orig, *mri_norm, *mri_bias ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mri_apply_bias.c,v 1.5 2011/03/02 00:04:13 nicks Exp $",
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
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  mri_orig = MRIread(argv[1]) ;
  if (mri_orig == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s",Progname, argv[1]) ;
  mri_bias = MRIread(argv[2]) ;
  if (mri_bias == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s",Progname, argv[2]) ;

  if (xform_fname)
  {
    TRANSFORM    *transform ;
    MRI          *mri ;

    transform = TransformRead(xform_fname) ;
    if (transform == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load transform %s", Progname, xform_fname) ;
    mri = MRIcloneDifferentType(mri_orig, MRI_FLOAT) ;
    TransformApplyInverse(transform, mri_bias, mri) ;
    MRIfree(&mri_bias) ; mri_bias = mri ;
    
  }
  mri_norm = apply_bias(mri_orig, NULL, mri_bias);

  fprintf(stderr, "writing to %s...\n", argv[3]) ;
  MRIwrite(mri_norm, argv[3]) ;
  MRIfree(&mri_norm) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "bias correction took %d minutes and %d seconds.\n",
          minutes, seconds) ;
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
  if (!stricmp(option, "sdir")) {}
  else if (!stricmp(option, "debug_voxel")) 
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else switch (toupper(*option)) {
  case 'T':
    xform_fname = argv[2] ;
    nargs = 1 ;
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
  printf("usage: %s [options] <input vol> <bias vol> <output volume>\n", Progname) ;
  exit(code) ;
}
static MRI *
apply_bias(MRI *mri_orig, MRI *mri_norm, MRI *mri_bias) {
  MATRIX   *m_vox2vox;
  VECTOR   *v1, *v2;
  int      x, y, z ;
  double   xd, yd, zd, bias, val_orig, val_norm ;

  if (mri_norm == NULL)
    mri_norm = MRIclone(mri_orig, NULL) ;

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_orig, mri_bias) ;
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0 ;
  VECTOR_ELT(v2, 4) = 1.0 ;

  for (x = 0 ; x < mri_orig->width ; x++) {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_orig->height ; y++) {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_orig->depth ; z++) {
        V3_Z(v1) = z ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val_orig = MRIgetVoxVal(mri_orig, x, y, z, 0) ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xd = V3_X(v2) ;
        yd = V3_Y(v2) ;
        zd = V3_Z(v2);
        MRIsampleVolume(mri_bias, xd, yd, zd, &bias) ;
        val_norm = val_orig * bias ;
        if (mri_norm->type == MRI_UCHAR) {
          if (val_norm > 255)
            val_norm = 255 ;
          else if (val_norm < 0)
            val_norm = 0 ;
        }
        MRIsetVoxVal(mri_norm, x, y, z, 0, val_norm) ;
      }
    }
  }

  MatrixFree(&m_vox2vox) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  return(mri_norm) ;
}
