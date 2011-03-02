/**
 * @file  mri_compute_structure_transforms.c
 * @brief compute optimal linear transform for each struct to the atlas
 *
 * use an SVG pseudo inverse to
 * compute optimal linear transform for each struct to the atlas in a
 * previously compute .m3z file.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.3 $
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
#include "gcamorph.h"
#include "transform.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs, x, y, z, n, i, label ;
  char         *m3z_fname, *out_dir, fname[STRLEN] ;
  int          msec, minutes, seconds, label_counts[MAX_CMA_LABELS] ;
  struct timeb start ;
  GCA_MORPH      *gcam ;
  TRANSFORM      _transform, *transform = &_transform ;
  MATRIX         *mX, *mY, *m, *mXinv ;
  LTA            *lta ;
  MRI            *mri_aseg ;
  float          xa, ya, za ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_compute_structure_transforms.c,v 1.3 2011/03/02 00:04:14 nicks Exp $", "$Name:  $");
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
   m3z_fname = argv[1] ;
  out_dir = argv[3] ;
  gcam = GCAMread(m3z_fname) ;
  if (gcam == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not ready .m3z from %s", Progname, m3z_fname);
  transform->xform = (void *)gcam ;
  transform->type = MORPH_3D_TYPE ;
  mri_aseg = MRIread(argv[2]) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_BADFILE, "%s: could not read aseg from %s\n", Progname, argv[2]) ;

  memset(label_counts, 0, sizeof(label_counts)) ;
  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
        label_counts[label]++ ;
      }

  TransformInvert(transform, mri_aseg) ;
  lta = LTAalloc(1, NULL) ;
  for (n = 0 ; n < MAX_CMA_LABELS ; n++)
  {
    if (label_counts[n] > 12)
    {
      mX = MatrixAlloc(4, label_counts[n], MATRIX_REAL) ;
      mY = MatrixAlloc(4, label_counts[n], MATRIX_REAL) ;
      printf("computing linear approximation for %s (%d)\n", cma_label_to_name(n), n) ;

      for (i = 1, x = 0 ; x < mri_aseg->width ; x++)
        for (y = 0 ; y < mri_aseg->height ; y++)
          for (z = 0 ; z < mri_aseg->depth ; z++)
          {
            label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
            if (label != n)
              continue ;
            TransformSample(transform,x, y, z, &xa, &ya, &za);

            *MATRIX_RELT(mX, 1, i) = x ;
            *MATRIX_RELT(mX, 2, i) = y ;
            *MATRIX_RELT(mX, 3, i) = z ;
            *MATRIX_RELT(mX, 4, i) = 1.0 ;

            *MATRIX_RELT(mY, 1, i) = xa ;
            *MATRIX_RELT(mY, 2, i) = ya ;
            *MATRIX_RELT(mY, 3, i) = za ;
            *MATRIX_RELT(mY, 4, i) = 1.0 ;
            i++ ;
          }

      mXinv = MatrixRightPseudoInverse(mX, NULL) ;
      m = MatrixMultiply(mY, mXinv, NULL) ;
      MatrixCopy(m, lta->xforms[0].m_L) ;
      copyVolGeom(&gcam->image, &lta->xforms[0].src) ;
      copyVolGeom(&gcam->atlas, &lta->xforms[0].dst) ;
      sprintf(fname, "%s/%s.lta", out_dir, cma_label_to_name(n)) ;
      sprintf(fname, "%s/label%d.lta", out_dir, n) ;
      LTAwrite(lta, fname) ;
      MatrixFree(&mX) ; MatrixFree(&mY) ; MatrixFree(&m) ; MatrixFree(&mXinv) ;
    }
  }
  MRIfree(&mri_aseg) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "structure transform calculation took %d minutes"
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
  switch (toupper(*option)) {
  case 'N':
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
  printf("usage: %s [options] <.m3z file> <output directory>\n",Progname) ;
  exit(code) ;
}





