/**
 * @file  mri_make_uchar.c
 * @brief change a volume to 8 bits/voxel
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2010/12/02 18:06:54 $
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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRIconvertToUchar(MRI *mri_in, LTA *tal_xform, MRI *mri_out) ;

char *Progname ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MRI       *mri_in, *mri_out ;
  LTA *tal_xform ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_make_uchar.c,v 1.1 2010/12/02 18:06:54 fischl Exp $", "$Name:  $");
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


  mri_in = MRIread(argv[1]) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname, argv[1]) ;
  tal_xform = LTAread(argv[2]) ;
  if (tal_xform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s", Progname, argv[2]) ;

  mri_out = MRIconvertToUchar(mri_in, tal_xform, NULL) ;
  MRIwrite(mri_out, argv[3]) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "type change took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  printf("usage: %s [options] <input volume> <talairach xform> <output volume>",Progname) ;
  exit(code) ;
}




#define FIRST_PERCENTILE    0.01
#define SECOND_PERCENTILE   0.99

MRI *
MRIconvertToUchar(MRI *mri_in, LTA *tal_xform, MRI *mri_out)
{
  int       x, y, z, i1, i2 ;
  double    xt, yt, zt, r, bin_size, x1, y1, x2, y2, m, b, val ;
  VECTOR    *v_X, *v_Y ;
  HISTOGRAM *h, *hcum ;

  v_X = MatrixAlloc(4, 1, MATRIX_REAL) ;
  v_Y = MatrixAlloc(4, 1, MATRIX_REAL) ;

  if (mri_out == NULL)
    mri_out = MRIcloneDifferentType(mri_in, MRI_UCHAR) ;

#define MAX_R 50

  VECTOR_ELT(v_X, 4) = 1.0 ;
  VECTOR_ELT(v_Y, 4) = 1.0 ;
  for (x = 0 ; x < mri_in->width ; x++)
  {
    V3_X(v_X) = x ;
    for (y = 0 ; y < mri_in->height ; y++)
    {
      V3_Y(v_X) = y ;
      for (z = 0 ; z < mri_in->depth ; z++)
      {
        V3_Z(v_X) = z ;
        MRIvoxelToTalairach(mri_in, x, y, z, &xt, &yt, &zt) ;
        r = sqrt(xt*xt + yt*yt + zt*zt) ;

        MRIvoxelToWorld(mri_in, x, y, z, &xt, &yt, &zt) ;
        V3_X(v_X) = xt ; V3_Y(v_X) = yt ; V3_Z(v_X) = zt ; 
        LTAtransformPoint(tal_xform, v_X, v_Y) ;
        xt = V3_X(v_Y) ; yt = V3_Y(v_Y) ; zt = V3_Z(v_Y) ;
        r = sqrt(xt*xt + yt*yt + zt*zt) ;
        if (r < MAX_R)
          MRIsetVoxVal(mri_out, x, y, z, 0, MRIgetVoxVal(mri_in, x, y, z, 0)) ;
      }
    }
  }
  MatrixFree(&v_X) ; MatrixFree(&v_Y) ;


  h = MRIhistogram(mri_out, 100) ;  // only in a radius around the center of the brain

  HISTOclearZeroBin(h) ;
  HISTOmakePDF(h, h) ;
  hcum = HISTOmakeCDF(h, NULL) ;

  bin_size = h->bins[2] - h->bins[1] ;
  i1 = HISTOfindBinWithCount(hcum, FIRST_PERCENTILE) ;
  i2 = HISTOfindBinWithCount(hcum, SECOND_PERCENTILE) ;

  x1 = h->bins[i1] ; x2 = h->bins[i2] ;
  y1 = FIRST_PERCENTILE*255 ;

  /*
    since the ball around tal (0,0,0) will contain almost only brain we
    want to map this intensity range into somewhere around the middle of
    the uchar range so that wm winds up around 128
  */
  y2 = (0.5*(FIRST_PERCENTILE+SECOND_PERCENTILE))*255 ;
  m = (y2 - y1) / (x2 - x1) ;
  b = y2 - m * x2 ;
  printf("mapping (%2.0f, %2.0f) to (%2.0f, %2.0f)\n", x1, x2, y1, y2) ;

  for (x = 0 ; x < mri_in->width ; x++)
  {
    for (y = 0 ; y < mri_in->height ; y++)
    {
      for (z = 0 ; z < mri_in->depth ; z++)
      {
        val = MRIgetVoxVal(mri_in, x, y, z, 0) ;
        val = val * m + b ;
        if (val < 0)
          val = 0 ;
        else if (val > 255)
          val = 255 ;
        MRIsetVoxVal(mri_out, x, y, z, 0, val) ;
      }
    }
  }

  return(mri_out) ;
}

