/**
 * @brief change a volume to 8 bits/voxel
 *
 * uses the Tal xform to find a ball of voxels that are mostly brain.
 * the top of the intensity histogram in this ball will then be white matter,
 * which allows us to center it at the desired value approximately (110)
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
#include "mrinorm.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRIconvertToUchar(MRI *mri_in, LTA *tal_xform, MRI *mri_out) ;

const char *Progname ;
static void usage_exit(int code) ;

double FIRST_PERCENTILE = 0.01;
double WM_PERCENTILE = 0.90;
double MAX_R = 50;
char *hcumfile = NULL;
char *radfile = NULL;


int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI       *mri_in, *mri_out ;
  LTA *tal_xform ;

  nargs = handleVersionOption(argc, argv, "mri_make_uchar");
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

  if (argc < 3)
    usage_exit(1) ;

  printf("FIRST_PERCENTILE %lf\n",FIRST_PERCENTILE);
  printf("WM_PERCENTILE    %lf\n",WM_PERCENTILE);
  printf("MAX_R %lf\n",MAX_R);


  mri_in = MRIread(argv[1]) ;
  if (mri_in == NULL)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read input volume %s", Progname, argv[1]) ;
  tal_xform = LTAread(argv[2]) ;
  if (tal_xform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s", Progname, argv[2]) ;

  mri_out = MRIconvertToUchar(mri_in, tal_xform, NULL) ;
  MRIwrite(mri_out, argv[3]) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr,
          "type change took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  switch (toupper(*option))
  {
  case 'N':
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  case 'F':
    sscanf(argv[2],"%lf",&FIRST_PERCENTILE);
    nargs = 1;
    break ;
  case 'W':
    sscanf(argv[2],"%lf",&WM_PERCENTILE);
    nargs = 1;
    break ;
  case 'R':
    sscanf(argv[2],"%lf",&MAX_R);
    nargs = 1;
    break ;
  case 'H':
    hcumfile = argv[2];
    nargs = 1;
    break ;
  case 'V':
    radfile = argv[2];
    nargs = 1;
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
  printf(
    "Uses the Tal xform to find a ball of voxels that are mostly brain.\n"
    "The top of the intensity histogram in this ball will then be white\n"
    "matter, which allows us to center it at the desired value,\n"
    "approximately (110).\n\n"
    "usage: %s [options] <input volume> <talairach xform> <output volume>\n",
    Progname) ;

  printf("-f FIRST_PERCENTILE (default %g)\n",FIRST_PERCENTILE);
  printf("-w WM_PERCENTILE (default %g)\n",WM_PERCENTILE);
  printf("-r MAX_R (default %g)\n",MAX_R);
  printf("-h cumulative histo file\n");
  printf("-v vradvol : volume with everthing outside of MAX_R set to 0\n");


  exit(code) ;
}

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
        V3_X(v_X) = xt ;
        V3_Y(v_X) = yt ;
        V3_Z(v_X) = zt ;
        LTAtransformPoint(tal_xform, v_X, v_Y) ;
        xt = V3_X(v_Y) ;
        yt = V3_Y(v_Y) ;
        zt = V3_Z(v_Y) ;
        r = sqrt(xt*xt + yt*yt + zt*zt) ;
        if (r < MAX_R)
          MRIsetVoxVal(mri_out, x, y, z, 0, MRIgetVoxVal(mri_in, x, y, z, 0)) ;
      }
    }
  }
  MatrixFree(&v_X) ;
  MatrixFree(&v_Y) ;
  if(radfile) MRIwrite(mri_out,radfile);

  // only in a radius around the center of the brain
  h = MRIhistogram(mri_out, 100) ;

  HISTOclearZeroBin(h) ;
  HISTOmakePDF(h, h) ;
  hcum = HISTOmakeCDF(h, NULL) ;

  //HISTOwriteTxt(h, "makeuchar.histo.txt") ;
  if(hcumfile) HISTOwriteTxt(hcum, hcumfile);

  bin_size = h->bins[2] - h->bins[1] ;
  i1 = HISTOfindBinWithCount(hcum, FIRST_PERCENTILE) ;
  i2 = HISTOfindBinWithCount(hcum, WM_PERCENTILE) ;
  printf("i1 = %d, i2 = %d\n",i1,i2);

  x1 = h->bins[i1] ;
  x2 = h->bins[i2] ;
  y1 = FIRST_PERCENTILE*255 ;

  /*
    since the ball around tal (0,0,0) will contain almost only brain we
    want to map this intensity range into somewhere around the middle of
    the uchar range so that wm winds up around 110
  */
  y2 = DEFAULT_DESIRED_WHITE_MATTER_VALUE ;
  m = (y2 - y1) / (x2 - x1) ;
  b = y2 - m * x2 ;

  long nzero=0, nsat=0;
  for (x = 0 ; x < mri_in->width ; x++)  {
    for (y = 0 ; y < mri_in->height ; y++)    {
      for (z = 0 ; z < mri_in->depth ; z++)      {
        val = MRIgetVoxVal(mri_in, x, y, z, 0) ;
        val = val * m + b ;
        if (val < 0){
          val = 0 ;
	  nzero ++;
	}
        else if (val > 255){
          val = 255 ;
	  nsat++;
	}
        MRIsetVoxVal(mri_out, x, y, z, 0, val) ;
      }
    }
  }

  printf("#mri_make_uchar# mapping %2.0f %2.0f to %2.0f %2.0f  :  b %g m %g : thresh %g maxsat %g : nzero %ld nsat %ld\n", 
	 x1, x2, y1, y2,b,m, -b/m, (255-b)/m, nzero, nsat) ;


  return(mri_out) ;
}

