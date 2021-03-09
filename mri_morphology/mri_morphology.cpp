/**
 * @brief applies morphological operations to a volume.
 *
 * Program for applying morphological operations open, close, dilate, and erode either to an
 * entire volume or to only a label (specified with -l <label>) within it.
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
#include "timer.h"
#include "version.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


const char *Progname ;

#define ERODE          1
#define DILATE         2
#define CLOSE          3
#define OPEN           4
#define DILATE_LABEL   5
#define MODE_FILTER    6
#define ERODE_THRESH   7
#define DILATE_THRESH  8
#define ERODE_BOTTOM   9
#define FILL_HOLES     10


static void usage_exit(int code) ;

static int label = -1 ;

MRI *MRIerodeBottom(MRI *mri_src, int label, MRI *mri_dst) ;

static MRI *mri_mask = NULL ;
int
main(int argc, char *argv[]) {
  char   *out_fname, **av ;
  int    ac, nargs, niter, operation, i ;
  MRI    *mri_src, *mri_dst = nullptr, *mri_saved_src = nullptr ;
  int          msec, minutes, seconds ;
  Timer start ;

  nargs = handleVersionOption(argc, argv, "mri_morphology");
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

  if (argc < 5)
    usage_exit(1) ;

  out_fname = argv[4] ;
  mri_src = MRIread(argv[1]) ;
  if (mri_src == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s\n",
              Progname,argv[1]);

  if (!stricmp(argv[2], "dilate"))
    operation = DILATE ;
  else  if (!stricmp(argv[2], "open"))
    operation = OPEN ;
  else  if (!stricmp(argv[2], "close"))
    operation = CLOSE ;
  else  if (!stricmp(argv[2], "erode"))
    operation = ERODE ;
  else  if (!stricmp(argv[2], "erode_thresh"))
    operation = ERODE_THRESH ;
  else  if (!stricmp(argv[2], "dilate_thresh"))
    operation = DILATE_THRESH ;
  else  if (!stricmp(argv[2], "mode"))
    operation = OPEN ;
  else  if (!stricmp(argv[2], "erode_bottom"))
    operation = ERODE_BOTTOM ;
  else  if (!stricmp(argv[2], "fill_holes"))
    operation = FILL_HOLES ;
  else {
    operation = 0 ;
    ErrorExit(ERROR_UNSUPPORTED, "morphological operation '%s'  is not supported", argv[2]) ;
  }

  if (label > 0)  // erase everything but label in src, and save orig vol
  {
    MRI *mri_tmp ;

    mri_saved_src = MRIcopy(mri_src, NULL) ;
    mri_tmp = MRIclone(mri_src, NULL);
    MRIcopyLabel(mri_src, mri_tmp,label) ;
    MRIfree(&mri_src) ;
    mri_src = mri_tmp ;
  }
    
  niter = atoi(argv[3]) ;
  switch (operation) {
  case MODE_FILTER: {
    MRI *mri_tmp ;
    if (mri_src->type != MRI_UCHAR) {
      mri_tmp = MRIchangeType(mri_src, MRI_UCHAR, 0, 1, 1) ;
      MRIfree(&mri_src) ;
      mri_src = mri_tmp ;
    }
    printf("applying mode filter %d times\n", niter) ;
    mri_dst = MRImodeFilter(mri_src, NULL, niter) ;
    break ;
  }
  case FILL_HOLES:
    {
      int x, y, z, nfilled = 0 ;
      mri_dst = MRIbinarize(mri_src, NULL, 1, 0, 1) ;
      for (x = 0 ; x < mri_src->width ; x++)
        for (y = 0 ; y < mri_src->height ; y++)
          for (z = 0 ; z < mri_src->depth ; z++)
          {
            if (x == Gx && y == Gy && z == Gz)
              DiagBreak() ;
            if ((MRIgetVoxVal(mri_src, x, y, z, 0) == 0) &&
                MRIneighborsOn3x3(mri_src, x, y, z, 0) >= niter)
	    {
	      nfilled++ ;
              MRIsetVoxVal(mri_dst, x, y, z, 0, 1) ;
	    }
          }
      printf("%d holes filled with more than %d neighbors on\n",nfilled,niter);
      break ;
    }
  case ERODE_BOTTOM:
    if (label < 0)
      ErrorExit(ERROR_BADPARM, "%s: must specify label with -l <label>", Progname) ;
    mri_dst = MRIerodeBottom(mri_src, label, NULL) ;
    break ;
  case DILATE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case CLOSE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case OPEN:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIdilate(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case ERODE:
    mri_dst = NULL ;
    for (i = 0 ; i < niter ; i++) {
      mri_dst = MRIerode(mri_src, mri_dst) ;
      MRIcopy(mri_dst, mri_src) ;
    }
    break ;
  case ERODE_THRESH:
    {
      double thresh = atof(argv[5]) ;
      char   *intensity_fname = argv[6] ;
      MRI    *mri_intensity ;

      if (argc < 7)
        ErrorExit(ERROR_BADPARM, "%s: for erode_thresh must specify <thresh> <intensity vol> on cmdline", Progname) ;
      mri_intensity = MRIread(intensity_fname) ;
      if (mri_intensity == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s",
                  Progname, intensity_fname) ;
      mri_dst = NULL ;
      printf("eroding %s with threshold %2.1f from volume %s\n",
             argv[1], thresh, intensity_fname) ;
      for (i = 0 ; i < niter ; i++) {
        mri_dst = MRIerodeThresh(mri_src, mri_intensity, thresh, mri_dst) ;
        MRIcopy(mri_dst, mri_src) ;
      }
    }
    break ;
  case DILATE_THRESH:
    {
      double thresh = atof(argv[5]) ;
      char   *intensity_fname = argv[6] ;
      MRI    *mri_intensity ;

      if (argc < 7)
        ErrorExit(ERROR_BADPARM, "%s: for erode_thresh must specify <thresh> <intensity vol> on cmdline", Progname) ;
      mri_intensity = MRIread(intensity_fname) ;
      if (mri_intensity == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s",
                  Progname, intensity_fname) ;
      mri_dst = NULL ;
      printf("dilating %s with threshold %2.1f from volume %s\n",
             argv[1], thresh, intensity_fname) ;
      for (i = 0 ; i < niter ; i++) {
        mri_dst = MRIdilate6Thresh(mri_src, mri_intensity, thresh, mri_dst) ;
        MRIcopy(mri_dst, mri_src) ;
      }
    }
    break ;
  default:
    break ;
  }
  if ((label > 0) && (operation != ERODE_BOTTOM))
  {
    MRIreplaceValues(mri_saved_src, mri_saved_src, label, 0) ;
    MRIcopyLabel(mri_dst, mri_saved_src, label) ;
    MRIcopy(mri_saved_src, mri_dst) ;
    MRIfree(&mri_saved_src) ;
  }

  if (mri_mask)
    MRImask(mri_dst, mri_mask, mri_dst, 0, 0) ;
  fprintf(stderr, "writing to %s...\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  MRIfree(&mri_dst) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphological processing took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "mask")) 
  {
    mri_mask = MRIread(argv[2]) ;
    if (mri_mask == NULL)
      exit(Gerror) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else switch (toupper(*option)) {
  case 'L':
    label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("only applying operations to label %d\n", label) ;
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
  printf("usage: %s [options] <volume> <operation> <# iter> <out volume>\n", Progname) ;
  printf("\twhere <operation> can be [open,close,dilate,erode,mode,fill_holes,erode_bottom,dilate_thresh,erode_thresh]\n");
  printf("\tvalid options are:\n") ;
  printf("\t-l <label>  only apply operations to <label> instead of all nonzero voxels\n") ;
  exit(code) ;
}

MRI *
MRIerodeBottom(MRI *mri_src, int label, MRI *mri_dst)
{
  int  neroded, x, y, z, olabel ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;


  do
  {
    neroded = 0 ;
    for (y = mri_dst->height-1 ; y  >= 0 ; y--)
      for (x = 0 ; x < mri_dst->width ; x++)
        for (z = 0 ; z < mri_dst->depth ; z++)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          olabel = nint(MRIgetVoxVal(mri_dst, x, y, z, 0)) ;
          if (olabel ==  label)
          {
            if (MRIgetVoxVal(mri_dst, x, y-1, z, 0) == label)
            {
              neroded++ ;
              MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
            }
          }
          else if (olabel > 0)
            MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
        }
  } while (neroded > 0) ;
  return(mri_dst) ;
}

