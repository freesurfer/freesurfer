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
#include "transform.h"
#include "cma.h"
#include "mri.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MRImarkTemporalWM(MRI *mri_seg, MRI *mri_dst) ;

const char *Progname ;
static void usage_exit(int code) ;

static const char *seg_dir = "seg" ;

static char subjects_dir[STRLEN] ;

int
main(int argc, char *argv[]) {
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI          *mri_seg, *mri_dst ;

  nargs = handleVersionOption(argc, argv, "mri_mark_temporal_lobe");
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

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
    if (argc < 3)
      usage_exit(1) ;
  }


  out_fname = argv[argc-1] ;

  subject_name = argv[1] ;
  sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, seg_dir) ;
  mri_seg = MRIread(fname) ;
  if (!mri_seg)
    ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
              Progname, fname) ;
  mri_dst = MRImarkTemporalWM(mri_seg, NULL) ;

  sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, out_fname) ;
  printf("writing labeled temporal lobe to %s...\n", fname) ;
  MRIwrite(mri_dst, fname) ;
  MRIfree(&mri_dst) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("temporal lobe marking took %d minutes"
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
  if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "PARC_DIR") || !stricmp(option, "SEG_DIR")) {
    seg_dir = argv[2] ;
    nargs = 1 ;
    printf("reading segmentation from subject's mri/%s directory\n",
           seg_dir) ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else switch (toupper(*option)) {
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
  printf("usage: %s [options] <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  printf(
    "\t-spacing  - spacing of classifiers in canonical space\n");
  printf("\t-gradient - use intensity gradient as input to classifier.\n") ;
  exit(code) ;
}

MRI *
MRImarkTemporalWM(MRI *mri_seg, MRI *mri_dst) {
  int   x, y, z, yi, found_other, found_gray, nchanged, label, i,
  width, height, depth, left ;
  MRI   *mri_tmp ;


  width = mri_seg->width ;
  height = mri_seg->height ;
  depth = mri_seg->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_seg, NULL) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z== Gz)
          DiagBreak() ;
        label = MRIvox(mri_seg, x, y, z) ;
        if (!IS_WM(label))
          continue ;
        left =
          (label == Left_Cerebral_White_Matter) ||
          (label == Left_Temporal_Cerebral_White_Matter) ;

        for (found_other = found_gray = 0, i = 1 ; i <= 3 ; i++) {
          yi = mri_seg->yi[y+i] ;
          label = MRIvox(mri_seg, x, yi, z) ;
          if (IS_CORTEX(label))
            found_gray = 1 ;   /* gray matter inferior */

          yi = mri_seg->yi[y-i] ;
          label = MRIvox(mri_seg, x, yi, z) ;
          if (IS_HIPPO(label) || IS_AMYGDALA(label))
            found_other = 1 ;
        }
        if (found_other && found_gray) {
          MRIvox(mri_dst, x, y, z) = left ?
                                     Left_Cerebral_White_Matter:
                                     Right_Cerebral_White_Matter ;
          MRIvox(mri_seg, x, y, z) = left ?
                                     Left_Temporal_Cerebral_White_Matter:
                                     Right_Temporal_Cerebral_White_Matter ;
        }
      }
    }
  }

  do {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == Gx && y == Gy && z== Gz)
            DiagBreak() ;
          if (MRIvox(mri_dst, x, y, z) > 0)
            continue ;
          label = MRIvox(mri_seg, x, y, z) ;
          if (!IS_WM(label))
            continue ;
          left =
            (label == Left_Cerebral_White_Matter) ||
            (label == Left_Temporal_Cerebral_White_Matter) ;

          for (found_other = found_gray = 0, i = 1 ; i <= 3 ; i++) {
            yi = mri_seg->yi[y+i] ;
            label = MRIvox(mri_seg, x, yi, z) ;
            if (IS_CORTEX(label))
              found_gray = 1 ;   /* gray matter inferior */

            yi = mri_seg->yi[y-i] ;  /* lateral ventricle superior */
            label = MRIvox(mri_seg, x, yi, z) ;
            if (IS_LAT_VENT(label))
              found_other = 1 ;
          }
          if (found_other && found_gray && MRIneighborsOn(mri_dst, x, y,z,1)>0) {
            nchanged++ ;
            MRIvox(mri_dst, x, y, z) = left ?
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter ;
            MRIvox(mri_seg, x, y, z) = left ?
                                       Left_Temporal_Cerebral_White_Matter :
                                       Right_Temporal_Cerebral_White_Matter ;
          }
        }
      }
    }
  } while (nchanged > 0) ;

  for (i = 0 ; i < 3 ; i++) {
    MRIcopy(mri_dst, mri_tmp) ;
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == Gx && y == Gy && z== Gz)
            DiagBreak() ;
          label = MRIvox(mri_seg, x, y, z) ;
          if (!IS_WM(label) || (MRIvox(mri_dst, x, y, z) > 0))
            continue ;
          left =
            (label == Left_Cerebral_White_Matter) ||
            (label == Left_Temporal_Cerebral_White_Matter) ;

          if (MRIneighborsOn(mri_dst, x, y,z, 1)>= 2) {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = left ?
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter ;
            MRIvox(mri_seg, x, y, z) = left ?
                                       Left_Temporal_Cerebral_White_Matter :
                                       Right_Temporal_Cerebral_White_Matter ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
  }
  /* now dilate twice to get stray voxels and a bit of main temporal lobe */
  for (i = 0 ; i < 2 ; i++) {
    MRIcopy(mri_dst, mri_tmp) ;
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == Gx && y == Gy && z== Gz)
            DiagBreak() ;
          label = MRIvox(mri_seg, x, y, z) ;
          if (!IS_WM(label) || (MRIvox(mri_dst, x, y, z) > 0))
            continue ;
          left =
            (label == Left_Cerebral_White_Matter) ||
            (label == Left_Temporal_Cerebral_White_Matter) ;

          if (MRIneighborsOn(mri_dst, x, y,z, 1)>= 1) {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = left ?
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter ;
            MRIvox(mri_seg, x, y, z) = left ?
                                       Left_Temporal_Cerebral_White_Matter :
                                       Right_Temporal_Cerebral_White_Matter ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
  }

  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
