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
#include "timer.h"
#include "proto.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "cma.h"
#include "version.h"

static int mle_label(MRI *mri_T1, MRI *mri_out_labeled, int x, int y,
                     int z, int wsize, int l1, int l2) ;
static int change_label(MRI *mri_T1, MRI *mri_out_labeled, int x, int y,
                        int z, int wsize, int left) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;

static MRI *edit_border_voxels(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) ;
static MRI *edit_ventricular_unknowns(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) ;
static int sagittalNeighbors(MRI *mri, int x, int y,int z,int whalf,int label);
static int neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label);
static int distance_to_label(MRI *mri_labeled, int label, int x,
                             int y, int z, int dx, int dy,
                             int dz, int max_dist) ;
static MRI *edit_hippocampus(MRI *mri_in_labeled, MRI *mri_T1,
                             MRI *mri_out_labeled) ;

static MRI *edit_amygdala(MRI *mri_in_labeled, MRI *mri_T1,
                          MRI *mri_out_labeled) ;
static MRI *edit_caudate(MRI *mri_in_labeled, MRI *mri_T1,
                         MRI *mri_out_labeled) ;
static MRI *edit_lateral_ventricles(MRI *mri_in_labeled, MRI *mri_T1,
                                    MRI *mri_out_labeled) ;
static MRI *edit_cortical_gray_matter(MRI *mri_in_labeled, MRI *mri_T1,
                                      MRI *mri_out_labeled) ;


const char *Progname ;

static int unknown_only = 0 ;
static int border_only = 0 ;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  MRI    *mri_in_labeled, *mri_T1, *mri_out_labeled = NULL ;
  char   *in_fname, *T1_fname, *out_fname ;
  int          msec, minutes, seconds ;
  Timer start ;

  nargs = handleVersionOption(argc, argv, "mri_edit_segmentation");
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

  if (argc < 4)
    usage_exit() ;

  in_fname = argv[1] ;
  T1_fname = argv[2] ;
  out_fname = argv[3] ;

  printf("reading from %s...\n", in_fname) ;
  mri_in_labeled = MRIread(in_fname) ;
  if (!mri_in_labeled)
    ErrorExit(ERROR_NO_FILE, "%s: could not open labeled file %s",
              Progname, in_fname) ;

  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NO_FILE, "%s: could not open T1 file %s",
              Progname, in_fname) ;

  start.reset() ;

  if (unknown_only == 0 && border_only == 0) {
    mri_out_labeled = edit_hippocampus(mri_in_labeled, mri_T1, NULL);
    edit_amygdala(mri_out_labeled, mri_T1, mri_out_labeled);
    edit_caudate(mri_out_labeled, mri_T1, mri_out_labeled);
    edit_lateral_ventricles(mri_out_labeled, mri_T1, mri_out_labeled);  /* must be after hippo */
    edit_cortical_gray_matter(mri_out_labeled, mri_T1, mri_out_labeled);
  } else
    mri_out_labeled = MRIcopy(mri_in_labeled, NULL) ;

  if (border_only == 0 || unknown_only != 0)
    edit_ventricular_unknowns(mri_out_labeled, mri_T1, mri_out_labeled) ;
  if (border_only != 0 || unknown_only == 0)
    edit_border_voxels(mri_out_labeled, mri_T1, mri_out_labeled) ;

  printf("writing output volume to %s...\n", out_fname) ;
  MRIwrite(mri_out_labeled, out_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("segmentation adjustment took %d minutes and %d seconds.\n",
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
  if (!stricmp(option, "no1d")) {
    printf("disabling 1d normalization...\n") ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "UNKNOWN")) {
    unknown_only = 1 ;
    printf("only correcting unknown voxels between ventricle and wm or near hypointensities\n") ;
  } else if (!stricmp(option, "BORDER")) {
    border_only = 1 ;
    printf("only correcting border voxels\n") ;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      printf("usage: %s [input directory] [output directory]\n", Progname) ;
      exit(1) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
static void
usage_exit(void) {
  printf("usage: %s [options] "
         "<input segmentation> <T1 volume> <output segmentation>\n",
         Progname) ;
  exit(0) ;
}

static int
distance_to_label(MRI *mri_labeled, int label, int x, int y, int z, int dx,
                  int dy, int dz, int max_dist) {
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++) {
    xi = x + d * dx ;
    yi = y + d * dy ;
    zi = z + d * dz ;
    xi = mri_labeled->xi[xi] ;
    yi = mri_labeled->yi[yi] ;
    zi = mri_labeled->zi[zi];
    if (MRIvox(mri_labeled, xi, yi, zi) == label)
      break ;
  }

  return(d) ;
}
static MRI *
edit_hippocampus(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int   width, height, depth, x, y, z, nchanged, dleft, label,
  dright, dpos, dant, dup, ddown, i, left, dgray, dhippo, dwhite, olabel ;
  MRI   *mri_tmp ;

  nchanged = 0 ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;

  /* change all gm within 2 mm of ventricle to wm */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri_out_labeled, x, y, z) ;
        if (IS_GM(label) == 0)
          continue ;
        left = label == Left_Cerebral_Cortex ;
        olabel = left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle;
        if (MRIneighborsInWindow(mri_out_labeled, x, y, z, 5, olabel) >= 1) {
          nchanged++ ;
          MRIvox(mri_out_labeled, x, y, z) = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
        }
      }
    }
  }

  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  /* change gray to hippocampus based on wm */
  for (i = 0 ; i < 3 ; i++) {
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;
          if (x == 160 && y == 127 && z == 118)
            DiagBreak() ;

          left = 0 ;
          switch (label) {
          case Left_Cerebral_White_Matter:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Cerebral_White_Matter:
            dgray = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_Cortex :
                                      Right_Cerebral_Cortex,
                                      x,y,z,0,1,0,3);

            dwhite = distance_to_label(mri_out_labeled,
                                       left ? Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,0,-1,0,1);

            dhippo = distance_to_label(mri_out_labeled,
                                       left ? Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,0,-1,0,3);
            /* change bright hippocampus that borders white matter to white matter */
            if (dgray <= 2 && dwhite <= 1 && dhippo <= 3) {
              mle_label(mri_T1, mri_tmp, x, y, z, 9,
                        left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex,
                        left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter) ;
            }
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Left_Hippocampus:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Hippocampus:
            dgray = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_Cortex :
                                      Right_Cerebral_Cortex,
                                      x,y,z,0,1,0,3);

            dwhite = distance_to_label(mri_out_labeled,
                                       left ? Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,0,1,0,1);

            dhippo = distance_to_label(mri_out_labeled,
                                       left ? Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,0,-1,0,1);
            /* change bright hippocampus that borders white matter to white matter */
            if (dgray <= 3 && dwhite <= 1 && dhippo <= 1) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
            }
            break ;
          case Unknown: {
            int dhippo, dl, dr, dgray, dwhite ;

            if (x == 160 && y == 129 && z == 121)
              DiagBreak() ;

            /*
             if the current label is unknown, and there is white matter above
             and hippocampus above and gray matter below, change to gray matter.
            */
            dl = distance_to_label(mri_out_labeled,
                                   Left_Hippocampus,
                                   x,y,z,0,-1,0,4);
            dr = distance_to_label(mri_out_labeled,
                                   Right_Hippocampus,
                                   x,y,z,0,-1,0,4);
            if (dl > 4 && dr > 4) /* could be inferior to inf-lat-ven also */
            {
              dl = distance_to_label(mri_out_labeled,
                                     Left_Inf_Lat_Vent,
                                     x,y,z,0,-1,0,4);
              dr = distance_to_label(mri_out_labeled,
                                     Right_Inf_Lat_Vent,
                                     x,y,z,0,-1,0,4);
            }

            if (dl < dr) {
              left = 1 ;
              dhippo = dl ;
              dgray = distance_to_label(mri_out_labeled,
                                        Left_Cerebral_Cortex,
                                        x,y,z,0,1,0,2);
              dwhite = distance_to_label(mri_out_labeled,
                                         Left_Cerebral_White_Matter,
                                         x,y,z,0,-1,0,2);
            } else {
              left = 0 ;
              dhippo = dr ;
              dwhite = distance_to_label(mri_out_labeled,
                                         Right_Cerebral_White_Matter,
                                         x,y,z,0,-1,0,2);
              dgray = distance_to_label(mri_out_labeled,
                                        Right_Cerebral_Cortex,
                                        x,y,z,0,1,0,2);
            }
            if (dhippo <= 4 && dwhite <= 2 && dgray <= 2) {
              MRIvox(mri_tmp, x, y, z) =
                left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
              nchanged++ ;
              continue ;
            }

            break ;
          }
          case Left_Cerebral_Cortex:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Cerebral_Cortex:
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 1) {
              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }

            /*
              if the current label is gray, and there is white matter below
              and hippocampus above, change to hippocampus.
            */
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            /*
              if the current label is gray, and there is white matter above
              and hippocampus below, change to hippocampus.
            */
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Cerebral_White_Matter :
                                    Right_Cerebral_White_Matter,
                                    x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,
                                    x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            dleft = distance_to_label(mri_out_labeled, left ?
                                      Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,-1,0,0,3);
            dright = distance_to_label(mri_out_labeled, left ?
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,1,0,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dleft <= 2 && dright <= 2 && dup <= 1) {
              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
            dright = distance_to_label(mri_out_labeled, left ?
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,1,0,0,3);
            dleft = distance_to_label(mri_out_labeled, left ?
                                      Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,-1,0,0,3);
            if (dleft <= 1 && dright <= 1) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            dleft = distance_to_label(mri_out_labeled, left ?
                                      Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,-1,0,0,3);
            dright = distance_to_label(mri_out_labeled, left ?
                                       Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,1,0,0,3);
            if (dleft <= 1 && dright <= 1) {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            dright = distance_to_label(mri_out_labeled, left ?
                                       Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,1,0,0,3);
            dleft = distance_to_label(mri_out_labeled, left ?
                                      Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,-1,0,0,3);
            if (dleft <= 1 && dright <= 1) {
              label = left ?
                      Left_Hippocampus : Right_Hippocampus;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }

            dpos = distance_to_label(mri_out_labeled,
                                     Right_Cerebral_White_Matter,x,y,z,0,0,-1,3);
            dant = distance_to_label(mri_out_labeled,
                                     Right_Cerebral_White_Matter,x,y,z,0,0,1,3);


            /* if the current label is gray and the is white matter directly above,
              and hippocampus within 3 mm, then change label to MLE of either gray or
              white
            */

            if (x == 156 && y == 128 && z == 115)
              DiagBreak() ;
            dgray = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_Cortex :
                                      Right_Cerebral_Cortex,
                                      x,y,z,0,1,0,1);

            dwhite = distance_to_label(mri_out_labeled,
                                       left ? Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,0,-1,0,1);

            dhippo = distance_to_label(mri_out_labeled,
                                       left ? Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,0,-1,0,3);
            /* change bright hippocampus that borders white matter to white matter */
            if (dgray <= 1 && dwhite <= 1 && dhippo <= 3) {
              mle_label(mri_T1, mri_tmp, x, y, z, 9,
                        left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex,
                        left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter) ;
            }
            break ;
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
  }

  /* make gray matter that is superior to hippocampus into hippocampus */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        int yi ;

        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri_out_labeled, x, y, z) ;
        if (!IS_HIPPO(label))
          continue ;
        left = label == Left_Hippocampus ;

        /* search for first non-hippocampal voxel */
        for (yi = y-1 ; yi >= 0 ; yi--) {
          label = MRIvox(mri_out_labeled, x, yi, z) ;
          if (!IS_HIPPO(label))
            break ;
        }
        i = 0 ;
        while (IS_GM(label) && yi >= 0) {
          nchanged++ ;
          MRIvox(mri_out_labeled, x, yi, z) = left ? Left_Hippocampus : Right_Hippocampus;
          yi-- ;
          label = MRIvox(mri_out_labeled, x, yi, z) ;
          if (++i >= 4)
            break ;  /* don't let it go too far */
        }

      }
    }
  }


  /* go through and change wm labels that have hippo both above and below them to hippo */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri_tmp, x, y, z) ;
        if (x == 160 && y == 127 && z == 118)
          DiagBreak() ;

        if (!IS_WM(label))
          continue ;
        left = label == Left_Cerebral_White_Matter ;
        if (IS_HIPPO(MRIvox(mri_out_labeled, x, mri_out_labeled->yi[y-1], z)) &&
            IS_HIPPO(MRIvox(mri_out_labeled, x, mri_out_labeled->yi[y+1], z))) {
          nchanged++ ;
          MRIvox(mri_out_labeled, x, y, z) = left ? Left_Hippocampus : Right_Hippocampus;
        }
      }
    }
  }

  MRIfree(&mri_tmp) ;
  printf("%d hippocampal voxels changed.\n", nchanged) ;
  return(mri_out_labeled) ;
}
static MRI *
edit_amygdala(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int   width, height, depth, x, y, z, nchanged, label, total_changed,
  dup, ddown, left ;
  MRI   *mri_tmp ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  total_changed = 0 ;
  do {
    nchanged = 0 ;

    /* change gray to wm if near amygdala */
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == 95 && y == 127 && z == 119)
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label) {
          case Left_Cerebral_Cortex:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Cerebral_Cortex:
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Amygdala :
                                    Right_Amygdala,x,y,z,0,-1,0,2);
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            if (dup <= 1 && ddown <= 1) {
              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
    total_changed += nchanged ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d amygdala voxels changed.\n", total_changed) ;
  return(mri_out_labeled) ;
}
static MRI *
edit_caudate(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int   width, height, depth, x, y, z, nchanged, label, total_changed,
  left, niter, change ;
  MRI   *mri_tmp ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  niter = total_changed = 0 ;
  do {
    nchanged = 0 ;

    /* change gray to wm if near amygdala */
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == 95 && y == 127 && z == 119)
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          change = left = 0 ;
          switch (label) {
          case Unknown:
            if (neighborLabel(mri_tmp, x, y, z, 1, Left_Caudate))
              change = left = 1 ;
            if (neighborLabel(mri_tmp, x, y, z, 1, Right_Caudate))
              change = 1 ;

            if (change) {
              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
    total_changed += nchanged ;
  } while ((nchanged > 0) && (niter++ < 3)) ;

  MRIfree(&mri_tmp) ;
  printf("%d unknown voxels bordering caudate changed to wm.\n",total_changed);
  return(mri_out_labeled) ;
}
float
label_mean(MRI *mri_T1, MRI *mri_labeled, int x, int y, int z, int wsize,
           int label) {
  int   xi, yi, zi, xk, yk, zk, whalf, nvox ;
  float mean ;


  whalf = (wsize-1)/2 ;

  for (mean = 0.0, nvox = 0, xk = -whalf ; xk <= whalf ; xk++) {
    xi = mri_T1->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++) {
      yi = mri_T1->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++) {
        zi = mri_T1->zi[z+zk] ;
        if (MRIvox(mri_labeled, xi, yi, zi) != label)
          continue ;
        nvox++ ;
        mean += (float)MRIvox(mri_T1, xi, yi, zi) ;
      }
    }
  }

  if (nvox > 0)
    mean /= (float)nvox ;
  else
    mean = 0.0f ;
  return(mean) ;
}
static int
mle_label(MRI *mri_T1, MRI *mri_labeled, int x, int y, int z, int wsize, int l1, int l2) {
  float  l1_mean, l2_mean, val ;
  int    label ;

  if (x == 95 && y == 127 && z == 119) /* dark wm (68) */
    DiagBreak() ;
  if (x == 94 && y == 126 && z == 119)  /* bright hippo (104) */
    DiagBreak() ;
  val = (float)MRIvox(mri_T1, x, y, z) ;
  l1_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize, l1) ;
  l2_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize, l2) ;
  if (fabs(l1_mean-val) < fabs(l2_mean-val))
    label = l1 ;
  else
    label = l2 ;

  MRIvox(mri_labeled, x, y, z) = label ;
  return(NO_ERROR) ;
}
static int
change_label(MRI *mri_T1,MRI *mri_labeled,int x,int y,int z,int wsize,int left) {
  float  wm_mean, hippo_mean, val ;
  int    label ;

  if (x == 95 && y == 127 && z == 119) /* dark wm (68) */
    DiagBreak() ;
  if (x == 94 && y == 126 && z == 119)  /* bright hippo (104) */
    DiagBreak() ;
  val = (float)MRIvox(mri_T1, x, y, z) ;
  wm_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize,
                       left ? Left_Cerebral_White_Matter :
                       Right_Cerebral_White_Matter) ;
  hippo_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize,
                          left ? Left_Hippocampus :
                          Right_Hippocampus) ;
  if (fabs(wm_mean-val) < fabs(hippo_mean-val))
    label = left ?
            Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
  else
    label = left ?
            Left_Hippocampus : Right_Hippocampus;

  MRIvox(mri_labeled, x, y, z) = label ;
  return(NO_ERROR) ;
}

static int
neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label) {
  int xi, yi, zi, xk, yk, zk ;

  for (zk = -whalf ; zk <= whalf ; zk++) {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++) {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++) {
        if (abs(xk)+abs(yk)+abs(zk) > 1) /* only 6-connected neighbors */
          continue ;
        xi = mri->xi[x+xk] ;
        if (MRIvox(mri, xi, yi, zi) == label)
          return(1) ;
      }
    }
  }
  return(0) ;
}
static int
sagittalNeighbors(MRI *mri, int x, int y, int z, int whalf, int label) {
  int nbrs = 0, yi, zi, yk, zk ;

  for (zk = -whalf ; zk <= whalf ; zk++) {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++) {
      yi = mri->yi[y+yk] ;
      if (MRIvox(mri, x, yi, zi) == label)
        nbrs++ ;
    }
  }
  return(nbrs) ;
}

static MRI *
edit_cortical_gray_matter(MRI *mri_in_labeled,MRI *mri_T1,MRI *mri_out_labeled) {
  int   width, height, depth, x, y, z, nchanged, label, total_changed,
  left, niter, change ;
  MRI   *mri_tmp ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  niter = total_changed = 0 ;
  do {
    nchanged = 0 ;

    /* change gray to wm if near amygdala */
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == 95 && y == 127 && z == 119)
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          change = left = 0 ;
          switch (label) {
          case Left_Cerebral_Cortex:
          case Unknown:
          case Right_Cerebral_Cortex:
            if (neighborLabel(mri_tmp,x,y,z,1,Left_Cerebral_White_Matter) &&
                neighborLabel(mri_tmp,x,y,z,1,Right_Cerebral_White_Matter)) {
              left =
                sagittalNeighbors(mri_tmp,x,y,z,3,Left_Cerebral_White_Matter)>
                sagittalNeighbors(mri_tmp,x,y,z,3,Right_Cerebral_White_Matter);

              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
    total_changed += nchanged ;
  } while ((nchanged > 0) && (niter++ < 3)) ;

  MRIfree(&mri_tmp) ;
  printf("%d midline voxels changed.\n", total_changed) ;
  return(mri_out_labeled) ;
}
static MRI *
edit_lateral_ventricles(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int   width, height, depth, x, y, z, nchanged, label, total_changed,
  left, niter, change, dvent, dwhite, olabel ;
  MRI   *mri_tmp ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ;
  height = mri_T1->height ;
  depth = mri_T1->depth ;

  niter = total_changed = 0 ;
  do {
    nchanged = 0 ;

    /* change gray to wm if near amygdala */
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          change = left = 0 ;
          switch (label) {
          case Left_Inf_Lat_Vent:
            left = 1 ;
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
          case Right_Inf_Lat_Vent:
            dwhite = distance_to_label(mri_out_labeled,
                                       left ? Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,0,1,0,3);
            dvent = distance_to_label(mri_out_labeled,
                                      left ? Left_Inf_Lat_Vent :
                                      Right_Inf_Lat_Vent,
                                      x,y,z,0,-1,0,3);
            if (dvent <= 1 && dwhite <= 1)  /* borders ventricle superior and wm inferior */
            {
              mle_label(mri_T1, mri_tmp, x, y, z, 9,
                        left ? Left_Inf_Lat_Vent : Right_Inf_Lat_Vent,
                        left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter) ;
            }
            break ;
          case Unknown:
            if (neighborLabel(mri_tmp, x, y, z, 1, Left_Lateral_Ventricle) &&
                neighborLabel(mri_tmp, x, y, z,1,Left_Cerebral_White_Matter))
              change = left = 1 ;
            if (neighborLabel(mri_tmp, x, y, z, 1, Right_Lateral_Ventricle) &&
                neighborLabel(mri_tmp, x, y, z,1,Right_Cerebral_White_Matter))
              change = 1 ;

            if (change) {
#if 0
              label = left ?
                      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              /*              MRIvox(mri_tmp, x, y, z) = label ;*/
#else
              mle_label(mri_T1, mri_tmp, x, y, z, 9,
                        left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle,
                        left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter) ;
#endif
              nchanged++ ;
              continue ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
    total_changed += nchanged ;
  } while ((nchanged > 0) && (++niter < 1)) ;

  /* change all gm within 2 mm of ventricle to wm */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri_out_labeled, x, y, z) ;
        if (IS_GM(label) == 0)
          continue ;
        left = label == Left_Cerebral_Cortex ;
        olabel = left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle;
        if (MRIneighborsInWindow(mri_out_labeled, x, y, z, 5, olabel) >= 1) {
          total_changed++ ;
          MRIvox(mri_out_labeled, x, y, z) = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
        }
      }
    }
  }

#if 0
  /* change all unknown within 1 mm of ventricle and 1mm of  wm  to one of them */
  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = MRIvox(mri_out_labeled, x, y, z) ;
        if (IS_UNKNOWN(label) == 0)
          continue ;
        left = label == Left_Cerebral_Cortex ;
        olabel = left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle;
        if (MRIneighborsInWindow(mri_out_labeled, x, y, z, 5, olabel) >= 1) {
          total_changed++ ;
          MRIvox(mri_out_labeled, x, y, z) = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
        }
      }
    }
  }
#endif


  MRIfree(&mri_tmp) ;
  printf("%d unknown voxels bordering ventricles changed to wm.\n",
         total_changed) ;
  return(mri_out_labeled) ;
}
static MRI *
edit_ventricular_unknowns(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int     label, x, y, z, xm1, xp1, ym1, yp1, zm1, zp1, lxp1, lxm1, lyp1, lym1, lzp1, lzm1, nchanged,
  l1 = 0, l2 = 0, change, i, total_changed, dwm, dhypo, dven, left ;
  MRI     *mri_in ;

  mri_in = MRIcopy(mri_in_labeled, NULL) ;
  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;

  total_changed = 0 ;

  /* first change unknowns near hypointensities to wm or hypo */
  do {
    nchanged = 0 ;
    for (x = 0 ; x < mri_T1->width ; x++) {
      for (y = 0 ; y < mri_T1->height ; y++) {
        for (z = 0 ; z < mri_T1->depth ; z++) {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          change = 0 ;
          left = 0 ;
          label = MRIvox(mri_in, x, y, z) ;
          if ((IS_GM(label) == 0) && (IS_UNKNOWN(label) == 0))
            continue ;

          dhypo = distance_to_label(mri_in, WM_hypointensities, x, y, z, 0, 1, 0, 2) ; /* look inferior */
          if (dhypo > 2)
            dhypo = distance_to_label(mri_in, WM_hypointensities, x, y, z, 0, 0, -1, 2) ;  /* look inferior */
          if (dhypo > 2)
            dhypo = distance_to_label(mri_in, WM_hypointensities, x, y, z, 0, 0, 1, 2) ;  /* look anterior */
          if (dhypo > 2)
            dhypo = distance_to_label(mri_in, WM_hypointensities, x, y, z, 0, 0, -1, 2) ;  /* look inferior */
          if (dhypo > 2)
            continue ;

#define WM_DIST   4

          dwm = distance_to_label(mri_in, Left_Cerebral_White_Matter, x, y, z, 0, -1, 0, WM_DIST) ; /* superior */
          if (dwm > WM_DIST) {
            dwm = distance_to_label(mri_in, Right_Cerebral_White_Matter, x, y, z, 0, -1, 0, WM_DIST) ;
          } else
            left = 1 ;

          if (dwm > WM_DIST) /* not superior, look anterior */
          {
            dwm = distance_to_label(mri_in, Right_Cerebral_White_Matter, x, y, z, 0, 0, 1, WM_DIST) ;
            if (dwm > WM_DIST) {
              dwm = distance_to_label(mri_in, Left_Cerebral_White_Matter, x, y, z, 0, 0, 1, WM_DIST) ;
              if (dwm <= WM_DIST)
                left = 1 ;
            }
          }

          if (dwm > WM_DIST) /* not superior, or anterior look posterior */
          {
            dwm = distance_to_label(mri_in, Right_Cerebral_White_Matter, x, y, z, 0, 0, -1, WM_DIST) ;
            if (dwm > WM_DIST) {
              dwm = distance_to_label(mri_in, Left_Cerebral_White_Matter, x, y, z, 0, 0, -1, WM_DIST) ;
              if (dwm <= WM_DIST)
                left = 1 ;
            }
          }

          if (dwm > WM_DIST)    /* couldn't find wm anywhere (didn't look inferior), check for region of hypos */
          {
            dwm = distance_to_label(mri_in, WM_hypointensities, x, y, z, 0, -1, 0, WM_DIST) ;
            if (distance_to_label(mri_in, Left_Cerebral_White_Matter, x, y, z, 0, -1, 0, 30)  <= 30)
              left = 1 ;
          }
          if (dwm > WM_DIST)
            continue ;
          dven = distance_to_label(mri_in, left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle, x, y, z, 0, 1, 0, 15) ;
          if (dven > 15) {
            dven = distance_to_label(mri_in, left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle, x, y, z, 0, 0, -1, 3) ;
            if (dven > 3)
              distance_to_label(mri_in, left ? Left_Lateral_Ventricle : Right_Lateral_Ventricle, x, y, z, 0, 0, 1, 3) ;
            if (dven > 3)
              continue ;
          }

          nchanged++ ;
          mle_label(mri_T1, mri_out_labeled, x, y, z, 15, left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter,
                    WM_hypointensities) ;
        }
      }
    }

    total_changed += nchanged ;
    MRIcopy(mri_out_labeled, mri_in) ;
  } while (nchanged > 0) ;

  for (i = 0 ; i < 2 ; i++) {
    for ( nchanged = x = 0 ; x < mri_T1->width ; x++) {
      for (y = 0 ; y < mri_T1->height ; y++) {
        for (z = 0 ; z < mri_T1->depth ; z++) {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          change = 0 ;
          label = MRIvox(mri_in, x, y, z) ;
          if (IS_UNKNOWN(label) == 0)
            continue ;
          xm1 = mri_T1->xi[x-1] ;
          xp1 = mri_T1->xi[x+1] ;
          ym1 = mri_T1->yi[y-1] ;
          yp1 = mri_T1->yi[y+1] ;
          zm1 = mri_T1->zi[z-1] ;
          zp1 = mri_T1->zi[z+1] ;
          lxp1 = MRIvox(mri_in, xp1, y, z) ;
          lxm1 = MRIvox(mri_in, xm1, y, z) ;
          lyp1 = MRIvox(mri_in, x, yp1, z) ;
          lym1 = MRIvox(mri_in, x, ym1, z) ;
          lzp1 = MRIvox(mri_in, x, y, zp1) ;
          lzm1 = MRIvox(mri_in, x, y, zm1) ;
          /* if it has ventricle on one side and wm on the other, change it to wm or ventricle
            these are caused by slight discrepencies between the FreeSurfer ribbon and the
            CMA subcortical segs.
          */
          if ((IS_WMH(lxm1) && IS_LAT_VENT(lxp1)) ||
              (IS_WMH(lxp1) && IS_LAT_VENT(lxm1))) {
            change = 1 ;
            l1 = lxm1 ;
            l2 = lxp1 ;
          }
          if ((IS_WMH(lym1) && IS_LAT_VENT(lyp1)) ||
              (IS_WMH(lyp1) && IS_LAT_VENT(lym1))) {
            change = 1 ;
            l1 = lym1 ;
            l2 = lyp1 ;
          }
          if ((IS_WMH(lzm1) && IS_LAT_VENT(lzp1)) ||
              (IS_WMH(lzp1) && IS_LAT_VENT(lzm1))) {
            change = 1 ;
            l1 = lzm1 ;
            l2 = lzp1 ;
          }
          if (change) {
            nchanged++ ;
            mle_label(mri_T1, mri_out_labeled, x, y, z, 15, l1, l2) ;
          } else /* could be two unknowns in a row */
          {
            int yi, olabel, wm, un, ven ;
#define WLEN 4
            ven = un = wm = 0 ;
            change = 0 ;
            for (yi = y-1 ; yi >= MAX(0,y-WLEN) ; yi--) {
              olabel = MRIvox(mri_in, x, yi, z) ;
              /* should be only white matter and unkowns above it */
              if (IS_WMH(olabel)) {
                if (IS_WM(olabel))
                  l1 = olabel ;
                wm++ ;
              } else if (IS_UNKNOWN(olabel))
                un++ ;
              else
                change = -1 ;  /* something else there - don't change */
            }
            /* if there is some wm above and no other labels */
            if ((wm >= WLEN/2) && ((wm+un) >= (WLEN-1)) && change >= 0) {
              un = 0 ;
              for (yi = y+1 ; yi <= MIN(mri_in->height-1,y+WLEN) ; yi++) {
                olabel = MRIvox(mri_in, x, yi, z) ;
                /* should be only ventricle and unkowns below it */
                if (IS_LAT_VENT(olabel)) {
                  ven++ ;
                  l2 = olabel ;
                } else if (IS_UNKNOWN(olabel))
                  un++ ;
                else
                  change = -1 ;
              }
              if (change >= 0 && ((ven+un) >= WLEN) && (ven >= WLEN/2))
                change = 1 ;
              if (change > 0) {
                nchanged++ ;
                mle_label(mri_T1, mri_out_labeled, x, y, z, 15, l1, l2) ;
              }
            }
            if (change <= 0)  /* look in posterior/anterior direction. If everthing is wm and ven change it */
            {
              int zi ;

              change = wm = ven = un = 0 ;
              for (zi = z-WLEN ; zi <= z+WLEN ; zi++) {
                if (zi < 0 || zi >= mri_in->depth)
                  continue ;
                olabel = MRIvox(mri_in, x, y, zi) ;
                /* should be only white matter and unkowns above it */
                if (IS_WMH(olabel)) {
                  if (IS_WM(olabel))
                    l1 = olabel ;
                  wm++ ;
                } else if (IS_UNKNOWN(olabel))
                  un++ ;
                else if (IS_LAT_VENT(olabel))
                  ven++ ;
                else
                  change = -1 ;  /* something else there - don't change */
              }
              if (change >= 0 && ((ven+wm) >= WLEN) && (ven >= WLEN/2) && (wm >= WLEN/2))
                change = 1 ;
              if (change > 0) {
                nchanged++ ;
                mle_label(mri_T1, mri_out_labeled, x, y, z, 15, l1, l2) ;
              }
            }
            if (change <= 0)  /* look in medial/lateral direction. If everthing is wm and ven change it */
            {
              int xi ;

              change = wm = ven = un = 0 ;
              for (xi = x-WLEN ; xi <= x+WLEN ; xi++) {
                if (xi < 0 || xi >= mri_in->width)
                  continue ;
                olabel = MRIvox(mri_in, xi, y, z) ;
                /* should be only white matter and unkowns above it */
                if (IS_WMH(olabel)) {
                  if (IS_WM(olabel))
                    l1 = olabel ;
                  wm++ ;
                } else if (IS_UNKNOWN(olabel))
                  un++ ;
                else if (IS_LAT_VENT(olabel)) {
                  ven++ ;
                  l2 = olabel ;
                } else
                  change = -1 ;  /* something else there - don't change */
              }
              if (change >= 0 && ((ven+wm) >= WLEN) && (ven >= WLEN/2) && (wm >= WLEN/2))
                change = 1 ;
              if (change > 0) {
                nchanged++ ;
                mle_label(mri_T1, mri_out_labeled, x, y, z, 15, l1, l2) ;
              }
            }
          }
        }
      }
    }

    /* now relabel cortex that is between WM and hypointensity */
    for (x = 0 ; x < mri_T1->width ; x++) {
      for (y = 0 ; y < mri_T1->height ; y++) {
        for (z = 0 ; z < mri_T1->depth ; z++) {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          change = 0 ;
          label = MRIvox(mri_in, x, y, z) ;
          if (IS_GM(label) == 0 && (IS_UNKNOWN(label) == 0))
            continue ;
          xm1 = mri_T1->xi[x-1] ;
          xp1 = mri_T1->xi[x+1] ;
          ym1 = mri_T1->yi[y-1] ;
          yp1 = mri_T1->yi[y+1] ;
          zm1 = mri_T1->zi[z-1] ;
          zp1 = mri_T1->zi[z+1] ;
          lxp1 = MRIvox(mri_in, xp1, y, z) ;
          lxm1 = MRIvox(mri_in, xm1, y, z) ;
          lyp1 = MRIvox(mri_in, x, yp1, z) ;
          lym1 = MRIvox(mri_in, x, ym1, z) ;
          lzp1 = MRIvox(mri_in, x, y, zp1) ;
          lzm1 = MRIvox(mri_in, x, y, zm1) ;
          if ((IS_WM(lxm1) && IS_HYPO(lxp1)) ||
              (IS_WM(lxp1) && IS_HYPO(lxm1))) {
            change = 1 ;
            l1 = lxm1 ;
            l2 = lxp1 ;
          }
          if ((IS_WM(lym1) && IS_HYPO(lyp1)) ||
              (IS_WM(lyp1) && IS_HYPO(lym1))) {
            change = 1 ;
            l1 = lym1 ;
            l2 = lyp1 ;
          }
          if ((IS_WM(lzm1) && IS_HYPO(lzp1)) ||
              (IS_WM(lzp1) && IS_HYPO(lzm1))) {
            change = 1 ;
            l1 = lzm1 ;
            l2 = lzp1 ;
          }


          if ((IS_HYPO(lxm1) && IS_HYPO(lxp1)) ||
              (IS_HYPO(lxp1) && IS_HYPO(lxm1))) {
            change = 1 ;
            l1 = lxm1 ;
            l2 = lxp1 ;
          }
          if ((IS_HYPO(lym1) && IS_HYPO(lyp1)) ||
              (IS_HYPO(lyp1) && IS_HYPO(lym1))) {
            change = 1 ;
            l1 = lym1 ;
            l2 = lyp1 ;
          }
          if ((IS_HYPO(lzm1) && IS_HYPO(lzp1)) ||
              (IS_HYPO(lzp1) && IS_HYPO(lzm1))) {
            change = 1 ;
            l1 = lzm1 ;
            l2 = lzp1 ;
          }
          if (change) {
            nchanged++ ;
            mle_label(mri_T1, mri_out_labeled, x, y, z, 15, l1, l2) ;
          }
        }
      }
    }
    MRIcopy(mri_out_labeled, mri_in) ;
    total_changed += nchanged ;
  }


  printf("%d unknown voxels changed to wm or ventricle\n", total_changed) ;
  MRIfree(&mri_in) ;
  return(mri_out_labeled) ;
}

static MRI *
edit_border_voxels(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled) {
  int     label, x, y, z, xm1, xp1, ym1, yp1, zm1, zp1, lxp1, lxm1, lyp1, lym1, lzp1, lzm1, nchanged,
  change, olabel ;
  float   means[MAX_CMA_LABELS], xp1d, xm1d, yp1d, ym1d, zp1d, zm1d, val, ld ;
  MRI     *mri_tmp ;

  if (mri_in_labeled == mri_out_labeled) {
    mri_tmp = MRIcopy(mri_in_labeled, NULL) ;
    mri_in_labeled = mri_tmp;
  } else
    mri_tmp = NULL ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  for (nchanged = x = 0 ; x < mri_T1->width ; x++) {
    for (y = 0 ; y < mri_T1->height ; y++) {
      for (z = 0 ; z < mri_T1->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        change = 0 ;
        olabel = label = MRIvox(mri_in_labeled, x, y, z) ;
        if (IS_UNKNOWN(label))
          continue ;
        xm1 = mri_T1->xi[x-1] ;
        xp1 = mri_T1->xi[x+1] ;
        ym1 = mri_T1->yi[y-1] ;
        yp1 = mri_T1->yi[y+1] ;
        zm1 = mri_T1->zi[z-1] ;
        zp1 = mri_T1->zi[z+1] ;
        lxp1 = MRIvox(mri_in_labeled, xp1, y, z) ;
        lxm1 = MRIvox(mri_in_labeled, xm1, y, z) ;
        lyp1 = MRIvox(mri_in_labeled, x, yp1, z) ;
        lym1 = MRIvox(mri_in_labeled, x, ym1, z) ;
        lzp1 = MRIvox(mri_in_labeled, x, y, zp1) ;
        lzm1 = MRIvox(mri_in_labeled, x, y, zm1) ;
        if (label == lxm1 &&
            label == lxp1 &&
            label == lym1 &&
            label == lyp1 &&
            label == lzm1 &&
            label == lzp1)
          continue ;  /* it's not a border label */
        if (!(IS_WM(label) || IS_HIPPO(label) || IS_GM(label) || IS_UNKNOWN(label) || IS_LAT_VENT(label)))
          continue ;


        val = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        means[label] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, label) ;
        ld = fabs(means[label] - val) ;

        if (ld/means[label] < 0.1)
          continue ;

        ld *= 0.5 ;  /* bias towards keeping same label */


        if (lxm1 != label)
          means[lxm1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lxm1) ;
        if (lxp1 != label)
          means[lxp1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lxp1) ;
        if (lym1 != label)
          means[lym1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lym1) ;
        if (lyp1 != label)
          means[lyp1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lyp1) ;
        if (lzm1 != label)
          means[lzm1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lzm1) ;
        if (lzp1 != label)
          means[lzp1] = label_mean(mri_T1, mri_in_labeled, x, y, z, 15, lzp1) ;
        xp1d = fabs(means[lxp1] - val) ;
        if (xp1d < ld) {
          ld = xp1d ;
          olabel = lxp1 ;
        }
        xm1d = fabs(means[lxm1] - val) ;
        if (xm1d < ld) {
          ld = xm1d ;
          olabel = lxm1 ;
        }
        yp1d = fabs(means[lyp1] - val) ;
        if (yp1d < ld) {
          olabel = lyp1 ;
          ld = yp1d ;
        }
        ym1d = fabs(means[lym1] - val) ;
        if (ym1d < ld) {
          olabel = lym1 ;
          ld = ym1d ;
        }
        zp1d = fabs(means[lzp1] - val) ;
        if (zp1d < ld) {
          olabel = lzp1 ;
          ld = zp1d ;
        }
        zm1d = fabs(means[lzm1] - val) ;
        if (zp1d < ld) {
          olabel = lzp1 ;
          ld = zp1d ;
        }

        /* only let certain labels change */
        if (!(IS_WM(olabel) || IS_HIPPO(olabel) || IS_GM(olabel) || IS_UNKNOWN(olabel) || IS_LAT_VENT(olabel)))
          continue ;
        if ((IS_GM(label) && IS_HIPPO(olabel)) ||
            (IS_GM(olabel) && IS_HIPPO(label)))
          continue ;  /* don't change hippo to gm based on intensity - too similar */
        if ((label != olabel) && MRIneighborsInWindow(mri_in_labeled, x, y, z, 3, olabel) > 5) {
          nchanged++ ;
          MRIvox(mri_out_labeled, x, y, z) = olabel ;
        }
      }
    }
  }

  if (mri_tmp)
    MRIfree(&mri_tmp) ;
  printf("%d border voxels changed\n", nchanged) ;
  return(mri_out_labeled) ;
}
