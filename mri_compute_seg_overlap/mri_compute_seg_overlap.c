/**
 * @file  mri_compute_seg_overlap.c
 * @brief Compute Dice coefficent comparing two segmentation volumes.
 *
 * This program compares two segmentation volumes and
 * computes the Dice and Jaccard Coefficients.
 * It considers only 9 major structures.
 */
/*
 * Original Authors: Xiao Han, Nick Schmansky 
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/11/11 22:41:09 $
 *    $Revision: 1.10 $
 *
 * Copyright (C) 2006-2010,
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
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"
#include "cma.h"

static void usage(int exit_val);
static int  get_option(int argc, char *argv[]) ;

char *Progname;

static char *log_fname = NULL ; //dice coeff of individual structures
static char *mlog_fname = NULL ; //mean of individual dice
static char *slog_fname = NULL ; //std of individual dice
static char *olog_fname = NULL ; //overall dice for subcortical structures

static const int num_labels       = 24;
static const int labels_of_interest[24] = {
      Left_Cerebral_White_Matter, Right_Cerebral_White_Matter,
      Left_Cerebral_Cortex, Right_Cerebral_Cortex,
      Left_Lateral_Ventricle, Right_Lateral_Ventricle,
      Left_Hippocampus, Right_Hippocampus,
      Left_Thalamus_Proper, Right_Thalamus_Proper,
      Left_Caudate, Right_Caudate,
      Left_Putamen, Right_Putamen,
      Left_Pallidum,Right_Pallidum,
      Left_Amygdala, Right_Amygdala,
      Left_Accumbens_area, Right_Accumbens_area,
      Third_Ventricle, Fourth_Ventricle,
      Left_Inf_Lat_Vent, Right_Inf_Lat_Vent
    };

/* Note: these are the labels included in the
     'overall subcortical Dice coefficient' calculations.
     It excludes:
     Left/Right-Cerebral-White-Matter (labels 2 and 41),
     Left/Right-Cerebral-Cortex (labels 3 and 42),
     Left/Right-Accumbens-area (labels 26 and 58) */
static const int num_labels_overall_Dice = 18;
static const int labels_overall_Dice[18] = {
      Left_Lateral_Ventricle, Right_Lateral_Ventricle,
      Left_Hippocampus, Right_Hippocampus,
      Left_Thalamus_Proper, Right_Thalamus_Proper,
      Left_Caudate, Right_Caudate,
      Left_Putamen, Right_Putamen,
      Left_Pallidum,Right_Pallidum,
      Left_Amygdala, Right_Amygdala,
      Third_Ventricle, Fourth_Ventricle,
      Left_Inf_Lat_Vent, Right_Inf_Lat_Vent
    };
// returns 1 if volVal is one of the labels to include in overall Dice calc
static int isOverallDiceLabel(int volVal) {
  int i;
  for (i=0; i < num_labels_overall_Dice; i++) {
    if (volVal == labels_overall_Dice[i]) return 1;
  }
  return 0;
}

/* maximum number of classes */
#define MAX_CLASSES 256
#define MAX_CLASS_NUM 255

int main(int argc, char *argv[]) {
  MRI *mri_seg1, *mri_seg2;
  int nargs, ac;
  char **av;
  int width, height, depth, x, y, z, f, nframes;
  int v1, v2;
  int i;
  FILE *log_fp;

  int Volume_union[MAX_CLASSES];
  int Volume_from1[MAX_CLASSES];
  int Volume_from2[MAX_CLASSES];
  int Volume_overlap[MAX_CLASSES];
  int subcorvolume1, subcorvolume2;
  int subcorvolume_overlap;

  double mean1, std1, mean2, std2;

  float correct_ratio[MAX_CLASSES];
  float correct_ratio2[MAX_CLASSES];

  Progname = argv[0];

  nargs = 
    handle_version_option
    (argc, argv,
     "$Id: mri_compute_seg_overlap.c,v 1.10 2010/11/11 22:41:09 mreuter Exp $",
     "$Name:  $");
  argc -= nargs ;
  if (1 == argc)
    exit (0);

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 3)
    usage(1);

  mri_seg1 = MRIread(argv[1]) ;
  if (!mri_seg1)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume1 %s",
              Progname, argv[1]) ;

  mri_seg2 = MRIread(argv[2]) ;
  if (!mri_seg2)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume2 %s",
              Progname, argv[2]) ;

  if ((mri_seg1->width != mri_seg2->width) ||
      (mri_seg1->height != mri_seg2->height) ||
      (mri_seg1->depth != mri_seg2->depth))
    ErrorExit(ERROR_BADPARM,
              "%s: two input label volumes have different sizes \n",
              Progname);

  for (i=0; i < MAX_CLASSES; i++) {
    Volume_union[i] = 0;
    Volume_overlap[i] = 0;
    Volume_from1[i] = 0;
    Volume_from2[i] = 0;
  }

  width = mri_seg1->width ;
  height = mri_seg1->height ;
  depth = mri_seg1->depth ;
  nframes = mri_seg1->nframes ;
  if (nframes == 0) nframes = 1;

  subcorvolume_overlap = 0;
  subcorvolume1 = 0;
  subcorvolume2 = 0;

  for (f = 0 ; f < nframes ; f++) {
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          v1 = (int) MRIgetVoxVal(mri_seg1,x,y,z,f);
          v2 = (int) MRIgetVoxVal(mri_seg2,x,y,z,f);

          if (v1 > MAX_CLASS_NUM ||
              v1 <= 0 ||
              v2 > MAX_CLASS_NUM ||
              v2 <= 0) continue;

          /* do not include these in the overall Dice coefficient calculations:
             Left/Right-Cerebral-White-Matter (labels 2 and 41),
             Left/Right-Cerebral-Cortex (labels 3 and 42),
             Left/Right-Accumbens-area (labels 26 and 58)
             Notice that these labels are not included in the 'if' checks: */

          if (v1 == v2) {
            if (isOverallDiceLabel(v1)) subcorvolume_overlap++;
          }

          if (isOverallDiceLabel(v1)) subcorvolume1++;
          if (isOverallDiceLabel(v2)) subcorvolume2++;

          Volume_from1[v1]++;
          Volume_from2[v2]++;

          if (v1 == v2) {
            Volume_overlap[v1]++;
            Volume_union[v1]++;
          } else {
            Volume_union[v1]++;
            Volume_union[v2]++;
          }
        }
      }
    }
  }

  for (i=0; i < MAX_CLASSES; i++) {
    correct_ratio[i] =
      (double)Volume_overlap[i]/((double)Volume_union[i] + 1e-10);
    correct_ratio2[i] =
      (double)Volume_overlap[i]*2.0/
      ((double)Volume_from1[i] + (double)Volume_from2[i] + 1e-10);
  }

  printf("Jaccard Coefficients:\n");
  mean1 = 0;
  std1 = 0;
  for (i=0; i < num_labels; i++) {
    printf("correct ratio for label %d = %g\n",
           labels_of_interest[i], correct_ratio[labels_of_interest[i]]);
    mean1 += correct_ratio[labels_of_interest[i]];
    std1 += correct_ratio[labels_of_interest[i]] *
            correct_ratio[labels_of_interest[i]];
  }
  mean1 /= num_labels;
  std1 /= num_labels;
  std1 = sqrt(std1 - mean1*mean1);
  printf("mean +/- std = %6.4f +/- %6.4f\n", mean1, std1);

  printf("Dice Coefficients:\n");
  // printf("ratio of overlap to volume of input1:\n");
  mean2 = 0;
  std2 = 0;
  for (i=0; i < num_labels; i++) {
    printf("label %d = %g\n",
           labels_of_interest[i], correct_ratio2[labels_of_interest[i]]);
    mean2 += correct_ratio2[labels_of_interest[i]];
    std2 += correct_ratio2[labels_of_interest[i]] *
            correct_ratio2[labels_of_interest[i]];
  }
  mean2 /= num_labels;
  std2 /= num_labels;
  std2 = sqrt(std2 - mean2*mean2);
  printf("mean +/- std = %6.4f +/- %6.4f \n", mean2, std2);

  if (log_fname != NULL) {
    log_fp = fopen(log_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
    for (i=0; i < num_labels; i++) {
      fprintf(log_fp, "%6.4f ", correct_ratio2[labels_of_interest[i]]);
    }
    fprintf(log_fp, "%6.4f ", mean2);
    fprintf(log_fp, "%6.4f ", std2);
    fprintf(log_fp, "%6.4f \n",
            subcorvolume_overlap*2.0/(float)(subcorvolume1 + subcorvolume2));
    fclose(log_fp);
  }

  if (mlog_fname != NULL) {
    log_fp = fopen(mlog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, mlog_fname) ;
    fprintf(log_fp, "%6.4f \n", mean2);
    fclose(log_fp);
  }

  if (slog_fname != NULL) {
    log_fp = fopen(slog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, slog_fname) ;
    fprintf(log_fp, "%6.4f \n", std2);
    fclose(log_fp);
  }

  if (olog_fname != NULL) {
    log_fp = fopen(olog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, olog_fname) ;
    fprintf(log_fp, "%6.4f \n",
            subcorvolume_overlap*2.0/(float)(subcorvolume1 + subcorvolume2));
    fclose(log_fp);
  }

  printf("Overall subcortical dice = %6.4f \n",
         subcorvolume_overlap*2.0/(float)(subcorvolume1 + subcorvolume2));

  exit(0);

}  /*  end main()  */


static void usage(int exit_val) {
  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <seg vol1> <seg vol2>\n", Progname);
  fprintf(fout, "This program compares two segmentation volumes and \n"
          "computes the Dice and Jaccard Coefficients. \n"
          "It considers only 9 major structures:\n"
          "  L/R Hippocampus\n"
          "  L/R Caudate\n"
          "  L/R Putamen\n"
          "  L/R Pallidum\n"
          "  L/R Amygdala\n"
          "  L/R Thalamus_Proper\n"
          "  L/R Lateral_Ventricle\n"
          "  Third and Fourth Ventricles\n"
          "  L/R Inf_Lat_Vent\n");
  fprintf(fout, "Options:\n");
  fprintf(fout, "   -log %%s   log_file for individual Dice \n");
  fprintf(fout, "   -mlog %%s  log_file for mean Dice \n");
  fprintf(fout, "   -slog %%s  log_file for std Dice \n");
  fprintf(fout, "   -olog %%s  log_file for overall Dice \n");
  exit(exit_val);

}  /*  end usage()  */


static int get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    usage(0) ;
  else if (!stricmp(option, "mlog")) {
    mlog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging mean Dice to %s\n", mlog_fname) ;
  } else if (!stricmp(option, "log") || !stricmp(option, "L")) {
    log_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging individual Dice to %s\n", log_fname) ;
  } else if (!stricmp(option, "slog")) {
    slog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging std Dice to %s\n", slog_fname) ;
  } else if (!stricmp(option, "olog")) {
    olog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging overall Dice to %s\n", olog_fname) ;
  } else {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(0) ;
    exit(1) ;
  }

  return(nargs) ;
}


/*  EOF  */
