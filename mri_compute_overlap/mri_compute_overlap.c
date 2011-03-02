/**
 * @file  mri_compute_overlap.c
 * @brief computing label overlap measures
 *
 * Computing three label overlap measures: volume difference, Dice and Jaccard.
 */
/*
 * Original Author: Nick S.?
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.17 $
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
#include "timer.h"
#include "version.h"
#include "gca.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static char *log_fname = NULL ;
static void usage_exit(int code) ;
static void print_help(void) ;

static int quiet = 0 ;
static int all_flag = 0 ;

static int in_label = -1 ;
static int out_label = -1 ;

static int isSeg  = 0;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, lno, nshared, nvox1, nvox2, total_nvox1, total_nvox2,
    total_nshared, nlabels, i ;
  int          msec, minutes, seconds/*, wrong, total, correct*/ ;
  struct timeb start ;
  MRI    *mri1, *mri2 ;
  FILE   *log_fp ;
  float  nvox_mean, nunion, total_nunion ;
  float tmp = 0.0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mri_compute_overlap.c,v 1.17 2011/03/02 00:04:14 nicks Exp $",
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

  if (argc < 3)
    usage_exit(1) ;

  mri1 = MRIread(argv[1]) ;
  if (!mri1)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[1]) ;
  mri2 = MRIread(argv[2]) ;
  if (!mri2)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname,
              argv[2]) ;

  if (in_label >= 0) {
    MRIreplaceValues(mri1, mri1, in_label, out_label) ;
    MRIreplaceValues(mri2, mri2, in_label, out_label) ;
  }

  if (log_fname) {
    //log_fp = fopen(log_fname, "a+") ;
    log_fp = fopen(log_fname, "w+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
  } else
    log_fp = NULL ;
  
  nlabels = total_nvox1 = total_nvox2 = total_nshared = 0;
  total_nunion = 0 ;
  if (all_flag) {
    MRI *mri1_label = NULL, *mri2_label = NULL ;
    int lnoLimit = 1000;

    mri1_label = MRIclone(mri1, NULL) ;
    mri2_label = MRIclone(mri2, NULL) ;

    if (isSeg)
      lnoLimit = MAX_CMA_LABEL;

    for (lno = 0 ; lno < lnoLimit ; lno++) {
#if 1
      nvox1 = MRIvoxelsInLabel(mri1, lno) ;
      nvox2 = MRIvoxelsInLabel(mri2, lno) ;
      if (!nvox1 && !nvox2)
        continue ;
      nvox_mean = (float)(nvox1+nvox2)/2.0f ;
      nshared = MRIlabelOverlap(mri1, mri2, lno) ;
      nunion  = (float)MRIlabelUnion(mri1, mri2, lno) ;
      if (nunion > 0.0) tmp = (float)nshared/nunion;
      else tmp = 0.0;

      if (!quiet) {
	if (isSeg) {
	  printf
	    ("label = %d (%s), volume diff = |(%d - %d)|"
	     " / %2.1f (%%)= %2.2f%%\n",
	     lno, cma_label_to_name(lno),
	     nvox1, nvox2, nvox_mean,
	     100.0f*(float)abs(nvox1-nvox2)/nvox_mean);
	  printf
	    ("label = %d (%s), volume overlap (Dice) = "
	     "%d / %2.1f (%%)= %2.2f %%\n",
	     lno, cma_label_to_name(lno), nshared, nvox_mean,
	     100.0f*(float)nshared/nvox_mean) ;
	  printf
	    ("label = %d (%s), volume overlap (Jaccard) = "
	     "%d / %2.1f (%%)= %2.2f%%\n",
	     lno, cma_label_to_name(lno), nshared, nunion,
	     100.0f*tmp) ;
	} else {
	  printf("volume diff = |(%d - %d)| / %2.1f (%%)= %2.2f %%\n",
		 nvox1, nvox2, nvox_mean,
		 100.0f*(float)abs(nvox1-nvox2)/nvox_mean);
	  printf("volume overlap (Dice) = %d / %2.1f (%%)= %2.2f%%\n",
		 nshared, nvox_mean, 100.0f*(float)nshared/nvox_mean) ;
	  printf("volume overlap (Jaccard) = %d / %2.2f (%%)= %2.2f%%\n",
		 nshared, nunion, 100.0f*tmp) ;
	}
      }
      if (log_fp) {
        fprintf(log_fp, "%d\t%2.2f\t%2.2f\t%2.2f\n", lno,
                100.0f*(float)abs(nvox1-nvox2)/nvox_mean,
                100.0f*(float)nshared/nvox_mean,
		100.0f*(float)nshared/nunion) ;
      }
#else
      nvox1 = MRIcopyLabel(mri1, mri1_label, lno) ;
      nvox2 = MRIcopyLabel(mri2, mri1_label, lno) ;
      if (!nvox1 && !nvox2)
        continue ;
      correct = MRIlabelOverlap(mri1, mri2, lno) ;
      total = (nvox1+nvox2)/2 ;
      wrong = total-correct ;
      printf("label %03d: %d of %d voxels correctly labeled - %2.2f%%\n",
             lno, correct, total, 100.0f*(float)correct/(float)total) ;
      if (log_fp) {
        fprintf(log_fp,"%d\t%d\t%d\t%2.4f\n", lno, correct, total,
                100.0f*(float)correct/(float)total) ;
        fclose(log_fp) ;
      }
#endif
      total_nvox1 += nvox1 ;
      total_nvox2 += nvox2 ;
      total_nshared += nshared ;
      total_nunion += nunion ;
      nlabels++ ;
    }
    if (log_fp)
      fclose(log_fp) ;
  } else {
    
    for (i = 3 ; i < argc ; i++) {
      float volume_overlap, volume_diff, volume_overlap_jacc ;

      lno = atoi(argv[i]) ;
      // only counts number of lno label
      nvox1 = MRIvoxelsInLabel(mri1, lno) ;
      nvox2 = MRIvoxelsInLabel(mri2, lno) ;
      nvox_mean = (float)(nvox1+nvox2)/2 ;
      // if both mri1 and mri2 has the same label, count it.
      nshared = MRIlabelOverlap(mri1, mri2, lno) ;
      nunion  = MRIlabelUnion(mri1, mri2, lno) ;
      volume_diff = 100.0f*(float)abs(nvox1-nvox2)/nvox_mean ;
      volume_overlap = 100.0f*(float)nshared/nvox_mean ;
      if(nunion>0.)
	volume_overlap_jacc = 100.0f*(float)nshared/nunion ;
      else
	volume_overlap_jacc = 0.0 ;
      if (!quiet) {
        printf("label %d: volume diff = |(%d - %d)| / %2.1f (%%)= %2.2f%%\n",
               lno,nvox1,nvox2,nvox_mean, volume_diff);
        printf("label %d: volume overlap (Dice) = %d / %2.1f (%%)= %2.2f%%\n",
               lno, nshared, nvox_mean, volume_overlap) ;
        printf("label %d: volume overlap (Jaccard) = %d / %2.1f (%%)= %2.2f%%\n",
               lno, nshared, nunion, volume_overlap_jacc) ;
      }
      if (log_fp) {
        fprintf(log_fp, "%2.2f\t%2.2f\t%2.2f\n", volume_diff, volume_overlap, volume_overlap_jacc) ;
        fclose(log_fp) ;
      }
      total_nvox1 += nvox1 ;
      total_nvox2 += nvox2 ;
      total_nshared += nshared ;
      total_nunion += nunion ;
      nlabels++ ;
    }
  }

  if (nlabels > 1) {
    nvox_mean = (total_nvox1 + total_nvox2) / 2 ;
    printf("total: volume diff = |(%d - %d)| / %2.1f (%%)= %2.2f%%\n",
	   total_nvox1,total_nvox2,nvox_mean,
	   100.0f*(float)abs(total_nvox1-total_nvox2)/nvox_mean);
    printf("total: volume overlap (Dice) = %d / %2.1f (%%) = %2.2f%%\n",
	   total_nshared, nvox_mean,
	   100.0f*(float)total_nshared/nvox_mean) ;
    if (total_nunion>0.0)
      tmp = (float)total_nshared/total_nunion;
    else tmp = 0.0;
    printf("total: volume overlap (Jaccard) = %d / %2.1f (%%)= %2.2f%%\n",
	   total_nshared, total_nunion, 100.0f*tmp) ;
    }
  
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "overlap calculation took %d minutes and %d seconds.\n",
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
  switch (toupper(*option)) {
  case 'T':
    in_label = atoi(argv[2]) ;
    out_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("translating label %d to label %d\n", in_label, out_label) ;
    break ;
  case 'Q':
    quiet = 1 ;
    break ;
  case 'A':
    all_flag = 1 ;
    fprintf(stderr, "print all labels\n");
    break ;
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "logging results to %s\n", log_fname) ;
    quiet = 1;
    break ;
  case 'S':
    isSeg = 1;
    fprintf(stderr, "show segmentation label names\n");
    break;
  case 'H':
    print_help();
    break;
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
  printf("usage: %s [options] <volume 1> <volume 2> [label numbers]\n",
         Progname) ;
  
  printf("\n") ;
  printf("Options:\n") ;
  printf("  -a compute overlap of all labels (if missing, labels of interest should be listed)\n");
  printf("  -s show label name for segmentation\n");
  printf("  -l <fname> filename to write results to\n");
  printf("  -q do not display results on std display (default is false; if -l is used, this option is set)\n");
  printf("  -t <l1>  <l2> translating label l1 to label l2\n");
  printf("  -h print help\n");

  exit(code) ;
}

static void print_help(void) {
  printf("\n");
  printf("mri_compute_overlap\n");
  printf("\n");
  printf("Computes three different types of overlap measures either for all the existing labels\n");
  printf("in the input volumes or a subset of them that the users lists in the command line. The\n");
  printf("three overlap measures are the following: \n");
  printf("(1) `volume difference` = 2*|A|-|B|/|A|+|B|\n");
  printf("(2) `volume overlap (Dice)`    = 2*|AB|/|A|+|B|; the intersection over the mean of the segmentation volumes\n");
  printf("(3) `volume overlap (Jaccard)` = |AB|/|A+B|; the intersection over the union of the segmentation volumes\n");
  printf("\n");
  printf("\n");
  usage_exit(1) ;
}
