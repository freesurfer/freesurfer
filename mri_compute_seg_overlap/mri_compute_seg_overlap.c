/* This is a simple program to compare two segmentation volumes to compute
 * the DICE and Jaccard coefficients. It only considers 12 major structures
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

#define CREAT_CSF 1
#define CHOICE 1

void usage(int exit_val);
MRI *fliplr(MRI *src);
static int  get_option(int argc, char *argv[]) ;

char *Progname;

static char *log_fname = NULL ; //dice coeff of individual structures
static char *mlog_fname = NULL ; //mean of individual dice
static char *slog_fname = NULL ; //std of individual dice
static char *olog_fname = NULL ; //overall dice for subcortical structures

static int num_labels  = 24;
static int labels_of_interest[24] = {2, 41, 3, 42,
                                     4, 43, 17, 53,
                                     10, 49, 11, 50,
                                     12, 51, 13, 52,
                                     18, 54, 26, 58,
                                     14, 15, 5, 44};
/* maximum number of classes */
#define MAX_CLASSES 256

int main(int argc, char *argv[])
{
  MRI *mri_seg1, *mri_seg2;
  int nargs, ac;
  char **av;
  int width, height, depth, x, y, z, f,nframes ;
  int v1, v2;
  int i;
  FILE *log_fp;

  int Volume_union[256];
  int Volume_from1[256];
  int Volume_from2[256];
  int Volume_overlap[256];
  int subcorvolume1, subcorvolume2;
  int subcorvolume_overlap;

  double mean1, std1, mean2, std2;

  float correct_ratio[256];
  float correct_ratio2[256];

  Progname = argv[0];

  nargs = handle_version_option
    (argc, argv,
     "$Id: mri_compute_seg_overlap.c,v 1.2 2006/08/24 18:22:01 nicks Exp $",
     "$Name:  $");
  argc -= nargs ;
  if (1 == argc)
    exit (0);

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
    {
      nargs = get_option(argc, argv) ;
      argc -= nargs ;
      argv += nargs ;
    }

  if(argc != 3)
    usage(1);

  mri_seg1 = MRIread(argv[1]) ;
  if (!mri_seg1)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume1 %s",
              Progname, argv[1]) ;

  mri_seg2 = MRIread(argv[2]) ;
  if (!mri_seg2)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume2 %s",
              Progname, argv[2]) ;

  if((mri_seg1->width != mri_seg2->width) ||
     (mri_seg1->height != mri_seg2->height) ||
     (mri_seg1->depth != mri_seg2->depth))
    ErrorExit(ERROR_BADPARM,
              "%s: two input label volumes have different sizes \n",
              Progname);

  for(i=0; i < MAX_CLASSES; i++){
    Volume_union[i] = 0;
    Volume_overlap[i] = 0;
    Volume_from1[i] = 0;
    Volume_from2[i] = 0;
  }

  width = mri_seg1->width ;
  height = mri_seg1->height ;
  depth = mri_seg1->depth ;
  nframes = mri_seg1->nframes ;
  if(nframes == 0) nframes = 1;

  subcorvolume_overlap = 0;
  subcorvolume1 = 0;
  subcorvolume2 = 0;

  for (f = 0 ; f < nframes ; f++){
    for (z = 0 ; z < depth ; z++){
      for (y = 0 ; y < height ; y++){
        for (x = 0 ; x < width ; x++){
          v1 = (int) MRIgetVoxVal(mri_seg1,x,y,z,f);
          v2 = (int) MRIgetVoxVal(mri_seg2,x,y,z,f);

          if(v1 > 255 || v1 <= 0 || v2 > 255 || v2 <= 0) continue;

          if(v1 == v2){
            if((v1>=4 && v1 <= 5) ||
               (v1 >= 10 && v1 <= 15) ||
               (v1 >= 17 && v1 <= 18) ||
               (v1>=43 && v1 <= 44) ||
               (v1 >= 49 && v1 <= 54 ))
              subcorvolume_overlap++;
          }

          if((v1>=4 && v1 <= 5) ||
             (v1 >= 10 && v1 <= 15) ||
             (v1 >= 17 && v1 <= 18) ||
             (v1>=43 && v1 <= 44) ||
             (v1 >= 49 && v1 <= 54 ))
            subcorvolume1++;

          if((v2>=4 && v2 <= 5) ||
             (v2 >= 10 && v2 <= 15) ||
             (v2 >= 17 && v2 <= 18) ||
             (v2>=43 && v2 <= 44) ||
             (v2 >= 49 && v2 <= 54 ) )
            subcorvolume2++;

          Volume_from1[v1]++;
          Volume_from2[v2]++;

          if(v1 == v2){
            Volume_overlap[v1]++;
            Volume_union[v1]++;
          }else{
            Volume_union[v1]++;
            Volume_union[v2]++;
          }
        }
      }
    }
  }

  for(i=0; i < MAX_CLASSES; i++){
    correct_ratio[i] =
      (double)Volume_overlap[i]/((double)Volume_union[i] + 1e-10);
    correct_ratio2[i] =
      (double)Volume_overlap[i]*2.0/
      ((double)Volume_from1[i] + (double)Volume_from2[i] + 1e-10);
  }

  printf("Jaccard Coefficients:\n");
  mean1 = 0; std1 = 0;
  for(i=0; i < num_labels; i++){
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
  mean2 = 0; std2 = 0;
  for(i=0; i < num_labels; i++){
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

  if(log_fname != NULL){
    log_fp = fopen(log_fname, "w+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
    for(i=0; i < num_labels; i++){
      fprintf(log_fp, "%6.4f ", correct_ratio2[labels_of_interest[i]]);
    }
    fprintf(log_fp, "%6.4f ", mean2);
    fprintf(log_fp, "%6.4f ", std2);
    fclose(log_fp);
  }

  if(mlog_fname != NULL){
    log_fp = fopen(mlog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, mlog_fname) ;
    fprintf(log_fp, "%6.4f \n", mean2);
    fclose(log_fp);
  }

  if(slog_fname != NULL){
    log_fp = fopen(slog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, slog_fname) ;
    fprintf(log_fp, "%6.4f \n", std2);
    fclose(log_fp);
  }

  if(olog_fname != NULL){
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


void usage(int exit_val)
{
  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <seg vol1> <seg vol2>\n", Progname);
  fprintf(fout, "This program compares two segmentation volumes and \n"
          "computes the DICE and Jaccard Coefficients. \n"
          "It considers only 12 major structures. \n") ;
  fprintf(fout, "Options:\n");
  fprintf(fout, "   -log %%s   log_file for individual DICE \n");
  fprintf(fout, "   -mlog %%s  log_file for mean DICE \n");
  fprintf(fout, "   -slog %%s  log_file for std DICE \n");
  fprintf(fout, "   -olog %%s  log_file for overall DICE \n");
  exit(exit_val);

}  /*  end usage()  */


static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    usage(0) ;
  else if (!stricmp(option, "mlog")
           )
    {
      mlog_fname = argv[2];
      nargs = 1;
      fprintf(stderr, "logging mean-dice to %s\n", mlog_fname) ;
    }
  else if (!stricmp(option, "log") ||
           !stricmp(option, "L")
           )
    {
      log_fname = argv[2];
      nargs = 1;
      fprintf(stderr, "logging individual dice to %s\n", log_fname) ;
    }
  else if (!stricmp(option, "slog")
           )
    {
      slog_fname = argv[2];
      nargs = 1;
      fprintf(stderr, "logging std-dice to %s\n", slog_fname) ;
    }
  else if (!stricmp(option, "olog")
           )
    {
      olog_fname = argv[2];
      nargs = 1;
      fprintf(stderr, "logging overall-dice to %s\n", olog_fname) ;
    }
  else{
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(0) ;
    exit(1) ;
  }

  return(nargs) ;
}


/*  EOF  */
