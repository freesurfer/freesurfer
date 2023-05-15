/**
 * @brief Compute Dice coefficent comparing two segmentation volumes.
 *
 * This program compares two segmentation volumes and
 * computes the Dice and Jaccard Coefficients.
 */
/*
 * Original Authors: Xiao Han, Nick Schmansky
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
#include "colortab.h"



static int do_cortex = 1 ;
static int do_wm = 1 ;

static void usage(int exit_val);
static int  get_option(int argc, char *argv[]) ;

const char *Progname;

static char *log_fname = NULL ; //dice coeff of individual structures
static char *mlog_fname = NULL ; //mean of individual dice
static char *slog_fname = NULL ; //std of individual dice
static char *olog_fname = NULL ; //overall dice for subcortical structures

static const int num_labels       = 24;
static const int labels_of_interest[24] =
{
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
static const int labels_overall_Dice[18] =
{
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
static int isOverallDiceLabel(int volVal)
{
  int i;
  for (i=0; i < num_labels_overall_Dice; i++)
  {
    if (volVal == labels_overall_Dice[i])
    {
      return 1;
    }
  }
  return 0;
}

/* maximum number of classes: 
   large enough to handle the maximum label value
   in FreeSurferColorLUT.txt
*/
#define MAX_CLASSES 15000
#define MAX_CLASS_NUM 14999

int all_labels_flag = FALSE;
int num_all_labels = 0;
int all_labels_of_interest[MAX_CLASSES];
COLOR_TABLE *ctab=NULL;
char *table_fname;
FILE *tablefp;

int main(int argc, char *argv[])
{
  MRI *mri_seg1, *mri_seg2;
  int nargs, ac;
  char **av;
  int width, height, depth, x, y, z, f, nframes;
  int v1, v2;
  int i, skipped;
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

  nargs = handleVersionOption(argc, argv, "mri_compute_seg_overlap");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 3)
  {
    usage(1);
  }

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

  for (i=0; i < MAX_CLASSES; i++)
  {
    Volume_union[i] = 0;
    Volume_overlap[i] = 0;
    Volume_from1[i] = 0;
    Volume_from2[i] = 0;
  }

  width = mri_seg1->width ;
  height = mri_seg1->height ;
  depth = mri_seg1->depth ;
  nframes = mri_seg1->nframes ;
  if (nframes == 0)
  {
    nframes = 1;
  }

  subcorvolume_overlap = 0;
  subcorvolume1 = 0;
  subcorvolume2 = 0;

  for (f = 0 ; f < nframes ; f++)  {
    for (z = 0 ; z < depth ; z++)    {
      for (y = 0 ; y < height ; y++)      {
        for (x = 0 ; x < width ; x++)        {
          v1 = (int) MRIgetVoxVal(mri_seg1,x,y,z,f);
          v2 = (int) MRIgetVoxVal(mri_seg2,x,y,z,f);

          if (v1 > MAX_CLASS_NUM || v1 < 0 || v2 > MAX_CLASS_NUM || v2 < 0) continue;

          /* do not include these in the overall Dice coefficient calculations:
             Left/Right-Cerebral-White-Matter (labels 2 and 41),
             Left/Right-Cerebral-Cortex (labels 3 and 42),
             Left/Right-Accumbens-area (labels 26 and 58)
             Notice that these labels are not included in the 'if' checks: */

          if (v1 == v2){
            if (all_labels_flag)             subcorvolume_overlap++;
            else if (isOverallDiceLabel(v1)) subcorvolume_overlap++;
          }

          if(all_labels_flag)              subcorvolume1++;
          else if (isOverallDiceLabel(v1)) subcorvolume1++;
          if (all_labels_flag)             subcorvolume2++;
          else if (isOverallDiceLabel(v2)) subcorvolume2++;

          Volume_from1[v1]++;
          Volume_from2[v2]++;

          if (v1 == v2) {
            Volume_overlap[v1]++;
            Volume_union[v1]++;
          }
          else {
            Volume_union[v1]++;
            Volume_union[v2]++;
          }
        }
      }
    }
  }

  for (i=0; i < MAX_CLASSES; i++)
  {
    correct_ratio[i] =
      (double)Volume_overlap[i]/((double)Volume_union[i] + 1e-10);
    correct_ratio2[i] =
      (double)Volume_overlap[i]*2.0/
      ((double)Volume_from1[i] + (double)Volume_from2[i] + 1e-10);
  }

  printf("Jaccard Coefficients:\n");
  mean1 = 0;
  std1 = 0;
  if (all_labels_flag)
  {
    num_all_labels = 0;
    for (i=0; i < MAX_CLASSES; i++)
    {
      if(correct_ratio[i] > 0.0) // This will include zero overlap areas as well (not just non-existing labels. If it is a problem, should flag existing labels....
      {
        printf("correct ratio for label %d = %g\n",i, correct_ratio[i]);
        mean1 += correct_ratio[i];
        std1 += correct_ratio[i] * correct_ratio[i];
        num_all_labels++;
      }
    }
    mean1 /= num_all_labels;
    std1 /= num_all_labels;
    std1 = sqrt(std1 - mean1*mean1);
    printf("mean +/- std = %6.4f +/- %6.4f\n", mean1, std1);
  }
  else
  {
    for (skipped = i=0; i < num_labels; i++)
    {
      if (do_cortex == 0 && IS_CORTEX(labels_of_interest[i]))
      {
	skipped++ ;
	continue ;
      }
      if (do_wm == 0 && IS_WHITE_CLASS(labels_of_interest[i]))
      {
	skipped++ ;
	continue ;
      }
      
      printf("correct ratio for label %d = %g\n",
             labels_of_interest[i], correct_ratio[labels_of_interest[i]]);
      mean1 += correct_ratio[labels_of_interest[i]];
      std1 += correct_ratio[labels_of_interest[i]] *
	correct_ratio[labels_of_interest[i]];
    }

    mean1 /= (num_labels-skipped);
    std1 /= (num_labels-skipped);
    std1 = sqrt(std1 - mean1*mean1);
    printf("mean +/- std = %6.4f +/- %6.4f\n", mean1, std1);
  }

  if(table_fname) tablefp = fopen(table_fname, "w") ;

  printf("Dice Coefficients:\n");
  // printf("ratio of overlap to volume of input1:\n");
  mean2 = 0;
  std2 = 0;
  if (all_labels_flag)
  {
    num_all_labels = 0;
    for (i=0; i < MAX_CLASSES; i++)
    {
      if(correct_ratio2[i] > 0.0) // This will include zero overlap areas as well (not just non-existing labels. If it is a problem, should flag existing labels....
      {
        all_labels_of_interest[num_all_labels] = i;
	if(ctab == NULL) printf("label %d = %g\n",i, correct_ratio2[i]);
	else             printf("%4d %s %8.6lf\n",i, ctab->entries[i]->name,correct_ratio2[i]);
        mean2 += correct_ratio2[i];
        std2 += correct_ratio2[i]*correct_ratio2[i];
        num_all_labels++;
      }
    }
    mean2 /= num_all_labels;
    std2 /= num_all_labels;
    std2 = sqrt(std2 - mean2*mean2);
    printf("mean +/- std = %6.4f +/- %6.4f \n", mean2, std2);
  }
  else  {
    for (skipped=i=0; i < num_labels; i++)    {
      if (do_cortex == 0 && IS_CORTEX(labels_of_interest[i]))      {
	skipped++ ;
	continue ;
      }
      if (do_wm == 0 && IS_WHITE_CLASS(labels_of_interest[i]))      {
	skipped++ ;
	continue ;
      }
      int j = labels_of_interest[i];
      double voldiff;
      voldiff = 100*(Volume_from1[j]-Volume_from2[j])/(0.5*((Volume_from1[j]+Volume_from2[j])));
      if(ctab == NULL) printf("label %d = %g\n",j, correct_ratio2[j]);
      else             printf("%4d  %-30s %8.6lf %8.6lf %6d %6d %6d %9.4lf\n",
			      j,ctab->entries[j]->name,correct_ratio[j],correct_ratio2[j],
			      Volume_from1[j],Volume_from2[j],Volume_overlap[j],voldiff);
      if(table_fname) fprintf(tablefp,"%4d  %-30s %8.6lf %8.6lf %6d %6d %6d %9.4lf\n",
			      j,ctab->entries[j]->name,correct_ratio[j],correct_ratio2[j],
			      Volume_from1[j],Volume_from2[j],Volume_overlap[j],voldiff);


      mean2 += correct_ratio2[labels_of_interest[i]];
      std2 += correct_ratio2[labels_of_interest[i]] * correct_ratio2[labels_of_interest[i]];
    }
    mean2 /= (num_labels-skipped);
    std2 /= (num_labels-skipped);
    std2 = sqrt(std2 - mean2*mean2);
    printf("mean +/- std = %6.4f +/- %6.4f \n", mean2, std2);
  }
  if(table_fname) fclose(tablefp);

  if (log_fname != NULL)
  {
    log_fp = fopen(log_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, log_fname) ;
    if (all_labels_flag)
    {
      for (i=0; i < num_all_labels; i++)
      {
        fprintf(log_fp, "%6.4f ", correct_ratio2[all_labels_of_interest[i]]);
      }
    }
    else
    {
      for (i=0; i < num_labels; i++)
      {
	if (do_cortex == 0 && IS_CORTEX(labels_of_interest[i]))
	  continue ;
	if (do_wm == 0 && IS_WHITE_CLASS(labels_of_interest[i]))
	  continue ;

        fprintf(log_fp, "%6.4f ", correct_ratio2[labels_of_interest[i]]);
      }
    }
    fprintf(log_fp, "%6.4f ", mean2);
    fprintf(log_fp, "%6.4f ", std2);
    fprintf(log_fp, "%6.4f \n",
            subcorvolume_overlap*2.0/(float)(subcorvolume1 + subcorvolume2));
    fclose(log_fp);
  }

  if (mlog_fname != NULL)
  {
    log_fp = fopen(mlog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, mlog_fname) ;
    fprintf(log_fp, "%6.4f \n", mean2);
    fclose(log_fp);
  }

  if (slog_fname != NULL)
  {
    log_fp = fopen(slog_fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, slog_fname) ;
    fprintf(log_fp, "%6.4f \n", std2);
    fclose(log_fp);
  }

  if (olog_fname != NULL)
  {
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


#include "mri_compute_seg_overlap.help.xml.h"
static void usage(int exit_val)
{
  outputHelpXml(mri_compute_seg_overlap_help_xml,
                mri_compute_seg_overlap_help_xml_len);
  printf("In the table, the rows are:\n");
  printf("1. Segmentation Index\n");
  printf("2. Segmentation Name (if avail)\n");
  printf("3. Not sure\n");
  printf("4. Dice\n");
  printf("5. nvoxels in vol 1\n");
  printf("6. nvoxels in vol 2\n");
  printf("7. nvoxels in intersection\n");
  printf("8. Percentage volume diff\n");

  exit(exit_val);
}  /*  end usage()  */


static int get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
  {
    usage(0) ;
  }
  else if (!stricmp(option, "-all_labels"))
  {
    all_labels_flag = TRUE;
    fprintf(stderr, "Computing overlap measures for all commonly existing labels. \n") ;
  }
  else if (!stricmp(option, "mlog"))
  {
    mlog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging mean Dice to %s\n", mlog_fname) ;
  }
  else if (!stricmp(option, "log") || !stricmp(option, "L"))
  {
    log_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging individual Dice to %s\n", log_fname) ;
  }
  else if (!stricmp(option, "default-ctab")){
    char tmpstr[2000];
    sprintf(tmpstr,"%s/FreeSurferColorLUT.txt",getenv("FREESURFER_HOME"));
    ctab = CTABreadASCII(tmpstr);
    printf("Using ctab %s\n",tmpstr);
  }
  else if (!stricmp(option, "table")){
    char tmpstr[2000];
    table_fname = argv[2];
    nargs = 1;
    sprintf(tmpstr,"%s/FreeSurferColorLUT.txt",getenv("FREESURFER_HOME"));
    ctab = CTABreadASCII(tmpstr);
    printf("Using ctab %s\n",tmpstr);
  }
  else if (!stricmp(option, "dice")){
    // Stand-alone method to compute dice. Has  more flexibility
    // 7: seg1 seg2 ctab ReportEmpty01 ExcludeId datfile tablefile
    // seg1 - first segmentation (for fdr and tdr, this is ground truth)
    // seg2 - second segmentation (for fdr and tdr, this is the test seg)
    // ctab - color table of segmentations to report on. can used "embedded"
    // ReportEmpty - 0=do not report segs that are empty in both, 1=report 
    // ExcludeId - exclude this seg (eg, 0 to exclude Unknown)
    // datfile - save the dice for each seg on a single line without anymore info
    // tablefile - save as a table with each seg on a row followed by 
    //   count1 count2 dice tpr fdr 
    SegDice sd;
    sd.seg1 = MRIread(argv[2]);
    if(sd.seg1==NULL) exit(1);
    sd.seg2 = MRIread(argv[3]);
    if(sd.seg2==NULL) exit(1);
    if(strcmp(argv[4],"embedded")!=0){
      sd.ctab = CTABreadASCII(argv[4]);
      if(sd.ctab==NULL) exit(1);
    } 
    else {
      printf("Using embedded color table\n");
      if(sd.seg1->ct != NULL) sd.ctab = sd.seg1->ct;
      else if(sd.seg2->ct != NULL) sd.ctab = sd.seg2->ct;
      else {
	printf("ERROR: neither seg1 nor seg2 have an embedded color table\n");
	exit(1);
      }
    }
    sscanf(argv[5],"%d",&sd.ReportEmpty);
    int ExcludeId;
    sscanf(argv[6],"%d",&ExcludeId);
    sd.excludelist.push_back(ExcludeId);
    int err = sd.ComputeDice();
    if(err) exit(err);
    sd.WriteDiceDat(argv[7]);
    sd.WriteDiceTable(argv[8]);
    exit(0);
  }
  else if (!stricmp(option, "slog"))
  {
    slog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging std Dice to %s\n", slog_fname) ;
  }
  else if (!stricmp(option, "olog"))
  {
    olog_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging overall Dice to %s\n", olog_fname) ;
  }
  else if (!stricmp(option, "wm"))
  {
    do_wm = atoi(argv[2]);
    nargs = 1;
    fprintf(stderr, "%sincluding cerebral white matter\n", do_wm ? "" : "NOT ") ;
  }
  else if (!stricmp(option, "cortex"))
  {
    do_cortex = atoi(argv[2]);
    nargs = 1;
    fprintf(stderr, "%sincluding cerebral cortex\n", do_cortex ? "" : "NOT ") ;
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(0) ;
    exit(1) ;
  }

  return(nargs) ;
}

