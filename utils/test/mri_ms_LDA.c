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


//
// mri_ms_LDA.c
//
// original author: Xiao Han
//
// Warning: Do not edit the following four lines.  CVS maintains them.
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "numerics.h"
#include "mri.h"
#include "matrix.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "mrinorm.h"
#define DEBUG 0
#define TEST 1
#define normflag 0


/* Ignore the off-diagonal elements of SW */
/* The synthesized image is more like synthW but poorer */
/* So, what's the use of it ?? */
/* But with pure LDA, the output image looks so noisy */
/*Maybe I should use the whole volume SW to compute LDA to suppress noise */
/* try it using mri_ms_LDA_whole volume */
#define CHOICE 0

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

static int USE_ONE = 0; /* use SW from WM only */
static int out_type = 3; /* MRI_FLOAT */
static double noise_threshold = 0.1; /* threshold for background noise */

MRI *MRInormalizeXH(MRI *mri_src, MRI *mri_dst, MRI *mri_mask);
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;

/* Compute LDA only in a local neighborhood of the
 * specified debug voxel
 */
static int debug_flag = 0;
static int window_flag = 0;
static int window_size = 30;

static float shift_value = -1;

static int ldaflag = 0;
static int conform = 0 ;
static int whole_volume = 0;
static int have_weight = 0; /* read in weights */
static int just_test = 0; /* if 1, try SW = Identity */
static int compute_m_distance = 0; /* if 1, compute distance
                                      in original space */

/* LDA is performed only within ROI */
static char *mask_fname = NULL; /* filename for ROI mask */
static char *label_fname = NULL; /* filename for segmentation */
static char *weight_fname = NULL; /* filename for storing LDA weights */
static char *synth_fname = NULL; /* filename for synthesized LDA volume */

static void usage_exit(int code) ;

void input_weights_to_file(float *weights, char *weight_fname, int VN);

void output_weights_to_file(float *weights, char *weight_fname, int VN);

void  update_LDAmeans(MRI **mri_flash,
                      MRI *mri_label,
                      MRI *mri_mask,
                      float *LDAmean1, float *LDAmean2,
                      int nvolumes_total, int classID1, int classID2);

void computeLDAweights(float *weights, MRI **mri_flash, MRI *mri_label,
                       MRI *mri_mask, float *LDAmean1, float *LDAmean2,
                       int nvolumes_total, int classID1, int classID2);

static int class1 = 0; /* to be used for LDA */
static int class2 = 0; /* to be used for LDA */

/* eps and lambda are used for covariance regularization */
static double eps = 1e-20;
static double lambda = 0.1;
static int regularize = 0;

#define MAX_IMAGES 200

int
main(int argc, char *argv[])
{
  char   **av, *in_fname;
  int    ac, nargs, i, j,  x, y, z, width, height, depth;
  MRI    *mri_flash[MAX_IMAGES], *mri_label, *mri_mask, *mri_tmp;
  int    msec, minutes, seconds, nvolumes, nvolumes_total ;
  Timer start ;
  float max_val, min_val, value;
  float *LDAmean1, *LDAmean2, *LDAweight;
  int label;
  double sum_white, sum_gray;
  int count_white, count_gray;

  nargs = handleVersionOption(argc, argv, "mri_ms_LDA");
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

  if (argc < 2)
    usage_exit(1) ;

  printf("command line parsing finished\n");

  if (have_weight == 0 && ldaflag == 0)
  {
    printf("Use -lda option to specify two class labels to optimize CNR on \n");
    usage_exit(0);
  }

  if (have_weight == 0 && label_fname == NULL)
  {
    printf("Use -label option to specify file for segmentation \n");
    usage_exit(0);
  }


  if (have_weight == 1 && weight_fname == NULL)
  {
    printf("Use -weight option to specify file for input LDA weights \n") ;
    usage_exit(0);
  }

  if (have_weight == 1 && synth_fname == NULL)
  {
    printf("Use -synth option to specify file for output synthesized volume \n") ;
    usage_exit(0);
  }

  //////////////////////////////////////////////////////////////////////////////////
  /*** Read in the input multi-echo volumes ***/
  nvolumes = 0 ;
  for (i = 1 ; i < argc; i++)
  {
    in_fname = argv[i] ;
    printf("reading %s...\n", in_fname) ;

    mri_flash[nvolumes] = MRIread(in_fname) ;
    if (mri_flash[nvolumes] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s",
                Progname, in_fname) ;
    /* conform will convert all data to UCHAR, which will reduce data resolution*/
    printf("%s read in. \n", in_fname) ;
    if (conform)
    {
      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      MRIfree(&mri_flash[nvolumes]);
      mri_flash[nvolumes] = mri_tmp ;
    }

    /* Change all volumes to float type for convenience */
    if (mri_flash[nvolumes]->type != MRI_FLOAT)
    {
      printf("Volume %d type is %d\n", nvolumes+1, mri_flash[nvolumes]->type);
      printf("Change data to float type \n");
      mri_tmp = MRIchangeType(mri_flash[nvolumes], MRI_FLOAT, 0, 1.0, 1);
      MRIfree(&mri_flash[nvolumes]);
      mri_flash[nvolumes] = mri_tmp; //swap
    }

    nvolumes++ ;
  }

  printf("All data read in\n");

  ///////////////////////////////////////////////////////////////////////////
  nvolumes_total = nvolumes ;   /* all volumes read in */

  for (i = 0 ; i < nvolumes ; i++)
  {
    for (j = i+1 ; j < nvolumes ; j++)
    {
      if ((mri_flash[i]->width != mri_flash[j]->width) ||
          (mri_flash[i]->height != mri_flash[j]->height) ||
          (mri_flash[i]->depth != mri_flash[j]->depth))
        ErrorExit(ERROR_BADPARM, "%s:\nvolumes %d (type %d) and %d (type %d) don't match (%d x %d x %d) vs (%d x %d x %d)\n",
                  Progname, i, mri_flash[i]->type, j, mri_flash[j]->type, mri_flash[i]->width,
                  mri_flash[i]->height, mri_flash[i]->depth,
                  mri_flash[j]->width, mri_flash[j]->height, mri_flash[j]->depth) ;
    }
  }

  width = mri_flash[0]->width;
  height = mri_flash[0]->height;
  depth = mri_flash[0]->depth;

  if (label_fname != NULL)
  {
    mri_label = MRIread(label_fname);
    if (!mri_label)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s\n",
                Progname, label_fname);

    if ((mri_label->width != mri_flash[0]->width) ||
        (mri_label->height != mri_flash[0]->height) ||
        (mri_label->depth != mri_flash[0]->depth))
      ErrorExit(ERROR_BADPARM, "%s: label volume size doesn't match data volumes\n", Progname);

    /* if(mri_label->type != MRI_UCHAR)
       ErrorExit(ERROR_BADPARM, "%s: label volume is not UCHAR type \n", Progname); */
  }

  if (mask_fname != NULL)
  {
    mri_mask = MRIread(mask_fname);
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s\n",
                Progname, mask_fname);

    if ((mri_mask->width != mri_flash[0]->width) ||
        (mri_mask->height != mri_flash[0]->height) ||
        (mri_mask->depth != mri_flash[0]->depth))
      ErrorExit(ERROR_BADPARM, "%s: mask volume size doesn't macth data volumes\n", Progname);

    if (mri_mask->type != MRI_UCHAR)
      ErrorExit(ERROR_BADPARM, "%s: mask volume is not UCHAR type \n", Progname);
  }
  else
  {

    if (have_weight == 1)
      noise_threshold = - 1e20;

    printf("Threshold input vol1 at %g to create mask \n", noise_threshold);
    printf("this threshold is useful to process skull-stripped data \n");

    mri_mask = MRIalloc(mri_flash[0]->width, mri_flash[0]->height, mri_flash[0]->depth, MRI_UCHAR);
    MRIcopyHeader(mri_flash[0], mri_mask);

    /* Simply set mask to be 1 everywhere */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {

          if ((float)MRIgetVoxVal(mri_flash[0], x, y,z,0) < noise_threshold)
            MRIvox(mri_mask, x, y,z) = 0;
          else
            MRIvox(mri_mask, x, y,z) = 1;
        }
  }

  /* Normalize input volumes */
  if (normflag)
  {
    printf("Normalize input volumes to zero mean, variance 1\n");
    for (i=0; i <nvolumes_total; i++)
    {
      mri_flash[i] = MRInormalizeXH(mri_flash[i], mri_flash[i], mri_mask);
    }
    printf("Normalization done.\n");
  }

  if (0)
  {
    printf("Using both hemi-sphere by changing rh-labels\n");
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          label = (int)MRIgetVoxVal(mri_label, x, y,z,0);
          if (label == 41) /* white matter */
            MRIsetVoxVal(mri_label, x, y, z, 0, 2);
          else if (label == 42) /* gm */
            MRIsetVoxVal(mri_label, x, y, z, 0, 3);
        }
  }


  if (debug_flag && window_flag)
  {
    /* Limit LDA to a local window */
    printf("Local window size = %d\n", window_size);
    window_size /= 2;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y,z) == 0) continue;

          if (z < (Gz - window_size) || z >(Gz + window_size)
              || y <(Gy - window_size) || y > (Gy + window_size)
              || x < (Gx - window_size) || x > (Gx + window_size))
            MRIvox(mri_mask, x, y,z) = 0;
        }

  }

  LDAmean1 = (float *)malloc(nvolumes_total*sizeof(float));
  LDAmean2 = (float *)malloc(nvolumes_total*sizeof(float));
  LDAweight = (float *)malloc(nvolumes_total*sizeof(float));

  if (have_weight)
  {
    printf("Read in LDA weights from weight-file\n");
    input_weights_to_file(LDAweight, weight_fname, nvolumes_total);
  }
  else
  { /* compute LDA weights */
    printf("Compute LDA weights to maximize CNR for region %d and region %d\n", class1, class2);
    /* Compute class means */
    update_LDAmeans(mri_flash, mri_label, mri_mask, LDAmean1, LDAmean2, nvolumes_total, class1, class2);
    printf("class means computed \n");

    /* Compute Fisher's LDA weights */
    computeLDAweights(LDAweight, mri_flash, mri_label, mri_mask, LDAmean1, LDAmean2, nvolumes_total, class1, class2);

    if (weight_fname != NULL)
    {
      output_weights_to_file(LDAweight, weight_fname, nvolumes_total);
    }
  }

  printf("LDA weights are: \n");
  for (i=0; i < nvolumes_total; i++)
  {
    printf("%g ", LDAweight[i]);
  }
  printf("\n");

  if (synth_fname != NULL)
  {
    /* linear projection of input volumes to a 1D volume */
    min_val = 10000.0;
    max_val = -10000.0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0) continue;

          value = 0.0;

          for (i=0; i < nvolumes_total; i++)
          {
            value += MRIFvox(mri_flash[i], x, y, z)*LDAweight[i];
          }

          //     if(value < 0) value = 0;

          if (max_val < value) max_val = value;
          if (min_val > value) min_val = value;

          /* Borrow mri_flash[0] to store the float values first */
          MRIFvox(mri_flash[0], x, y, z) = value;
        }

    printf("max_val = %g, min_val = %g \n", max_val, min_val);

    /* Check to make sure class1 has higher intensity than class2 */
    if (have_weight == 0)
    {
      sum_white =0;
      count_white = 0;
      sum_gray = 0;
      count_gray = 0;
      for (z=0; z < depth; z++)
      {
        if (count_white > 300 && count_gray > 300) break;
        for (y=0; y< height; y++)
        {
          for (x=0; x < width; x++)
          {

            if ((int)MRIgetVoxVal(mri_label, x, y,z,0) == class1)
            {
              sum_white += MRIFvox(mri_flash[0], x, y, z);
              count_white += 1;
            }
            else if ((int)MRIgetVoxVal(mri_label, x, y,z,0) == class2)
            {
              sum_gray += MRIFvox(mri_flash[0], x, y, z);
              count_gray += 1;
            }
          }
        }
      }

      if (count_white > 1 && count_gray > 1)
      {
        if (sum_white *count_gray < sum_gray*count_white)
        {
          for (z=0; z < depth; z++)
            for (y=0; y< height; y++)
              for (x=0; x < width; x++)
              {
                if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0) continue;
                value = MRIFvox(mri_flash[0], x, y, z);
                MRIFvox(mri_flash[0], x, y, z) = max_val - value;
              }
          max_val = max_val - min_val;
          min_val = 0;
        }
      }
    }

    /* The following is copied to be consistent with mri_synthesize */
    /* Don't know why add min_val, minus should make more sense */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0)
          {
            MRIFvox(mri_flash[0], x, y, z) = 0; /*background always set to 0 */
            continue;
          }
          /* Borrow mri_flash[0] to store the float values first */
          if (shift_value > 0)
          {
            value = MRIFvox(mri_flash[0], x, y, z) + shift_value;
            if (value < 0) value = 0;
            MRIFvox(mri_flash[0], x, y, z) = value;
          }
          else if (mask_fname != NULL)
            MRIFvox(mri_flash[0], x, y, z) -= min_val;
        }


    MRIfree(&mri_mask);
    if (mri_flash[0]->type == out_type)
    {
      mri_mask = MRIcopy(mri_flash[0], mri_mask);
    }
    else
    {
      mri_mask = MRIchangeType(mri_flash[0], out_type, 0.1, 0.99, 0);
    }

    /* Scale output to [0, 255] */
    if (0)
    {
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0) continue;

            value = (MRIFvox(mri_flash[0], x, y, z) - min_val)*255.0/(max_val - min_val) + 0.5; /* +0.5 for round-off */

            if (value > 255.0) value = 255.0;
            if (value < 0) value = 0;

            /* Borrow mri_flash[0] to store the float values first */
            MRIvox(mri_mask, x, y, z) = (BUFTYPE) value;
          }
    }

    /* Output synthesized volume */
    MRIwrite(mri_mask, synth_fname);

  }

  free(LDAmean1);
  free(LDAmean2);
  free(LDAweight);

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("LDA took %d minutes and %d seconds.\n", minutes, seconds) ;

  MRIfree(&mri_mask);

  if (label_fname)
    MRIfree(&mri_label);


  for (i=0; i < nvolumes_total; i++)
  {
    MRIfree(&mri_flash[i]);
  }
  exit(0);
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
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    debug_flag = 1;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "window"))
  {
    window_flag = 1 ;
    window_size = atoi(argv[2]) ;
    nargs = 1;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "shift"))
  {
    shift_value  = atof(argv[2]) ;
    nargs = 1;
    printf("shift output by %g before truncating at zero \n", shift_value) ;
  }
  else if (!stricmp(option, "out_type"))
  {
    out_type = atoi(argv[2]);
    nargs = 1;
    printf("Output type is %d\n", out_type);
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "use_one"))
  {
    USE_ONE = 1 ;
    printf("Using only the covariance matrix of first class as SW \n") ;
  }
  else if (!stricmp(option, "whole_volume"))
  {
    whole_volume = 1 ;
    printf("Synthesize background region too (if LDA)\n") ;
  }
  else if (!stricmp(option, "lda"))
  {
    ldaflag = 1;
    class1 = atoi(argv[2]) ;
    class2 = atoi(argv[3]) ;
    nargs = 2;
    printf("Using LDA method to generate synthesized volume (%d, %d) \n", class1, class2);
  }
  else if (!stricmp(option, "mask"))
  {
    mask_fname = argv[2];
    printf("using %s as mask for regions of interest \n", mask_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "label"))
  {
    label_fname = argv[2];
    printf("using %s as segmentation volume \n", label_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "synth"))
  {
    synth_fname = argv[2];
    printf("using %s as output for synthesized volume \n", synth_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "weight"))
  {
    weight_fname = argv[2];
    printf("using %s as input for LDA weights \n", weight_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
  }
  else if (!stricmp(option, "test"))
  {
    just_test = 1 ;
    printf("Test: set Sw to identity matrix.\n") ;
  }
  else if (!stricmp(option, "regularize"))
  {
    regularize = 1 ;
    lambda = atof(argv[2]);
    nargs = 1;
    printf("regularization for covarinace matrix, lambda = %g\n", lambda) ;
  }
  else if (!stricmp(option, "distance"))
  {
    compute_m_distance = 1 ;
    printf("Compute M distance between cluster centers.\n") ;
  }
  else switch (toupper(*option))
    {
    case 'T':
      noise_threshold = atof(argv[2]) ;
      printf("Threshold for background noise = %g\n", noise_threshold) ;
      printf("The threshold is applied on first input volume.\n") ;
      nargs = 1 ;
      break ;
    case 'W':
      have_weight = 1;
      break;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
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
  printf("usage: %s [options] <volume1> ... <volume N>\n", Progname) ;
  printf("This program takes an arbitrary # of FLASH images as input,\n"
         "and performs LDA dimension reduction on the multidimensional\n"
         "intensity space.\n");
  printf("Options includes:\n");
  printf("\t -lda %%d %%d to set the two class labels to optimize for\n");
  printf("\t -mask fname to set the brain mask volume\n");
  printf("\t -label fname to set the brain mask volume\n");
  printf("\t -weight fname for input LDA weights \n");
  printf("\t -shift # shift all values equal to -# to zero \n");
  printf("\t -synth fname for output synthesized volume \n");
  printf("\t -conform to conform input volumes (brain mask typically already conformed) \n");
  printf("\t -W to indicate weights are available from weight_fname\n");

  exit(code) ;

}


void  update_LDAmeans(MRI **mri_flash, MRI *mri_label, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2)
{

  int m, x, y, z, depth, height, width;
  float numer1, numer2; /* just to be consistent with LDA_all */
  double denom1, denom2;
  int label;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  for (m=0; m < nvolumes_total; m++)
  {
    numer1 = 0;
    denom1 = 0;
    numer2 = 0;
    denom2 = 0;

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) > 0)
          {
            label = (int) MRIgetVoxVal(mri_label, x, y, z,0);
            if (label == classID1)
            {
              numer1 += MRIFvox(mri_flash[m], x, y, z);
              denom1 += 1.0;
            }

            if (label == classID2)
            {
              numer2 += MRIFvox(mri_flash[m], x, y, z);
              denom2 += 1.0;
            }

          }

        }

    if (DEBUG)
    {
      printf("Class2 has size = %g, class 3 has size = %g\n", denom1, denom2);
    }
    LDAmean1[m] = numer1/(denom1 + 1e-15);
    LDAmean2[m] = numer2/(denom2 + 1e-15);

  }

  return;

}

void computeLDAweights(float *weights, MRI **mri_flash, MRI *mri_label, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2)
{
  /* To make it consistent with later CNR computation */

  int m1, m2, x, y, z, depth, height, width;
  double denom, denom1, denom2, sumw;
  float data1, data2;
  int label;
  double Mdistance;

  MATRIX *InvSW, *SW1, *SW2;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  SW1 = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  SW2 = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      SW1->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
      SW2->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
    }
  }

  /* printf("SW matrix initialized \n"); */
  denom1 = 0.0;
  denom2 = 0.0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;

        label = (int) MRIgetVoxVal(mri_label, x, y, z,0);

        if (label != classID1 &&
            label != classID2)
          continue;

        if (label == classID1)
        {
          denom1 += 1.0;
          for (m1=0; m1 < nvolumes_total; m1++)
          {
            data1 = MRIFvox(mri_flash[m1], x, y, z) - LDAmean1[m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              data2 = MRIFvox(mri_flash[m2], x, y, z) - LDAmean1[m2];
              SW1->rptr[m1+1][m2+1] += data1*data2;
            }
          }
        }
        else if (label == classID2)
        {
          denom2 += 1.0;
          for (m1=0; m1 < nvolumes_total; m1++)
          {
            data1 = MRIFvox(mri_flash[m1], x, y, z) - LDAmean2[m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              data2 = MRIFvox(mri_flash[m2], x, y, z) - LDAmean2[m2];
              SW2->rptr[m1+1][m2+1] += data1*data2;
            }
          }
        }

      } /* for all data points */

  if (denom1 <= 0.0 || denom2 <= 0)
    ErrorExit(ERROR_BADPARM, "%s: one or two classes is empty. \n", Progname);

  if (DEBUG)
    printf("brain size = %g\n", denom);

  if (USE_ONE)
  {
    printf("ONLY use SW from first class\n");
    printf("Seems reducing background noise\n");
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=m1; m2 <= nvolumes_total; m2++)
      {
        SW1->rptr[m1][m2] = SW1->rptr[m1][m2]/denom1;
        SW1->rptr[m2][m1] = SW1->rptr[m1][m2];
      }
    } /* for m1, m2 */

  }
  else
  {
    /* The following matches HBM2005 abstract's CNR definition */
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=m1; m2 <= nvolumes_total; m2++)
      {
        SW1->rptr[m1][m2] = SW1->rptr[m1][m2]/denom1 + SW2->rptr[m1][m2]/denom2;
        SW1->rptr[m2][m1] = SW1->rptr[m1][m2];
      }
    } /* for m1, m2 */

    if (regularize)
    {
      printf("regularization of the covariance estimate\n");
      for (m1=1; m1 <= nvolumes_total; m1++)
        SW1->rptr[m1][m1] += eps;  /* prevent SW1 to be singular */

      /* Borrow SW2 to store its inverse */
      SW2 = MatrixInverse(SW1, SW2);
      if (SW2 == NULL)
      {
        printf("Inverse matrix is NULL. Exit. \n");
        exit(1);
      }

      /* (1-lambda)* inv(SW + eps I) + labmda*I */
      for (m1=1; m1 <= nvolumes_total; m1++)
      {
        for (m2=m1; m2 <= nvolumes_total; m2++)
        {
          SW2->rptr[m1][m2] = (1.0 - lambda)*SW2->rptr[m1][m2];
          SW2->rptr[m2][m1] = SW2->rptr[m1][m2];
        }
        SW2->rptr[m1][m1] += lambda;
      }

      SW1 = MatrixInverse(SW2, SW1); // this inverse is quite redundant, since it will be inverted back again later

    }
  }

  if (0)
  {
    printf("SW is:\n");
    MatrixPrint(stdout, SW1);
  }

#if 0
  /* The following approach is equivalent to use -regularize; i.e., regularizing is equivalent to set SW to indentity */
  /* Compute inverse of SW */
  if (just_test == 0)
    InvSW = MatrixInverse(SW1, NULL);
  else
  {
    InvSW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=1; m2 <= nvolumes_total; m2++)
      {
        InvSW->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
      }
      InvSW->rptr[m1][m1] = 1.0;
    }
  }
#else
  /* Here, we try to ignore the covariance term */
  if (just_test)
  {
    for (m1=1; m1 < nvolumes_total; m1++)
    {
      for (m2=m1+1; m2 <= nvolumes_total; m2++)
      {
        SW1->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
        SW1->rptr[m2][m1] = 0.0; /* index starts from 1 for matrix */
      }
    }
  }

  InvSW = MatrixInverse(SW1, NULL);

#endif

  if (InvSW == NULL)
  { /* inverse doesn't exist */
    ErrorExit(ERROR_BADPARM, "%s: singular fuzzy covariance matrix.\n", Progname);
  }


  if (0)
  {
    printf("Inverse SW is:\n");
    MatrixPrint(stdout, InvSW);
  }

  if (0)
  {
    printf("Means for class 2 is \n");
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      printf("%g ", LDAmean1[m1-1]);
    }
    printf("\n");
    printf("Means for class 3 is \n");
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      printf("%g ", LDAmean2[m1-1]);
    }
  }
  /* Compute weights */
  denom = 0.0;
  sumw = 0.0;
  if (CHOICE == 1)
  {
    /* Do not use invSW, assume SW is diagonal */
    printf("Ignore off-diagonal of SW\n");
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      weights[m1-1]= (LDAmean1[m1-1] - LDAmean2[m1-1])/(SW1->rptr[m1][m1] + 1e-15);

      sumw += weights[m1-1];
      denom += weights[m1-1]*weights[m1-1];
    }
  }
  else
  {
    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      weights[m1-1]= 0.0;
      for (m2=1; m2 <= nvolumes_total; m2++)
      {
        weights[m1-1] += InvSW->rptr[m1][m2] *(LDAmean1[m2-1] - LDAmean2[m2-1]);
      }
      sumw += weights[m1-1];
      denom += weights[m1-1]*weights[m1-1];
    }
  }

  if (compute_m_distance)
  {
    Mdistance = 0;
    for (m1=0; m1 < nvolumes_total; m1++)
    {
      Mdistance += weights[m1]*(LDAmean1[m1] - LDAmean2[m1]);
    }
    printf("Mdistance = %g \n", Mdistance);
  }

  denom = sqrt(denom + 0.0000001);
  /* Normalized weights to have norm 1 */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    if (sumw > 0)
      weights[m1-1] /= denom;
    else
      weights[m1-1] /= -denom;
  }

  MatrixFree(&InvSW);

  MatrixFree(&SW1);
  MatrixFree(&SW2);

  return;
}

void input_weights_to_file(float *weights, char *weight_fname, int VN)
{
  int i;
  FILE *infile = NULL;

  infile = fopen(weight_fname, "r");

  if (!infile)
  {
    ErrorExit(ERROR_BADPARM, "%s: unable to open file %s for LDA weights .\n", Progname, weight_fname);
  }

  for (i=0; i < VN; i++)
  {
    if (feof(infile)!=0)
    {
      ErrorExit(ERROR_BADPARM, "%s: insufficient file length for LDA weights .\n", Progname);
    }
    fscanf(infile, "%f\n", &(weights[i]));
  }

  fclose(infile);

  return;
}

void output_weights_to_file(float *weights, char *weight_fname, int VN)
{
  int i;
  FILE *fp = NULL;

  fp = fopen(weight_fname, "w");
  if (!fp)
  {
    ErrorExit(ERROR_BADPARM, "%s: unable to open file %s for output LDA weights .\n", Progname, weight_fname);
  }

  for (i=0; i < VN; i++)
    fprintf(fp,"%g\n", weights[i]);

  fclose(fp);
  return;
}


MRI *MRInormalizeXH(MRI *mri_src, MRI *mri_dst, MRI *mri_mask)
{
  /* Normalize the source volume to be zero mean and variance 1*/
  /* mri_dst and mri_src can be the same */

  int width, height, depth, x, y, z;
  float mean, variance, total, tmpval;

  if (mri_src->type != MRI_FLOAT)
  {
    printf("Normalization is only applied for float-typed volume \n");
    mri_dst = MRIcopy(mri_src, mri_dst);
    return (mri_dst);
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  /* compute mean */
  mean = 0.0;
  total = 0.0;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        if (!mri_mask)
        {
          mean += MRIFvox(mri_src, x, y, z);
          total += 1;
        }
        else
        {
          if (MRIvox(mri_mask, x, y, z) >0)
          {
            mean += MRIFvox(mri_src, x, y, z);
            total += 1;
          }
        }
      }

  if (total > 0.0)
    mean = mean/total;

  /* compute variance */
  variance = 0.0;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        if (!mri_mask)
        {
          tmpval = MRIFvox(mri_src, x, y, z) - mean;
          variance += tmpval*tmpval;
        }
        else
        {
          if (MRIvox(mri_mask, x, y, z) >0)
          {
            tmpval = MRIFvox(mri_src, x, y, z) - mean;
            variance += tmpval*tmpval;
          }
        }
      }

  if (total > 0)
    variance = sqrt(variance/total);
  else
    variance = 1;

  /* normalization: invert variance first to save time */
  variance = 1.0/variance;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        tmpval = MRIFvox(mri_src, x, y, z) - mean;
        MRIFvox(mri_dst, x, y, z) = tmpval*variance;
      }

  return (mri_dst);
}
