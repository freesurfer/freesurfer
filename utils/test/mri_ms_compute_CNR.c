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
// mri_ms_compute_CNR.c
// Compute pair-wise CNR and Mahalanobis distances from meflash data
// Actually, optimal CNR and M-dist is the same!
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

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define MAX_IMAGES 200
#define DEBUG 0

static int xoff[6] =
  {
    1, -1, 0, 0, 0, 0
  };
static int yoff[6] =
  {
    0,  0, 1, -1, 0, 0
  };
static int zoff[6] =
  {
    0,  0, 0,  0, 1, -1
  };

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

/*
static int CNR_pairs = 34;
static int ilist[34] = {2, 17, 17, 17, 17, 17, 18, 18, 18, 18, 10, 10, 10, 10, 41, 53, 53, 53, 53, 53, 54, 54, 54, 54, 49, 49, 49, 49,  2,   3, 11, 41, 42, 50};
static int jlist[34] = {3, 2,  3,  18, 4,  5,  2,  3,  4,  5,  11, 4,  12, 2,  42, 41, 42, 54, 43, 44, 41, 42, 43, 44, 50, 43, 51, 41, 77,  77, 77, 77, 77, 77};
*/

/*
static int CNR_pairs = 32;

static int ilist[32] = {2, 10, 10, 10, 11, 11, 12, 12, 12, 13, 13, 17, 17, 17, 18, 18, 41, 49, 49, 49, 50, 50, 51, 51, 51, 52, 52, 53, 53, 53, 54, 54};
static int jlist[32] = {3,  2,  11, 4,  2,  4,  2, 3,  13,  2,  3, 18,  2,  3,  2,  3, 42, 41, 50, 43, 41, 43, 41, 42, 52, 41, 42, 54, 41, 42, 41, 42};
*/

static int CNR_pairs = 33;
static int ilist[33] =
  {
    219, 10, 10, 10, 11, 11, 12, 12, 12, 13, 13, 17, 17, 17, 18, 18, 219, 49, 49, 49, 50, 50, 51, 51, 51, 52, 52, 53, 53, 53, 54, 54, 2
  };
static int jlist[33] =
  {
    220,  219,  11, 4,  219,  4,  219, 220,  13,  219,  220, 18,  219,  220,  219,  220, 220, 219, 50, 43, 219, 43, 220, 220, 52, 219, 220, 54, 219, 220, 219, 220, 3
  };


const char *Progname ;

static int MINLABEL = 2;
static int MAXLABEL = 250;

/* Compute LDA only in a local neighborhood of the
 * specified debug voxel
 */
static int debug_flag = 0;
static int window_flag = 0;
static int window_size = 30;

static int conform = 0 ;
static int whole_volume = 0;
static int have_weight = 0; /* read in weights */
static int just_test = 0; /* if 1, try SW = Identity */
static int compute_m_distance = 0; /* if 1, compute distance in original space */

/* LDA is performed only within ROI */
static char *mask_fname = NULL; /* filename for ROI mask */
static char *label_fname = NULL; /* filename for segmentation */
static char *weight_fname = NULL; /* filename for storing LDA weights */
static char *synth_fname = NULL; /* filename for synthesized LDA volume */

static char *fname = NULL; /* filename to record CNR values */

static void usage_exit(int code) ;

void input_weights_to_file(float *weights, char *weight_fname, int VN);

void output_weights_to_file(float *weights, char *weight_fname, int VN);

void computeClassStats(float *LDAmean, MATRIX *SW, float *classsize, MRI **mri_flash, MRI *mri_label, MRI *mri_mask, int nvolumes_total, int classID);
MATRIX *ComputeAdjMatrix(MRI *mri_label, MRI *mri_mask, int minlabel, int maxlabel);

float computePairCNR(float *LDAmean1, float *LDAmean2, MATRIX *SW1, MATRIX *SW2, float classSize1, float classSize2, int nvolumes_total, float *weights, int flag);

static int class1 = 0;
static int class2 = 0;
static int ldaflag = 0;

int
main(int argc, char *argv[])
{
  char   **av, *in_fname;
  int    ac, nargs, i, j,  x, y, z, width, height, depth;
  MRI    *mri_flash[MAX_IMAGES], *mri_label, *mri_mask;
  int index;
  int    msec, minutes, seconds, nvolumes, nvolumes_total ;
  Timer start ;
  float max_val, min_val, value;
  float *LDAweight = NULL;
  float **LDAmeans = NULL; /* Centroid for each considered class */
  float *classSize =NULL; /* relative size of each class */
  MATRIX **SWs; /* Within class scatter-matrix for each considered class */
  MATRIX *AdjMatrix; /* Adjacency matrix of all classes */

  FILE *fp;

  int num_classes;
  double cnr;


  nargs = handleVersionOption(argc, argv, "mri_ms_compute_CNR");
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

  if (label_fname == NULL)
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


  if (ldaflag)
  {
    MINLABEL = MIN(class1, class2);
    MAXLABEL = MAX(class1, class2);
  }

  num_classes = MAXLABEL - MINLABEL + 1;

  printf("Total of %d classes considered in LDA training\n", num_classes);

  if (num_classes <= 1)
  {
    printf("Need to specify at least two classes to evaluate CNR\n");
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
      MRI *mri_tmp ;

      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      mri_flash[nvolumes] = mri_tmp ;
    }

    /* Change all volumes to float type for convenience */
    if (mri_flash[nvolumes]->type != MRI_FLOAT)
    {
      printf("Volume %d type is %d\n", nvolumes+1, mri_flash[nvolumes]->type);
      MRI *mri_tmp;
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
    mri_mask = MRIalloc(mri_flash[0]->width, mri_flash[0]->height, mri_flash[0]->depth, MRI_UCHAR);
    MRIcopyHeader(mri_flash[0], mri_mask);

    /* Simply set mask to be 1 everywhere */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          MRIvox(mri_mask, x, y,z) = 1;
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

  LDAweight = (float *)calloc(nvolumes_total, sizeof(float));

  /* Allocate memory */
  LDAmeans = (float **)malloc(num_classes*sizeof(float *));
  SWs = (MATRIX **)malloc(num_classes*sizeof(MATRIX *));
  classSize = (float *)malloc(num_classes*sizeof(float));
  for (i=0; i< num_classes; i++)
  {
    LDAmeans[i] = (float *)malloc(nvolumes_total*sizeof(float));
    SWs[i] = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
    if (SWs[i] == NULL || LDAmeans[i] == NULL)
      ErrorExit(ERROR_BADPARM, "%s: unable to allocate required memory \n", Progname);
  }

  if (ldaflag)
  {
    AdjMatrix = (MATRIX *)MatrixAlloc(num_classes, num_classes, MATRIX_REAL);

    /* The diagnoal entries of AdjMatrix is set to zero initially */
    for (i=1; i <= num_classes;i++)
      for (j=i; j <= num_classes; j++)
      {
        AdjMatrix->rptr[i][j] = 0.0;
        AdjMatrix->rptr[j][i] = 0.0;
      }

    AdjMatrix->rptr[class1-MINLABEL +1][class2-MINLABEL+1] = 1.0;
    AdjMatrix->rptr[class1-MINLABEL +1][class1-MINLABEL+1] = 1.0;
    AdjMatrix->rptr[class2-MINLABEL +1][class2-MINLABEL+1] = 1.0;
    AdjMatrix->rptr[class2-MINLABEL +1][class1-MINLABEL+1] = 1.0;
  }
  else if (MINLABEL <=2 && MAXLABEL >= 76)
  {
    printf("Manually set adjacent matrix \n");

    AdjMatrix = (MATRIX *)MatrixAlloc(num_classes, num_classes, MATRIX_REAL);

    /* The diagnoal entries of AdjMatrix is set to zero initially */
    for (i=1; i <= num_classes;i++)
      for (j=i; j <= num_classes; j++)
      {
        AdjMatrix->rptr[i][j] = 0.0;
        AdjMatrix->rptr[j][i] = 0.0;
      }

    for (index = 0; index < CNR_pairs; index++)
    {
      i = ilist[index] - MINLABEL;
      j = jlist[index] - MINLABEL;

      AdjMatrix->rptr[i+1][j+1] = 1.0;
      AdjMatrix->rptr[j+1][i+1] = 1.0;
    }

    /* left-hemisphere */
    /*
    AdjMatrix->rptr[2+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[17+1-MINLABEL][18+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[5+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[4+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[18+1-MINLABEL][2+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[18+1-MINLABEL][3+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[18+1-MINLABEL][5+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[18+1-MINLABEL][4+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][4+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][12+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][2+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][3+1-MINLABEL] = 1.0;
    */

    /* right-hemisphere */
    /*
    AdjMatrix->rptr[53+1-MINLABEL][41+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[53+1-MINLABEL][54+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[53+1-MINLABEL][42+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[53+1-MINLABEL][44+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[53+1-MINLABEL][43+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[54+1-MINLABEL][41+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[54+1-MINLABEL][42+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[54+1-MINLABEL][44+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[54+1-MINLABEL][43+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][50+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][43+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][51+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][41+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][42+1-MINLABEL] = 1.0;

    for(i=1; i < num_classes;i++)
      for(j=i+1; j <= num_classes; j++){
    if(AdjMatrix->rptr[i][j] > 0.5)
    AdjMatrix->rptr[j][i] = AdjMatrix->rptr[i][j];
    else
    AdjMatrix->rptr[i][j] = AdjMatrix->rptr[j][i];
      }
    */

    /*    AdjMatrix->rptr[2+1-MINLABEL][3+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][10+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][12+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][13+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][18+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][12+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][18+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[4+1-MINLABEL][10+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[4+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][13+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[12+1-MINLABEL][13+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[12+1-MINLABEL][26+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[17+1-MINLABEL][18+1-MINLABEL] = 1.0;
    */
    /* right-hemisphere */
    /* AdjMatrix->rptr[41+1-MINLABEL][42+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][49+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][50+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][51+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][52+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][53+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[41+1-MINLABEL][54+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[42+1-MINLABEL][51+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[42+1-MINLABEL][53+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[42+1-MINLABEL][54+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[43+1-MINLABEL][49+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[43+1-MINLABEL][50+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][50+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[49+1-MINLABEL][52+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[51+1-MINLABEL][52+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[51+1-MINLABEL][58+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[53+1-MINLABEL][54+1-MINLABEL] = 1.0;
    */
  }
  else
    AdjMatrix = ComputeAdjMatrix(mri_label, mri_mask, MINLABEL, MAXLABEL);
  /* AdjMatrix may need manual adjusted to avoid meaningless comparisons
   * such as computing CNR between left WM and right WM
   */


  for (i=1; i <= num_classes;i++)
    for (j=1; j <= num_classes; j++)
    {
      if (j==i) continue;
      /* the diagonal term will indicate whether the class is useful or not */
      AdjMatrix->rptr[i][i] +=  AdjMatrix->rptr[i][j] + AdjMatrix->rptr[j][i];
    }


  printf("Compute individual class statistics\n");
  /* Compute class means and covariance matrix */
  /* Note that here SWs will be covaraince matrix, not scatter matrix */
  for (i=0; i < num_classes; i++)
  {
    if (AdjMatrix->rptr[i+1][i+1] < 0.5) continue;
    computeClassStats(LDAmeans[i], SWs[i], &classSize[i], mri_flash, mri_label, mri_mask, nvolumes_total, MINLABEL + i);
  }
  printf("class statistics computed \n");

  if (fname != NULL)
    fp = fopen(fname, "w");
  else
    fp = 0;

  printf("compute pair-wise CNR/Mahalanobis distances \n");
  if (ldaflag)
  {
    for (i=0; i <num_classes-1;i++)
      for (j=i+1; j < num_classes; j++)
      {
        if (AdjMatrix->rptr[i+1][j+1] < 0.5) continue;
        cnr = computePairCNR(LDAmeans[i], LDAmeans[j], SWs[i], SWs[j], classSize[i],  classSize[j], nvolumes_total, LDAweight, 1);
        if (fp)
          fprintf(fp, "%9.4f ", (float)cnr);

        printf("CNR of class %d and class %d is %g\n", i+MINLABEL, j+MINLABEL, cnr);

      }

    if (fp)
    {
      fprintf(fp, "\n");
      fclose(fp);
    }

  }
  else
  {

    for (index = 0; index < CNR_pairs; index++)
    {
      i = ilist[index] - MINLABEL;
      j = jlist[index] - MINLABEL;

      if (AdjMatrix->rptr[i+1][j+1] < 0.5) continue;
      if (i== (2-MINLABEL) && j == (3-MINLABEL) && nvolumes_total > 1)
        cnr = computePairCNR(LDAmeans[i], LDAmeans[j], SWs[i], SWs[j], classSize[i],  classSize[j], nvolumes_total, LDAweight, 1);
      else
        cnr = computePairCNR(LDAmeans[i], LDAmeans[j], SWs[i], SWs[j], classSize[i], classSize[j], nvolumes_total, 0, 0);

      if (fp)
        fprintf(fp, "%9.4f ", (float)cnr);

      printf("CNR of class %d and class %d is %g\n", i+MINLABEL, j+MINLABEL, cnr);

    }

    if (fp)
    {
      fprintf(fp, "\n");
      fclose(fp);
    }
  }

  /* output weights for optimize CNR for class 2 and class 3 */
  if (weight_fname != NULL)
  {
    output_weights_to_file(LDAweight, weight_fname, nvolumes_total);
  }

  free(classSize);

  for (i=0; i< num_classes; i++)
  {
    free(LDAmeans[i]);
    MatrixFree(&SWs[i]);
  }
  free(LDAmeans);
  free(SWs);
  MatrixFree(&AdjMatrix);

  if (nvolumes_total > 1 && ((MINLABEL <=2 && MAXLABEL >= 3) || ldaflag))
  {
    if (ldaflag)
    {
      printf("LDA weights for %d-%d are: \n", class1, class2);
    }
    else
    {
      printf("LDA weights for 2-3 are: \n");
    }
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



            if (max_val < value) max_val = value;
            if (min_val > value) min_val = value;

            /* Borrow mri_flash[0] to store the float values first */
            MRIFvox(mri_flash[0], x, y, z) = value;
          }

      printf("max_val = %g, min_val = %g \n", max_val, min_val);

      /* Scale output to [0, 255] */
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

      /* Output synthesized volume */
      MRIwrite(mri_mask, synth_fname);
    }
  }

  free(LDAweight);

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("LDA took %d minutes and %d seconds.\n", minutes, seconds) ;

  MRIfree(&mri_mask);

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
    printf("compute CNR withing a window\n") ;
  }
  else if (!stricmp(option, "lda"))
  {
    ldaflag = 1;
    class1 = atoi(argv[2]) ;
    class2 = atoi(argv[3]) ;
    nargs = 2;
    printf("Evaluate CNR for two classes (%d, %d) \n", class1, class2);
  }
  else if (!stricmp(option, "min_label"))
  {
    MINLABEL = atoi(argv[2]);
    nargs = 1;
    printf("Label range starts at %d\n", MINLABEL);
  }
  else if (!stricmp(option, "max_label"))
  {
    MAXLABEL = atoi(argv[2]);
    nargs = 1;
    printf("Label range ends at %d\n", MAXLABEL);
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "whole_volume"))
  {
    whole_volume = 1 ;
    printf("Synthesize background region too (if LDA)\n") ;
  }
  else if (!stricmp(option, "mask"))
  {
    mask_fname = argv[2];
    printf("using %s as mask for regions of interest \n", mask_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "f")
           || !stricmp(option, "fname"))
  {
    fname = argv[2];
    nargs = 1;
    printf("Output CNR values to the file %s\n", fname);
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
  else if (!stricmp(option, "distance"))
  {
    compute_m_distance = 1 ;
    printf("Compute M distance between cluster centers.\n") ;
  }
  else switch (toupper(*option))
    {
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
         "and a label volume, then compute pair-wise CNR and M-distances\n");
  printf("Options includes:\n");
  printf("\t -mask fname to set the brain mask volume\n");
  printf("\t -label fname to set the brain mask volume\n");
  printf("\t -weight fname for input LDA weights \n");
  printf("\t -synth fname for output synthesized volume \n");
  printf("\t -conform to conform input volumes (brain mask typically already conformed) \n");
  printf("\t -W to indicate weights are available from weight_fname\n");
  printf("\t -f filename to output CNR values to the specified file\n");
  printf("\t -window %%d, used together with -debug_voxel allows computing CNR and synth within a local window with size %%d in each direction\n");
  printf("\t -debug_voxel (x,y,z) gives coordinates of a voxel\n");
  printf("\t -whole_volume: apply weights computed locally to whole volume \n");
  exit(code) ;

}

float computePairCNR(float *LDAmean1, float *LDAmean2, MATRIX *SW1, MATRIX *SW2, float classSize1, float classSize2, int nvolumes_total, float *weights, int flag)
{
  /* if flag == 1, output the optimal weights to weights */
  float cnr;
  int m1, m2;
  MATRIX *SW, *SB, *InvSW;
  float *weight;
  double value1, value2;

  if (classSize1 < 0.5 || classSize2 < 0.5)
  {
    printf("One of the two classes is empty\n");
    return 0;
  }

  SW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  SB = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
  weight = (float *)malloc(nvolumes_total*sizeof(float));

  /* initialize all */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] = 0.5*(SW1->rptr[m1][m2] + SW2->rptr[m1][m2]); /* index starts from 1 for matrix */
      SB->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
    }
    /* weight[m1-1] = 0.0; */ /*initialize later */
  }

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    value1 = LDAmean1[m1-1]  - LDAmean2[m1-1];
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      value2 = LDAmean1[m2-1]  - LDAmean2[m2-1];
      SB->rptr[m1][m2] += value1*value2;
      SB->rptr[m2][m1]  = SB->rptr[m1][m2];
    }
  }

  InvSW = MatrixInverse(SW, NULL);

  /* compute the optimal weight that will give largest CNR for the two
   * classes! Note that the weight differs depending on which formula
   * to use to evaluate final CNR
   * Here, the weight is consistent with the following definition:
   * CNR = |\mu_1 - \mu_2|\sqrt((\sigma_1^2 + \sigma_2^2)*0.5)
   */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    weight[m1-1]= 0.0;
    for (m2=1; m2 <= nvolumes_total; m2++)
    {
      weight[m1-1] += InvSW->rptr[m1][m2] *(LDAmean1[m2-1] - LDAmean2[m2-1]);
    }
  }

  /* compute CNR */
  cnr = 0 ;
  value1 = 0.0;
  value2 = 0.0; /* borrowed */
  for (m1=0; m1 < nvolumes_total; m1++)
  {
    cnr += weight[m1]*(LDAmean1[m1] - LDAmean2[m1]);
    if (flag)
    {
      value1 += weight[m1];
      value2 += weight[m1]*weight[m1];
    }
  }

  if (flag && weights != NULL)
  {
    value2 = sqrt(value2);
    for (m1=0; m1 < nvolumes_total; m1++)
    {
      if (value1 > 0)
        weights[m1] = weight[m1]/(value2 + 1e-15);
      else
        weights[m1] = -weight[m1]/(value2 + 1e-15);
    }
  }

  cnr = sqrt(cnr);

  MatrixFree(&SW);
  MatrixFree(&InvSW);
  MatrixFree(&SB);
  free(weight);

  return cnr;
}

void computeClassStats(float *LDAmean, MATRIX *SW, float *classsize, MRI **mri_flash, MRI *mri_label, MRI *mri_mask, int nvolumes_total, int classID)
{
  /* Compute LDAmean and SW from given data and label */

  int m1, m2, x, y, z, depth, height, width;
  int denom;
  double data1, data2;
  int label;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  if (!LDAmean)
    LDAmean = (float *)malloc(nvolumes_total*sizeof(float));

  if (!SW)
    SW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    LDAmean[m1-1] = 0.0;
    for (m2=1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
    }
  }

  for (m1=0; m1 < nvolumes_total; m1++)
    LDAmean[m1] = 0;

  /* Compute Class mean */
  denom = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;

        label = (int) MRIgetVoxVal(mri_label, x, y, z,0);

        if (label != classID)
          continue;

        denom += 1;
        for (m1=0; m1 < nvolumes_total; m1++)
          LDAmean[m1] += MRIFvox(mri_flash[m1], x, y, z);


      } /* for all data points */

  *classsize = denom;

  if (denom < 1) return;

  for (m1=0; m1 < nvolumes_total; m1++)
  {
    LDAmean[m1] /= (float) denom;
  } /* for m1, m2 */

  /* Compute scatter matrix; if divided by N (or N-1), will give covariance
     matrix */
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;
        label = (int) MRIgetVoxVal(mri_label, x, y, z,0);

        if (label != classID)
          continue;

        for (m1=0; m1 < nvolumes_total; m1++)
        {
          data1 = MRIFvox(mri_flash[m1], x, y, z) - LDAmean[m1];
          for (m2=m1; m2 < nvolumes_total; m2++)
          {
            data2 = MRIFvox(mri_flash[m2], x, y, z) - LDAmean[m2];
            SW->rptr[m1+1][m2+1] += data1*data2;
          }
        }
      } /* for all data points */

  /* Do not normalize */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] /= (float)denom;
      SW->rptr[m2][m1] = SW->rptr[m1][m2];
    }
  } /* for m1, m2 */

  return;
}

MATRIX *ComputeAdjMatrix(MRI *mri_label, MRI *mri_mask, int minlabel, int maxlabel)
{
  MATRIX *AdjMatrix;
  int i, j, label1, label2, offset, numLabels;
  int depth, width, height;
  int x,y,z,cx,cy,cz;

  numLabels = maxlabel - minlabel + 1;

  depth = mri_label->depth;
  width = mri_label->width;
  height = mri_label->height;

  AdjMatrix = (MATRIX *)MatrixAlloc(numLabels, numLabels, MATRIX_REAL);

  if (!AdjMatrix)
    ErrorExit(ERROR_BADPARM, "%s: unable to allowcate memory.\n", Progname);
  /* The diagnoal entries of AdjMatrix is set to zero and remain zero */
  for (i=1; i <= numLabels;i++)
    for (j=i; j <= numLabels; j++)
    {
      AdjMatrix->rptr[i][j] = 0.0;
      AdjMatrix->rptr[j][i] = 0.0;
    }


  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;
        label1 = (int) MRIgetVoxVal(mri_label, x, y, z,0);
        if (label1 < minlabel || label1 > maxlabel) continue;

        /* Find all 6-neighbor with different label */
        for (offset = 0; offset < 6; offset++)
        {
          cx = x + xoff[offset];
          cy = y + yoff[offset];
          cz = z + zoff[offset];
          if (cx < 0 || cx >= width || cy < 0 || cy >= height
              || cz < 0 || cz >= depth) continue;

          label2 = (int) MRIgetVoxVal(mri_label, cx, cy, cz,0);
          if (label2 < minlabel || label2 > maxlabel || label2 == label1)
            continue;

          AdjMatrix->rptr[label1-minlabel+1][label2-minlabel+1] = 1.0;
          AdjMatrix->rptr[label2-minlabel+1][label1-minlabel+1] = 1.0;
        } /* for_offset */
      }

  return (AdjMatrix);
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
