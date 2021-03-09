/**
 * @brief performs EM-segmentation on the multidimensional intensity space.
 *
 * This program takes an arbitrary # of FLASH images as input,
 * and performs EM-segmentation on the multidimensional intensity space.
 *
 * intended for gca guided
 * 1.Find xn, yn, zn for each voxel once and store them,
 * instead of have to constantly do the transform
 *
 * Things to do:
 * 1. Need to postprocessing with fixed class stats but ignoring
 * atlas prior, since atlas is really way-off at some places (or
 * should I ignore voxel-wise prior from the beginning, just fix
 * class-prior from atlas??
 * 2. May need MRF prior to make sure segmentation is smooth
 * 3. PVE model may be needed; better try on mri_ms_EM first
 */
/*
 * Original Author: Xiao Han
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
#include "gca.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "mrinorm.h"

#define TEST 1
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

/* maximum number of classes */
#define MAX_CLASSES 20

/* If set at 0.7, then Mahalanobis distance is 26.x; if set to 0.5, it's only 14.6. Of cause, lower threshold includes more partial volumed voxels, and increasecovariance */
/* For fair comparison, need to generate a hard segmentation! */

#define MEMTHRESHOLD 0.7

static int debug_flag = 0;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;

static int normflag = 0; /* normalize input volume */

/* Values lower than threshold are considered as background, and not counted
 * in clustering
 */
static double noise_threshold = 1.0; /* threshold for background noise */

static double tolerance = 0.01; /* Convergence criterion */

// static int num_classes = 3; /* default number of classes */

static int max_iters = 100; /* defualt maximal iteration number */

static int hard_segmentation = 0;

static char *label_fname = NULL; /* filename for initial segmentation */

static char *gca_fname = NULL; /* filename for GCA */

static char *xform_fname = NULL; /* filename for xform to GCA */

static int conform = 0 ;

static int mask_subcortical = 0;

static int HARD_ONLY = 0;

static int rescale = 0;

static int label_of_interest = 0; /* the label to be refined */

/* Clustering is performed only within ROI */
static char *mask_fname = NULL; /* filename for ROI mask */

static void usage_exit(int code) ;

MRI *MRInormalizeXH(MRI *mri_src, MRI *mri_dst, MRI *mri_mask);

void update_centroids(MRI **mri_flash, MRI **mri_mem, MRI **mri_likelihood, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes);

void update_F(MATRIX **F, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes);

void compute_detF(MATRIX **F, double *detF, int num_classes);

void compute_inverseF(MATRIX **F, int num_classes);

double distancems(MRI **mri_flash, float *centroids, MATRIX *F, double detF, int x, int y, int z);

static int InterpMethod = SAMPLE_TRILINEAR;  /*E* prev default behavior */
static int sinchalfwindow = 3;

static int xoff[6] =
  {
    1, 0, -1, 0, 0, 0
  };
static int yoff[6] =
  {
    0, 1, 0, -1, 0, 0
  };
static int zoff[6] =
  {
    0, 0, 0, 0, 1, -1
  };

static int fix_class_size = 1;

//kappa only used in commented code, treating warnings as errors
//so removing the definition here
//static double kappa = 0.001; /* 0.01 is too big */

/* eps and lambda are used for covariance regularization */
static double eps = 1e-30;
static double lambda = 0.2;
static int regularize = 0;

#define MAX_IMAGES 200

int
main(int argc, char *argv[])
{
  char   **av, *in_fname, fname[100];
  int    ac, nargs, i, j,  x, y, z, width, height, depth, iter, c;
  int num_classes = 0;
  MRI    *mri_flash[MAX_IMAGES], *mri_mem[MAX_CLASSES], *mri_mask, *mri_label, *mri_tmp;
  char   *out_prefx ;
  int    msec, minutes, seconds, nvolumes, nvolumes_total, label ;
  Timer start ;
  MATRIX *F[MAX_CLASSES];
  float *centroids[MAX_CLASSES];
  float value;
  float max_change = 1000.0;
  double detF[MAX_CLASSES]; /* determinant of the covariance matrix */
  double distance2, sum_of_distance; /* actually the inverse distance */
  double oldMems[MAX_CLASSES];
  float classSize[MAX_CLASSES]; // wouldn't be used here
  int brainsize;
  GCA *gca;
  GCA_PRIOR *gcap;
  TRANSFORM *transform;
  int xn, yn, zn, n;
  int label_to_index[256];
  int index_to_label[256];
  int total_label;
  int MLE_label;
  double max_prior;

  nargs = handleVersionOption(argc, argv, "mri_ms_gca_EM");
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

  out_prefx = argv[argc-1] ;

  printf("command line parsing finished\n");

  if (label_of_interest == 0)
  {
    printf("Use -label option to specify the structure to be refined. Exit.\n");
    exit(0);
  }


  if (gca_fname == NULL)
  {
    printf("Use -gca option to specify GCA filename. Exit. \n");
    exit(1);
  }

  if (xform_fname == NULL)
  {
    printf("use -xform option to specify gca-transform filename. Exit. \n");
    exit(1);
  }

  //////////////////////////////////////////////////////////////////////////////////
  fprintf(stderr, "reading gca from %s...\n", gca_fname);
  gca = GCAread(gca_fname);
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s:could not read classifier array from %s", Progname, gca_fname);

  if (gca->flags & GCA_NO_MRF)
    ErrorExit(ERROR_BADPARM, "%s: gca %s built without Markov priors", Progname, gca_fname);


  /*** Read in the input multi-echo volumes ***/
  nvolumes = 0 ;
  for (i = 1 ; i < argc-1 ; i++)
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

  fprintf(stderr, "reading transform from %s ...\n", xform_fname);
  transform = TransformRead(xform_fname);
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s:could not open transform %s", Progname, xform_fname);
  printf("Transform read in, now computing its inverse...\n");
  TransformInvert(transform, mri_flash[0]);
  printf("Transform for GCA is inverted \n");

  printf("Need to read in initial segmentation from aseg result\n");

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


  if (label_fname == NULL)
  {
    printf("use GCA prior to generate initial segmentation\n");
    mri_label = MRIalloc(mri_flash[0]->width, mri_flash[0]->height, mri_flash[0]->depth, MRI_SHORT);
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          MRIsetVoxVal(mri_label, x, y, z, 0, 0);

          if (!GCAsourceVoxelToNode(gca, mri_flash[0], transform, x, y, z, &xn, &yn, &zn))
          {
            gcap = getGCAP(gca, mri_flash[0], transform, x, y, z) ;

            if (gcap == NULL)
            {
              printf("Failed to get prior for voxel (%d,%d,%d)\n",x,y,z);
              continue;
            }

            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("gcap->nlabels = %d\n", gcap->nlabels);
              printf("(xn,yn,zn) = (%d, %d, %d)\n", xn, yn, zn);
            }
            if (gcap->nlabels < 1) continue;
            max_prior = gcap->priors[0];
            MLE_label = gcap->labels[0];
            // going through gcap labels
            for (n = 1 ; n < gcap->nlabels ; n++)
            {
              if (debug_flag && x == Gx && y == Gy && z == Gz)
              {
                printf("gcap->labels[%d]=%d\n", n, gcap->labels[n]);
              }
              if (gcap->priors[n] > max_prior)
              {
                max_prior = gcap->priors[n];
                MLE_label = gcap->labels[n];
              }
            }

            MRIsetVoxVal(mri_label, x, y, z, 0, MLE_label);
          }


        }

    MRIwrite(mri_label, "testGCA.mgz");

    // MRIfree(&mri_label);

    // for(i=0; i < nvolumes_total; i++){
    //  MRIfree(&mri_flash[i]);
    //}
    //  exit(0);

  }

  if (label_fname != NULL)
  {
    printf("Read in initial segmentation...\n");
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

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          MRIvox(mri_mask, x, y,z) = 1;
        }
#if 0
    if (label_fname != NULL)
    {
      mri_tmp = MRIclone(mri_mask, NULL);
      printf("Use label volume to define regions to be segmented. \n");
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            label = (int)MRIgetVoxVal(mri_label, x, y,z,0);
            /* if((label >=9 && label <= 12) || (label >= 48 && label <= 51))
            MRIvox(mri_mask, x, y,z) = 0; */
            if (label <= 0)
              MRIvox(mri_tmp, x, y,z) = 0;
            else
              MRIvox(mri_tmp, x, y,z) = 1;
          }
      /* dilate the region with label by 2 */
      mri_mask = MRIdilateLabelUchar(mri_tmp, mri_mask, 1, 2);
      MRIfree(&mri_tmp);
    }
#endif
    /* Compute the mask from the first input volume */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIFvox(mri_flash[0], x, y, z) < (float) noise_threshold)
            MRIvox(mri_mask, x, y,z) = 0;
        }

  }

  /* Now find all relevant labels, i.e., structures adjacent to structure of interest */
  for (i=0; i<256;i++)
    label_to_index[i] = -1;

  label_to_index[label_of_interest] = 0;
  index_to_label[0] = label_of_interest;
  total_label = 1;

  for (z=1; z < (depth-1); z++)
    for (y=1; y< (height-1); y++)
      for (x=1; x < (width-1); x++)
      {
        if ((int)MRIgetVoxVal(mri_label, x, y, z,0) != label_of_interest)
          continue;

        for (i=0; i < 6; i++)
        {
          xn = x + xoff[6];
          yn = y + yoff[6];
          zn = z + zoff[6];
          label = (int)MRIgetVoxVal(mri_label, xn, yn, zn,0);
          if (label <= 0) continue;
          if (label_of_interest == 13 && label == 3) continue;
          if (label_to_index[label] < 0)
          {
            label_to_index[label] = total_label;
            index_to_label[total_label] = label;
            total_label++;
          }
        }
      }

  num_classes = total_label;
  printf("total of classes involved = %d, they are: \n", num_classes);

  for (c=0; c < num_classes; c++)
    printf("%d ", index_to_label[c]);
  printf("\n");

  //  printf("find the region of interest, so as to reduce volume sizes\n");
  printf("Be careful in later use of GCA if truncating input volumes to smaller size here\n");

  brainsize = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;
        label = (int)MRIgetVoxVal(mri_label, x, y, z,0);
        if (label_to_index[label] < 0) //an irrelevant label
          MRIvox(mri_mask,x, y, z) = 0;
        else
          brainsize++;
      }

  printf("Brain size = %d\n", brainsize);
  if (brainsize < 1)
  {
    printf("region to be segment has zero size. exit.\n");
    exit(1);
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

  /* Allocate memory for membership volumes */
  for (i=0; i < num_classes; i++)
  {
    mri_mem[i] = MRIclone(mri_flash[0], NULL);
    centroids[i] = (float *)malloc(nvolumes_total*sizeof(float));
    F[i] = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
    //    classSize[i] = 1.0/(float)num_classes;
  }


  /* Initialize class means and covariance matrices */
  /* actually only need to initialize membership function */
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;
        label = (int)MRIgetVoxVal(mri_label, x, y, z,0);
        i = label_to_index[label];

        MRIFvox(mri_mem[i],x,y,z)  = 1;
      }

  /* Use the membership function obtained above to get initial centroids
   * for the multi-dimensional clustering
   */

  /* Perform iterative fuzzy-clustering until convergence */
  iter = 0;
  max_change = 10000.0;
  while (max_change > tolerance && iter < max_iters)
  {
    iter++;
    max_change = 0.0;

    printf("iteration # %d:\n", iter);
    /* Update centroids */
    update_centroids(mri_flash, mri_mem, NULL, mri_mask, centroids, nvolumes_total, num_classes);

    for (c=0; c < num_classes; c++)
    {
      printf("Centroids for class %d:", c+1);
      for (i=0; i < nvolumes_total; i++)
        printf(" %g ", centroids[c][i]);
      printf("\n");
    }

    /* Update covariance matrix and prior class-probability */
    /* Compute fuzzy covariance matrix F, one for each class */
    update_F(F, mri_flash, mri_mem, NULL, mri_mask, centroids, nvolumes_total, num_classes);

    if (fix_class_size ==1)
    {
      for (c=0; c < num_classes; c++)
        classSize[c] = 1;
    }
    else
    {
      /* update class size */
      for (c=0; c < num_classes; c++)
        classSize[c] = 0;
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            if (MRIvox(mri_mask, x, y, z) == 0) continue;
            for (c=0; c < num_classes; c++)
            {
              classSize[c] += MRIFvox(mri_mem[c], x, y, z);
            }
          }

      for (c=0; c < num_classes; c++)
        classSize[c] /= (float)brainsize;
    }

    if (nvolumes_total == 1)
    {
      for (c=0; c < num_classes; c++)
      {

        F[c]->rptr[1][1] = 1.0/F[c]->rptr[1][1];
        detF[c] = F[c]->rptr[1][1];
      }
    }
    else
    {

      /* Compute the inverse, and store in the original place */
      compute_inverseF(F, num_classes);

      /* compute determinant of inverseF */
      compute_detF(F, detF, num_classes);

    }

    for (c=0; c < num_classes; c++)
      printf("classSize[%d] = %g, DetF[%d] = %g\n", c, classSize[c], c, detF[c]);

    /* Update membership function */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0) continue;

          gcap = getGCAP(gca, mri_flash[0], transform, x, y, z);
          if (gcap == NULL)
          {
            MRIvox(mri_mask, x, y, z) = 0; // clear this voxel
            continue;
          }

          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            for (c=0; c < num_classes; c++)
            {
              printf("Prior for label %d = %g\n", index_to_label[c], getPrior(gcap, index_to_label[c]));
            }
          }

          sum_of_distance = 0;
          /* Compute distance */
          for (c=0; c < num_classes; c++)
          {
            /* record old membership values */
            oldMems[c] = MRIFvox(mri_mem[c], x, y, z);
            distance2 = distancems(mri_flash, centroids[c], F[c], detF[c], x, y, z);

            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("distance2 =%g \n",distance2);
            }

            /* detF is already the inverse */
            distance2 = exp(-0.5*distance2)*getPrior(gcap, index_to_label[c])*sqrt(detF[c]);
            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("prob =%g \n",distance2);
            }

            /* printf("distance2 = %g\n", distance2); */
            MRIFvox(mri_mem[c], x, y, z) = distance2;
            sum_of_distance += distance2;
          }

          if (sum_of_distance <= 1e-10)
          { /* Outlier */
            if (iter > 15)
              MRIvox(mri_mask, x, y, z) = 0;

            sum_of_distance = 1.0;
            //printf("(x,y,z) = (%d,%d,%d), sum_of_distance = %g\n",x,y,z, sum_of_distance);
            //ErrorExit(ERROR_BADPARM, "%s: overflow in computing membership function.\n", Progname);
          }

          for (c=0; c < num_classes; c++)
          {
            /* borrow distance2 here */
            distance2 = MRIFvox(mri_mem[c], x, y, z)/sum_of_distance;
            MRIFvox(mri_mem[c], x, y, z) = distance2;

            distance2 -= oldMems[c];
            if (distance2 < 0) distance2 = -distance2;

            /* if(distance2 > 0.99 && iter > 1){
               printf("oldmem = %g, newmem =%g, sum_of_distance = %g\n", oldMems[c], MRIFvox(mri_mem[c], x, y, z), sum_of_distance);
               } */

            if (max_change < distance2) max_change = distance2;

          }


        } /* end of all data points */

    printf("maxchange = %g\n", max_change);

  } /* end of while-loop */

  TransformFree(&transform) ;
  GCAfree(&gca);


  if (hard_segmentation)
  {
    MRI *mri_seg = NULL;

    mri_seg = MRIcopy(mri_mask, mri_seg);

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0) continue;
          i = 0;
          value = MRIFvox(mri_mem[0], x, y, z);
          for (j=1; j < num_classes; j++)
          {
            if (value < MRIFvox(mri_mem[j], x, y, z))
            {
              i = j;
              value = MRIFvox(mri_mem[j], x, y, z);
            }
          }

          MRIvox(mri_seg, x, y, z) = index_to_label[i];

        }

    sprintf(fname,"%sEM_hseg.mgz", out_prefx);
    MRIwrite(mri_seg, fname);
    MRIfree(&mri_seg);
  }

  if (HARD_ONLY == 0)
  {
    /* output membership functions */
    printf("Convert membership function to BYTE and write out\n");
    mri_tmp = MRIclone(mri_mask, NULL);

    for (c=0; c < num_classes; c++)
    {

      /*
      printf("Output membership volume for class %d\n",c) ;
      mri_tmp = MRIconform(mri_mem[c]) ;
      */
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            if (MRIvox(mri_mask, x, y, z) == 0)
            {
              MRIvox(mri_tmp,x,y,z) = 0;
              continue;
            }
            value =  255.0*MRIFvox(mri_mem[c], x, y, z) + 0.5;
            if (value > 255.0) value = 255;
            if (value < 0.0) value = 0;

            MRIvox(mri_tmp,x,y,z) = (unsigned char)value;
          }

      sprintf(fname,"%sEM_mem%d.mgz", out_prefx, c+1);
      MRIwrite(mri_tmp, fname);

    }

    MRIfree(&mri_tmp);
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("parameter estimation took %d minutes and %d seconds.\n", minutes, seconds) ;

  MRIfree(&mri_mask);

  MRIfree(&mri_label);

  for (i=0; i < num_classes; i++)
  {
    MRIfree(&mri_mem[i]);
    free(centroids[i]);
    MatrixFree(&F[i]);
  }

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
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "tissue"))
  {
    label_of_interest = atoi(argv[2]);
    nargs = 1;
    printf("refine the segmentation of structure with label %d\n",label_of_interest);
  }
  else if (!stricmp(option, "hardonly"))
  {
    HARD_ONLY = 1;
    printf("Do not output membership functions\n") ;
  }
  else if (!stricmp(option, "norm"))
  {
    normflag = 1;
    printf("Normalize input volumes to N(0,1)\n");
  }
  else if (!stricmp(option, "mask"))
  {
    mask_fname = argv[2];
    printf("using %s as mask for regions of interest \n", mask_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "rescale"))
  {
    rescale = 1;
    printf("Rescale the membership function to improve contrast.\n") ;
  }
  else if (!stricmp(option, "hard_seg"))
  {
    hard_segmentation = 1;
    printf("Output a hard segmentation to out_pre.hseg \n") ;
  }
  else if (!stricmp(option, "window"))
  {
    printf("window option not implemented\n");
    /*E* window_flag = 1 ; */
  }

  /*E* Interpolation method.  Default is trilinear, other options are
    nearest, cubic, sinc.  You can say -foo or -interp foo.  For sinc,
    you can say -interp sinc 3 or -interp sinc -hw 3 or -sinc 3 or
    -sinc -hw 3.  Maybe -hw 3 should imply sinc, but right now it
    doesn't.  */

  else if (!stricmp(option, "st") ||
           !stricmp(option, "sample") ||
           !stricmp(option, "sample_type") ||
           !stricmp(option, "interp"))
  {
    InterpMethod = MRIinterpCode(argv[2]) ;
    nargs = 1;
    if (InterpMethod==SAMPLE_SINC)
    {
      if ((argc<4) || !strncmp(argv[3],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
      {
        printf("using sinc interpolation (default windowwidth is 6)\n");
      }
      else
      {
        sinchalfwindow = atoi(argv[3]);
        nargs = 2;
        printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
      }
    }
  }
  else if (!stricmp(option, "sinc"))
  {
    InterpMethod = SAMPLE_SINC;
    if ((argc<3) || !strncmp(argv[2],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
    {
      printf("using sinc interpolation (default windowwidth is 6)\n");
    }
    else
    {
      sinchalfwindow = atoi(argv[2]);
      nargs = 1;
      printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
    }
  }
  else if (!stricmp(option, "sinchalfwindow") ||
           !stricmp(option, "hw"))
  {
    /*E* InterpMethod = SAMPLE_SINC; //? */
    sinchalfwindow = atoi(argv[2]);
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
  }
  else if (!stricmp(option, "gca"))
  {
    gca_fname = argv[2];
    printf("using %s as GCA \n", gca_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "xform"))
  {
    xform_fname = argv[2];
    printf("using %s as transformation to GCA \n", xform_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "label"))
  {
    label_fname = argv[2];
    printf("using %s as segmentation volume \n", label_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "mask_subcortical"))
  {
    mask_subcortical = 1;
    printf("mask subcortical GM region too. \n");
  }
  else if (!stricmp(option, "trilinear"))
  {
    InterpMethod = SAMPLE_TRILINEAR;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "cubic"))
  {
    InterpMethod = SAMPLE_CUBIC;
    printf("using cubic interpolation\n");
  }
  else if (!stricmp(option, "nearest"))
  {
    InterpMethod = SAMPLE_NEAREST;
    printf("using nearest-neighbor interpolation\n");
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
  }
  else if (!stricmp(option, "regularize"))
  {
    regularize = 1;
    lambda = atof(argv[2]);
    nargs = 1;
    printf("Regularize the covariance matrix with lambda = %g\n", lambda);
  }
  else switch (toupper(*option))
    {
      //    case 'M':
      //num_classes = atoi(argv[2]) ;
      //nargs = 1 ;
      //printf("Number of classes=%d\n", num_classes) ;
      //if(num_classes > MAX_CLASSES)
      // ErrorExit(ERROR_BADPARM, "%s: too many desired classes.\n", Progname);
      // break ;
    case 'T':
      noise_threshold = atof(argv[2]) ;
      printf("Threshold for background noise = %g\n", noise_threshold) ;
      printf("The threshold is applied on first input volume.\n") ;
      nargs = 1 ;
      break ;
    case 'E':
      tolerance = atof(argv[2]) ;
      printf("End tolerance = %g\n", tolerance) ;
      nargs = 1 ;
      break ;
    case 'R':
      max_iters = atoi(argv[2]) ;
      printf("Maximum number of iterations = %d\n", max_iters) ;
      nargs = 1 ;
      break ;
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
  printf("usage: %s [options] <volume> ... <output prefix>\n", Progname) ;
  printf("This program takes an arbitrary # of FLASH images as input,\n"
         "and performs EM-segmentation on the multidimensional\n"
         "intensity space.\n");
  printf("Options includes:\n");
  printf("\t -mask fname to set the brain mask volume\n");
  //  printf("\t -M # to set the number of classes \n");
  printf("\t -E # to set the convergence tolerance \n");
  printf("\t -R # to set the max number of iterations \n");
  printf("\t -T # to set the threshold for background \n");
  printf("\t -regularize # to set covariance matrix regularization parameter\n");
  printf("\t -norm to normalize input volumes before clustering \n");
  printf("\t -conform to conform input volumes (brain mask typically already conformed) \n");
  printf("\t -noconform to prevent conforming to COR \n");
  printf("\t -tissue #: tissue to be refined \n");
  printf("\t -label aseg_file: initial aseg segmentation \n");
  printf("\t -gca $GCA the atlas to be used \n");
  printf("\t -xform xfm.m3z xform for mapping to GCA \n");
  printf("\t -synthonly to prevent output membership functions \n");
  exit(code) ;

}

void update_centroids(MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes)
{
  int m, c, x, y, z, depth, height, width;
  double numer, denom;
  double mem, data;
  double scale ;
  //  double kappa = exp(-4.5);

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  /* Put the voxel iteration inside makes it easier to
   * program, but may take longer time. Otherwise, need to
   * declare numer and denom as matrices!
   */
  for (c = 0; c < num_classes; c++)
  {
    for (m=0; m < nvolumes_total; m++)
    {
      numer = 0;
      denom = 0;

      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            if (MRIvox(mri_mask, x, y, z) > 0)
            {
              // if(mri_lihood != NULL){
              //  scale =  MRIFvox(mri_lihood[c], x, y, z);
              // scale = scale/(scale + kappa);
              // }else
              scale = 1;
              mem = MRIFvox(mri_mem[c], x, y, z)*scale;
              data = MRIFvox(mri_flash[m], x, y, z);
              numer += mem*data;
              denom += mem;
            }

          }
      if (denom != 0.0)
        centroids[c][m] = numer/denom;
      else
      {
        ErrorExit(ERROR_BADPARM, "%s: overflow in computing centroids.\n", Progname);
      }
    }
  }


  return;
}

double distancems(MRI **mri_flash, float *centroids, MATRIX *F, double detF, int x, int y, int z)
{
  /* F would be the inverse of the covariance matrix of the class */
  int i, j, rows;
  double data1, data2;

  double mydistance = 0;

  rows = F->rows; /* actually the data dimensions */

  for (i=0; i < rows; i++)
  {
    if (debug_flag && x == Gx && y == Gy && z == Gz)
      printf("(mri_flash[%d] = %g, centroid[%d] = %g\n", i, MRIFvox(mri_flash[i], x, y, z), i, centroids[i]);
    data1 =  MRIFvox(mri_flash[i], x, y, z) - centroids[i];
    for (j=0; j< rows; j++)
    {
      data2 =  MRIFvox(mri_flash[j], x, y, z) - centroids[j];
      if (debug_flag && x == Gx && y == Gy && z == Gz)
        printf("(mri_flash[%d] = %g, centroid[%d] = %g\n", i, MRIFvox(mri_flash[i], x, y, z), i, centroids[i]);
      if (debug_flag && x == Gx && y == Gy && z == Gz)
        printf("data1 = %g, data2 = %g, F[%d][%d] = %g\n", data1, data2, i,j, F->rptr[i+1][j+1]);

      mydistance += data1*data2*F->rptr[i+1][j+1];
    }

  }

  /* No multiply by detF in the EM algorithm*/
  /* mydistance *= detF; */

  return mydistance;
}

void update_F(MATRIX **F, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes)
{

  int m1, m2, c, x, y, z, depth, height, width;
  double denom;
  double mem, data1, data2;
  double scale ;
  //  double kappa = exp(-4.5);
  MATRIX * tmpM = 0; /*used for covariance regularization */

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  for (c=0; c < num_classes; c++)
  {

    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=m1; m2 <= nvolumes_total; m2++)
      {
        F[c]->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
      }
    }

    denom = 0.0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0) continue;
          //   if(mri_lihood != NULL){
          //scale =  MRIFvox(mri_lihood[c], x, y, z);
          //  scale = scale/(scale + kappa);
          //}else
          scale = 1;
          mem = MRIFvox(mri_mem[c], x, y, z)*scale;
          /* mem = mem*mem; */ /* here differs from FCM */
          denom +=  mem;

          for (m1=0; m1 < nvolumes_total; m1++)
          {
            data1 = MRIFvox(mri_flash[m1], x, y, z) - centroids[c][m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              data2 = MRIFvox(mri_flash[m2], x, y, z) - centroids[c][m2];
              F[c]->rptr[m1+1][m2+1] += data1*data2*mem;
            }
          }

        } /* for all data points */

    if (denom <= 0.0)
      ErrorExit(ERROR_BADPARM, "%s: overflow in computing fuzzy covariance matrix.\n", Progname);

    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=m1; m2 <= nvolumes_total; m2++)
      {
        F[c]->rptr[m1][m2] /= denom;
        F[c]->rptr[m2][m1] = F[c]->rptr[m1][m2];
      }

      /* F[c]->rptr[m1][m1] += 1e-10; */ /* prevent F to be singular */
    } /* for m1, m2 */

    if (regularize)
    {
      for (m1=1; m1 <= nvolumes_total; m1++)
        F[c]->rptr[m1][m1] += eps;  /* prevent F to be singular */

      tmpM = MatrixInverse(F[c], tmpM);
      if (tmpM == NULL) continue;

      /* (1-lambda)* inv(F + eps I) + labmda I */
      for (m1=1; m1 <= nvolumes_total; m1++)
      {
        for (m2=m1; m2 <= nvolumes_total; m2++)
        {
          tmpM->rptr[m1][m2] = (1.0 - lambda)*tmpM->rptr[m1][m2];
          tmpM->rptr[m2][m1] = tmpM->rptr[m1][m2];
        }
        tmpM->rptr[m1][m1] += lambda;
      }

      F[c] = MatrixInverse(tmpM, F[c]);
    }

  } /* for(c== 0) */

  if (tmpM)
    MatrixFree(&tmpM);

  return;
}

void compute_detF(MATRIX **F, double *detF, int num_classes)
{

  int c, i, n;
  float tmpv;

  n = F[0]->rows;

  for (c=0; c < num_classes; c++)
  {

    tmpv = MatrixDeterminant(F[c]);

    if (tmpv < 0.000001)
    { /*singular */
      for (i=1; i <= F[c]->rows; i++)
        F[c]->rptr[i][i] += 0.01; /* try to recover it */

      tmpv = MatrixDeterminant(F[c]);
    }

    detF[c] = tmpv;

    /* detF[c] = powf(tmpv, 1.0/n); */ /* Already took the 1-nth power */
  }

  return;
}

void compute_inverseF(MATRIX **F, int num_classes)
{
  MATRIX *mTmp;
  int c, row, rows, cols;

  rows =  F[0]->cols;
  cols =  F[0]->cols;

  for (c=0; c < num_classes; c++)
  {
    mTmp = NULL;

    mTmp = MatrixInverse(F[c], mTmp);

    if (mTmp == NULL)
    { /* inverse doesn't exist */
      ErrorExit(ERROR_BADPARM, "%s: singular fuzzy covariance matrix.\n", Progname);
    }

    for (row=1; row <= rows; row++)
    {
      memmove((char *)(F[c]->rptr[row]), (char *)mTmp->rptr[row],
             (cols+1)*sizeof(float)) ;
    }

    MatrixFree(&mTmp);
    mTmp = NULL;
  }


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
