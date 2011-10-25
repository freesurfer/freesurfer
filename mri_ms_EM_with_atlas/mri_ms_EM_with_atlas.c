/**
 * @file  mri_ms_EM_with_atlas.c
 * @brief takes an arbitrary # of FLASH images, performs EM-seg
 *
 * Trying to use more heuristics for MEF segmentation
 * assume using two averages, first is 30, second is 5
 * first use their ratio to find CSF! Then use CSF to refine brain mask
 * e.g., using an initial WM segmentation or aseg to initialize
 * may better use an atlas (why not)
 * Bias correction may need be done earlier (yes, or refine it here too)
 */
/*
 * Original Author: Xiao Han
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/10/25 14:13:02 $
 *    $Revision: 1.10 $
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
#include "PoissonSolver.h"

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

#define TEST_BRAIN_MASK 0
#define FIX_DURA 1

/* maximum number of classes */
#define MAX_CLASSES 20
#define MEMTHRESHOLD 0.7

#define weight30   0.990247
#define weight5  -0.139292

//static double detF_lowbound = 5e-7;

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

/* If set at 0.7, then Mahalanobis distance is 26.x; if set to 0.5, it's only 14.6. Of course, lower threshold includes more partial volumed voxels, and increasecovariance */
/* For fair comparison, need to generate a hard segmentation! */

//static int  tests = 0; //doesn't help much

//static int fix_class_size = 1;

static int no_INU = 0;
static double kappa = 1e-8; /* 0.01 is too big */ //this number seems need to be reduced if using 16 echos, but need be this large for just using average; so critical, not good. How to make it robust?

static int debug_flag = 0;

static char *gca_fname = NULL; /* filename for GCA */

static char *xform_fname = NULL; /* filename for xform to GCA */

static int T1_PD = 0; // if 1, input 1 is T1, and input 2 is PD
static int invert = 0; // if 1, invert output image contrast

static int remove_cerebellum = 0; // if 1, remove cerebellum in final synthesized image

void mri_invert_contrast_and_scale(MRI *mri, MRI *mri_mask, double scale);

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static int normflag = 0; /* normalize input volume */

/* Values lower than threshold are considered as background, and not counted
 * in clustering
 */
static double noise_threshold = 1.0; /* threshold for background noise */

static double tolerance = 0.01; /* Convergence criterion */

static double lap_weight = 2.0; /* Weight for bias-field smoothness */

static int num_classes = 4; /* default number of classes */

static int clear_dura = 0;

static int max_iters = 200; /* defualt maximal iteration number */

static int whole_volume = 0; /* when LDA, project whole volume */
/* Is this necessary? it may affect the scaling
 * and reduce dynamic range of useful regions
* But without it, how do I measure background
* noise; but who cares about background noise!
*/

static int hard_segmentation = 0;

static char *label_fname = NULL; /* filename for segmentation */

static int conform = 1 ;

static int mask_subcortical = 0;

static int SYNTH_ONLY = 0;

static int rescale = 0;

static int fuzzy_lda = 0; /* whether use Fuzzy cov and fuzzy centroid for LDA */

static double mybeta = 0; /*weight for MRF */

/* eps and lambda are used for covariance regularization */
static double eps = 1e-20; //such a big change seems having no effect
static double lambda = 0.1; //0.2;
static int regularize = 0;
static int noprior = 0;

/* Clustering is performed only within ROI */
static char *mask_fname = NULL; /* filename for ROI mask */

void indexx(unsigned long n, float arr[], unsigned long indx[]);

static void usage_exit(int code) ;

MRI *MRInormalizeXH(MRI *mri_src, MRI *mri_dst, MRI *mri_mask);

double NbhdLikelihood(MRI **mri_mem, MRI *mri_mask, int x, int y, int z, int label, int num_classes);

void fuzzyLDAweights(float *weights, MATRIX **F, float **centroids, int classID1, int classID2, float size1, float size2,  int nvolumes_total);

void update_centroids1D(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int num_classes);

void MRI_FCM(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int  num_classes);

void MRI_EM(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int  num_classes);

void compute_bias(MRI **mri_mult, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI *mri_mask, float **centroids, MATRIX **F, int nvolumes_total, int num_classes);

void update_centroids(MRI **mri_flash, MRI **mri_mem, MRI **mri_likelihood, MRI **, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes);

void update_F(MATRIX **F, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI **, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes);

void compute_detF(MATRIX **F, double *detF, int num_classes);

void compute_inverseF(MATRIX **F, int num_classes);

double distancems(MRI **mri_flash, MRI **, float *centroids, MATRIX *F, double detF, int x, int y, int z);
double distancems_noF(MRI **mri_flash, MRI**mri_mult, float *centroids, MATRIX *F, double detF, int x, int y, int z); //ignore off-diagnal terms

//this one ignores covariance matrix
double distancems2(MRI **mri_flash, MRI**mri_mult, float *centroids, MATRIX *F, double detF, int x, int y, int z, int channel);
void  update_LDAmeans(MRI **mri_flash, MRI **mri_mem, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2, float threshold);

void computeLDAweights(float *weights, MRI **mri_flash, MRI **mri_mem, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2, float threshold);

void  set_equilavent_classes(int *equivalent_classes);

static int InterpMethod = SAMPLE_TRILINEAR;
static int sinchalfwindow = 3;
static int ldaflag = 0; /* Using LDA methods to synthesize volume */
static int class1 = 0; /* to be used for LDA */
static int class2 = 0; /* to be used for LDA */

static char *hseg_fname  = NULL;
static char *synth_fname  = NULL;

#define MAX_IMAGES 200

int
main(int argc, char *argv[])
{
  char   **av, *in_fname, fname[100];
  int    ac, nargs, i, j,  x, y, z, width, height, depth, iter, c;
  MRI    *mri_flash[MAX_IMAGES], *mri_mem[MAX_CLASSES], *mri_mask, *mri_label;
  MRI *mri_flash_backup[MAX_IMAGES];
  MRI  *mri_mult[MAX_IMAGES];
  MRI *cerebellum_mask;
  MRI  *mri_prior[MAX_IMAGES]; //prior from atlas; set to equal for all classes at background? Where is prior for dura?
  MRI *mri_mem_old[MAX_CLASSES];
  MRI *mri_tmp = NULL;
  MRI *mri_tmp2 = NULL;
  char   *out_prefx ;
  int    msec, minutes, seconds, nvolumes, nvolumes_total, label, clabel ;
  struct timeb start ;
  MATRIX *F[MAX_CLASSES];
  float *centroids[MAX_CLASSES];
  float centroids1D[MAX_CLASSES];
  float max_val, min_val, value;
  float max_change = 1000.0;
  double detF[MAX_CLASSES]; /* determinant of the covariance matrix */
  double distance2, sum_of_distance; /* actually the inverse distance */
  double oldMems[MAX_CLASSES];
  float classSize[MAX_CLASSES];
  int brainsize;
  int num_outlier;
  int nx=0, ny=0, nz=0, nc;
  int n;
  int xn, yn, zn, index;
  float mean1, std1, mean2, std2, denom1, denom2, threshold;
  float scale;
  int MLE_label;
  double max_prior;
  int equiv_class[MAX_CLASSES];
  double nbhdP;
  double inv_sqrt_2piM; /* 1/sqrt((2 pi)^M) */ //normalization factor
  float total_prior;
  float *LDAmean1, *LDAmean2, *LDAweight;
  int class_dura;
  GCA *gca;
  GCA_PRIOR *gcap;
  TRANSFORM *transform;
  float val1;

  /* Used for sorting using Numerical recipe */
  float NRarray[MAX_CLASSES+1];
  unsigned long NRindex[MAX_CLASSES+1];
  int indexmap[MAX_CLASSES + 1];

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_ms_EM_with_atlas.c,v 1.10 2011/10/25 14:13:02 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (ldaflag)
  {
    if (class1 == class2 || class1 < 1 || class1 > num_classes
        || class2 < 1 || class2 > num_classes)
    {
      ErrorExit(ERROR_NOFILE, "%s: The two class IDs must be within [1, %d]",
                Progname, num_classes) ;
    }

  }

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
  {
    usage_exit(1) ;
  }

  out_prefx = argv[argc-1] ;

  printf("command line parsing finished\n");

  if (gca_fname == NULL || xform_fname == NULL)
  {
    printf("This version assumes that an atlas and xformation is needed! Exit\n");
    exit(0);
  }

  if (num_classes != 3)
  {
    printf("Currently assuming num_classes  is 3.\n");
    num_classes = 3;
  }

  //////////////////////////////////////////////////////////////////////////////////
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
      mri_tmp = 0;
    }



    /* Change all volumes to float type for convenience */
    if (mri_flash[nvolumes]->type != MRI_FLOAT)
    {
      printf("Volume %d type is %d\n", nvolumes+1, mri_flash[nvolumes]->type);
      printf("Change data to float type \n");
      mri_tmp = MRIchangeType(mri_flash[nvolumes], MRI_FLOAT, 0, 1.0, 1);
      MRIfree(&mri_flash[nvolumes]);
      mri_flash[nvolumes] = mri_tmp;
      mri_tmp = 0; //swap
    }

    mri_flash_backup[nvolumes] = MRIcopy(mri_flash[nvolumes], NULL);

    MRImakePositive(mri_flash[nvolumes], mri_flash[nvolumes]);
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

  printf("-----I need to fix mri_likelihood for dura to make sure it doesn't got exchanged with CSF ---\n");
  printf("---- Maybe unnecessary, Dura seems mainly compete with GM ---\n");

  printf("How to make dura a subclass of GM?? ---\n");

  printf("After restricting the detF, this method seems working for all pfizer dat: still some questions are pounding: 1. should I reduce the detF threshold; 2. should I normalize the volume of all classes, like in G-K scheme; 3. weirdly, for some subjects, the cross-correlation, i.e., off-diagnal of F[1] (dura) is positive, some are negative! 4. Should make sure the detF threshold independent of input signal scales; 5. why when I didn't make all images to positive, convergence become worse? Due to gain field? Seems so. 6. Relative class size seems alright, i.e., CSF = 3.6, dura = 1.3, GM = 4.8, WM = 4.4, but from the GM seems too much from the output segmentation. If looking at the synthesize image, GM is much less, somehow EM makes many CSF classified as GM; 7. should I tie the centroids of dura to CSF and GM, as I did now?? well, I guess it still work if I didn't. 8. It's now clear that every algorithm requires fine-tuning to be working! 9. Can I get rid of dependency on aseg? May be using histogram to segment WM, or using FReesurfer on flash30.  \n");

  printf("May need to simply normalize all covariance matrices, and only rely on classSize to control the relative size of each class, like in G-H scheme\n");

  //This brain-masking procedure works perfect for 1001, but for 1006, eye-sockets are included, while brain.mgz doesn't! So, still no clear advantage over mri_watershed
  // after I add an additional threshold, i.e., flash30 cannot be higher than flash5*1.5, then eye-socks are gone

  printf("read in atlas and registration file\n");
  fprintf(stderr, "reading gca from %s...\n", gca_fname);
  gca = GCAread(gca_fname);
  if (!gca)
  {
    ErrorExit(ERROR_NOFILE, "%s:could not read classifier array from %s", Progname, gca_fname);
  }


#if 0
  mri_tmp = GCAbuildMostLikelyVolumeFrame(gca, NULL, 0);
  MRIwrite(mri_tmp, "gca_channel0.mgz");
  MRIfree(&mri_tmp);
  mri_tmp = GCAbuildMostLikelyVolumeFrame(gca, NULL, 1);
  MRIwrite(mri_tmp, "gca_channel1.mgz");
  MRIfree(&mri_tmp);
  exit(0);
#endif

  fprintf(stderr, "reading transform from %s ...\n", xform_fname);
  transform = TransformRead(xform_fname);
  if (!transform)
  {
    ErrorExit(ERROR_NOFILE, "%s:could not open transform %s", Progname, xform_fname);
  }
  printf("Transform read in, now computing its inverse...\n");
  TransformInvert(transform, mri_flash[0]);
  printf("Transform for GCA is inverted \n");

  printf("Build atlas prior for each class \n");

  /* Allocate memory for membership volumes */
  for (i=0; i < num_classes; i++)
  {
    mri_prior[i] = MRIclone(mri_flash[0], NULL);
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
    {
      ErrorExit(ERROR_BADPARM, "%s: mask volume size doesn't macth data volumes\n", Progname);
    }

    if (mri_mask->type != MRI_UCHAR)
    {
      ErrorExit(ERROR_BADPARM, "%s: mask volume is not UCHAR type \n", Progname);
    }
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

    /* Compute the mask from the first input volume */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIFvox(mri_flash[0], x, y, z) < (float) noise_threshold)
          {
            MRIvox(mri_mask, x, y,z) = 0;
          }
        }

  }

#if 0
  if (T1_PD == 1)
  {
    printf("further threshold with T1 map \n");
    /* Compute the mask from the first input volume */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIFvox(mri_flash[0], x, y, z) > 3000)
          {
            MRIvox(mri_mask, x, y,z) = 0;
          }
        }

  }
#endif

  set_equilavent_classes(equiv_class);


  mri_label = MRIclone(mri_mask, NULL);
  cerebellum_mask = MRIclone(mri_mask, NULL);

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        MRIvox(cerebellum_mask, x, y, z) = 0;
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }
        for (i=0; i < num_classes; i++)
        {
          MRIsetVoxVal(mri_prior[i], x, y, z, 0, 0.01);
        }

        if (!GCAsourceVoxelToNode(gca, mri_flash[0], transform, x, y, z, &xn, &yn, &zn))
        {
          gcap = getGCAP(gca, mri_flash[0], transform, x, y, z) ;

          if (gcap == NULL)
          {
            printf("Failed to get prior for voxel (%d,%d,%d)\n",x,y,z);
            MRIvox(mri_mask, x, y, z) = 0;
            continue;
          }

          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("gcap->nlabels = %d\n", gcap->nlabels);
            printf("(xn,yn,zn) = (%d, %d, %d)\n", xn, yn, zn);
          }

          if (gcap->nlabels < 1)
          {
            MRIvox(mri_mask, x, y, z) = 0;
            continue;
          }
          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("gcap->labels[%d]=%d\n", 0, gcap->labels[0]);
          }

          max_prior = gcap->priors[0];
          MLE_label =  gcap->labels[0];
          // going through gcap labels
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("gcap->labels[%d]=%d\n", n, gcap->labels[n]);
            }

            clabel = equiv_class[gcap->labels[n]];
            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("equiv_label=%d, prior = %g\n", clabel,gcap->priors[n] );
            }
            if (clabel > 0 && clabel <= num_classes)
            {
              MRIFvox(mri_prior[clabel-1],x, y,z) += gcap->priors[n];
            }

            if (gcap->priors[n] > max_prior)
            {
              max_prior = gcap->priors[n];
              MLE_label = gcap->labels[n];
            }
          }

          if (MLE_label == 7 || MLE_label == 8 || MLE_label == 16 || MLE_label == 46 || MLE_label == 47)
          {
            MRIvox(cerebellum_mask, x, y, z) = 1;  //for later removal of cerebellum and brain stem and maybe sth else if needed
          }


          MLE_label = equiv_class[MLE_label];
          if (MLE_label == 0)
          {
            MRIvox(mri_label, x, y, z) = (unsigned char) MLE_label;
          }
          else if (MLE_label == 1)
          {
            MRIvox(mri_label, x, y, z) = 50;
          }
          else if (MLE_label == 2 || MLE_label == 4)
          {
            MRIvox(mri_label, x, y, z) = 41;
          }
          else if (MLE_label == 3)
          {
            MRIvox(mri_label, x, y, z) = 42;
          }
        }

      }

  TransformFree(&transform) ;
  GCAfree(&gca);

#if 0
  //  MRIwrite(mri_prior[0],"prior_CSF.mgz");
  // MRIwrite(mri_prior[1],"prior_WM.mgz");
  // MRIwrite(mri_prior[2],"prior_GM.mgz");

  MRIwrite(mri_label, "MLE_label.mgz");
  exit(0);
#endif

  /* Allocate memory for likelihood volumes */
  for (i=0; i < num_classes; i++)
  {
    mri_mem_old[i] = MRIclone(mri_flash[0], NULL);
  }

  /* Allocate memory for bias-field term */
  for (i=0; i <nvolumes_total; i++)
  {
    mri_mult[i] = MRIclone(mri_flash[i], NULL);
  }

  brainsize = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        for (i=0; i <nvolumes_total; i++)
        {
          MRIFvox(mri_mult[i], x, y, z) = 1.0;
        }

        if (MRIvox(mri_mask, x, y, z) <= 0)
        {
          continue;
        }

        /* since mri_mem_old is now likelihood, it need be initialized */
        total_prior = 0.0;
        for (i=0; i < num_classes; i++)
        {
          total_prior += MRIFvox(mri_prior[i],x,y,z);
        }
        for (i=0; i < num_classes; i++)
        {
          if (total_prior < 0.7)
          {
            MRIFvox(mri_mem_old[i],x,y,z) = 0;
          }
          else
          {
            MRIFvox(mri_mem_old[i],x,y,z) = 1.0;  //or just set all to total_prior??
          }
        }

        brainsize++;
      }

  printf("Brain size = %d\n", brainsize);
  if (brainsize < 1)
  {
    printf("region to be segment has zero size. exit.\n");
    exit(1);
  }

  /* Allocate memory for membership volumes */
  for (i=0; i < num_classes; i++)
  {
    mri_mem[i] = MRIclone(mri_flash[0], NULL);
    //    mri_mem_old[i] = MRIclone(mri_flash[0], NULL);
    centroids[i] = (float *)malloc(nvolumes_total*sizeof(float));
    F[i] = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);
    classSize[i] = 1.0/(float)num_classes;
  }

  /* initialize F's to identity matrices */
  for (c=0; c < num_classes; c++)
  {
    for (i=1; i <= nvolumes_total; i++)
    {
      for (j=i+1; j <= nvolumes_total; j++)
      {
        F[c]->rptr[i][j] = 0.0;
        F[c]->rptr[j][i] = 0.0;
      }
      F[c]->rptr[i][i] = 1.0;
    }
    detF[c] = 1.0;
  }

  printf("nvolumes_total = %d\n", nvolumes_total);

  if (nvolumes_total != 2)
  {
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          for (c=0; c < num_classes; c++)
          {
            MRIFvox(mri_mem[c],x,y,z) = 0;
          }

          if (MRIvox(mri_mask,x,y,z) == 0)
          {
            continue;
          }
          //threshold on prior to get initial centroids estimation
          for (c=0; c < num_classes; c++)
          {
            if (MRIFvox(mri_prior[c], x, y, z) > 0.5)
            {
              MRIFvox(mri_mem[c],x,y,z) = MRIFvox(mri_prior[c], x, y, z);
            }
            else
            {
              MRIFvox(mri_mem[c], x, y, z) = 0;
            }
          }

        }
  }
  else
  {
    if (T1_PD == 0)
    {
      printf("Using ratio of volume1 and volume 2 to find CSF, assuming volume 1 is flash30 and volume 2 is flash 5. Of cause should be able to decide from fa of each volume\n");

      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            for (c=0; c < num_classes; c++)
            {
              MRIFvox(mri_mem[c],x,y,z) = 0;
            }
            if (MRIvox(mri_mask,x,y,z) == 0)
            {
              continue;
            }

            //threshold on prior to get initial centroids estimation
            for (c=1; c < num_classes; c++)
            {
              if (MRIFvox(mri_prior[c], x, y, z) > 0.5)
              {
                MRIFvox(mri_mem[c],x,y,z) = MRIFvox(mri_prior[c], x, y, z);
              }
              else
              {
                MRIFvox(mri_mem[c], x, y, z) = 0;
              }
            }

            value = MRIFvox(mri_flash[0],x,y,z);

            if (MRIFvox(mri_flash[1], x, y, z) > 2.5*value && MRIFvox(mri_prior[0],x, y, z) > 0.3)
            {
              MRIFvox(mri_mem[0],x,y,z) = 1;
            }
            else
            {
              MRIFvox(mri_mem[0],x,y,z) = MRIFvox(mri_prior[0],x,y,z);
            }
          }
    }
    else
    {
      printf("Threshold T1 to get CSF map \n");

      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            for (c=0; c < num_classes; c++)
            {
              MRIFvox(mri_mem[c],x,y,z) = 0;
            }
            if (MRIvox(mri_mask,x,y,z) == 0)
            {
              continue;
            }

            //threshold on prior to get initial centroids estimation
            for (c=1; c < num_classes; c++)
            {
              if (MRIFvox(mri_prior[c], x, y, z) > 0.5)
              {
                MRIFvox(mri_mem[c],x,y,z) = MRIFvox(mri_prior[c], x, y, z);
              }
              else
              {
                MRIFvox(mri_mem[c], x, y, z) = 0;
              }
            }

            value = MRIFvox(mri_flash[0],x,y,z);

            if (MRIFvox(mri_flash[0], x, y, z) > 2200 && MRIFvox(mri_prior[0],x, y, z) > 0.3)
            {
              MRIFvox(mri_mem[0],x,y,z) = 1;
            }
            else
            {
              MRIFvox(mri_mem[0],x,y,z) = MRIFvox(mri_prior[0],x,y,z);
            }
          }

    }

    //    MRIwrite(mri_mem[0],"mem_CSF_init.mgz");
    // MRIwrite(mri_mem[1],"mem_WM_init.mgz");
    // MRIwrite(mri_mem[2],"mem_GM_init.mgz");
  }

  printf("Note that the above step was performed on original inputs, i.e., unnormalized!!\n");

  if (T1_PD == 1)
  {
    printf("Invert image constrast so that WM is brightest\n");
    if (nvolumes_total == 2)
    {
      mri_invert_contrast_and_scale(mri_flash[0], mri_mask, 30);
      mri_invert_contrast_and_scale(mri_flash[1], mri_mask, 100);
    }
    else
    {
      for (i=0; i <nvolumes_total; i++)
      {
        mri_invert_contrast_and_scale(mri_flash[i], mri_mask, 1.0);
      }
    }
  }
  /* Normalize input volumes */
  if (normflag)
  {
    printf("Normalize input volumes to variance 1 (not zero-mean though!)\n");
    for (i=0; i <nvolumes_total; i++)
    {
      mri_flash[i] = MRInormalizeXH(mri_flash[i], mri_flash[i], mri_mask);
    }
    printf("Normalization done.\n");
  }

  /* Use the membership function obtained above to get initial centroids
   * for the multi-dimensional clustering
   */

  inv_sqrt_2piM = 1.0/pow(6.28318531, 0.5*nvolumes_total);

  printf("Compute centroids for first time, and use it to adjust prior (due to bad MEF atlas, maybe I should use MPRAGE atlas)\n");

  update_centroids(mri_flash, mri_mem, mri_mem_old, mri_mult, mri_mask, centroids, nvolumes_total, num_classes);
  update_F(F, mri_flash, mri_mem, mri_mem_old, mri_mult, mri_mask, centroids, nvolumes_total, num_classes);
  mean1 =  centroids[0][0];
  std1 = sqrt(F[0]->rptr[1][1]);

  if (T1_PD == 1)
  {
    MRIwrite(mri_flash[0],"flash_0.mgz");
    MRIwrite(mri_flash[1],"flash_1.mgz");
  }

  printf("mean1 =%g, std1 = %g\n", mean1, std1);

  if (nvolumes <= 2)
  {
#if 0 //old version
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask,x,y,z) == 0)
          {
            continue;
          }
          if (MRIFvox(mri_prior[1], x, y, z) > 0.1)
          {
            continue;
          }
          val1 = MRIFvox(mri_flash[0], x, y, z);
          if (val1 < mean1 + std1)
          {
            val2 =  0.5*(MRIFvox(mri_prior[0], x, y, z) + MRIFvox(mri_prior[2], x, y, z));
            MRIFvox(mri_prior[0], x, y, z) =  val2; //CSF
            MRIFvox(mri_prior[2], x, y, z) =  val2; //GM
          }
        }

# else
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask,x,y,z) == 0)
          {
            continue;
          }
          if (MRIFvox(mri_prior[2], x, y, z) < 0.5)
          {
            continue;
          }
          val1 = MRIFvox(mri_flash[0], x, y, z);
          //The following doesn't work for BIRN data
          // seems I need to use true CSF to estimate initial threshold; or just finish the atlas!!
          // note this ad hoc fix of atlas works only for INU corrected flash30
          if (val1 < mean1 + std1)
          {
            //   if(val1 < mean1 + 0.5*std1){
            //     val2 =  0.5*(MRIFvox(mri_prior[0], x, y, z) + MRIFvox(mri_prior[2], x, y, z));
            MRIFvox(mri_prior[0], x, y, z) = 0.8; // val2; //CSF
            MRIFvox(mri_prior[2], x, y, z) = 0.2; // val2; //GM
          }
        }

#endif

    // MRIwrite(mri_prior[0],"prior_CSF.mgz");
    // MRIwrite(mri_prior[1],"prior_WM.mgz");
    // MRIwrite(mri_prior[2],"prior_GM.mgz");
    // exit(0);
  }

  /* Perform iterative fuzzy-clustering until convergence */
  iter = 0;
  max_change = 10000.0;
  while (max_change > tolerance && iter < max_iters)
  {
    iter++;
    max_change = 0.0;

    printf("iteration # %d:\n", iter);
    /* Update centroids */
    update_centroids(mri_flash, mri_mem, mri_mem_old, mri_mult, mri_mask, centroids, nvolumes_total, num_classes);

    /* Update covariance matrix and prior class-probability */
    /* Compute fuzzy covariance matrix F, one for each class */
    update_F(F, mri_flash, mri_mem, mri_mem_old, mri_mult, mri_mask, centroids, nvolumes_total, num_classes);

    if (T1_PD == 1)
    {
      printf("set F os WM to be equal to GM \n");
      for (i=1; i <= nvolumes_total; i++)
      {
        for (j=1; j <= nvolumes_total; j++)
        {
          F[1]->rptr[i][j] = 0.1*F[1]->rptr[i][j] + 0.9* F[2]->rptr[i][j];
        }

      }
    }

    for (c=0; c < num_classes; c++)
    {
      printf("Centroids for class %d:", c+1);
      for (i=0; i < nvolumes_total; i++)
      {
        printf(" %g ", centroids[c][i]);
      }
      printf("\n");
    }

    MatrixPrint(stdout, F[0]);
    MatrixPrint(stdout, F[1]);
    MatrixPrint(stdout, F[2]);

    if (nvolumes_total == 1)
    {
      for (c=0; c < num_classes; c++)
      {

        F[c]->rptr[1][1] = 1.0/F[c]->rptr[1][1];
        detF[c] = F[c]->rptr[1][1]; /* detF is actually the det of inverseF */
      }
    }
    else
    {
      /* Compute the inverse, and store in the original place */
      compute_inverseF(F, num_classes);

      /* compute determinant of inversed F */
      compute_detF(F, detF, num_classes);

      /* Use identity */
      /* This will lead to every voxel belongs equally to all three classes */
      /* for(c=0; c < num_classes; c++){
      for(i=1; i <= nvolumes_total; i++){
      for(j=i+1; j <= nvolumes_total; j++){
      F[c]->rptr[i][j] = 0.0;
      F[c]->rptr[j][i] = 0.0;
      }
      F[c]->rptr[i][i] = 1.0;
      }
      detF[c] = 1.0;
      }
      */
    }

    for (c=0; c < num_classes; c++)
    {
      printf("DetF[%d] = %g\n", c, detF[c]);
    }

    /* Update membership function */
    /* first copy old_ones */
    /* I will use mri_mem_old to store likelihood */
    /*
      for(c=0; c < num_classes; c++){
      mri_mem_old[c] = MRIcopy(mri_mem[c], mri_mem_old[c]);
      }
    */

    /* Need to change this to Gauss-Seidel, otherwise it's too slow */
    num_outlier = 0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          sum_of_distance = 1e-20;
          /* Compute distance */
          for (c=0; c < num_classes; c++)
          {
            /* record old membership values */
            oldMems[c] = MRIFvox(mri_mem[c], x, y, z);
            distance2 = distancems(mri_flash, mri_mult, centroids[c], F[c], detF[c], x, y, z);
            distance2 = exp(-0.5*distance2);
            /* Distance2 now is the likelihood, store it */
            MRIFvox(mri_mem_old[c], x, y, z) = distance2;

            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("exp(-0.5*distance(%d)) =%g \n", c, distance2);
            }

            //     distance2 *= (classSize[c]*sqrt(detF[c])*inv_sqrt_2piM);
            //distance2 *= (classSize[c]*sqrt(detF[c]));

            //just for test
            if (noprior)
            {
              distance2 *= sqrt(detF[c]);
            }
            else
            {
              distance2 *= (sqrt(detF[c])*MRIFvox(mri_prior[c],x,y,z));
            }


            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("LE(%d) = %g, prior  =%g \n", c, distance2,MRIFvox(mri_prior[c],x,y,z));
            }

            /* Compute neighborhood PDF */
            /* Need to find the ML estimate for each neighbor and then sum up
               the prior. This will be very time-consuming; easier if use atlas
               To make it faster, have to store previous membership functions
            */

            //     nbhdP = 100* NbhdLikelihood(mri_mem_old, mri_mask, x, y, z, c, num_classes);
            /* Use current membership */
            /* Gives better convergence rate */
            /* Maybe a Gauss-seidel iteration is even better */
            if (!FZERO(mybeta))
            {
              //     if(0){
              nbhdP = NbhdLikelihood(mri_mem, mri_mask, x, y, z, c, num_classes);
              if (debug_flag && x == Gx && y == Gy && z == Gz)
              {
                printf("c= %d, nbhdP =%g \n", c, nbhdP);
              }
              distance2 *= nbhdP;
            }

            /* printf("distance2 = %g\n", distance2); */
            MRIFvox(mri_mem[c], x, y, z) = distance2;
            sum_of_distance += distance2;
          }

          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("sum_of_distance= %g \n", sum_of_distance);
          }

          if (sum_of_distance <= 0)   /* Outlier */
          {
            sum_of_distance = 1.0; //this hard truncation seems a bad idea! but sounds good! Anyway, using it lead to slow convergence
            num_outlier++; //the number will keep growing
            // printf("(x,y,z) = (%d,%d,%d), sum_of_distance = %g\n",x,y,z, sum_of_distance);
            // ErrorExit(ERROR_BADPARM, "%s: overflow in computing membership function.\n", Progname);
          }

          for (c=0; c < num_classes; c++)
          {
            /* borrow distance2 here */
            distance2 = MRIFvox(mri_mem[c], x, y, z)/sum_of_distance;
            MRIFvox(mri_mem[c], x, y, z) = distance2;

            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("memship(%d) =%g \n", c, distance2);
            }

            distance2 -= oldMems[c];
            if (distance2 < 0)
            {
              distance2 = -distance2;
            }


            //     if(distance2 > 0.99 && iter > 1){
            // printf("oldmem = %g, newmem =%g, sum_of_distance = %g\n", oldMems[c], MRIFvox(mri_mem[c], x, y, z), sum_of_distance);
            //}

            if (max_change < distance2)
            {
              max_change = distance2;
              nx = x, ny = y, nz = z;
              nc = c;
            }
          }


        } /* end of all data points */

    printf("num_outlier = %d\n", num_outlier);
    printf("maxchange = %g (at (%d,%d,%d))\n", max_change, nx, ny, nz);
    /* printf("oldmem(%d) = %g, new = %g\n", nc, MRIFvox(mri_mem_old[nc], nx, ny, nz),MRIFvox(mri_mem[nc], nx, ny, nz)); */ // no longer valid

    printf("compute bias field for each channel and correct INU...\n");

    //    if(iter == 9){
    //  MRIwrite(mri_mem[0],"mem_CSF_iter7.mgz");
    //  MRIwrite(mri_mem[1],"mem_WM_iter7.mgz");
    // MRIwrite(mri_mem[2],"mem_GM_iter7.mgz");
    //}

    //    if(!no_INU &&  max_change > 0.05)
    if (!no_INU &&  max_change > 0.15)
      //    if(!no_INU &&  iter == 1)
    {
      compute_bias(mri_mult, mri_flash, mri_mem, mri_mem_old, mri_mask, centroids, F, nvolumes_total, num_classes);
    }

    // if(iter == 9){
    //  MRIwrite(mri_mult[0], "mult_0_iter7.mgz");
    //  MRIwrite(mri_mult[1], "mult_1_iter7.mgz");
    // }

  } /* end of while-loop */

  scale = (centroids[1][0] - centroids[2][0])/(centroids[2][0] - centroids[0][0]);
  if (scale > 1.0)
  {
    scale = 1.0;
  }

  printf("Compute final membership function without atlas prior\n");
  //  printf("scale distance to CSF by %g\n", scale);
#if 0
  printf("enhance CSF by moving centroids up\n");
  if (num_classes == 3)
  {
    centroids[0][0] = (1-scale)*(2*centroids[2][0] - centroids[1][0]) + scale*centroids[0][0];
    printf("centroid for CSF is reset to %g\n", centroids[0][0]);
  }
#endif

  printf("How to use prior in the fuzzification??? Direct apply not working well\n");
#if 1
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }

        sum_of_distance = 1e-20;
        /* Compute distance */
        for (c=0; c < num_classes; c++)
        {
          // distance2 = distancems(mri_flash, mri_mult, centroids[c], F[c], detF[c], x, y, z);

          //this distance ignores off-diagonal terms of covariance matrix
          // CSF is more smeared!
          distance2 = distancems_noF(mri_flash, mri_mult, centroids[c], F[c], detF[c], x, y, z);
          //   distance2 = exp(-0.5*distance2);

          //  if(c == 0) distance2 *= scale;

          //10-25-05: try to add in prior; this indeed removed PVE voxels, but outside is noisy
          //the atlas need to be improved for outside! Maybe I need to use MPRAGE atlas

          //SQUARE THE DISTANCE SOMEHOW MAKE the combined image look more smooth (not more fuzzy though), more like hard segmentation but look nicer, as if prior was used
          if (noprior)
          {
            distance2 = 1.0/(distance2*distance2 + 1e-30);
          }
          else
          {
            distance2 = 1.0/(distance2*distance2 + 1e-30);
          }
          //using atlas will be even smoother
          //   distance2 = sqrt(MRIFvox(mri_prior[c],x,y,z))/(distance2 + 1e-30);

          //10-25-05: try to add in MRF term
          if (!FZERO(mybeta))
          {
            nbhdP = NbhdLikelihood(mri_mem, mri_mask, x, y, z, c, num_classes);
            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("c= %d, nbhdP =%g \n", c, nbhdP);
            }
            distance2 /= nbhdP;
          }


          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("1.0/distance(%d) =%g \n", c, distance2);
          }

          //   distance2 *= sqrt(detF[c]);
          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("LE(%d) =%g \n", c, distance2);
          }

          /* printf("distance2 = %g\n", distance2); */
          MRIFvox(mri_mem[c], x, y, z) = distance2;
          sum_of_distance += distance2;
        }

        if (debug_flag && x == Gx && y == Gy && z == Gz)
        {
          printf("sum_of_distance= %g \n", sum_of_distance);
        }

        for (c=0; c < num_classes; c++)
        {
          /* borrow distance2 here */
          distance2 = MRIFvox(mri_mem[c], x, y, z)/sum_of_distance;
          MRIFvox(mri_mem[c], x, y, z) = distance2;

          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("memship(%d) =%g \n", c, distance2);
          }

        }
      }
#endif

  printf("Edit membership to enhance CSF \n");

  mean1 = centroids[0][0]; //CSF
  std1 = 0.5*F[0]->rptr[1][1]; //scale variance by 2
  mean2 = centroids[2][0]; //GM
  std2 = F[2]->rptr[1][1];
  printf("CSF mean = %g, 1/std^2 = %g\n", mean1, std1);
  printf("GM mean = %g, 1/std^2 = %g\n", mean2, std2);

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }
        if (MRIFvox(mri_mem[0], x, y, z) > MRIFvox(mri_mem[2], x, y, z))
        {
          continue;  //already CSF
        }

        value = MRIFvox(mri_flash[0], x, y, z);
        if (value < mean1 || value > mean2)
        {
          continue;
        }
        //The value is between CSF and GM
        denom1 = (value - mean1);
        denom1 = 1.0/(denom1*denom1*std1 + 1e-30);

        denom2 = (value - mean2);
        denom2 = 1.0/(denom2*denom2*std2 + 1e-30);

        //membership for CSF
        scale = denom1/(denom2 + denom1);

        if (scale < 0.5)
        {
          continue;
        }

        scale = (scale - 0.5)/0.5;

        MRIFvox(mri_mem[2], x, y, z) *= exp( -scale);
      }


  if (!no_INU)
  {
    /* Output INU-corrected volume */
    for (i=0; i <nvolumes_total; i++)
    {
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            MRIFvox(mri_flash_backup[i],x,y,z) /= (MRIFvox(mri_mult[i],x,y,z) + 1e-20);
          }

      sprintf(fname,"%scorrected_input_%d.mgz", out_prefx, i);
      MRIwrite(mri_flash_backup[i], fname);

      sprintf(fname,"%schannel_%d_gain.mgz", out_prefx, i);
      MRIwrite(mri_mult[i], fname);
    }
  }

  /* How to map the multi-dimensional data to 1D using the class centers and
   * distance may worth some further investigation.
   * May need to optimize the location of the 1D centroids so as
   * to minimize the distance distortion. Or can I set the 1D centroid,
   * and only concerns about the mapping?
   */
  /* Will the problem becomes easier if my aim is to generate two volumes
   * one for depicting GM/WM contrast, one for GM/CSF contrast? should be!
   * Note that the centroid has no ordering in the nD space, but in 1D, they
   * clearly have an ordering!
   */
  /* Maybe the mapping is not important. Important is that the clustering is
   * good. If it is, then most voxels will be close to one centroid! So the
   * mapping seems indeed not important! If the clustering is bad, most
   * pixel will be fuzzy (membership close to 0.5), then no matter how to
   * choose the map. The synthesized image will have poor contrast!
   */


  /* Update centroids using first image and use it in corting */
  update_centroids1D(mri_flash[0], mri_mem, mri_mask, centroids1D, num_classes);

  for (i=0; i < num_classes; i++)
  {
    NRarray[i+1] = centroids1D[i];
    printf("centroids1D[%d] = %g\n",i, centroids1D[i]);
  }

  /* NR sort in ascending order */
  memset(NRindex, 0, MAX_CLASSES+1);
  indexx(num_classes, NRarray, NRindex);

  printf("Sorted centroids\n");
  for (i=1; i <= num_classes; i++)
  {
    printf("NRindex[%d] = %ld\n", i, NRindex[i]);
    indexmap[NRindex[i]-1] = i;
  }
  for (i=0; i < num_classes; i++)
  {
    printf("indexmap[%d] = %d\n", i, indexmap[i]);
  }

  if (num_classes == 4)
  {
    /* compute size as number of voxels adjacent to background */
    /* Will use that size to tell which is dura */
    for (c=0; c < num_classes; c++)
    {
      classSize[c] = 0;
    }
    for (z=1; z < (depth-1); z++)
      for (y=1; y< (height-1); y++)
        for (x=1; x < (width-1); x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }
          /* find hard-seg */
          label = 0;
          for (c=1; c < num_classes; c++)
          {
            if (MRIFvox(mri_mem[c], x, y, z) > MRIFvox(mri_mem[label],x,y,z))
            {
              label = c;
            }
          }

          for (index = 0; index < 6; index++)
          {
            xn = x + xoff[index];
            yn = y + yoff[index];
            zn = z + zoff[index];
            if (MRIvox(mri_mask, xn, yn, zn) == 0)
            {
              classSize[label] += 1.0;
            }
          }
        }

    if (classSize[NRindex[2]] > classSize[NRindex[1]])
    {
      printf("Switch to make sure dura has more background neighbors, and as final label 1 \n");
      i = NRindex[1] -1;
      j = NRindex[2] - 1;
      c = NRindex[1];
      NRindex[1] = NRindex[2];
      NRindex[2] = c;
      indexmap[i] = 2;
      indexmap[j] = 1;
    }
  }

  if (hard_segmentation || (hseg_fname != NULL))
  {
    MRI *mri_seg = NULL;

    mri_seg = MRIcopy(mri_mask, mri_seg);

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }
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

          MRIvox(mri_seg, x, y, z) = indexmap[i];

        }

    if (hseg_fname == NULL)
    {
      sprintf(fname,"%sEM_hseg.mgz", out_prefx);
      MRIwrite(mri_seg, fname);
    }
    else
    {
      MRIwrite(mri_seg, hseg_fname);
    }

    MRIfree(&mri_seg);
  }

  if (1)
  {
    MRI *mri_seg = NULL;

    mri_seg = MRIcopy(mri_mask, mri_seg);

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }
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

          j = indexmap[i];
          if (j == num_classes)
          {
            MRIvox(mri_seg, x, y, z) = 110;
          }
          else
          {
            MRIvox(mri_seg, x, y, z) = 0;
          }

        }

    sprintf(fname,"%sWM_seg.mgz", out_prefx);
    MRIwrite(mri_seg, fname);
    MRIfree(&mri_seg);
  }

  /* Use the INU_corrected volume for synthesization */
  if (nvolumes_total == 2 && num_classes == 4)
  {
    mri_tmp = MRIclone(mri_flash_backup[0], NULL);
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          //   if(MRIvox(mri_brainmask, x, y, z) == 0 || MRIvox(mri_mask,x,y,z) == 0){
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            MRIsetVoxVal(mri_tmp, x, y, z, 0, 0);
          }
          else
          {
            scale = MRIFvox(mri_flash_backup[1],x,y,z)/(MRIFvox(mri_flash_backup[0],x,y,z) + 1e-10);
            /*     if(scale < 1.2)
             value =  0.9527*MRIFvox(mri_flash_backup[0],x,y,z) - 0.3039*(1 + 0.5*MRIFvox(mri_mem[1], x, y, z))* MRIFvox(mri_flash_backup[1],x,y,z);
             else
             value =  0.9527*MRIFvox(mri_flash_backup[0],x,y,z) - 0.3039*MRIFvox(mri_flash_backup[1],x,y,z);
            */

            value =  weight30*MRIFvox(mri_flash_backup[0],x,y,z) + weight5 * MRIFvox(mri_flash_backup[1],x,y,z);
            MRIsetVoxVal(mri_tmp, x, y, z, 0, value);
          }
        }

    /* Now scale mri_tmp to [0, 255], and with WM around 160 */
    mean1 = 0;
    std1 = 0;
    denom1 = 0;
    mean2 = 0;
    std2 = 0;
    denom2 = 0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask,x,y,z) == 0)
          {
            continue;
          }
          value = MRIFvox(mri_mem[3],x,y,z);
          mean1 += MRIgetVoxVal(mri_tmp,x,y,z,0)*value;
          denom1 += value;

          value = MRIFvox(mri_mem[0],x,y,z); //csf
          mean2 += MRIgetVoxVal(mri_tmp,x,y,z,0)*value;
          std2 += MRIgetVoxVal(mri_tmp,x,y,z,0)*MRIgetVoxVal(mri_tmp,x,y,z,0)*value;
          denom2 += value;
        }
    mean1 /= (denom1 + 1e-30);
    printf("WM mean in the synthesized volume before scaling is %g\n", mean1);

    mean2 /= (denom2 + 1e-30);
    std2 /= (denom2 + 1e-30);

    std2 = sqrt(std2 - mean2*mean2);
    printf("CSF mean and std in the synthesized volume before scaling are %g and %g\n", mean2, std2);

    threshold = mean2 - 2.0*std2;
    if (mean1 < threshold)
    {
      printf("Something is wrong: WM mean is less than lowest CSF in synthesized image\n");
    }

    mri_tmp2 =  MRIalloc(width, height, depth, MRI_UCHAR);
    MRIcopyHeader(mri_flash_backup[0], mri_tmp2);

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask,x,y,z) == 0)
          {
            MRIvox(mri_tmp2, x, y, z) = 0;
            continue;
          }
          value =  MRIgetVoxVal(mri_tmp,x,y,z,0);
          value = (value -threshold)*160/(mean1 - threshold);
          if (value < -0.4 )
          {
            value = -0.4;  //this is wrong! Many places are negative
          }
          if (value > 254.4)
          {
            value = 254.4;
          }
          MRIvox(mri_tmp2, x, y, z) = (value + 0.5);
        }

    if (synth_fname == NULL)
    {
      MRIwrite(mri_tmp2, "synth.mgz");
    }
    else
    {
      MRIwrite(mri_tmp2, synth_fname);
    }

    MRIfree(&mri_tmp);
    MRIfree(&mri_tmp2);
  }

  if (!ldaflag)
  {
    /* synthesize the 1D image */
    if (num_classes == 3)
    {
      j = NRindex[0+1] - 1;
      if (invert)
      {
        centroids1D[j] = 120;
      }
      else
      {
        centroids1D[j] = 30;
      }


      j = NRindex[1+1] - 1;
      centroids1D[j] = 75;
      j = NRindex[2+1] - 1;
      if (invert)
      {
        centroids1D[j] = 30;
      }
      else
      {
        centroids1D[j] = 120;
      }
    }
    else
    {
      for (i=0; i < num_classes; i++)
      {
        j = NRindex[i+1] - 1;
        centroids1D[j] = 15.0 + i*(225.0 - 15.0)/(num_classes-1.0);
      }
    }

    if (rescale)
    {
      /* Nonlinear scale of the membership functions */
      for (z=0; z < depth; z++)
        for (y=0; y< height; y++)
          for (x=0; x < width; x++)
          {
            if (MRIvox(mri_mask, x, y, z) == 0)
            {
              continue;
            }

            /* Scale by a sin-function */
            for (i=0; i < num_classes; i++)
            {
              value = MRIFvox(mri_mem[i], x, y, z);
              MRIFvox(mri_mem[i], x, y, z) = 0.5*(sin((value-0.5)*3.14159265) + 1.0);
            }

            /* Renormalize to sum up to 1 */
            value = 0.0;
            for (i=0; i < num_classes; i++)
            {
              value += MRIFvox(mri_mem[i], x, y, z);
            }

            for (i=0; i < num_classes; i++)
            {
              MRIFvox(mri_mem[i], x, y, z) /= (value + 0.0000000001);
            }
          }
    }

    /* borrow mri_mask for the output image */
    /* Linear interpolation will blur things, better performed in
     * float image. So the conformal switch should be one!
     */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          if (remove_cerebellum)
          {
            if (MRIvox(cerebellum_mask, x, y,z) == 1)
            {
              MRIvox(mri_mask, x, y, z) = (BUFTYPE) 0;
              continue;
            }
          }

          value = 0.0;
          for (i=0; i < num_classes; i++)
          {
            value += MRIFvox(mri_mem[i], x, y, z)*centroids1D[i];
          }

          MRIvox(mri_mask, x, y, z) = (BUFTYPE) value;
        }

    /* Output synthesized volume */
    sprintf(fname,"%sEM_combined.mgz", out_prefx);
    MRIwrite(mri_mask, fname);
  }
  else
  {
    /* Using LDA method for synthesize image */
    printf("Start LDA processing ...\n");
    /* map the output class ID to the true class ID */
    class1 = NRindex[class1] - 1;
    class2 = NRindex[class2] - 1;
    if (clear_dura)
    {
      class_dura = NRindex[1] - 1;
    }
    else
    {
      class_dura = -1;
    }

    printf("class1 = %d, class2 = %d, class_dura = %d\n", class1, class2, class_dura);

    LDAmean1 = (float *)malloc(nvolumes_total*sizeof(float));
    LDAmean2 = (float *)malloc(nvolumes_total*sizeof(float));
    LDAweight = (float *)malloc(nvolumes_total*sizeof(float));

    if (fuzzy_lda == 0)
    {
      /* Compute class means */
      update_LDAmeans(mri_flash, mri_mem, mri_mask, LDAmean1, LDAmean2, nvolumes_total, class1, class2, MEMTHRESHOLD); /* 0.7 is a threshold, as to what voxels will be counted as in class1 */
      printf("class means computed \n");

      /* Compute Fisher's LDA weights */
      computeLDAweights(LDAweight, mri_flash, mri_mem, mri_mask, LDAmean1, LDAmean2, nvolumes_total, class1, class2, MEMTHRESHOLD);
    }
    else
    {
      /* Use fuzzy covariance matrix and centroids to compute LDA weights */
      update_F(F, mri_flash, mri_mem, mri_mem_old, mri_mult, mri_mask, centroids, nvolumes_total, num_classes);

      fuzzyLDAweights(LDAweight, F, centroids, class1, class2, 1, 1, nvolumes_total);

    }

    printf("LDA weights are: \n");
    for (i=0; i < nvolumes_total; i++)
    {
      printf("%g ", LDAweight[i]);
    }
    printf("\n");
    /* linear projection of input volumes to a 1D volume */
    min_val = 10000.0;
    max_val = -10000.0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          value = 0.0;

          if (clear_dura == 1 &&  MRIFvox(mri_mem[class_dura], x, y, z) > 0.25)
          {
            /* value = 0.0; */ /* Note that 0 is not minimum value!
              So this won't work */
            MRIvox(mri_mask, x , y, z) = 0;
            continue; /* This agrees with later processing */
          }
          else
          {
            for (i=0; i < nvolumes_total; i++)
            {
              value += MRIFvox(mri_flash[i], x, y, z)*LDAweight[i];
            }

          }

          if (max_val < value)
          {
            max_val = value;
          }
          if (min_val > value)
          {
            min_val = value;
          }

          /* Borrow mri_flash[0] to store the float values first */
          MRIFvox(mri_flash[0], x, y, z) = value;
        }

    printf("max_val = %g, min_val = %g \n", max_val, min_val);

    /* Scale output to [0, 255] */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (whole_volume == 0 && MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          value = (MRIFvox(mri_flash[0], x, y, z) - min_val)*255.0/(max_val - min_val) + 0.5; /* +0.5 for round-off */

          if (value > 255.0)
          {
            value = 255.0;
          }
          if (value < 0)
          {
            value = 0;
          }

          /* Borrow mri_flash[0] to store the float values first */
          MRIvox(mri_mask, x, y, z) = (BUFTYPE) value;
        }

    /* Output synthesized volume */
    sprintf(fname,"%sLDA_combined.mgz", out_prefx);
    MRIwrite(mri_mask, fname);

    free(LDAmean1);
    free(LDAmean2);
    free(LDAweight);
  }

  if (SYNTH_ONLY == 0)
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
            if (value > 255.0)
            {
              value = 255;
            }
            if (value < 0.0)
            {
              value = 0;
            }

            MRIvox(mri_tmp,x,y,z) = (unsigned char)value;
          }

      sprintf(fname,"%sEM_mem%d.mgz", out_prefx, indexmap[c]);
      MRIwrite(mri_tmp, fname);

    }

    MRIfree(&mri_tmp);
  }

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("parameter estimation took %d minutes and %d seconds.\n", minutes, seconds) ;

  MRIfree(&mri_mask);

  if (mri_label != NULL)
  {
    MRIfree(&mri_label);
  }

  for (i=0; i < num_classes; i++)
  {
    MRIfree(&mri_mem[i]);
    MRIfree(&mri_prior[i]);
    MRIfree(&mri_mem_old[i]);
    free(centroids[i]);
    MatrixFree(&F[i]);
  }

  for (i=0; i < nvolumes_total; i++)
  {
    MRIfree(&mri_flash[i]);
    MRIfree(&mri_mult[i]);
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
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "no_INU"))
  {
    no_INU = 1 ;
    printf("Do not perform INU correction\n") ;
  }
  else if (!stricmp(option, "fuzzy_lda"))
  {
    fuzzy_lda = 1 ;
    printf("Using fuzzy LDA weighting scheme\n") ;
  }
  else if (!stricmp(option, "remove_cerebellum"))
  {
    remove_cerebellum = 1 ;
    printf("Will remove cerebellum in the final synthesized image \n") ;
  }
  else if (!stricmp(option, "clear_dura"))
  {
    clear_dura = 1 ;
    printf("Remove voxels belongs to second class\n") ;
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
  else if (!stricmp(option, "synthonly"))
  {
    SYNTH_ONLY = 1;
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
  else if (!stricmp(option, "noprior"))
  {
    noprior = 1;
    printf("Don't use prior in computing memberships.\n") ;
  }
  else if (!stricmp(option, "synth_f"))
  {
    synth_fname = argv[2];
    nargs = 1;
    printf("Output synthesized image to file %s\n", synth_fname) ;
  }
  else if (!stricmp(option, "hseg_f"))
  {
    hseg_fname =  argv[2];
    nargs = 1;
    printf("Output hard segmentation to file %s\n", hseg_fname) ;
  }
  else if (!stricmp(option, "T1_PD"))
  {
    T1_PD = 1;
    printf("Em in parameter space\n") ;
  }
  else  if (!stricmp(option, "invert"))
  {
    invert = 1;
    printf("invert combined image contrast\n");
  }
  else if (!stricmp(option, "kappa"))
  {
    kappa = atof(argv[2]);
    nargs =1;
    printf("Typicality parameter kappa = %g.\n", kappa) ;
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
  else if (!stricmp(option, "beta"))
  {
    mybeta = atof(argv[2]);
    printf("weight for MRF = %g\n", mybeta);
    nargs = 1;
  }
  else if (!stricmp(option, "regularize"))
  {
    regularize = 1;
    lambda = atof(argv[2]);
    nargs =1;
    printf("Regularize the covariance matrix, lambda = %g\n", lambda);
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
  else switch (toupper(*option))
    {
    case 'M':
      num_classes = atoi(argv[2]) ;
      nargs = 1 ;
      printf("Number of classes=%d\n", num_classes) ;
      if (num_classes > MAX_CLASSES)
      {
        ErrorExit(ERROR_BADPARM, "%s: too many desired classes.\n", Progname);
      }
      break ;
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
  printf("\t -gca fname to specify the atlas name\n");
  printf("\t -xform fname to specify the atlas registration file name\n");
  printf("\t -M # to set the number of classes \n");
  printf("\t -E # to set the convergence tolerance \n");
  printf("\t -R # to set the max number of iterations \n");
  printf("\t -T # to set the threshold for background \n");
  printf("\t -norm to normalize input volumes before clustering \n");
  printf("\t -conform to conform input volumes (brain mask typically already conformed) \n");
  printf("\t -noconform to prevent conforming to COR \n");
  printf("\t -synthonly to prevent output membership functions \n");
  printf("\t -synth_f %%s to output synthesized volume to a file \n");
  printf("\t -hseg_f %%s to output hard segmentation for a file \n");
  printf("\t -striponly %%s to output computed brain mask to a file \n");
  printf("\t -noprior don't use prior in computing memberships \n");
  printf("\t -remove_cerebellum will remove cerebellum in the synth output \n");
  exit(code) ;

}

double NbhdLikelihood(MRI **mri_mem, MRI *mri_mask, int x, int y, int z, int label, int num_classes)
{
  int c, index, nx, ny, nz, nlabel;
  int width, height, depth;

  double p, maxp, tmpp, currentp;
  int total, diff;

  width = mri_mask->width;
  height = mri_mask->height;
  depth = mri_mask->depth;

  currentp = MRIFvox(mri_mem[label], x, y, z);

  total = 0;
  diff = 0;
  p = 0;
  for (index = 0; index < 6; index++)
  {
    nx = x + xoff[index];
    ny = y + yoff[index];
    nz = z + zoff[index];

    if (nx < 0 || nx >= width || ny < 0 || ny >= height || nz < 0 || nz >=depth)
    {
      continue;
    }
    if (MRIvox(mri_mask, nx, ny, nz) == 0)
    {
      continue;
    }

    total++;

    if (1)
    {
      /* find the hard-label for this neighbor */
      nlabel = 0;
      maxp = MRIFvox(mri_mem[0], nx, ny, nz);
      for (c = 1; c < num_classes; c++)
      {
        tmpp = MRIFvox(mri_mem[c], nx, ny, nz);
        if (tmpp > maxp)
        {
          maxp = tmpp;
          nlabel = c;
        }
      }

      if (nlabel != label)
      {
        diff++;
      }
    }
    else     /* soft penalty, will tie more closely to a smoothing of probability itself*/
    {
      /* This will make all p equal, and the hardsegmentation is more like noise */
      tmpp = currentp -  MRIFvox(mri_mem[label], nx, ny, nz);
      if (tmpp < 0)
      {
        p -= tmpp;
      }
      else
      {
        p += tmpp;
      }
    }
  }

  p = exp(-mybeta*diff/(total + 1e-10));

  return p;
}


void  update_LDAmeans(MRI **mri_flash, MRI **mri_mem, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2, float threshold)
{
  /* maybe I should design a Fuzzy LDA !! */

  int m, x, y, z, depth, height, width;
  double numer1, numer2, denom1, denom2;
  double mem;

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
            mem = MRIFvox(mri_mem[classID1], x, y, z);
            if (mem >= threshold)
            {
              numer1 += MRIFvox(mri_flash[m], x, y, z);
              denom1 += 1.0;
            }

            mem = MRIFvox(mri_mem[classID2], x, y, z);
            if (mem >= threshold)
            {
              numer2 += MRIFvox(mri_flash[m], x, y, z);
              denom2 += 1.0;
            }

          }

        }

    LDAmean1[m] = numer1/(denom1 + 0.0001);
    LDAmean2[m] = numer2/(denom2 + 0.0001);

  }

  return;

}

void fuzzyLDAweights(float *weights, MATRIX **F, float **centroids, int classID1, int classID2, float size1, float size2,  int nvolumes_total)
{


  int m1, m2;
  double denom, sumw;

  MATRIX *InvSW, *SW;

  SW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] = F[classID1]->rptr[m1][m2] *size1 + F[classID2]->rptr[m1][m2] * size2;    /* index starts from 1 for matrix */
    }

    /* SW->rptr[m1][m1] += 0.00001; */ /* prevent SW to be singular */
  }

  /* Compute inverse of SW */
  InvSW = MatrixInverse(SW, NULL);

  if (InvSW == NULL)   /* inverse doesn't exist */
  {
    ErrorExit(ERROR_BADPARM, "%s: singular fuzzy covariance matrix.\n", Progname);
  }

  /* Compute weights */
  denom = 0.0;
  sumw = 0.0;
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    weights[m1-1]= 0.0;
    for (m2=1; m2 <= nvolumes_total; m2++)
    {
      weights[m1-1] += InvSW->rptr[m1][m2] *(centroids[classID1][m2-1] - centroids[classID2][m2-1]);
    }
    sumw += weights[m1-1];
    denom += weights[m1-1]*weights[m1-1];
  }

  denom = sqrt(denom + 0.0000001);
  /* Normalized weights to have norm 1 */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    if (sumw > 0)
    {
      weights[m1-1] /= denom;
    }
    else
    {
      weights[m1-1] /= -denom;
    }
  }

  MatrixFree(&InvSW);

  MatrixFree(&SW);


  return;

}


void computeLDAweights(float *weights, MRI **mri_flash, MRI **mri_mem, MRI *mri_mask, float *LDAmean1, float *LDAmean2, int nvolumes_total, int classID1, int classID2, float threshold)
{

  int m1, m2, x, y, z, depth, height, width;
  double denom, sumw;
  double data1, data2;

  double Mdistance;

  MATRIX *InvSW, *SW;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  SW = (MATRIX *)MatrixAlloc(nvolumes_total, nvolumes_total, MATRIX_REAL);

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] = 0.0; /* index starts from 1 for matrix */
    }
  }

  /* printf("SW matrix initialized \n"); */
  denom = 0.0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }

        if (MRIFvox(mri_mem[classID1], x, y, z) < threshold &&
            MRIFvox(mri_mem[classID2], x, y, z) < threshold)
        {
          continue;
        }

        denom +=  1.0;

        if (MRIFvox(mri_mem[classID1], x, y, z) >= threshold)
        {
          for (m1=0; m1 < nvolumes_total; m1++)
          {
            data1 = MRIFvox(mri_flash[m1], x, y, z) - LDAmean1[m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              data2 = MRIFvox(mri_flash[m2], x, y, z) - LDAmean1[m2];
              SW->rptr[m1+1][m2+1] += data1*data2;
            }
          }
        }
        else
        {
          for (m1=0; m1 < nvolumes_total; m1++)
          {
            data1 = MRIFvox(mri_flash[m1], x, y, z) - LDAmean2[m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              data2 = MRIFvox(mri_flash[m2], x, y, z) - LDAmean2[m2];
              SW->rptr[m1+1][m2+1] += data1*data2;
            }
          }
        }

      } /* for all data points */

  if (denom <= 0.0)
  {
    ErrorExit(ERROR_BADPARM, "%s: overflow in computing fuzzy covariance matrix.\n", Progname);
  }

  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    for (m2=m1; m2 <= nvolumes_total; m2++)
    {
      SW->rptr[m1][m2] /= denom;
      SW->rptr[m2][m1] = SW->rptr[m1][m2];
    }

    SW->rptr[m1][m1] += 0.00001; /* prevent SW to be singular */
  } /* for m1, m2 */


  /* Compute inverse of SW */
  InvSW = MatrixInverse(SW, NULL);

  if (InvSW == NULL)   /* inverse doesn't exist */
  {
    ErrorExit(ERROR_BADPARM, "%s: singular fuzzy covariance matrix.\n", Progname);
  }

  /* Compute weights */
  denom = 0.0;
  sumw = 0.0;
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

  if (1)
  {
    Mdistance = 0.0;
    for (m1=0; m1 < nvolumes_total; m1++)
    {
      Mdistance += weights[m1]*(LDAmean1[m1] - LDAmean2[m1]);
    }

    printf("Mahalanobis Distance between the two classes in the original feature space is %g\n", Mdistance);

  }

  denom = sqrt(denom + 0.0000001);
  /* Normalized weights to have norm 1 */
  for (m1=1; m1 <= nvolumes_total; m1++)
  {
    if (sumw > 0)
    {
      weights[m1-1] /= denom;
    }
    else
    {
      weights[m1-1] /= -denom;
    }
  }

  MatrixFree(&InvSW);

  MatrixFree(&SW);

  return;
}

void update_centroids(MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI **mri_mult, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes)
{
  /* This step stays the same as FCM, just that the membership function is now the probability function. No, not the same! No more power of membership
   */
  int m, c, x, y, z, depth, height, width;
  double numer, denom;
  double mem, data;
  double scale =1.0;
  //  double kappa = exp(-4.5);
  double gain_term;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;

  //  printf("kappa = %g \n", kappa);

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
              /* Use this "typicallity scale leads to poor results */
              scale =  MRIFvox(mri_lihood[c], x, y, z);
              scale = scale/(scale + kappa);
              mem = MRIFvox(mri_mem[c], x, y, z)*scale;
              gain_term = MRIFvox(mri_mult[m], x, y, z);
              data = MRIFvox(mri_flash[m], x, y, z);
              numer += mem*data*gain_term;
              denom += mem*gain_term*gain_term;
            }

          }
      if (denom != 0.0)
      {
        centroids[c][m] = numer/denom;
      }
      else
      {
        ErrorExit(ERROR_BADPARM, "%s: overflow in computing centroids.\n", Progname);
      }
    }
  }


  return;
}

double distancems(MRI **mri_flash, MRI**mri_mult, float *centroids, MATRIX *F, double detF, int x, int y, int z)
{
  /* F would be the inverse of the covariance matrix of the class */
  int i, j, rows;
  double data1, data2;

  double mydistance = 0;

  rows = F->rows; /* actually the data dimensions */

  for (i=0; i < rows; i++)
  {
    data1 =  MRIFvox(mri_flash[i], x, y, z) -  MRIFvox(mri_mult[i], x, y, z)*centroids[i];
    for (j=0; j< rows; j++)
    {
      data2 =  MRIFvox(mri_flash[j], x, y, z) -  MRIFvox(mri_mult[j], x, y, z)*centroids[j];

      mydistance += data1*data2*F->rptr[i+1][j+1];
    }

  }

  /* No multiply by detF in the EM algorithm*/
  /* mydistance *= detF; */

  return mydistance;
}


double distancems_noF(MRI **mri_flash, MRI**mri_mult, float *centroids, MATRIX *F, double detF, int x, int y, int z)
{
  /* F would be the inverse of the covariance matrix of the class */
  int i, rows;
  double data1;

  double mydistance = 0;

  rows = F->rows; /* actually the data dimensions */

  for (i=0; i < rows; i++)
  {
    data1 =  MRIFvox(mri_flash[i], x, y, z) -  MRIFvox(mri_mult[i], x, y, z)*centroids[i];
    mydistance += data1*data1*F->rptr[i+1][i+1];
  }


  /* No multiply by detF in the EM algorithm*/
  /* mydistance *= detF; */

  return mydistance;
}

double distancems2(MRI **mri_flash, MRI**mri_mult, float *centroids, MATRIX *F, double detF, int x, int y, int z, int channel)
{
  //This is actually a 1D distance, using only one channel
  /* F would be the inverse of the covariance matrix of the class */
  double data1;

  double mydistance;

  data1 =  MRIFvox(mri_flash[channel], x, y, z) -  MRIFvox(mri_mult[channel], x, y, z)*centroids[channel];

  mydistance = data1*data1*F->rptr[channel+1][channel+1];

  return mydistance;
}


void compute_bias(MRI **mri_mult, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI *mri_mask, float **centroids, MATRIX **F, int nvolumes_total, int num_classes)
{
  /* F should be the inverse of the covariance matrix */
  MRI *g_term, *h_term, *f_term;
  int c, i, x, y, z, depth, height, width;
  float tmp_h, tmp_f;
  float max_b, min_b;
  double sum_h, sum_f;
  double scaler, scale;

  depth = mri_flash[0]->depth;
  width = mri_flash[0]->width;
  height = mri_flash[0]->height;


  g_term = MRIclone(mri_flash[0], NULL);
  h_term = MRIclone(mri_flash[0], NULL);
  f_term = MRIclone(mri_flash[0], NULL);

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        MRIFvox(g_term, x, y, z) = lap_weight;
      }

  for (i=0; i < nvolumes_total; i++)
  {
    /* construct the bias field PDE coefficients */
    sum_h = 0;
    sum_f = 0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) <= 0)
          {
            MRIFvox(h_term, x, y, z) = 0;
            MRIFvox(f_term, x, y, z) = 0;
          }
          else
          {
            tmp_h = 0;
            tmp_f = 0;

            for (c=0; c < num_classes; c++)
            {
              scale =  MRIFvox(mri_lihood[c], x, y, z);
              scale = scale/(scale + kappa);
              if (T1_PD == 0)
              {
                tmp_h += centroids[c][i] * centroids[c][i]*F[c]->rptr[i+1][i+1] * MRIFvox(mri_mem[c], x, y, z)*scale;
              }
              else
              {
                tmp_h += centroids[c][i] * centroids[c][i] * MRIFvox(mri_mem[c], x, y, z)*scale;
              }
              if (T1_PD == 0)
              {
                tmp_f += MRIFvox(mri_flash[i], x, y, z) * centroids[c][i]*F[c]->rptr[i+1][i+1] * MRIFvox(mri_mem[c], x, y, z)*scale;
              }
              else
              {
                tmp_f += MRIFvox(mri_flash[i], x, y, z) * centroids[c][i] * MRIFvox(mri_mem[c], x, y, z)*scale;
              }

            }

            /* If I swicth the sign below, it doesn't converge anymore */
            MRIFvox(h_term, x, y, z) =  tmp_h*0.0001;
            MRIFvox(f_term, x, y, z) = -tmp_f*0.0001;

            sum_h += tmp_h*0.0001;
            sum_f += tmp_f*0.0001;
            if (debug_flag && x == Gx && y == Gy && z == Gz)
            {
              printf("h= %g, f =%g \n", tmp_h*0.1, -tmp_f*0.1);
            }
          }
        }

    printf("sum_f = %g, sum_h = %g\n", sum_f, sum_h);
    // scale so that mean gain field is 1
    // this normalization seems help faster convergence
    // why fantasm doesn't need this?? Maybe my weights for smoothing term is
    // still too small!!
#if 1
    scaler = sum_f/(sum_h + 1e-20);
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) <= 0)
          {
            continue;
          }
          MRIFvox(h_term, x, y, z) *= scaler;
        }
#endif

    //    if( i== 0){
    if (0)
    {
      MRIwrite(h_term, "h_term.mgz");
      MRIwrite(f_term, "f_term.mgz");
      MRIwrite(mri_mult[0], "init_mult.mgz");
    }
    printf("start multigrid Poisson solver ...\n");
    //    multigrid(mri_mult[i], f_term, g_term, h_term, width, height, depth);
    multigrid(mri_mult[i], f_term, h_term, width, height, depth);
    printf("multigrid solver converged\n");
    if (0)
    {
      MRIwrite(mri_mult[0], "res_mult.mgz");
      exit(0);
    }

    /* Correct for image intensity */
    max_b = -1000;
    min_b = 1000;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          tmp_h = MRIFvox(mri_mult[i],x,y,z);
          if (debug_flag && x == Gx && y == Gy && z == Gz)
          {
            printf("gain= %g \n", tmp_h);
          }

          if (max_b < tmp_h)
          {
            max_b = tmp_h;
          }
          if (min_b > tmp_h)
          {
            min_b = tmp_h;
          }
          //  MRIFvox(mri_flash[i],x,y,z) *= (1e-20 + tmp_h);
        }
    printf("max_gain = %g, min_gain = %g\n", max_b, min_b);
  }



  MRIfree(&g_term);
  MRIfree(&h_term);
  MRIfree(&f_term);
  return;
}

void update_F(MATRIX **F, MRI **mri_flash, MRI **mri_mem, MRI **mri_lihood, MRI **mri_mult, MRI *mri_mask, float **centroids, int nvolumes_total, int num_classes)
{

  int m1, m2, c, x, y, z, depth, height, width;
  double denom;
  double mem, data1, data2;
  double gain1, gain2;
  double scale = 1.0;
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

    //    printf("kappa = %g\n", kappa);
    denom = 0.0;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }
          /* This scale is necessary, but kappa cannot be too big!
           *  Setting scale to constant one will not converge
           */
          scale =  MRIFvox(mri_lihood[c], x, y, z);
          scale = scale/(scale + kappa);
          mem = MRIFvox(mri_mem[c], x, y, z)*scale;
          /* mem = mem*mem; */ /* here differs from FCM */
          denom +=  mem;

          for (m1=0; m1 < nvolumes_total; m1++)
          {
            gain1 = MRIFvox(mri_mult[m1], x, y, z);
            data1 = MRIFvox(mri_flash[m1], x, y, z) - gain1*centroids[c][m1];
            for (m2=m1; m2 < nvolumes_total; m2++)
            {
              gain2 = MRIFvox(mri_mult[m2], x, y, z);
              data2 = MRIFvox(mri_flash[m2], x, y, z) - gain2*centroids[c][m2];
              F[c]->rptr[m1+1][m2+1] += data1*data2*mem;
            }
          }

        } /* for all data points */

    if (denom <= 0.0)
    {
      ErrorExit(ERROR_BADPARM, "%s: overflow in computing fuzzy covariance matrix.\n", Progname);
    }

    for (m1=1; m1 <= nvolumes_total; m1++)
    {
      for (m2=m1; m2 <= nvolumes_total; m2++)
      {
        F[c]->rptr[m1][m2] /= denom;
        F[c]->rptr[m2][m1] = F[c]->rptr[m1][m2];
      }

      F[c]->rptr[m1][m1] += eps; /* prevent F to be singular */
    } /* for m1, m2 */

    if (regularize)
    {
      //for(m1=1; m1 <= nvolumes_total; m1++)
      //F[c]->rptr[m1][m1] += eps;  /* prevent F to be singular */

      tmpM = MatrixInverse(F[c], tmpM);
      if (tmpM == NULL)
      {
        continue;
      }

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
  {
    MatrixFree(&tmpM);
  }

  return;
}

void compute_detF(MATRIX **F, double *detF, int num_classes)
{

  int c, n;
  float tmpv;

  n = F[0]->rows;

  for (c=0; c < num_classes; c++)
  {

    tmpv = MatrixDeterminant(F[c]);

#if 0
    /* the following may lead to unexpected performance */
    if (tmpv < 0.000001)   /*singular */
    {
      for (i=1; i <= F[c]->rows; i++)
      {
        F[c]->rptr[i][i] += 0.01;  /* try to recover it */
      }

      tmpv = MatrixDeterminant(F[c]);
    }
#endif

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

    if (mTmp == NULL)   /* inverse doesn't exist */
    {
      ErrorExit(ERROR_BADPARM, "%s: singular fuzzy covariance matrix.\n", Progname);
    }

    for (row=1; row <= rows; row++)
    {
      memcpy((char *)(F[c]->rptr[row]), (char *)mTmp->rptr[row],
             (cols+1)*sizeof(float)) ;
    }

    MatrixFree(&mTmp);
    mTmp = NULL;
  }


  return;
}

void update_centroids1D(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int num_classes)
{

  int  c, x, y, z, depth, height, width;
  double numer, denom;
  double mem, data;

  depth = mri_flash->depth;
  width = mri_flash->width;
  height = mri_flash->height;

  /* Put the voxel iteration inside makes it easier to
   * program, but may take longer time. Otherwise, need to
   * declare numer and denom as matrices!
   */
  for (c = 0; c < num_classes; c++)
  {
    numer = 0;
    denom = 0;

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) > 0)
          {
            mem = MRIFvox(mri_mem[c], x, y, z);
            data = MRIFvox(mri_flash, x, y, z);
            /* mem = mem*mem; */ /* or powf(mem, q)  for q != 2 */
            numer += mem*data;
            denom += mem;
          }

        }
    if (denom != 0.0)
    {
      centroids1D[c] = numer/denom;
    }
    else
    {
      ErrorExit(ERROR_BADPARM, "%s: overflow in computing centroids.\n", Progname);
    }

  }


  return;
}

void MRI_EM(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int  num_classes)
{
  /* Simplified, assuming all class have equal size and variance 1 */
  /* 1D EM */
  int depth, width, height, x, y, z, c;
  float max_change = 100.0;
  float distance2, sum_of_distance;
  float *oldMems;
  float *sigmaI; /* variance */
  float *piI; /* class size */
  int total_num;


  depth = mri_flash->depth;
  width =  mri_flash->width;
  height =  mri_flash->height;

  oldMems = (float *)malloc(sizeof(float)*num_classes);
  sigmaI = (float *)malloc(sizeof(float)*num_classes);
  piI = (float *)malloc(sizeof(float)*num_classes);

  /* Initialize membership values from the given centroids. Just random */
  total_num = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }
        total_num++;

        for (c=0; c < num_classes; c++)
        {
          MRIFvox(mri_mem[c], x, y, z) = 1.0/num_classes;
        }
      }

  printf("total_num = %d\n", total_num);

  for (c=0; c < num_classes; c++)
  {
    sigmaI[c] = 1;
    piI[c] = 1.0/num_classes;
  }

  max_change = 100.0;
  while (max_change > 0.01)
  {
    max_change = 0.0;

    /* Update membership function */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          sum_of_distance = 1e-20;
          /* Compute distance */
          for (c=0; c < num_classes; c++)
          {
            /* record old membership values */
            oldMems[c] = MRIFvox(mri_mem[c], x, y, z);

            distance2 = MRIFvox(mri_flash, x, y, z) - centroids1D[c];

            distance2 = piI[c]*exp(-0.5*distance2*distance2/sigmaI[c])/(sqrt(sigmaI[c]) + 1e-30);

            MRIFvox(mri_mem[c], x, y, z) = distance2;

            sum_of_distance += distance2;
          }

          if (sum_of_distance <= 0.0)
          {
            ErrorExit(ERROR_BADPARM, "%s: overflow in computing membership function.\n", Progname);
          }

          for (c=0; c < num_classes; c++)
          {
            /* borrow distance2 here */
            distance2 = MRIFvox(mri_mem[c], x, y, z)/sum_of_distance;
            MRIFvox(mri_mem[c], x, y, z) = distance2;

            distance2 -= oldMems[c];
            if (distance2 < 0)
            {
              distance2 = -distance2;
            }

            if (max_change < distance2)
            {
              max_change = distance2;
            }

          }

        } /* end of all data points */


    printf("max_change = %g\n", max_change);

    /* Update centroids */
    update_centroids1D(mri_flash, mri_mem, mri_mask, centroids1D, num_classes);
    printf("Centroids: ");
    for (c=0; c < num_classes; c++)
    {
      printf(" %g,",centroids1D[c]);
    }
    printf("\n");

    /* Update class variance and priors */
    for (c=0; c < num_classes; c++)
    {
      sigmaI[c] = 0;
      piI[c] = 0;
    }

    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }
          for (c=0; c < num_classes; c++)
          {
            distance2 = MRIFvox(mri_flash, x, y, z) - centroids1D[c];
            sigmaI[c] += distance2*distance2*MRIFvox(mri_mem[c], x, y, z);
            piI[c] += MRIFvox(mri_mem[c], x, y, z);
          }
        }
    /* Update class-priors */
    for (c=0; c < num_classes; c++)
    {
      sigmaI[c] /= (piI[c] + 1e-30);
      piI[c] /= (float)(total_num + 1e-30);
    }

    printf("variance and class Size: \n");
    for (c=0; c < num_classes; c++)
    {
      printf("std[%d] = %g, piI[%d]=%g\n",c, sigmaI[c], c, piI[c]);
    }
    printf("\n");
  }/* end of while-loop */

  free(oldMems);
  free(sigmaI);
  free(piI);

  return;
}

void mri_invert_contrast_and_scale(MRI *mri, MRI *mri_mask, double scale)
{
  double max_v, min_v, value;
  int depth, width, height, x, y, z;

  depth = mri->depth;
  width =  mri->width;
  height =  mri->height;

  max_v =  -1e30;
  min_v = 1e30;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }
        value = MRIgetVoxVal(mri, x, y, z, 0);
        if (max_v < value)
        {
          max_v = value;
        }
      }

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          MRIsetVoxVal(mri, x, y, z, 0, 0);
          continue;
        }
        value = (max_v - MRIgetVoxVal(mri, x, y, z, 0))/scale;
        MRIsetVoxVal(mri, x, y, z, 0, value);
      }


  return;

}
void MRI_FCM(MRI *mri_flash, MRI **mri_mem, MRI *mri_mask, float *centroids1D, int  num_classes)
{
  /* 1D FCM */
  int depth, width, height, x, y, z, c;
  float max_change = 100.0;
  float distance2, sum_of_distance;
  float *oldMems;

  depth = mri_flash->depth;
  width =  mri_flash->width;
  height =  mri_flash->height;

  oldMems = (float *)malloc(sizeof(float)*num_classes);

  /* Initialize membership values from the given centroids. Just random */
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0)
        {
          continue;
        }

        for (c=0; c < num_classes; c++)
        {
          MRIFvox(mri_mem[c], x, y, z) = 1.0/num_classes;
        }
      }

  max_change = 100.0;
  while (max_change > 0.01)
  {
    max_change = 0.0;

    /* Update membership function */
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {
          if (MRIvox(mri_mask, x, y, z) == 0)
          {
            continue;
          }

          sum_of_distance = 0.0;
          /* Compute distance */
          for (c=0; c < num_classes; c++)
          {
            /* record old membership values */
            oldMems[c] = MRIFvox(mri_mem[c], x, y, z);

            distance2 = MRIFvox(mri_flash, x, y, z) - centroids1D[c];
            distance2 = distance2*distance2;

            if (distance2 == 0.0)
            {
              MRIFvox(mri_mem[c], x, y, z) = 100000000.0;
            }
            else
            {
              MRIFvox(mri_mem[c], x, y, z) = 1.0/distance2;
            }

            sum_of_distance += MRIFvox(mri_mem[c], x, y, z);
          }

          if (sum_of_distance <= 0.0)
          {
            ErrorExit(ERROR_BADPARM, "%s: overflow in computing membership function.\n", Progname);
          }

          for (c=0; c < num_classes; c++)
          {
            /* borrow distance2 here */
            distance2 = MRIFvox(mri_mem[c], x, y, z)/sum_of_distance;
            MRIFvox(mri_mem[c], x, y, z) = distance2;

            distance2 -= oldMems[c];
            if (distance2 < 0)
            {
              distance2 = -distance2;
            }

            if (max_change < distance2)
            {
              max_change = distance2;
            }

          }

        } /* end of all data points */

    printf("Centroids: ");
    for (c=0; c < num_classes; c++)
    {
      printf(" %g,",centroids1D[c]);
    }
    printf("\n");

    /* Update centroids */
    update_centroids1D(mri_flash, mri_mem, mri_mask, centroids1D, num_classes);

  }/* end of while-loop */

  free(oldMems);

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
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

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
  {
    mean = mean/total;
  }

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
  {
    variance = sqrt(variance/total);
  }
  else
  {
    variance = 1;
  }

  /* normalization: invert variance first to save time */
  variance = 1.0/variance;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++)
      {
        // tmpval = MRIFvox(mri_src, x, y, z) - mean;
        tmpval = MRIFvox(mri_src, x, y, z);
        MRIFvox(mri_dst, x, y, z) = tmpval*variance;
      }

  return (mri_dst);
}

/*
 * Indexes an array arr[], outputs the array indx[] such that arr[indx[j]]
 * is in ascending order for j = 1, 2, ... N.  The input quantities n and arr
 * are not changed.
 * That is, it produces a second array (indx[]) that contains pointers to the
 * elements of the original array (arr[]) in the order of their size.
 *
 * This routine is from Numerical Recipes for C.
 *
 * !!! NJS note: the guts of this routine have been removed, as we
 * cannot distribute NRC routines.  A suitable replacement will need to be
 * created by someone to get this code working again!
 */

void indexx(unsigned long n, float arr[], unsigned long indx[])
{
  ErrorExit(
    ERROR_BADPARM,
    "%s: ERROR: NRC routine indexx has been removed from the source.\n"
    "%s is not usable until a suitable replacement is found!\n",
    Progname, Progname);
  return;
}


void  set_equilavent_classes(int *equivalent_classes)
{
  int i;

  for (i=0; i < MAX_CLASSES; i++)
  {
    equivalent_classes[i] = i;
    // Need some clever string comparison to automatically build equivalent classes; for now I will just manually set
    //for(j=i+1; j < MAX_CLASSES; j++){

    //    }
  }

  //CSF 1, GM 3, WM 2
  equivalent_classes[0] = 1;
  equivalent_classes[1] = 1;
  equivalent_classes[2] = 2;
  equivalent_classes[3] = 3;
  equivalent_classes[4] = 1;
  equivalent_classes[5] = 1;
  equivalent_classes[6] = 0;
  equivalent_classes[7] = 2;
  equivalent_classes[8] = 3;
  equivalent_classes[9] = 3;
  equivalent_classes[10] = 3;
  equivalent_classes[11] = 3;
  equivalent_classes[12] = 3;
  equivalent_classes[13] = 3;
  equivalent_classes[14] = 1;
  equivalent_classes[15] = 1;
  equivalent_classes[16] = 2; //this is brain stem, really WM??
  equivalent_classes[17] = 3;
  equivalent_classes[18] = 3;
  equivalent_classes[19] = 0; //Left_Insula
  equivalent_classes[20] = 0; //Left_Operculum
  equivalent_classes[21] = 0; //Line_1
  equivalent_classes[22] = 0; //Line_2
  equivalent_classes[23] = 0; //Line_3
  equivalent_classes[24] = 1;
  equivalent_classes[25] = 0; //Left_lesion
  equivalent_classes[26] = 0; //left_accumbens_area
  equivalent_classes[27] = 0; //Left_Substancia_Nigra
  equivalent_classes[28] = 0; //Left_VentralDC;
  equivalent_classes[29] = 0; //left_undetermined
  equivalent_classes[30] = 0; // Left_vessel
  equivalent_classes[31] = 0; //left_choroid_plexus
  equivalent_classes[32] = 0; //lEFT_f3ORB
  equivalent_classes[33] = 0;
  equivalent_classes[34] = 0;
  equivalent_classes[35] = 0;
  equivalent_classes[36] = 0;
  equivalent_classes[37] = 0;
  equivalent_classes[38] = 0;
  equivalent_classes[39] = 0;
  equivalent_classes[40] = 1;
  equivalent_classes[41] =equivalent_classes[2];
  equivalent_classes[42] =equivalent_classes[3];
  equivalent_classes[43] =equivalent_classes[4];
  equivalent_classes[44] =equivalent_classes[5];
  equivalent_classes[45] =equivalent_classes[6];
  equivalent_classes[46] =equivalent_classes[7];
  equivalent_classes[47] =equivalent_classes[8];
  equivalent_classes[48] =equivalent_classes[9];
  equivalent_classes[49] =equivalent_classes[10];
  equivalent_classes[50] =equivalent_classes[11];
  equivalent_classes[51] =equivalent_classes[12];
  equivalent_classes[52] =equivalent_classes[13];
  equivalent_classes[53] =equivalent_classes[17];
  equivalent_classes[54] =equivalent_classes[18];
  equivalent_classes[55] =equivalent_classes[19];
  equivalent_classes[56] =equivalent_classes[20];
  equivalent_classes[57] =equivalent_classes[25];
  equivalent_classes[58] =equivalent_classes[26];
  equivalent_classes[59] =equivalent_classes[27];
  equivalent_classes[60] =equivalent_classes[28];
  equivalent_classes[61] =equivalent_classes[29];
  equivalent_classes[62] =equivalent_classes[30];
  equivalent_classes[63] =equivalent_classes[31]; //choroid_plexus
  equivalent_classes[64] =equivalent_classes[32];
  equivalent_classes[65] =equivalent_classes[33];
  equivalent_classes[66] =equivalent_classes[34];
  equivalent_classes[67] =equivalent_classes[35];
  equivalent_classes[68] =equivalent_classes[36];
  equivalent_classes[69] =equivalent_classes[37];
  equivalent_classes[70] =equivalent_classes[38];
  equivalent_classes[71] =equivalent_classes[39];
  equivalent_classes[74] =equivalent_classes[73];
  equivalent_classes[76] =equivalent_classes[75];
  equivalent_classes[79] =equivalent_classes[78];
  equivalent_classes[82] =equivalent_classes[81];
  equivalent_classes[84] =equivalent_classes[83];
  equivalent_classes[187] =equivalent_classes[186];

  return;
}
