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


/* Compute the mean and std of a labeled region */
/* manually chose which class-pairs to compute */

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

void usage(int exit_val);
MRI *fliplr(MRI *src);

const char *Progname;

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

static int MINLABEL = 2;
static int MAXLABEL = 39;

/* Compute CNR inside a window only */

static int debug_flag = 0;
static int window_flag = 0;
static int window_size = 30;

static char *fname = NULL;

void computeClassStats(float *LDAmean, float *variance, float *classsize, MRI *mri_src, MRI *mri_label, MRI *mri_mask, int classID);
MATRIX *ComputeAdjMatrix(MRI *mri_label, MRI *mri_mask, int minlabel, int maxlabel);
static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{

  char **av;
  MRI *mri_src, *mri_label, *mri_mask;
  int ac, nargs;
  int width, height, depth, x, y, z, i,j;
  FILE *fp;

  float *LDAmeans = NULL;
  float *classSize = NULL;
  MATRIX *AdjMatrix;
  float *Variances = NULL;
  int num_classes, label;

  double cnr;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_compute_all_CNR");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage(1);

  mri_src = MRIread(argv[1]) ;
  if (!mri_src)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;
  mri_label = MRIread(argv[2]) ;
  if (!mri_label)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume %s",
              Progname, argv[2]) ;


  if ((mri_src->width != mri_label->width) ||
      (mri_src->height != mri_label->height) ||
      (mri_src->depth != mri_label->depth))
    ErrorExit(ERROR_BADPARM, "%s: source (type %d) and label (type %d) volumes don't match (%d x %d x %d) vs (%d x %d x %d)\n",
              Progname, mri_src->type, mri_label->type, mri_src->width,
              mri_src->height, mri_src->depth,
              mri_label->width, mri_label->height, mri_label->depth) ;

  if (argc == 4)
    MINLABEL = atoi(argv[3]);
  if (argc == 5)
  {
    MAXLABEL = atoi(argv[4]);
  }


  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  printf("Clear cerebellum related labels  6, 7, 8 and 16, 24, 28\n");
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        label = (int)MRIgetVoxVal(mri_label, x, y,z,0);
        if ((label >= 6 && label <= 8) || label == 16 || label == 24 || label == 28)
          MRIsetVoxVal(mri_label, x, y, z, 0, 0);

      }


  mri_mask = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_UCHAR);
  MRIcopyHeader(mri_src, mri_mask);

  /* Simply set mask to be 1 everywhere */
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        MRIvox(mri_mask, x, y,z) = 1;
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

          if (z < (Gz - window_size) || z >(Gz + window_size)
              || y <(Gy - window_size) || y > (Gy + window_size)
              || x < (Gx - window_size) || x > (Gx + window_size))
            MRIvox(mri_mask, x, y,z) = 0;
        }

  }

  num_classes = MAXLABEL - MINLABEL + 1;
  printf("Total of %d classes considered in LDA training\n", num_classes);

  /* Allocate memory */
  LDAmeans = (float *)malloc(num_classes*sizeof(float));
  Variances = (float *)malloc(num_classes*sizeof(float));
  classSize = (float *)malloc(num_classes*sizeof(float));

  /* Note that the scatter matrix is just class-size*covariance matrix */
  for (i=0; i < num_classes; i++)
  {
    computeClassStats(&LDAmeans[i], &Variances[i], &classSize[i], mri_src, mri_label, mri_mask, MINLABEL + i);
  }

  if (MINLABEL <=2 && MAXLABEL >= 20)
  {

    AdjMatrix = (MATRIX *)MatrixAlloc(num_classes, num_classes, MATRIX_REAL);

    if (!AdjMatrix)
      ErrorExit(ERROR_BADPARM, "%s: unable to allowcate memory.\n", Progname);
    /* The diagnoal entries of AdjMatrix is set to zero and remain zero */
    for (i=1; i <= num_classes;i++)
      for (j=i; j <= num_classes; j++)
      {
        AdjMatrix->rptr[i][j] = 0.0;
        AdjMatrix->rptr[j][i] = 0.0;
      }

    AdjMatrix->rptr[2+1-MINLABEL][3+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][10+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][12+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][13+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][17+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[2+1-MINLABEL][18+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][12+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[3+1-MINLABEL][18+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[10+1-MINLABEL][11+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[12+1-MINLABEL][13+1-MINLABEL] = 1.0;
    AdjMatrix->rptr[17+1-MINLABEL][18+1-MINLABEL] = 1.0;
  }
  else
  {
    printf("Compute Co-occurance matrix\n");
    AdjMatrix = ComputeAdjMatrix(mri_label, mri_mask, MINLABEL, MAXLABEL);
  }
  if (fname != NULL)
  {
    fp = fopen(fname, "w");
  }
  else
    fp = 0;

  for (i=0; i <num_classes-1;i++)
    for (j=i+1; j < num_classes; j++)
    {
      if (AdjMatrix->rptr[i+1][j+1] < 0.5) continue;
      cnr = LDAmeans[i] - LDAmeans[j];
      cnr /= sqrt(Variances[i] + Variances[j] + 1e-15);
      if (cnr < 0) cnr = -cnr;
      if (fp)
        fprintf(fp,"%9.4f ",(float)cnr);

      printf("CNR of class %d and class %d is %g\n", i+MINLABEL, j+MINLABEL, cnr);
    }

  if (fp)
    fprintf(fp,"\n");
  if (fp)
  {
    fclose(fp);
  }

  MRIfree(&mri_src);
  MRIfree(&mri_label);

  free(LDAmeans);
  free(classSize);
  free(Variances);

  exit(0);

}  /*  end main()  */
void computeClassStats(float *LDAmean, float *variance, float *classsize, MRI *mri_src, MRI *mri_label, MRI *mri_mask, int classID)
{
  /* Compute LDAmean and SW from given data and label */

  int x, y, z, depth, height, width;
  double numer, denom, meanV;
  double data1;
  int label;

  depth = mri_src->depth;
  width = mri_src->width;
  height = mri_src->height;

  /* Compute Class mean */
  denom = 0.0;
  numer = 0.0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;

        label = (int) MRIgetVoxVal(mri_label, x, y, z,0);

        if (label != classID)
          continue;

        denom += 1;
        numer += MRIgetVoxVal(mri_src, x, y, z, 0);

      } /* for all data points */

  *classsize = denom;

  meanV = numer/(denom + 1e-15);

  *LDAmean = meanV;

  /* compute class variance */
  numer = 0.0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_mask, x, y, z) == 0) continue;

        label = (int) MRIgetVoxVal(mri_label, x, y, z,0);

        if (label != classID)
          continue;

        data1 =  MRIgetVoxVal(mri_src,x,y,z,0) - meanV;
        numer += data1*data1;
      } /* for all data points */

  *variance = numer/(denom + 1e-15);

  //  printf("mean=%g, variance = %g\n", meanV, *variance);
  return;
}


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
  else if (!stricmp(option, "f")
           || !stricmp(option, "fname"))
  {
    fname = argv[2];
    nargs = 1;
    printf("Output CNR values to the file %s\n", fname);
  }

  return(nargs) ;
}

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <in vol> <label vol> label [label2] \n", Progname);
  fprintf(fout, "this program computes the mean and std of a labeled region. \n") ;
  fprintf(fout, "If label2 is also specified, the program computes the CNR of both regions. \n") ;
  fprintf(fout, "Using '-f fname' to specify the file to store CNR values \n") ;

  exit(exit_val);

}  /*  end usage()  */

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

/*  EOF  */

