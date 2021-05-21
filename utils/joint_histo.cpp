/**
 * @brief 1d histogram utilities.
 *
 * Utilities for computing and analyzing 2 dimensional histograms.
 */

/**
 Not exactly corresponding to histo.h and histo.c!!
**/

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "histo.h"
#include "macros.h"
#include "mri.h"
#include "proto.h"

#include "joint_histo.h"

static int min(int a, int b) { return (a < b) ? a : b; }
static int max(int a, int b) { return (a > b) ? a : b; }
static float flmin(float a, float b) { return (a < b) ? a : b; }
static float flmax(float a, float b) { return (a > b) ? a : b; }

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int JHISTOfree(JOINT_HISTOGRAM **pjhisto)
{
  JOINT_HISTOGRAM *jhisto;

  jhisto = *pjhisto;
  *pjhisto = NULL;
  if (jhisto) {
    if (jhisto->counts) {
      free(jhisto->counts);
      jhisto->counts = NULL;
    }
    else
      DiagBreak();
    free(jhisto);
  }
  else
    DiagBreak();

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int JHISTOdump(JOINT_HISTOGRAM *jhisto, FILE *fp)
{
  int bin_no_1, bin_no_2;

  if (!jhisto)
    fprintf(stderr, "NULL joint histogram");
  else {
    fprintf(
        fp, "nbins_1 = %d, nbins_2 = %d, sample_count = %d\n", jhisto->nbins_1, jhisto->nbins_2, jhisto->sample_count);
    for (bin_no_1 = 0; bin_no_1 < jhisto->nbins_1; bin_no_1++)
      for (bin_no_2 = 0; bin_no_2 < jhisto->nbins_2; bin_no_2++)
        fprintf(fp, "%f ", jhisto->counts[bin_no_1 * jhisto->nbins_2 + bin_no_2]);
    fprintf(fp, "\n");
  }
  return (NO_ERROR);
}

int JHISTOwriteInto(JOINT_HISTOGRAM *h, FILE *fp)
{
  int b1, b2;

  fwriteInt(h->nbins_1, fp);
  fwriteInt(h->nbins_2, fp);

  // fwriteInt(h->total_bins, fp) ;
  fwriteInt(h->sample_count, fp);

  for (b1 = 0; b1 < h->nbins_1; b1++)
    for (b2 = 0; b2 < h->nbins_2; b2++) fwriteInt(h->counts[b1 * h->nbins_2 + b2], fp);

  return (NO_ERROR);
}

JOINT_HISTOGRAM *JHISTOreadFrom(FILE *fp)
{
  int b1, b2, nbins_1, nbins_2;  //, sample_count;
  JOINT_HISTOGRAM *jh;

  nbins_1 = freadInt(fp);
  nbins_2 = freadInt(fp);
  jh = JHISTOalloc(nbins_1, nbins_2);

  // total_bins = freadInt(fp) ;
  // sample_count =
  freadInt(fp);

  for (b1 = 0; b1 < jh->nbins_1; b1++)
    for (b2 = 0; b2 < jh->nbins_2; b2++) jh->counts[b1 * jh->nbins_2 + b2] = freadInt(fp);

  return (jh);
}

JOINT_HISTOGRAM *JHISTOalloc(int nbins_1, int nbins_2)
{
  JOINT_HISTOGRAM *jhisto;

  jhisto = (JOINT_HISTOGRAM *)calloc(1, sizeof(JOINT_HISTOGRAM));
  if (!jhisto) ErrorExit(ERROR_NO_MEMORY, "JHISTOalloc(%d, %d): allocation failed", nbins_1, nbins_2);

  // histo->bins = (float *)calloc(nbins, sizeof(float)) ;

  jhisto->counts = (float *)calloc(nbins_1 * nbins_2, sizeof(float));
  // fprintf(stderr, "histo->bins and ->counts allocated %d bins\n", nbins);

  jhisto->nbins_1 = nbins_1;
  jhisto->nbins_2 = nbins_2;

  return (jhisto);
}

JOINT_HISTOGRAM *JHISTOrealloc(JOINT_HISTOGRAM *jhisto, int nbins_1, int nbins_2)
{
  if (jhisto == NULL) return (JHISTOalloc(nbins_1, nbins_2));

  if (jhisto->counts) free(jhisto->counts);
  jhisto->counts = (float *)calloc(nbins_1 * nbins_2, sizeof(float));

  jhisto->nbins_1 = nbins_1;
  jhisto->nbins_2 = nbins_2;

  return (jhisto);
}

//////////////////
//////////////////

int JHISTOfindBin(JOINT_HISTOGRAM *jhisto, double val1, double val2)
{
  int bins1 = jhisto->nbins_1;
  int bins2 = jhisto->nbins_2;

  int i = (int)((val1 * bins1) / (jhisto->max - jhisto->min));
  int j = (int)((val2 * bins2) / (jhisto->max - jhisto->min));
  // insure they are in range
  i = max(min(i, bins1 - 1), 0);
  j = max(min(j, bins2 - 1), 0);
  // do index arithmetic
  return (i * bins2 + j);
}

void JHISTOfill(MRI *mri1, MRI *mri2, JOINT_HISTOGRAM *jhisto)
{
  // MRIcheckVolDims(mri1, mri2);

  int width = mri1->width;
  int height = mri1->height;
  int depth = mri1->depth;
  int frame = mri1->nframes;
  int f, x, y, z, index;
  double val1, val2;
  int count = 0;

  //
  float min1, min2, max1, max2;
  MRIvalRange(mri1, &min1, &max1);
  MRIvalRange(mri2, &min2, &max2);
  jhisto->min = flmin(min1, min2);
  jhisto->max = flmax(max1, max2);
  // printf("Min max of the joint histogram: (%f, %f)\n", jhisto->min, jhisto->max);
  //

  for (f = 0; f < frame; f++)
    for (z = 0; z < depth; z++)
      for (y = 0; y < height; y++)
        for (x = 0; x < width; x++) {
          MRIsampleVolumeFrame(mri1, x, y, z, f, &val1);
          MRIsampleVolumeFrame(mri2, x, y, z, f, &val2);
          if (!(val1 == 0 && val2 == 0)) {
            index = JHISTOfindBin(jhisto, val1, val2);
            jhisto->counts[index]++;
            count++;
          }
        }
  jhisto->sample_count = count;
}

double JHISTOgetEntropy(JOINT_HISTOGRAM *jhisto)
{
  double result = 0, p;
  int i, count;
  int total_bucket_count = jhisto->nbins_1 * jhisto->nbins_2;
  int sample_count = jhisto->sample_count;

  for (i = 0; i < total_bucket_count; i++) {
    count = jhisto->counts[i];
    if (count != 0) {
      p = ((double)count) / sample_count;
      result -= p * log(p);
    }
  }
  return result;
}

/*double
MRIcomputeMi(MRI *mri1, MRI *mri2, int bins1, int bins2)
{
  float min1, min2, max1, max2;
  MRIvalRange(mri1, &min1, &max1) ;
  MRIvalRange(mri2, &min2, &max2) ;

  JOINT_HISTOGRAM* jhisto = JHISTOalloc(bins1, bins2);
  JHISTOfill(mri1, mri2, jhisto); // TODO: want to make histo and jhisto be closer?

  //HISTOGRAM* histo1 = HISTOalloc(bins1);
  HISTOGRAM* histo1 = HISTObins(bins1, min1, max1);
  HISTOfill(mri1, histo1); // use HISTOcount

  //HISTOGRAM* histo2 = HISTOalloc(bins2);
  HISTOGRAM* histo2 = HISTObins(bins2, min2, max2);
  HISTOfill(mri2, histo2); // use HISTOcount

  double marginalentropy1 =   HISTOgetEntropy(histo1);
  //printf("Marginal entropy 1 = %f \n", marginalentropy1);
  double marginalentropy2 =   HISTOgetEntropy(histo2);
  //printf("Marginal entropy 2 = %f \n", marginalentropy2);
  double jointentropy     =   JHISTOgetEntropy(jhisto);
  //printf("Joint entropy = %f \n", jointentropy);

  double mi_score = marginalentropy1 + marginalentropy2 - jointentropy;

  return mi_score;
}
*/
