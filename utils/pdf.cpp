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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "proto.h"
#include "utils.h"

#include "pdf.h"

double round(double x);

/*----------------------------------------------------------------
  PDFtodSeed() - generate a random seed based on the time-of-day.
  The seed will be different from microsecond-to-microsecond. This
  can be used as input to srand48().
  -------------------------------------------------------------------*/
unsigned long PDFtodSeed(void)
{
  struct timeval tv;
  unsigned long seed;
  gettimeofday(&tv, NULL);
  // seed = ((unsigned long) 1000000)*tv.tv_sec + tv.tv_usec;
  seed = (unsigned long)(tv.tv_sec + tv.tv_usec);
  // printf("%ld %ld\n",tv.tv_sec,tv.tv_usec);
  return (seed);
}

/*********************************************************
 * Name:    PDFgaussian(void)
 * Purpose: generates random numbers that obey a gaussian
 *          distribution with zero mean and std dev of 1:
 *              pdf(x) = e^(x^2/2)/sqrt(2pi)
 ************************************************************/
double PDFgaussian(void)
{
  double v1, v2, r2;

  do {
    v1 = 2.0 * drand48() - 1.0;
    v2 = 2.0 * drand48() - 1.0;
    r2 = v1 * v1 + v2 * v2;
  } while (r2 > 1.0);

  return (v1 * sqrt(-2.0 * log(r2) / r2));
}
/*********************************************************
 * Name:    PDFerlang(order)
 * Purpose: generates random numbers that obey an erlang
 *          distribution with mean=1 and stddev = mean/sqrt(order):
 *   pdf(x) = r*((r*(x-avg+1))^(r-1)) * exp(-r*(x-avg+1)) / (r-1)!
 * when order=1, an exponential distribution is generated.
 ************************************************************/
double PDFerlang(int order)
{
  double v, n;

  v = 0;
  for (n = 0; n < order; n++) v = v + -log(drand48());
  v /= order;
  return (v);
}

/*----------------------------------------------------------------
  PDFsampleCDF() - sample a value from the given CDF. The resulting
  data will be distributed with the PDF that created the CDF.  cdf[n]
  is the probability that the random number will be <= xcdf[n].  See
  also PDFloadCDF().
  -------------------------------------------------------------------*/
double PDFsampleCDF(double *xcdf, double *cdf, int ncdf)
{
  double u;
  int n;

  u = drand48();
  n = PDFsearchOrderedTable(u, cdf, ncdf);
  return (xcdf[n]);

#if 0
  // This is the old brute-force method
  // This can be done much more efficiently by searching
  // an ordered table.
  dmin = fabs(u-cdf[0]);
  x = xcdf[0];
  for (n=1; n < ncdf; n++)
  {
    d = fabs(u-cdf[n]);
    if (dmin > d)
    {
      dmin = d;
      x = xcdf[n];
    }
  }
  return(x);
#endif
}
/*----------------------------------------------------------------
  PDFsearchOrderedTable() - returns the index in y such that y(index)
  is closest to u. Assumes that y is sorted from lowest to highest.
  ----------------------------------------------------------------*/
int PDFsearchOrderedTable(double u, double *y, int ny)
{
  int n1, n2, n3;

  n1 = 0;
  n2 = (int)round(ny / 2);
  n3 = ny - 1;
  while (n1 != n2 && n2 != n3) {
    // printf("n2 = %d, cdf[n2] = %g\n",n2,y[n2]);
    if (y[n2] <= u) {
      n1 = n2;
      n2 = (int)round((n2 + n3) / 2);
    }
    else {
      n3 = n2;
      n2 = (int)round((n1 + n2) / 2);
    }
  }
  // printf("n2 = %d, cdf[n2] = %g\n",n2,y[n2]);
  if (n2 + 1 < ny && fabs(y[n2] - u) > fabs(y[n2 + 1] - u)) n2 = n2 + 1;
  // printf("n2 = %d, cdf[n2] = %g\n",n2,y[n2]);
  return (n2);
}

/*----------------------------------------------------------------
  PDFloadCDF() - read in a CDF. The file format is that each row has
  two columns. The first column is the x at which the cdf is sampled,
  the second column is the value of the cdf. See also PDFsampleCDF().
  ----------------------------------------------------------------*/
int PDFloadCDF(char *fname, double **xcdf, double **cdf, int *ncdf)
{
  FILE *fp;
  int n;
  char tmpstring[1000];

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", fname);
    return (1);
  }

  // Count the number of rows
  *ncdf = 0;
  while (fgets(tmpstring, 1000, fp) != NULL) (*ncdf)++;
  fclose(fp);
  fp = fopen(fname, "r");
  // printf("ncdf = %d\n",*ncdf);

  *xcdf = (double *)calloc(*ncdf, sizeof(double));
  *cdf = (double *)calloc(*ncdf, sizeof(double));

  for (n = 0; n < *ncdf; n++) {
    if (fscanf(fp, "%lf %lf", (*xcdf + n), (*cdf + n)) != 2) {
      printf("ERROR: cannot read parameters\n");
    }
  }
  fclose(fp);

  return (0);
}
