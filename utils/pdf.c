#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdf.h"
#include "proto.h"
#include "utils.h"
#include <sys/time.h>


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
  //seed = ((unsigned long) 1000000)*tv.tv_sec + tv.tv_usec;
  seed = (unsigned long)( tv.tv_sec + tv.tv_usec);
  //printf("%ld %ld\n",tv.tv_sec,tv.tv_usec);
  return(seed);
}


/*********************************************************
 * Name:    PDFgaussian(void)
 * Purpose: generates random numbers that obey a gaussian 
 *          distribution with zero mean and std dev of 1:
 *              pdf(x) = e^(x^2/2)/sqrt(2pi)
 ************************************************************/
double PDFgaussian(void)
{
  double v1,v2,r2;

  do
    {
      v1 = 2.0 * drand48() - 1.0;
      v2 = 2.0 * drand48() - 1.0;
      r2 = v1*v1 + v2*v2;
    } 
  while( r2 > 1.0);

  return( v1 * sqrt( -2.0 * log(r2)/r2 ));
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
  for(n=0; n < order; n++)
    v = v + -log(drand48());
  v /= order;
  return(v);
}

/*----------------------------------------------------------------
  PDFsampleCDF() - sample a value from the given CDF. The resulting
  data will be distributed with the PDF that created the CDF. This is
  a fairly inefficient implementation. See also PDFloadCDF().
  -------------------------------------------------------------------*/
double PDFsampleCDF(double *xcdf, double *cdf, int ncdf)
{
  double u, x, d, dmin;
  int n;

  u = drand48();

  // This can be done much more efficiently by searching
  // an ordered table.
  dmin = 10;
  for(n=0; n < ncdf; n++){
    d = fabs(u-cdf[n]);
    if(dmin > d){
      dmin = d;
      x = xcdf[n];
    }
  }
  return(x);
}

/*----------------------------------------------------------------
  PDFloadCDF() - read in a CDF. The file format is that the first
  line has the number of rows. Each row has two columns. The first
  column is the x at which the cdf is sampled, the second column
  is the value of the cdf. See also PDFsampleCDF().
  ----------------------------------------------------------------*/
int PDFloadCDF(char *fname, double **xcdf, double **cdf, int *ncdf)
{
  FILE *fp;
  int n;

  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",fname);
    return(1);
  }

  fscanf(fp,"%d",ncdf);
  //printf("ncdf = %d\n",*ncdf);

  *xcdf = (double *) calloc(*ncdf,sizeof(double));
  *cdf  = (double *) calloc(*ncdf,sizeof(double));

  for(n=0; n < *ncdf; n++)
    fscanf(fp,"%lf %lf",(*xcdf+n),(*cdf+n));

  fclose(fp);

  return(0);
}
