#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdf.h"
#include "proto.h"
#include "utils.h"

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
 *   pdf(x) = r*((r*(x-avg))^(r-1)) * exp(-r*(x-avg)) / (r-1)!
 * when order=1, this generates an exponential distribution.
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

