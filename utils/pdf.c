#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pdf.h"

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

