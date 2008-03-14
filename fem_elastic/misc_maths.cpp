#include <iostream>
#include <map>

#include "misc_maths.h"

using namespace std;

vector<double> gauss_kernel(double dSigma, int max_len)
{
  int len = (int)floor( 8.0 * dSigma + 0.5 ) + 1;

  if ( max_len>0 && max_len<len )
    len = max_len;

  if ( (len>>1)<<1 == len )
    len++;

  int half_len = len>>1;

  double norm = 0.0;
  vector<double> kernel;
  double k;

  for ( int x=0; x<len; x++)
  {
    double dx = (x-half_len);
    if ( abs(dx) <= 2.0 * dSigma )
      k = exp( -dx*dx / (2.0*dSigma*dSigma) );
    else if ( 2.0*dSigma<abs(dx) && abs(dx)<=4.0*dSigma )
      k= 1.0 / (16.0 * M_E * M_E ) * pow(4.0 - abs(dx)/dSigma, 4.0);
    else
      k = 0.0;

    kernel.push_back(k);
    norm += k;
  }

  for (vector<double>::iterator it=kernel.begin();
       it != kernel.end(); it++)
    *it /= norm;

  return kernel;
}


std::vector<double> smooth_kernel(std::vector<double> vin,
                                  std::vector<double> vker)
{

  map<int, int> lut; // to treat borders

  int alpha = (vker.size()-1) >> 1;

  // warn user if kernel is even sized
  if ( vker.size()%2 == 0)
    cout << " smooth_kernel -> even kernel size -> discarding last value\n";

  // build LUT to handle borders
  for ( int i=-alpha; i<0; i++)
    lut[i]=0;

  for ( int i=0; i< int(vin.size()); i++)
    lut[i]=i;

  for ( int i=vin.size(); i<int(vin.size())+alpha; i++)
    lut[i]=vin.size()-1;

  // allocate ret vector
  vector<double> vres(vin.size());
  fill(vres.begin(), vres.end(), 0.0);

  for (int i=0; i< int(vin.size()); i++)
  {
    for ( int j=-alpha; j<=alpha; j++)
      vres[i] += vin[lut[i+j]] * vker[j+alpha];
  }

  return vres;
}



