#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h> // exit

#include "transformUtils.h"

float* read_transform(const char* fname)
{
  std::ifstream ifs(fname);

  if ( !ifs )
  {
    std::cerr << " Error opening file " << fname
    << " to read transform\n";
    return NULL;
  }

  float* fp = new float[12];
  unsigned int ui=0;
  for ( float *pbuf = fp; ui<12 && ifs; ++ui, ++pbuf)
    ifs >> *pbuf;


  if ( ui<12 )
  {
    delete[] fp;
    return NULL;
  }

  return fp;
}

void
write_transform(float* transform, const char* fname)
{
  std::ofstream ofs(fname);
  if ( !ofs )
  {
    std::cerr << " Error opening file " << fname
    << " for writing transform\n";
    return;
  }

  for (unsigned int ui=0; ui<12; ++ui)
    ofs << *(transform+ui) << " ";
  ofs << std::endl;

  ofs.close();
}

void
inv_transform(float* t,
              float** pinv_t)
{
  // compute the inverse of a 3x3 matrix - use the cofactor method

  float fdet = t[0] * ( t[4]*t[8] - t[5]*t[7] )
               - t[3] * ( t[1]*t[8] - t[2]*t[7] )
               + t[6] * ( t[1]*t[5] - t[2]*t[4] );

  if ( fabs(fdet) < 1.0e-5 )
  {
    std::cerr << " inv_transform -> null det \n";
    exit(1);
  }

  float *inv_t = new float[12];

  inv_t[0] = ( t[4]*t[8] - t[5]*t[7] ) / fdet;
  inv_t[1] = ( t[2]*t[7] - t[1]*t[8] ) / fdet;
  inv_t[2] = ( t[1]*t[5] - t[2]*t[4] ) / fdet;

  inv_t[3] = ( t[5]*t[6] - t[3]*t[8] ) / fdet;
  inv_t[4] = ( t[0]*t[8] - t[2]*t[6] ) / fdet;
  inv_t[5] = ( t[2]*t[3] - t[0]*t[5] ) / fdet;

  inv_t[6] = ( t[3]*t[7] - t[4]*t[6] ) / fdet;
  inv_t[7] = ( t[1]*t[6] - t[0]*t[7] ) / fdet;
  inv_t[8] = ( t[0]*t[4] - t[1]*t[3] ) / fdet;

  inv_t[9] = - t[9]* inv_t[0] - t[10]* inv_t[3] - t[11]* inv_t[6];
  inv_t[10]= - t[9]* inv_t[1] - t[10]* inv_t[4] - t[11]* inv_t[7];
  inv_t[11]= - t[9]* inv_t[2] - t[10]* inv_t[5] - t[11]* inv_t[8];

  *pinv_t = inv_t;
}
