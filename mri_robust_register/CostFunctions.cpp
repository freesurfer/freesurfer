/**
 * @file  CostFunctions.cpp
 * @brief A class that makes available many different cost functions for images
 *   and to combine multiple volumes by mean or median
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.16 $
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
 
#include "CostFunctions.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include "RobustGaussian.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "error.h"
#include "macros.h"
#include "mrimorph.h"
#include "matrix.h"

#ifdef __cplusplus
}
#endif
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

using namespace std;

float CostFunctions::mean(MRI *i)
{
  int count = 0;
  double d= 0.0;
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    d += (*it1);
    count++;
  }
  return (float)(d/count);
}

float CostFunctions::var(MRI *i)
{
  double m = mean(i);
  double d= 0.0;
  double dd;
  int count = 0;
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    dd = (*it1)-m;
    d += dd * dd;
    count++;
  }
  return (float)(d/count);
}

float CostFunctions::median(MRI *i)
{
  int n = i->width * i->height * i->depth;
  float* t = (float *)calloc(n, sizeof(float));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  float qs = RobustGaussian<float>::median(t,n);
  free(t);
  return qs;
}

float CostFunctions::mad(MRI *i,float d)
{
  int n = i->width * i->height * i->depth;
  float* t = (float *)calloc(n, sizeof(float));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  float qs = RobustGaussian<float>::mad(t,n,d);
  free(t);
  return qs;
}


float CostFunctions::leastSquares(MRI * i1, MRI * i2)
{
  assert(i1 != NULL);

  if (i2)
  {
    assert(i1->width  == i2->width);
    assert(i1->height == i2->height);
    assert(i1->depth  == i2->depth);
  }

  double d = 0;
  //cout << "CostFunctions::leastSquares chunk data" <<endl;
  double dd;
  if (i2 == NULL)
  {
    MRIiterator it1(i1);
    for (it1.begin();!it1.isEnd(); it1++)
      d += (*it1) * (*it1);
  }
  else
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    it2.begin();
    //assert(i1->type == i2->type);
    for (it1.begin();!it1.isEnd(); it1++)
    {
      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      dd = (double)(*it1) - (double)(*it2);
      d += dd * dd;
      it2++;
    }

  }
  //cout << " d: " << d << endl;
  return (float)d;
}


double CostFunctions::tukeyBiweight(MRI * i1, MRI * i2,double sat)
{
  assert(i1 != NULL);

  if (i2)
  {
    assert(i1->width  == i2->width);
    assert(i1->height == i2->height);
    assert(i1->depth  == i2->depth);
  }

  int n = i1->width * i1->height * i1->depth;
  float* diff = (float *)calloc(n, sizeof(float));

  int cc = 0;
  if (i2 == NULL)
  {
    MRIiterator it1(i1);
    for (it1.begin();!it1.isEnd(); it1++)
      diff[cc] = (*it1);
    cc++;
  }
  else
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    it2.begin();
    //assert(i1->type == i2->type);
    for (it1.begin();!it1.isEnd(); it1++)
    {
      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
      diff[cc] = (*it1) - (*it2);
      cc++;
      it2++;
    }

  }


  float sigma = RobustGaussian<float>::mad(diff,n);

  double d = 0;
  for (int i=0;i<n;i++)
    d += rhoTukeyBiweight(diff[i]/sigma,sat);

  free(diff);
  //cout << " d: " << d << endl;
  return d;
}


float CostFunctions::normalizedCorrelation(MRI * i1, MRI * i2)
{

  assert(i1->width  == i2->width);
  assert(i1->height == i2->height);
  assert(i1->depth  == i2->depth);

  double d   = 0;
  float d1  = 0;
  float d2  = 0;
  double dd1 = 0;
  double dd2 = 0;

  if (i1->ischunked && i2->ischunked)
  {
    MRIiterator it1(i1);
    MRIiterator it2(i2);
    for (it1.begin();!it1.isEnd(); it1++)
    {
      d1 = (*it1);
      d2 = (*it2);
      d += d1 * d2;
      dd1 += d1 *d1;
      dd2 += d2 *d2;
      it2++;
    }
  }
  else
  {
    for (int z = 0;z<i1->depth;z++)
      for (int y = 0;y<i1->height;y++)
        for (int x = 0;x<i1->width;x++)
        {
          d1 = MRIgetVoxVal(i1,x,y,z,0);
          d2 = MRIgetVoxVal(i2,x,y,z,0);
          d += d1 * d2;
          dd1 += d1 *d1;
          dd2 += d2 *d2;
        }
  }
  return (float)(d/(sqrt(dd1)*sqrt(dd2)));
}

double CostFunctions::moment(MRI *i, int x, int y, int z)
{
  double dd= 0.0;
  for (int d = 0 ; d<i->depth; d++)
    for (int h = 0 ; h<i->height; h++)
      for (int w = 0 ; w<i->width; w++)
      {
        dd += pow(d+1.0,z) * pow(h+1.0,y) * pow(w+1.0,x) * 
              MRIgetVoxVal(i, (int)w,(int)h,(int)d,0);
      }
  return dd;
}

std::vector < double > CostFunctions::centroid(MRI *i)
// M_100/M_000 , M_010/M_000 , M_001 / M_000
{
  //cout << "CostFunctions::centroid" << endl;
  std::vector < double > dd(3,0.0);
  double n = 0;
	double val;
  for (int d = 0 ; d<i->depth; d++)
    for (int h = 0 ; h<i->height; h++)
      for (int w = 0 ; w<i->width; w++)
      {
			  val = MRIgetVoxVal(i,w,h,d,0);
        n += val;
        dd[0] += (w+1) * val;
        dd[1] += (h+1) * val;
        dd[2] += (d+1) * val;				
      }
  dd[0] = dd[0]/n;
  dd[1] = dd[1]/n;
  dd[2] = dd[2]/n;
	if (isnan(dd[0] +dd[1] +dd[2]))
	{
	   cerr << "CostFunctions::centroid is NAN (empty image? n = " << n << " )"<<endl;
		 exit(1);
	}
	
  return dd;
}

vnl_matrix_fixed < double ,3 ,3 > CostFunctions::orientation(MRI *i)
// M_ijk below is the moment(i,j,k)
{
  // compute mean
  std::vector < double > dd(3,0.0);
  double n = 0;
  float val;
  int wp1, hp1, dp1;

//  MATRIX* cov = MatrixAlloc(3,3,MATRIX_REAL);
//  cov = MatrixZero(3,3,cov);
	vnl_matrix_fixed < double, 3,3> cov(0.0);

  for (int d = 0 ; d<i->depth; d++)
  {
    dp1 = d + 1;
    for (int h = 0 ; h<i->height; h++)
    {
      hp1 = h+1;
      for (int w = 0 ; w<i->width; w++)
      {
        wp1 = w+1;
        val = MRIgetVoxVal(i, w,h,d,0);

        n += val; // M_000
        dd[0] += (wp1) * val; // M_100
        dd[1] += (hp1) * val; // M_010
        dd[2] += (dp1) * val; // M_001

        cov[0][0] += wp1*wp1*val; // M_200
        cov[1][1] += hp1*hp1*val; // M_020
        cov[2][2] += dp1*dp1*val; // M_002
        cov[0][1] += wp1*hp1*val; // M_110
        cov[1][2] += hp1*dp1*val; // M_011
        cov[0][2] += wp1*dp1*val; // M_101


      }
    }
  }
  // compute means:
  dd[0] = dd[0]/n;
  dd[1] = dd[1]/n;
  dd[2] = dd[2]/n;

  // compute cov matrix
  cov[0][0] = (cov[0][0]/n)- dd[0]*dd[0]; // M_200/M_000 - mean_x^2
  cov[1][1] = (cov[1][1]/n)- dd[1]*dd[1];
  cov[2][2] = (cov[2][2]/n)- dd[2]*dd[2];
  cov[0][1] = (cov[0][1]/n)- dd[0]*dd[1]; // M_110/M_000 - mean_x * mean_y
  cov[1][2] = (cov[1][2]/n)- dd[1]*dd[2];
  cov[0][2] = (cov[0][2]/n)- dd[0]*dd[2];

  // mirror off diagonal elements: cov is symmetric
  cov[1][0] = cov[0][1];
  cov[2][1] = cov[1][2];
  cov[2][0] = cov[0][2];

  // compute Eigenvectors
  vnl_symmetric_eigensystem < double > SymEig(cov);
  // sort:
  unsigned int smallest = 0;
  unsigned int largest  = 0;
  for (uint i = 1; i<3;i++)
  {
    if (SymEig.D[largest] < SymEig.D[i]) largest = i;
    if (SymEig.D[smallest] > SymEig.D[i]) smallest = i;
  }
  unsigned int largetosmall[3];
  largetosmall[0] = largest;
  largetosmall[1] = 3-smallest-largest;
  largetosmall[2] = smallest;	
  // should be sorted when using this lookup:
  assert(SymEig.D[largetosmall[0]] >= SymEig.D[largetosmall[1]]);
  assert(SymEig.D[largetosmall[1]] >= SymEig.D[largetosmall[2]]);
  vnl_matrix_fixed < double,3 ,3 > evec;
  for (uint i =0;i<3;i++)
    evec.set_column(i,SymEig.V.get_column(largetosmall[i]));
	
//   float eval[3];
//   MATRIX* evec = MatrixEigenSystem(cov, eval ,NULL) ;
//   // should be sorted:
//   assert(eval[0] >= eval[1]);
//   assert(eval[1] >= eval[2]);

  // make det positive:
  //double d = MatrixDeterminant(evec);
  double d = SymEig.determinant();
  if (d<0)
  {
    //cout << "Orientation: neg. determinant ..  fixing" << endl;
    for (int r=0;r<3;r++)
      evec[r][0] = - evec[r][0];
  }


  //cout << " evals: " << eval[0] << " " << eval[1] << " " << eval[2] << endl;
  //cout << " evecs: " << endl;
  // MatrixPrintFmt(stdout,"% 2.8f",evec);

  //MatrixFree(&cov);
  return evec;

}
