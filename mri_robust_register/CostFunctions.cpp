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

using namespace std;

double CostFunctions::mean(MRI *i)
{
  int count = 0;
  double d= 0.0;
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    d += (*it1);
    count++;
  }
  return d/count;
}

double CostFunctions::var(MRI *i)
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
  return d/count;
}

double CostFunctions::median(MRI *i)
{
  int n = i->width * i->height * i->depth;
  double* t = (double *)calloc(n, sizeof(double));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  double qs = RobustGaussian::median(t,n);
  free(t);
  return qs;
}

double CostFunctions::mad(MRI *i,double d)
{
  int n = i->width * i->height * i->depth;
  double* t = (double *)calloc(n, sizeof(double));

  int cc = 0;
  // copy array to t
  MRIiterator it1(i);
  for (it1.begin();!it1.isEnd(); it1++)
  {
    t[cc] = (*it1);
    cc++;
  }

  double qs = RobustGaussian::mad(t,n,d);
  free(t);
  return qs;
}


double CostFunctions::leastSquares(MRI * i1, MRI * i2)
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
  return d;
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
  double* diff = (double *)calloc(n, sizeof(double));

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
      diff[cc] = (double)(*it1) - (double)(*it2);
      cc++;
      it2++;
    }

  }


  double sigma = RobustGaussian::mad(diff,n);

  double d = 0;
  for (int i=0;i<n;i++)
    d += rhoTukeyBiweight(diff[i]/sigma,sat);

  free(diff);
  //cout << " d: " << d << endl;
  return d;
}


double CostFunctions::normalizedCorrelation(MRI * i1, MRI * i2)
{

  assert(i1->width  == i2->width);
  assert(i1->height == i2->height);
  assert(i1->depth  == i2->depth);

  double d   = 0;
  double d1  = 0;
  double d2  = 0;
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
  return d/(sqrt(dd1)*sqrt(dd2));
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

std::vector <double> CostFunctions::centroid(MRI *i)
// M_100/M_000 , M_010/M_000 , M_001 / M_000
{
  std::vector < double > dd(3,0.0);
  double n = 0;
  for (int d = 0 ; d<i->depth; d++)
    for (int h = 0 ; h<i->height; h++)
      for (int w = 0 ; w<i->width; w++)
      {
        n += MRIgetVoxVal(i, w,h,d,0);
        dd[0] += (w+1) * MRIgetVoxVal(i, w,h,d,0);
        dd[1] += (h+1) * MRIgetVoxVal(i, w,h,d,0);
        dd[2] += (d+1) * MRIgetVoxVal(i, w,h,d,0);
      }
  dd[0] = dd[0]/n;
  dd[1] = dd[1]/n;
  dd[2] = dd[2]/n;
  return dd;
}

MATRIX* CostFunctions::orientation(MRI *i)
// M_ijk is the moment(i,j,k)
{
  // compute mean
  std::vector < double > dd(3,0.0);
  double n = 0;
  double val;
  int wp1, hp1, dp1;

  MATRIX* cov = MatrixAlloc(3,3,MATRIX_REAL);
  cov = MatrixZero(3,3,cov);

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

        cov->rptr[1][1] += wp1*wp1*val; // M_200
        cov->rptr[2][2] += hp1*hp1*val; // M_020
        cov->rptr[3][3] += dp1*dp1*val; // M_002
        cov->rptr[1][2] += wp1*hp1*val; // M_110
        cov->rptr[2][3] += hp1*dp1*val; // M_011
        cov->rptr[1][3] += wp1*dp1*val; // M_101


      }
    }
  }
  // compute means:
  dd[0] = dd[0]/n;
  dd[1] = dd[1]/n;
  dd[2] = dd[2]/n;

  // compute cov matrix
  cov->rptr[1][1] = (cov->rptr[1][1]/n)- dd[0]*dd[0]; // M_200/M_000 - mean_x^2
  cov->rptr[2][2] = (cov->rptr[2][2]/n)- dd[1]*dd[1];
  cov->rptr[3][3] = (cov->rptr[3][3]/n)- dd[2]*dd[2];
  cov->rptr[1][2] = (cov->rptr[1][2]/n)- dd[0]*dd[1]; // M_110/M_000 - mean_x * mean_y
  cov->rptr[2][3] = (cov->rptr[2][3]/n)- dd[1]*dd[2];
  cov->rptr[1][3] = (cov->rptr[1][3]/n)- dd[0]*dd[2];

  // mirror off diagonal elements
  cov->rptr[2][1] = cov->rptr[1][2];
  cov->rptr[3][2] = cov->rptr[2][3];
  cov->rptr[3][1] = cov->rptr[1][3];

  // compute Eigenvectors
  float eval[3];
  MATRIX* evec = MatrixEigenSystem(cov, eval ,NULL) ;
  // should be sorted:
  assert(eval[0] >= eval[1]);
  assert(eval[1] >= eval[2]);

  // make det positive:
  double d = MatrixDeterminant(evec);
  if (d<0)
  {
    //cout << "Orientation: neg. determinant ..  fixing" << endl;
    for (int r=1;r<4;r++)
      evec->rptr[r][1] = - evec->rptr[r][1];
  }


  //cout << " evals: " << eval[0] << " " << eval[1] << " " << eval[2] << endl;
  //cout << " evecs: " << endl;
  // MatrixPrintFmt(stdout,"% 2.8f",evec);

  MatrixFree(&cov);
  return evec;

}
