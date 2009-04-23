#include "CostFunctions.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include "RobustGaussian.h"

#ifdef __cplusplus
extern "C" {
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
