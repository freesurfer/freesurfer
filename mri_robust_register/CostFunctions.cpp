#include "CostFunctions.h"

#include <cassert>
#include <fstream>
#include <sstream>

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
	if (!i->ischunked)
	{
	   cerr<< "CostFunctions::mean need chunk MRI, set environment: FS_USE_MRI_CHUNK" << endl;
	   exit(1);
	}

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
	if (!i->ischunked)
	{
	   cerr<< "CostFunctions::mean need chunk MRI, set environment: FS_USE_MRI_CHUNK" << endl;
	   exit(1);
	}

	double m = mean(i);
	double d= 0.0;
	double dd;
	int count = 0;
	MRIiterator it1(i); 
	for (it1.begin();!it1.isEnd(); it1++)
	{
	   dd = (*it1)-m;
	   d += dd;
	   count++;
	}
       return d/count;
}


double CostFunctions::leastSquares(MRI * i1, MRI * i2)
{
   assert(i1 != NULL);
 
   assert(i1->width  == i2->width);
   assert(i1->height == i2->height);
   assert(i1->depth  == i2->depth);

   double d = 0;
   if(i1->ischunked && (i2==NULL || i2->ischunked))
   {
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
	   //assert(i1->type == i2->type);
	   for (it1.begin();!it1.isEnd(); it1++)
	   {
	      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
	      dd = (double)(*it1) - (double)(*it2);
	      d += dd * dd;
	      it2++;
	   }
	
	}
    }
    else
    {
       //cout << "CostFunctions::LeastSquares standart data" << endl;
	double dd;
	if (i2 == NULL)
	{
	   for (int z = 0;z<i1->depth;z++)
	   for (int y = 0;y<i1->height;y++)
	   for (int x = 0;x<i1->width;x++)
	   {
	      dd = (double)MRIgetVoxVal(i1,x,y,z,0) ;
	      d += dd*dd ;
	   }
	   
	}
	else
	{
	   for (int z = 0;z<i1->depth;z++)
	   for (int y = 0;y<i1->height;y++)
	   for (int x = 0;x<i1->width;x++)
	   {
	      dd = (double)MRIgetVoxVal(i1,x,y,z,0) - (double)MRIgetVoxVal(i2,x,y,z,0);
	      d += dd*dd ;
	   }
	   
	
	}
       
    
    }
    cout << " d: " << d << endl;
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
