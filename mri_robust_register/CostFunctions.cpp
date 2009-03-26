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
//   assert(i2 != NULL);
	if (!i1->ischunked || (i2!=NULL && ! i2->ischunked))
	{
	   cerr<< "CostFunctions::leastSquares need chunk MRI, set environment: FS_USE_MRI_CHUNK" << endl;
	   exit(1);
	}

	if (i2!=NULL && i1->bytes_total != i2->bytes_total)
	{
	   cerr<< "CostFunctions::leastSquares byte sizes of images differ" << endl;
	   exit(1);
	}


	double d = 0;
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
	   assert(i1->type == i2->type);
	   for (it1.begin();!it1.isEnd(); it1++)
	   {
	      //cout << "it1: " << *it1 << " it2: " << *it2 << endl;
	      dd = (*it1) - (*it2);
	      d += dd * dd;
	      it2++;
	   }
	
	}
	return d;
}

double CostFunctions::normalizedCorrelation(MRI * i1, MRI * i2)
{
	if (!i1->ischunked || ! i2->ischunked)
	{
	   cerr<< "CostFunctions::leastSquares need chunk MRI, set environment: FS_USE_MRI_CHUNK" << endl;
	   exit(1);
	}

	if (i1->bytes_total != i2->bytes_total)
	{
	   cerr<< "CostFunctions::leastSquares byte sizes of images differ" << endl;
	   exit(1);
	}


	double d = 0;
	double d1 = 0;
	double d2 = 0;
 	double dd1 =0;
	double dd2 = 0;

 	MRIiterator it1(i1); 
 	MRIiterator it2(i2); 
 	assert(i1->type == i2->type);
 	for (it1.begin();!it1.isEnd(); it1++)
 	{
	   d1 = (*it1);
	   d2 = (*it2);
 	   d += d1 * d2;
	   dd1 += d1 *d1;
	   dd2 += d2 *d2;
 	   it2++;
 	}
 	
	return d/(sqrt(d1)*sqrt(d2));
}
