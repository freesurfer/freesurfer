#ifndef __HausdorffMembershipFunction_txx
#define __HausdorffMembershipFunction_txx

#include <iostream>

#include <fstream>
#include "HausdorffMembershipFunction.h"


template < class TVector >
	HausdorffMembershipFunction< TVector >
::HausdorffMembershipFunction():Superclass()
{
}
template < class TVector >
double 
HausdorffMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const
{
//	std::cout << " hausdorff " << std::endl;
	typedef typename MeasurementVectorType::CellType CellType;
	const std::vector<CellType>* labels1 =m1->GetLabels();
	//const std::vector<CellType>* labels2 =m2->GetLabels();
	double max1=0.0, max2=0.0;	
	int numPoints =  labels1->size()-1;
	for(int i=0;i<numPoints;i++)
	{
		double min1=std::numeric_limits<double>::max();
		double min2=std::numeric_limits<double>::max();	
		for(int j=0;j<numPoints;j++)
		{
			double euclid=0; 
			double euclid2=0; 
			for(int k=0;k<3;k++)
			{	
				euclid+=pow((*m1)[i*3+k]-(*m2)[j*3+k],2);
				euclid2+=pow((*m2)[i*3+k]-(*m1)[j*3+k],2);
			}
			min1= std::min( euclid,min1);
			min2= std::min(euclid2, min2);	
		}	
		max1= std::max(max1,min1);
		max2= std::max(max2,min2);
	}

	return 1/(std::max(max1, max2)+1);
}

#endif
