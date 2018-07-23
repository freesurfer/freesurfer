#ifndef __itkLabelsHistogramMembershipFunction_txx
#define __itkLabelsHistogramMembershipFunction_txx

#include <iostream>

#include <fstream>
#include "itkLabelsHistogramMembershipFunction.h"
#include <limits>

#include <set>
template < class TVector >
double 
LabelsHistogramMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const
{
	double dist=0.0;	
	double dist_inv=0.0;	
	typedef typename MeasurementVectorType::CellType CellType;
	const std::vector<CellType> labels1 =*m1->GetLabels();
	const std::vector<CellType> labels2 =*m2->GetLabels();
	std::set<int> set1, set2;

	int validLabels=0;	
	for(int i=0;i<labels2.size();i++)
	{
		double max_d=0.0;
		for(int j=0;j<labels1.size();j++)
		{	
			//	int j=i;
			double d=0.0;
			for(int k=0;k<labels2[i].size();k++)
			{
				int label2 = labels2[i][k];
				int label1 = labels1[j][k];
				if( label2==label1)
				{
					set1.insert(label1);
					set2.insert(label2);

					dist++;
				}
			}
		}

	}
	return dist*(set1.size()+set2.size())/2;
}



#endif
