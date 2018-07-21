#ifndef __LabelsPointToPointMembershipFunction_txx
#define __LabelsPointToPointMembershipFunction_txx

#include <iostream>

#include <fstream>
#include "LabelsPointToPointMembershipFunction.h"
#include <limits>

#include <set>
template < class TVector >
double 
LabelsPointToPointMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const 
{
	double dist=0.0;	
	//double dist_inv=0.0;	
	typedef typename MeasurementVectorType::CellType CellType;
	const std::vector<CellType> labels1 =*m1->GetLabels();
	const std::vector<CellType> labels2 =*m2->GetLabels();
	//std::cout << " labels coutn " << this->m_labelsCount << std::endl;
	//int validLabels=0;	
	for(int i=0;i<labels2.size();i++)
	{
		for(int k=0;k<labels2[i].size() && k < this->m_labelsCount;k++)
		{
			int label2 = labels2[i][k];
			int label1 = labels1[i][k];
			if( label2==label1)
			{
				dist++;
			}
		}

	}
	return dist;
}


#endif
