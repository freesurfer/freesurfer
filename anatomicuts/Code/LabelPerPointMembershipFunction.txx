#include <iostream>

#include <fstream>
#include "LabelPerPointMembershipFunction.h"
#include <limits>

#include <set>
template < class TVector >
	LabelPerPointMembershipFunction< TVector >
::LabelPerPointMembershipFunction()
{
	this->m_Variance = 0 ;
	this->m_withCosine = false;
	this->m_withEuclid = false;

}
template < class TVector >
double 
LabelPerPointMembershipFunction< TVector >
::EvaluateNO2(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const
{
	std::cout << " hola evaluate " << std::endl;
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
	//return max(dist,dist_inv);
	return dist*(set1.size()+set2.size())/2;

}


template < class TVector >
void  
	LabelPerPointMembershipFunction< TVector >
::AddChild(const MeasurementVectorType* child) 
{
	this->m_Variance += this->Evaluate(child);

	childs.push_back(child);
}
template < class TVector >
void  
	LabelPerPointMembershipFunction< TVector >
::RecalculateCentroid() 
{	
	this->m_Variance =0;
	int numPoints = this->GetCentroid()->GetLabels()->size()-1;
	MeasurementVectorType averageMv(numPoints*3);

	for(int i=0;i<this->childs.size();i++)
	{
		double totalEuclid=0, totalEuclid_inv=0;
		for(int j=0;j<numPoints;j++)
		{
			double euclid=0, euclid_inv=0;	
			for(int k=0;k<3;k++)
			{
				euclid+=pow((*this->GetCentroid())[j*3+k]-(*this->childs[i])[j*3+k],2);
				euclid_inv+=pow((*this->GetCentroid())[j*3+k]-(*this->childs[i])[numPoints*3-j*3+k-3],2);
			}
			totalEuclid+= sqrt(euclid);
			totalEuclid_inv+= sqrt(euclid_inv);
		}

		if(totalEuclid < totalEuclid_inv)
		{
			for(int k=0; k< numPoints*3;k++)
				averageMv[k]=(*this->childs[i])[k]/this->childs.size();
		}
		else
		{
			for(int j=0;j<numPoints;j++)
			{	
				for(int k=0; k<3; k++)
				{	
					averageMv[j*3+k]=(*this->childs[i])[numPoints*3-j*3+k-3]/this->childs.size();
				}
			}
		}
	}
	double minDist = std::numeric_limits<double>::max();
	for(int i=0;i<this->childs.size();i++)
	{
		double totalEuclid=0, totalEuclid_inv=0;
		for(int j=0;j<numPoints;j++)
		{
			double euclid=0, euclid_inv=0;
			for(int k=0;k<3;k++)
			{
				euclid+=pow(averageMv[j*3+k]-(*this->childs[i])[j*3+k],2);
				euclid_inv+=pow((averageMv)[j*3+k]-(*this->childs[i])[numPoints*3-j*3+k-3],2);
			}

			totalEuclid+= sqrt(euclid);
			totalEuclid_inv+= sqrt(euclid_inv);
		}

		if(totalEuclid < minDist || totalEuclid_inv < minDist)
		{
			minDist = std::min(totalEuclid,totalEuclid_inv);
			this->SetCentroid(this->childs[i]);
		}

	}
}

template < class TVector >
void  
LabelPerPointMembershipFunction< TVector >
::PrintSelf(std::ostream& os, itk::Indent indent) const
{

	typedef typename MeasurementVectorType::CellType CellType;
	const std::vector<CellType>* labels1 =this->GetCentroid()->GetLabels();
	for(int i=0;i<labels1->size();i++)
	{
		for(int j=0;j<(*labels1)[i].size();j++)
			std::cout <<  (*labels1)[i][j] << " ";
	}
	std::cout << std::endl;
}
/*template < class TVector >
void  
LabelPerPointMembershipFunction< TVector >
::SetCentroid( const MeasurementVectorType* c)
		{
			this->m_Centroid = c; 
		}
template < class TVector >
		const typedef MeasurementVectorType* 
LabelPerPointMembershipFunction< TVector >::GetCentroid() const {
			return this->m_Centroid; 
		}

template < class TVector >
		std::vector<const MeasurementVectorType*>  
LabelPerPointMembershipFunction< TVector >::GetChilds()
		{
			return this->childs;
		}
template < class TVector >
		double 
LabelPerPointMembershipFunction< TVector >::GetVariance(){ return this->m_Variance/this->childs.size();}

template < class TVector >
		double
LabelPerPointMembershipFunction< TVector >::Evaluate(const MeasurementVectorType *measurement) const{return this->Evaluate(this->GetCentroid(), measurement);}

template < class TVector >
		double 
LabelPerPointMembershipFunction< TVector >::Evaluate(const MeasurementVectorType &measurement) const{ std::cout << "not implemented " << std::endl;return -1;}

template < class TVector >
		void 

LabelPerPointMembershipFunction< TVector >::WithEuclid(bool on)
		{
			this->m_withEuclid = on;
		}
template < class TVector >
		void 
LabelPerPointMembershipFunction< TVector >::WithCosine(bool on)
		{
			this->m_withCosine = on;
		}
		
template < class TVector >
		void 
LabelPerPointMembershipFunction< TVector >::ClearChilds(){ this->childs.clear();}

template < class TVector >
		int 
LabelPerPointMembershipFunction< TVector >::GetNumberOfChilds(){return this->childs.size();}
template < class TVector >
		void 
LabelPerPointMembershipFunction< TVector >::AddDirectionalNeighbors(vnl_matrix<int>* neighbors)
		{
			this->m_directionalNeighbors.push_back(neighbors);
		}
template < class TVector >
		void 
LabelPerPointMembershipFunction< TVector >::ClearDirectionalNeighbors()
		{
			this->m_Variance = 0;
			this->m_directionalNeighbors.clear();
		}
*/
