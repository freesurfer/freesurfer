#ifndef __EuclideanMembershipFunction_txx
#define __EuclideanMembershipFunction_txx

#include <iostream>

#include <fstream>
#include "EuclideanMembershipFunction.h"


template < class TVector >
	EuclideanMembershipFunction< TVector >
::EuclideanMembershipFunction():Superclass()
{
	this->m_Variance = 0 ;
	this->m_withCosine =false;
}
template < class TVector >
double 
EuclideanMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *measurement) const
{
	return this->Evaluate(this->GetCentroid(), measurement);
}
template < class TVector >
double 
EuclideanMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const
{
	typedef typename MeasurementVectorType::CellType CellType;
	const std::vector<CellType>* labels1 =m1->GetLabels();
	const std::vector<CellType>* labels2 =m2->GetLabels();
	double dist=0.0, dist_inv=0.0;	
	//double cos=1, cos_inv=1;
	//std::cout << labels1->size() << std::endl;	

	for(int i=0;i<labels1->size()-1;i++)
	{
		double euclid=0, euclid_inv=0;
		double cos=1, cos_inv=1;
		for(int k=0;k<3;k++)
		{
			euclid+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
			euclid_inv+=pow((*m1)[i*3+k]-(*m2)[(labels2->size()-1)*3-i*3+k-3],2);
		}
		if(this->m_withCosine )
		{
			double norm=0, norm2=0, norm2_inv=0;
			double vec[3],  vec2[3],vec2_inv[3];
			cos=0,cos_inv=0;
			int ii = (i==labels1->size()-1)?i-1:i;
			for(int k=0;k<3;k++)
			{
				vec[k]= (*m1)[ii*3+k] - (*m1)[(ii+1)*3+k];
				norm+= pow(vec[k],2);

				vec2[k]= (*m2)[ii*3+k] - (*m2)[(ii+1)*3+k];
				norm2+= pow(vec2[k],2);

				vec2_inv[k]= (*m2)[labels2->size()*3 - ii*3+k-3] - (*m2)[labels2->size()*3- (ii+1)*3+k-3];
				norm2_inv+= pow(vec2_inv[k],2);

				cos += vec[k]*vec2[k];
				cos_inv += vec[k]*vec2_inv[k];

			}
			norm2=  (sqrt(norm)*sqrt(norm2));
			norm2_inv =(sqrt(norm)*sqrt(norm2_inv));
			cos= abs(cos/norm2);
			cos_inv=abs(cos_inv/norm2_inv);	
		}
		//dist+= (1/sqrt(euclid+1.0))*cos;
		//dist_inv+= (1/sqrt(euclid_inv+1.0))*cos_inv;
		dist+= sqrt(euclid);
		dist_inv+= sqrt(euclid_inv);
		
	}

	//dist =	   max( dist, dist_inv); ///(labels1->size()*7); 
	dist = std::min(dist, dist_inv)/labels1->size();
	dist = 100000.0/(dist+1.0);
//	dist = exp(-dist/25);
//	std::cout << dist << std::endl ;
	return dist;
}


template < class TVector >
void  
	EuclideanMembershipFunction< TVector >
::AddChild(const MeasurementVectorType* child) 
{
	this->m_Variance += this->Evaluate(child);

	childs.push_back(child);
}
template < class TVector >
void  
	EuclideanMembershipFunction< TVector >
::RecalculateCentroid() 
{
	this->m_Variance =0;
	int numPoints = this->GetCentroid()->GetLabels()->size();
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
	/*double var = this->GetVariance();
	std::cout << "childs size" << childs.size() << std::endl;
	if(this->childs.size() > 0)
	{
		double minDist =  numeric_limits<double>::min();
		int offset = this->childs.size()/min(50,(int)this->childs.size());
		for(unsigned int i=0; i< this->childs.size(); i++)
		{
			double min_i = 0;
			//			int offset2= this->childs.size();
			for(unsigned int j=0; j< this->childs.size() ; j++)
			{
				double D = this->Evaluate(childs[i], childs[j]);
				min_i += D*D;

			}

			min_i /= this->childs.size();

			//	min_i =pow(var - min_i,2);		
			if(  min_i > minDist )
			{
				minDist = min_i;
				this->SetCentroid(this->childs[i]);
			}
			i+=offset;
		}

		this->m_Variance =  minDist/this->childs.size();
	
	}
	*/
}

template < class TVector >
void  
EuclideanMembershipFunction< TVector >
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
	Superclass::PrintSelf(os,indent);
}


#endif
