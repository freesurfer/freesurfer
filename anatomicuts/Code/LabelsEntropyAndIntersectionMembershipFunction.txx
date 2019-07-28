#include <iostream>

#include <fstream>
#include "LabelsEntropyAndIntersectionMembershipFunction.h"
#include <limits>

#include <set>
template < class TVector >
double 
LabelsEntropyAndIntersectionMembershipFunction< TVector >
::Evaluate(const MeasurementVectorType *m1, const MeasurementVectorType *m2 ) const
{
	typedef typename MeasurementVectorType::LabelsDirectionType LabelsDirectionsType;
	typedef typename MeasurementVectorType::LabelsMapType LabelsMapType;
	const  LabelsDirectionsType labels1 =m1->GetLabelsPerDirection();
	const  LabelsDirectionsType labels2 =m2->GetLabelsPerDirection();
	int numberOfPoints =  m1->GetNumberOfPoints();
	
//	std::cout << "number of points "<< numberOfPoints << std::endl;
	double intersection=0.0;	
	//double entropy=0.0;
	double normalizing=numberOfPoints ; //pow( pow(numberOfPoints,2)*labels1.size(),2);	
	//double dice=0.0;
	std::set<int> adhocEntropy;
	double affinity=1.0;
	if(this->m_labelsAndEuclid)
	{
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if(labels2[i].count( it->first )!= 0)
				{
					intersection+=labels2[i].find(it->first)->second *  it->second ; 
					adhocEntropy.insert(it->first);
				}
			}

		}
		affinity = intersection*adhocEntropy.size()/(pow(normalizing,2)*labels1.size()*5);

		typedef typename MeasurementVectorType::CellType CellType;
		const std::vector<CellType>* labels1_ =m1->GetLabels();
		const std::vector<CellType>* labels2_ =m2->GetLabels();
		double dist=0.0, dist_inv=0.0;	

		for(int i=0;i<labels1_->size();i++)
		{
			double euclid=0, euclid_inv=0;
			//double cos=1, cos_inv=1;
			for(int k=0;k<3;k++)
			{
				euclid+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
				euclid_inv+=pow((*m1)[i*3+k]-(*m2)[labels2_->size()*3-i*3+k-3],2);
			}
			dist+=1.0/sqrt(euclid+1);
			dist_inv+=1.0/sqrt(euclid_inv+1);
		}

		dist= std::max(dist, dist_inv)/labels1_->size();
		//std::cout << affinity << " " << dist << std::endl;
		affinity += dist;
		return affinity;
	
	}
	else if(this->m_labels)
	{

//		std::cout << "labels"<< this->m_labels <<  std::endl;	
//		std::cout << " labels1,siuze" << labels1[0].size()<<  "      " << labels2[0].size()<<std::endl;	
		if(labels1.size()==0 || labels2.size()==0)
		{
			return 0;
		}
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if(labels2[i].count( it->first )!= 0)
				{
					intersection+=labels2[i].find(it->first)->second *  it->second ; 
					adhocEntropy.insert(it->first);
				}
			}

		}
		return intersection*adhocEntropy.size()/normalizing;


	}
	else if(this->m_euclidean)
	{
		double euclidean=0.0;
		
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first) >0)
					euclidean += pow(it->second - labels2[i].find(it->first)->second,2);
				else
					euclidean += pow(it->second,2);
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				if( labels1[i].count(it->first)==0)
					euclidean += pow(it->second,2);
			}
		}
		affinity = 1/sqrt(euclidean+1);
		return affinity;		
	}
	else if(this->m_meanEuclidean)
	{
		//std::cout << "mean euc" <<std::endl;
		typedef typename MeasurementVectorType::CellType CellType;
		const std::vector<CellType>* labels1 =m1->GetLabels();
		const std::vector<CellType>* labels2 =m2->GetLabels();
		double dist=0.0, dist_inv=0.0;	
		//std::cout << " num points " << labels1->size() << std::endl;
		for(int i=0;i<labels1->size();i++)
		{
			double euclid=0, euclid_inv=0;
			//double cos=1, cos_inv=1;
			for(int k=0;k<3;k++)
			{
				euclid+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
				euclid_inv+=pow((*m1)[i*3+k]-(*m2)[labels2->size()*3-i*3+k-3],2);
			}
			dist+=sqrt(euclid);
			dist_inv+=sqrt(euclid_inv);
		}

		dist=std::min(dist/labels1->size(), dist_inv/labels1->size());
		dist =  1.0/(dist+1);
		//std::cout << dist << std::endl;
		return dist;
	}
	else if(this->m_gaussian)
	{
		typedef typename MeasurementVectorType::CellType CellType;
		const std::vector<CellType>* labels1 =m1->GetLabels();
		const std::vector<CellType>* labels2 =m2->GetLabels();
		double dist=0.0, dist_inv=0.0;	

		for(int i=0;i<labels1->size();i++)
		{
			double euclid=0, euclid_inv=0;
			//double cos=1, cos_inv=1;
			for(int k=0;k<3;k++)
			{
				euclid+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
				euclid_inv+=pow((*m1)[i*3+k]-(*m2)[labels2->size()*3-i*3+k-3],2);
			}
			dist+=sqrt(euclid);
			dist_inv+=sqrt(euclid_inv);
		}

		dist= std::min(dist/labels1->size(), dist_inv/labels1->size());
		return exp(-dist/800.0);
	}else if(this->m_meanClosestPointInvert)
	{
//		std::cout << "meann clos" << std::endl;
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
			max1+=sqrt(min1);
			max2+=sqrt(min2);
		}

		double dist = (max1+ max2)/(2*numPoints);
		return 1.0/(dist+1.0);
	}else if(this->m_meanClosestPointGaussian)
	{
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
			max1+=sqrt(min1);
			max2+=sqrt(min2);
		}

		double aff = (max1+ max2)/(2*numPoints);
		return exp(-aff/800.0);
	}
	else if(this->m_meanAndCovGaussian)
	{
		int numPoints=3;
		double dist=0.0;	

		for(int i=0;i<numPoints;i++)
		{
			for(int k=0;k<3;k++)
			{
				dist+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
			}
		}
		dist=sqrt(dist);
//		std::cout << " dist " <<dist/800.0 <<" " <<   exp(-(dist/800.0)) << std::endl;
		return exp(-dist/800.0);

	}
	else if(this->m_meanAndCovInvert)
	{
		int numPoints=3;
		double dist=0.0;	

		for(int i=0;i<numPoints;i++)
		{
			for(int k=0;k<3;k++)
			{
				dist+=pow((*m1)[i*3+k]-(*m2)[i*3+k],2);
			}
		}
		dist=sqrt(dist);
	//	std::cout << " dist " <<dist/800.0 <<" " <<   exp(-(dist/800.0)) << std::endl;
		return 1.0/(dist+1.0);

	}
	std::cout << " ups "<< std::endl;
	return -1;
/*	else if(this->m_entropy )
	{
//		std::cout << "labels"<< this->m_labels <<  std::endl;	
//		std::cout << " labels1,siuze" << labels1[0].size()<<  "      " << labels2[0].size()<<std::endl;	
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if(labels2[i].count( it->first )!= 0)
				{
					double p = labels2[i].find(it->first)->second *  it->second /normalizing;
					entropy+= (-1)*p*log(p);
					intersection+=p; 
					adhocEntropy.insert(it->first);
				}
			}

		}
	}
	else if(this->m_intersection)
	{
		double min=0;
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first)!=0)
					min+=  std::min(it->second, labels2[i].find(it->first)->second);
			}
		}
		affinity=min;;
	
	}
	else if(this->m_euclidean)
	{
		double euclidean=0.0;
		
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first) >0)
					euclidean += pow(it->second - labels2[i].find(it->first)->second,2);
				else
					euclidean += pow(it->second,2);
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				if( labels1[i].count(it->first)==0)
					euclidean += pow(it->second,2);
			}
		}
		affinity = 1/sqrt(euclidean+1);
		
	}
	else if(this->m_dice)
	{
		double dice=0;
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				dice+= pow(it->second,2);
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				dice+= pow(it->second,2);
				if( labels1[i].count(it->first)!=0)
					intersection+= labels1[i].find(it->first)->second * it->second;
			}
		}
		affinity = 2*intersection/dice;

	}
	else if(this->m_ruzicka)
	{
		double min=0.0, max=0.0;
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first)!=0)
				{
					min+=std::min( labels2[i].find(it->first)->second , it->second);
					max+= std::max(labels2[i].find(it->first)->second ,  it->second);
				}
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				if( labels1[i].count(it->first)==0)
				{
					max+=  it->second;
				}
			}
	}
		affinity+=min/(max);
//		std::cout << "affinity " << affinity <<  " min " << min << " max " << max << std::endl;
	}	else if(this->m_kulczynskis)
	{
		double min=0.0, diff=0.0;
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first)!=0)
				{
					min+=std::min( labels2[i].find(it->first)->second , it->second);
					diff+= std::abs(labels2[i].find(it->first)->second -  it->second);
				}
				else
				{
					diff+= it->second;
				}
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				if( labels1[i].count(it->first)==0)
				{
					diff+= it->second;
				}
			}}
		affinity*=min/(diff+1);
//		std::cout << "affinity " << affinity <<  " min " << min << " diff " << diff << std::endl;
	}
	else if(this->m_jensenShannon)
	{
		double dist=0;
		for( int i=0;i<labels1.size(); i++)
		{
			for(typename LabelsMapType::const_iterator it = labels1[i].cbegin();it != labels1[i].cend(); it++)
			{
				if( labels2[i].count(it->first)!=0)
					dist+= it->second * std::log(2*it->second/ (it->second + labels2[i].find(it->first)->second));
			        else
					dist+= it->second * std::log(2);
			}
			for(typename LabelsMapType::const_iterator it = labels2[i].cbegin();it != labels2[i].cend(); it++)
			{
				if( labels1[i].count(it->first)!=0)
					dist+= it->second * std::log(2*it->second/ (it->second + labels1[i].find(it->first)->second));
			        else
					dist+= it->second * std::log(2);
			}
		}
		affinity*=dist*0.5;
	}


	if (this->m_entropy)
		return entropy;/// labels1.size();
	else if(this->m_labels)
	else
		return affinity;
*/
}



