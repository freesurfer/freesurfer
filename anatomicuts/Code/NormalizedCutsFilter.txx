#ifndef _NormalizedCutsFilter_txx
#define _NormalizedCutsFilter_txx


#include "NormalizedCutsFilter.h"
#include "itkDecisionRule.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkDistanceToCentroidMembershipFunction.h"
//#include "itkCurrentsToCentroidMembershipFunction.h"
//#include "itkEuclideanToCentroidMembershipFunction.h"
//#include "itkKMeansClassifierFilter.h"
#include "itkNormalVariateGenerator.h"
#include <iostream>
#include <limits>
//#include <utility>
#include <algorithm>
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_math.h"
#include "math.h"
#include <stdlib.h> 
#include <vnl/vnl_sparse_matrix.h>
#include <vnl/algo/vnl_sparse_symmetric_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <set>
#include "ThreadedMembershipFunction.h"

template< class T>
class PriorityNode {
	typedef typename T::Pointer TPointer;
	public:
	TPointer _thing;
	int _priority ;
	std::string _id;
	PriorityNode(int priority, std::string id, TPointer thing) { _priority=priority; _id=id; _thing=thing; }
	bool operator<(const PriorityNode n) const { return n._priority > this->_priority; }
};

template< class TMesh,class  TMembershipFunctionType>
	void
NormalizedCutsFilter < TMesh,TMembershipFunctionType >::Update()
{
	if (!this->GetInput())
		itkExceptionMacro (<< "Mesh not set");

	int numberOfPoints  =  this->GetInput()->GetCells()->Begin().Value()->GetNumberOfPoints();
	std::cout << " number of points " << numberOfPoints << std::endl;

	std::vector<std::string> labels;
	labels.resize(this->GetInput()->GetNumberOfCells());

	typename SampleType::Pointer sample = SampleType::New();
	//itk 3
	for(unsigned int i =0; i<this->GetInput()->GetNumberOfCells();i++)
	{
		MeasurementVectorType mv(numberOfPoints*3);	      
		mv.SetCell(this->GetInput(), i) ; //inputCellIt.Value());
		sample->PushBack(mv);
	}
	std::string lastLabel="1";
	std::priority_queue<PriorityNode<SampleType>> queue;	
	queue.push(PriorityNode<SampleType>(sample->Size(),lastLabel,sample));
	while(queue.size()<this->GetNumberOfClusters())
	{
		PriorityNode<SampleType> node = queue.top();
		sample = node._thing;
		lastLabel=node._id;
		queue.pop();
		std::vector<std::pair<int,int>> centroidIndeces =  this->SelectCentroidsParallel( sample,(*this->GetMembershipFunctionVector())[0]);	

		typename SampleType::Pointer samplePositives = SampleType::New();
		typename SampleType::Pointer sampleNegatives = SampleType::New();
		
		if(sample->Size() > this->GetNumberOfFibersForEigenDecomposition())
		{
			//Multi-thread
			std::vector<std::pair<int, int>> inIndeces;
			std::vector<std::pair<int, int>> outIndeces;

			for( int i=0; i< centroidIndeces.size();i++)
			{	
				for(int j=0; j< sample->Size();j++)
				{
					inIndeces.push_back(std::pair<int,int>(j,centroidIndeces[i].second));
					outIndeces.push_back(std::pair<int,int>(j,i));

				}
			}
			typename ThreadedMembershipFunctionType::Pointer threadedMembershipFunction = ThreadedMembershipFunctionType::New();
			typename ThreadedMembershipFunctionType::DomainType domain;
			domain[0]=0;
			domain[1]= inIndeces.size()-1;
			typename MembershipFunctionType::Pointer hola = (*this->GetMembershipFunctionVector())[0];
			threadedMembershipFunction->SetStuff(sample,inIndeces, outIndeces,hola,sample->Size());
			threadedMembershipFunction->Execute(hola ,domain);
			//vnl_sparse_matrix<double>* ms= threadedMembershipFunction->GetResults();
			//std::cout << " finding maximum start " << std::endl;	
			std::vector<int> maxvals = threadedMembershipFunction->GetMaxIndeces();
			for(int j=0; j< sample->Size();j++)
			{
				int argmax =0; //, maxVal=0;
				argmax=maxvals[j];
				labels[sample->GetMeasurementVector(j).GetCellId()]=lastLabel +std::to_string(centroidIndeces[argmax].first) ;

				if(centroidIndeces[argmax].first==0)
				{
					samplePositives->PushBack(sample->GetMeasurementVector(j));
				}
				else
				{
					sampleNegatives->PushBack(sample->GetMeasurementVector(j));
				}
			}
			//std::cout << " finding maximum end " << std::endl;	
		}
		else
		{
			typename SampleType::ConstIterator iter = sample->Begin();
			int i=0;
			while ( iter != sample->End() )
			{
				labels[iter.GetMeasurementVector().GetCellId()]=lastLabel + std::to_string(centroidIndeces[i].first);

				if(centroidIndeces[i].first==0)
				{
					samplePositives->PushBack(iter.GetMeasurementVector());
				}
				else
				{
					sampleNegatives->PushBack(iter.GetMeasurementVector());
				}
				++iter;
				i++;
			}
		}
		
		if(samplePositives->Size() > 0 && sampleNegatives->Size()>0)
		{
			this->m_clusterIdHierarchy.push_back(std::make_pair(node._id,lastLabel+"0"));
			this->m_clusterIdHierarchy.push_back(std::make_pair(node._id,lastLabel+"1"));

			queue.push(PriorityNode<SampleType>(samplePositives->Size(),lastLabel+"0",samplePositives));
			queue.push(PriorityNode<SampleType>(sampleNegatives->Size(),lastLabel+"1",sampleNegatives));
			std::cout << "queue size " << queue.size() <<  " positives " << samplePositives->Size() << " negatives " << sampleNegatives->Size() << std::endl;
		}
		else
		{
			std::cout <<  " queue emplace " << samplePositives->Size() << std::endl;
			std::cout <<  " queue emplace " << sampleNegatives->Size() << std::endl;
			std::cin.get();
			throw std::runtime_error("Gak!");
		}
		//lastLabel+=2;

	}
	this->SetLabels(labels);
}

template< class TMesh,class  TMembershipFunctionType>
	std::vector<std::pair<int,int>>	
NormalizedCutsFilter < TMesh ,TMembershipFunctionType>::SelectCentroidsParallel(typename SampleType::Pointer samples, const typename MembershipFunctionType::Pointer membershipFunction )
{

	std::vector<std::pair<int,int>> indices;
	std::vector<int> selected;

	const unsigned int n =std::min(this->GetNumberOfFibersForEigenDecomposition(), (int)samples->Size());
	// =new  vnl_sparse_matrix<double>(n,n);
	vnl_sparse_matrix<double> identity(n,n);
	vnl_sparse_matrix<double> diagonal(n,n);
	
	int offset =(samples->Size()>n)? samples->Size()/n:1;
	for (unsigned i=0; i<n; i++) 
	{
		selected.push_back(i*offset);
	}
	//int zeros =0;
	std::vector<std::pair<int, int>> inIndeces;
	std::vector<std::pair<int, int>> outIndeces;
	for (unsigned i=0; i<n; i++) 
	{
		for (unsigned j=i; j<n; j++) 
		{
			inIndeces.push_back(std::pair<int,int>(selected[i],selected[j]));
			outIndeces.push_back(std::pair<int,int>(i,j));
		}
	}

	typename ThreadedMembershipFunctionType::Pointer threadedMembershipFunction = ThreadedMembershipFunctionType::New();
	typename ThreadedMembershipFunctionType::DomainType domain;
	domain[0]=0;
	domain[1]= inIndeces.size()-1;
	//std::cout <<  "domain "<< domain[1] << std::endl;
	typename MembershipFunctionType::Pointer hola = (*this->GetMembershipFunctionVector())[0];
	threadedMembershipFunction->SetStuff(samples,inIndeces, outIndeces,hola,n);
	threadedMembershipFunction->Execute(hola ,domain);
	vnl_sparse_matrix<double>* ms= threadedMembershipFunction->GetResults();
	for (unsigned i=0; i<n; i++) 
	{
		
		//for (unsigned j=i; j<n; j++) 
		//{
		//	double val = membershipFunction->Evaluate(&samples->GetMeasurementVector(selected[i]), &samples->GetMeasurementVector(selected[j]));
		//	ms(i,j) = ms(j,i) =  val;
		//}
		diagonal(i,i) =ms->sum_row(i);	
		identity(i,i)=1;
	}
	vnl_sparse_matrix<double> prod(n,n);
	diagonal.subtract(*ms,prod);

	vnl_sparse_symmetric_eigensystem es;
	int res = es.CalculateNPairs(prod, diagonal, n-1, 0.0000001,0,true, true,1000000,-1);//this->GetNumberOfClusters());
	if(res<0)
		std::cout << " ERROR " <<std::endl;

	vnl_vector< double > vector ;
	std::cout <<"e0 " <<  es.get_eigenvalue(0) << "e1 " << es.get_eigenvalue(1) <<std::endl;
	if(es.get_eigenvalue(0)>0.1e-10)
		vector  = es.get_eigenvector(0);
	else
		vector  = es.get_eigenvector(1);
	//double in=0,out=0,maximum=0;
	//int k_i=0;
	int positivos=0, negativos=0;
	for(int i=0;i<n;i++)
	{
		if(vector(i)> 0) //vector(k_i))
		{
			positivos++;
			indices.push_back(  std::pair<int, int>(0,selected[i]));
		}
		else
		{
			negativos++;
			indices.push_back(  std::pair<int, int >(1, selected[i]));
		}	
	}
	std::cout << " positivos " << positivos << " negativos " << negativos << std::endl;
	delete ms;
	return indices;
}
template< class TMesh,class  TMembershipFunctionType>
	std::vector<std::pair<int,int>>	
NormalizedCutsFilter < TMesh ,TMembershipFunctionType>::SelectCentroids(typename SampleType::Pointer samples, const typename MembershipFunctionType::Pointer membershipFunction )
{

	std::vector<std::pair<int,int>> indices;
	std::vector<int> selected;

	const unsigned int n =min(this->GetNumberOfFibersForEigenDecomposition(), (int)samples->Size()); //(int)samples->Size(); //this->GetNumberOfClusters()*10;
	vnl_sparse_matrix<double> ms(n,n);
	vnl_sparse_matrix<double> identity(n,n);
	vnl_sparse_matrix<double> diagonal(n,n);
	bool random=false;//true;
	int offset =(samples->Size()>n)? samples->Size()/n:1;
	for (unsigned i=0; i<n; i++) 
	{
		if(random)
			selected.push_back(rand()%n);
		else
			selected.push_back(i*offset);
	}
	int zeros =0;
	for (unsigned i=0; i<n; i++) 
	{
		
		//double diag =  membershipFunction->Evaluate(&samples->GetMeasurementVector(selected[i]), &samples->GetMeasurementVector(selected[i]));
	
		for (unsigned j=i; j<n; j++) 
		{
			double val = membershipFunction->Evaluate(&samples->GetMeasurementVector(selected[i]), &samples->GetMeasurementVector(selected[j]));
			ms(i,j) = ms(j,i) =  val;
		}
		diagonal(i,i) =ms.sum_row(i);	
		identity(i,i)=1;
	}
	//	std::cout << "zeros " << zeros << std::endl;
	vnl_sparse_matrix<double> prod(n,n);
	diagonal.subtract(ms,prod);
	//	diagonal_inv.mult(prod, ms);
	//	ms.mult(diagonal_inv, prod);
	//	prod.mult(diagonal, prod);
	//	identity.subtract(prod,ms);

	//	ms = prod;	
	//	ms = diagonal - ms ; //identity - diagonal * ms * diagonal;

	//	std::cout << std::endl;
	//	vnl_symmetric_eigensystem<double> ed(md);
	vnl_sparse_symmetric_eigensystem es;
	//	std::cout << " before CalculateNPairs " <<  std::endl;
	//int res = es.CalculateNPairs(prod, diagonal, 2, 0,0,true,true,1000,0.1);//this->GetNumberOfClusters());
	int res = es.CalculateNPairs(prod, diagonal, n-1, 0.0000001,0,true, true,1000000,-1);//this->GetNumberOfClusters());
	//	std::cout << " CalculateNPairs " << res<<  std::endl;
	if(res<0)
		std::cout << " ERROR " <<std::endl;

	//	TEST("vnl_sparse_symmetric_eigensystem::CalculateNPairs() succeeded",			es.CalculateNPairs(ms,nvals), 0);
/*	   if( es.get_eigenvalue(0) > 0.3 && es.get_eigenvalue(1) -es.get_eigenvalue(0) <0.01 )// 0.001	
	   {
	   std::cout << es.get_eigenvalue(1) - es.get_eigenvalue(0) << std::endl;
	   return indices;	
	   }
	   */
	vnl_vector< double > vector ;
	std::cout <<"e0 " <<  es.get_eigenvalue(0) << "e1 " << es.get_eigenvalue(1) <<std::endl;
	if(es.get_eigenvalue(0)>0.1e-10)
		vector  = es.get_eigenvector(0);
	else
		vector  = es.get_eigenvector(1);
	double in=0,out=0,maximum=0;
	int k_i=0;
	int positivos=0, negativos=0;
	for(int i=0;i<n;i++)
	{
//	std::cout << " - " <<vector(i);;
		if(vector(i)> 0) //vector(k_i))
		{
			positivos++;
			indices.push_back(  std::pair<int, int>(0,selected[i]));
			//indices.push_back(  std::pair<int, int>(0, i*offset));
			//indices.push_back(  std::pair<int, int>(0, all[i*offset]));
		}
		else
		{
			negativos++;
			indices.push_back(  std::pair<int, int >(1, selected[i]));
			//indices.push_back(  std::pair<int, int >(1, i*offset));
			//indices.push_back(  std::pair<int, int >(1, all[i*offset]));
		}	
	}
	std::cout << " positivos " << positivos << " negativos " << negativos << std::endl;
	return indices;
	/*
	// Report 'em.
	std::set<int> ind_rep;
	for (unsigned i=0; i<this->GetNumberOfClusters() && indices.size()<this->GetNumberOfClusters(); i++)
	{
	vnl_vector< double > vector  = es.get_eigenvector(i);
	std::cout <<  vector <<  std::endl;
	double min_dist = std::numeric_limits<double>::max();
	int min_index=0;

	for(int j=0;j<n;j++)
	{
	double dist = 0;
	for(int k=0;k<n;k++)
	{
	dist += pow(vector(k)- ms(j,k),2);
	}
	if(dist < min_dist)
	{ 
	std::cout << " dist " << dist << " " << "min_dist " << min_dist << " "   << min_index << std::endl;
	min_dist = dist;
	min_index=  j*offset;

	}
	}
	if(ind_rep.count(min_index)==0)
	indices.push_back(min_index);

	}
	return indices;*/

}
#endif
