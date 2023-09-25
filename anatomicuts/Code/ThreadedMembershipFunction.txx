#ifndef _ThreadedMembershipFunction_txx
#define _ThreadedMembershipFunction_txx


#include "ThreadedMembershipFunction.h"
#include <iostream>
#include "vnl/vnl_matrix.h"
#include <stdlib.h> 

template< class  TMembershipFunctionType> void
ThreadedMembershipFunction< TMembershipFunctionType >::BeforeThreadedExecution() 
{
#if ITK_VERSION_MAJOR >= 5
	const itk::ThreadIdType numberOfThreads = this->GetNumberOfWorkUnitsUsed();
#else
	const itk::ThreadIdType numberOfThreads = this->GetNumberOfThreadsUsed();
#endif	
	//std::cout << " number of threads "<< numberOfThreads << std::endl;
	//this->m_results.resize( numberOfThreads );
	this->m_maxIndex.resize( numberOfThreads );
	this->m_maxValue.resize( numberOfThreads );
	for( itk::ThreadIdType ii = 0; ii < numberOfThreads; ++ii )
	{
		this->m_maxIndex[ii].resize(m_matrixDim,0);
		this->m_maxValue[ii].resize(m_matrixDim,0);
//		this->m_results[ii] = new vnl_sparse_matrix<double>(m_matrixDim, m_matrixDim);
	}
	this->m_results2 = new int[m_indeces.size()];

}
template< class  TMembershipFunctionType> void
ThreadedMembershipFunction< TMembershipFunctionType >::ThreadedExecution(const DomainType& subDomain, const itk::ThreadIdType threadId) 
{
	for( itk::IndexValueType ii = subDomain[0]; ii <= subDomain[1]; ++ii )
	{	
		int i = m_indeces[ii].first;// [0];
		int j= m_indeces[ii].second; //[1];
		double val = m_membershipFunction->Evaluate(&m_samples->GetMeasurementVector(i), &m_samples->GetMeasurementVector(j));
		i = m_outIndeces[ii].first;// [0];
		j= m_outIndeces[ii].second; //[1];
		//(*m_results[threadId])(i,j)=(*m_results[threadId])(j,i)= val;
		m_results2[ii]=val;
		if( val > m_maxValue[threadId][i])
		{
			m_maxValue[threadId][i]=val;
			m_maxIndex[threadId][i]=j;
		}
		//m_results2[j * m_matrixDim +i]= val;
	}
}

template< class  TMembershipFunctionType> void
ThreadedMembershipFunction< TMembershipFunctionType >::AfterThreadedExecution() 
{

}
template< class  TMembershipFunctionType> std::vector<int>
ThreadedMembershipFunction< TMembershipFunctionType >::GetMaxIndeces() 
{
#if ITK_VERSION_MAJOR >= 5
	const itk::ThreadIdType numberOfThreads = this->GetNumberOfWorkUnitsUsed();
#else
	const itk::ThreadIdType numberOfThreads = this->GetNumberOfThreadsUsed();
#endif	
	std::vector<int> indeces(m_matrixDim,0);
	std::vector<double> values(m_matrixDim,0);
	for (int i =0; i<m_matrixDim;i++)
	{
		for( itk::ThreadIdType ii = 0; ii < numberOfThreads; ++ii )
		{
			if(m_maxValue[ii][i] > values[i])
			{
				values[i]= m_maxValue[ii][i];
				indeces[i]=m_maxIndex[ii][i];
			}
		}
	}
	return indeces;
}
template< class  TMembershipFunctionType> 
vnl_sparse_matrix<double>* ThreadedMembershipFunction< TMembershipFunctionType >::GetResults()
{
	//std::cout << " get results start " << std::endl;
	vnl_sparse_matrix<double>*res = new vnl_sparse_matrix<double>(m_matrixDim,m_matrixDim);
	//const itk::ThreadIdType numberOfThreads = this->GetNumberOfThreadsUsed();
	//this->m_results.resize( numberOfThreads );
	/*for( itk::ThreadIdType ii = 0; ii < numberOfThreads; ++ii )
	{
		(*res)+= (*m_results[ii]);
		delete m_results[ii];
	}*/

	for(int k=0;k<m_indeces.size();k++)
	{
		int i= m_outIndeces[k].first;
		int j= m_outIndeces[k].second;
		(*res)(i,j)=(*res)(j,i)=this->m_results2[k];
	}
	//std::cout << " get results end" << std::endl;
	return res;

}

#endif
