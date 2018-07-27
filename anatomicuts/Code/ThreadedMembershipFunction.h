#ifndef _ThreadedMembershipFunction_h
#define _ThreadedMembershipFunction_h

#include "itkDomainThreader.h"
#include "vnl/vnl_sparse_matrix.h"
#include "itkThreadedIndexedContainerPartitioner.h"

template<class TMembershipFunctionType> 
class ThreadedMembershipFunction :  public itk::DomainThreader<itk::ThreadedIndexedContainerPartitioner, TMembershipFunctionType>
{
	public :
		using Self = ThreadedMembershipFunction;
		using Superclass =  itk::DomainThreader<itk::ThreadedIndexedContainerPartitioner,TMembershipFunctionType>;
		using Pointer =  itk::SmartPointer<Self>;
		using ConstPointer = itk::SmartPointer<const Self>;

		using DomainType = typename Superclass::DomainType;
		itkNewMacro(Self);
		typedef TMembershipFunctionType MembershipFunctionType;
		typedef typename TMembershipFunctionType::MeasurementVectorType MeasurementVectorType; 
		typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;

		void SetStuff(typename SampleType::Pointer samples, std::vector<std::pair<int,int>> indeces, std::vector<std::pair<int,int>> outIndeces, typename MembershipFunctionType::Pointer msf,int n) 
		{
			m_samples = samples;
			m_indeces = indeces;
			m_outIndeces = outIndeces;
			m_membershipFunction =  msf;
			m_matrixDim = n;
		}
		vnl_sparse_matrix<double>* GetResults();
		std::vector<int> GetMaxIndeces(); //{return this->m_maxIndex;}

	protected:
		ThreadedMembershipFunction(){}
		~ThreadedMembershipFunction(){}

	private:
		int m_matrixDim;
		typename SampleType::Pointer  m_samples;
		std::vector<std::pair<int,int>> m_indeces; 
		std::vector<std::pair<int,int>> m_outIndeces; 
		//std::vector<vnl_sparse_matrix<double>*> m_results;
		std::vector<std::vector<int>> m_maxIndex;
		std::vector<std::vector<double>> m_maxValue;
		int* m_results2;
		typename MembershipFunctionType::Pointer m_membershipFunction;
		void BeforeThreadedExecution();
		void ThreadedExecution(const DomainType&, const itk::ThreadIdType);
		void AfterThreadedExecution();

};
#include "ThreadedMembershipFunction.txx"
#endif

