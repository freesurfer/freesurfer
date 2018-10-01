#ifndef __NormalizedCutsFilter_h
#define __NormalizedCutsFilter_h

#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkListSample.h"
#include "itkArray.h"
#include "itkVector.h"
//#include "itkImageKmeansModelEstimator.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkKMeansClassifierFilter.h"
#include "itkVariableLengthVector.h"
#include "itkPolylineCell.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkMeshToMeshFilter.h"
#include "ThreadedMembershipFunction.h"
#if ITK_VERSION_MAJOR < 4
#include "itkMaximumDecisionRule2.h"
#else
#include "itkMaximumDecisionRule.h"
#endif
using namespace itk;
using namespace Statistics;
template <class TMesh,class  TMembershipFunctionType>
class NormalizedCutsFilter :
	public ProcessObject
	//  public MeshToMeshFilter<TInputMesh, TOutputMesh>
{
	public: 
		typedef NormalizedCutsFilter       Self;
		typedef SmartPointer<Self>                              Pointer;
		typedef SmartPointer<const Self>                        ConstPointer;

		itkNewMacro  (Self);
		itkTypeMacro (NormalizedCutsFilter,ProcessObject);

		typedef TMesh                         MeshType;
		typedef typename MeshType::MeshTraits MeshTraits;
		typedef typename MeshType::PointType  PointType;
		typedef typename MeshType::PixelType  PixelType;
		typedef typename MeshType::CellType   CellType;
		typedef typename MeshType::Pointer   MeshPointerType;

		typedef typename std::vector<MeshPointerType>  ListOfOutputMeshTypePointer;     

		typedef typename MeshType::CellTraits      MeshCellTraits;
		typedef typename MeshType::CellIdentifier  MeshCellIdentifier;
		typedef typename MeshType::CellType        MeshCellType;

		typedef typename MeshType::CellAutoPointer MeshCellAutoPointer;
		typedef typename MeshType::PointIdentifier MeshPointIdentifier;
		typedef typename MeshCellTraits::PointIdIterator   MeshPointIdIterator;

		typedef typename MeshType::PointsContainerPointer       PointsContainerPointer;

		typedef typename MeshType::PointsContainer       PointsContainer;

		typedef typename MeshType::CellsContainer      CellsContainer;

		typedef typename MeshType::CellsContainerPointer      CellsContainerPointer;
		typedef PolylineCell<MeshCellType>                  PolylineCellType;


		typedef TMembershipFunctionType MembershipFunctionType;
		typedef typename MembershipFunctionType::Pointer                      MembershipFunctionPointer;
		typedef typename MembershipFunctionType::MeasurementVectorType MeasurementVectorType;

		typedef ThreadedMembershipFunction<MembershipFunctionType> ThreadedMembershipFunctionType;
	

		typedef ListSample< MeasurementVectorType > SampleType;
		typedef WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
		typedef typename TreeGeneratorType::KdTreeType TreeType;
		typedef KdTreeBasedKmeansEstimator<TreeType> EstimatorType;


#if ITK_VERSION_MAJOR < 4
		typedef MaximumDecisionRule2 DecisionRuleType;
#else
		typedef MaximumDecisionRule DecisionRuleType;
#endif

	//	typedef KMeansClassifierFilter<SampleType,MembershipFunctionType> ClassifierType;
//		typedef typename ClassifierType::ClassLabelVectorType                     ClassLabelVectorType;
	//	typedef typename ClassifierType::MembershipFunctionVectorType             MembershipFunctionVectorType;
	        typedef typename MembershipFunctionType::Pointer MembershipFunctionTypePointer;	
		typedef typename std::vector<MembershipFunctionTypePointer>             MembershipFunctionVectorType;



		ListOfOutputMeshTypePointer GetOutput()
		{
			return this->m_Output;
		}

		void SetNumberOfClusters(int n)
		{ this->numberOfClusters = n;}

		int GetNumberOfClusters()
		{
			return this->numberOfClusters;
		}
		void SetNumberOfFibersForEigenDecomposition(int e)
		{
			this->m_numberOfFibersForEigenDecomposition = e;
		}
		int GetNumberOfFibersForEigenDecomposition()
		{
			return m_numberOfFibersForEigenDecomposition;
		}

		std::vector<std::string> GetLabels()
		{ return this->labels;}
		void SetLabels(std::vector<std::string> labels)
		{ this->labels = labels; }

		MeshPointerType GetInput()
		{ return this->input; }
		void SetInput(MeshPointerType input)
		{	this->input = input ; }
		//void Update();
		virtual void Update(void);    

		itkGetMacro( SigmaCurrents, int );
		//  itkSetMacro( SigmaCurrents, int );
		void SetSigmaCurrents(int s){	this->m_SigmaCurrents =s; }

		itkSetMacro( NumberOfIterations, unsigned int );
		itkGetMacro( NumberOfIterations, unsigned int );

		MembershipFunctionVectorType* GetMembershipFunctionVector()
		{
			return this->m_membershipFunctions;
		}
		void SetMembershipFunctionVector(MembershipFunctionVectorType *functions)
		{
			this->m_membershipFunctions = functions;
		}

		std::vector<std::pair<std::string,std::string>> GetClusterIdHierarchy()
		{
			return this->m_clusterIdHierarchy;
		}
	protected:
		

		std::vector<std::pair<int,int>> SelectCentroids(typename SampleType::Pointer samples, const typename MembershipFunctionType::Pointer);
		std::vector<std::pair<int,int>> SelectCentroidsParallel(typename SampleType::Pointer samples, const typename MembershipFunctionType::Pointer);
		MeshPointerType input;
		std::vector<std::string> labels;
		ListOfOutputMeshTypePointer m_Output;
		int numberOfClusters;
		NormalizedCutsFilter() {}
		~NormalizedCutsFilter() {}

		//    virtual void GenerateData (void);
	private:
		unsigned int 		 		m_NumberOfIterations; 
		NormalizedCutsFilter (const Self&);
		std::vector<std::pair<std::string,std::string>> m_clusterIdHierarchy;
		void operator=(const Self&);    
		int m_SigmaCurrents;
		int m_numberOfFibersForEigenDecomposition;
//		void SaveClustersInMeshes(MembershipFunctionVectorType mfv);
		MembershipFunctionVectorType *m_membershipFunctions; 
};  
#include "NormalizedCutsFilter.txx"
#endif
