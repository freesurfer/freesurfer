#ifndef _HierarchicalClusteringPruner_h_
#define _HierarchicalClusteringPruner_h_

#include "itkMeshToMeshFilter.h"
#include "AppendBundleFilter.h"
enum FiberFormat {VTK=0, TRK=1};

template < class TOutputMesh,class TImageType>
class HierarchicalClusteringPruner: public MeshSource<TOutputMesh> 
{
	public:
		typedef HierarchicalClusteringPruner Self;
		typedef MeshSource<TOutputMesh>       Superclass;
		typedef SmartPointer<Self>                              Pointer;
		typedef SmartPointer<const Self>                        ConstPointer;

		itkNewMacro  (Self);
		itkTypeMacro (HierarchicalClusteringPruner, MeshSource);

		typedef TImageType				ImageType;
		typedef TOutputMesh                         OutputMeshType;
		typedef typename OutputMeshType::MeshTraits OutputMeshTraits;
		typedef typename OutputMeshType::PointType  OutputPointType;
		typedef typename OutputMeshType::PixelType  OutputPixelType;


		/** Some convenient typedefs. */
		typedef typename OutputMeshType::Pointer         OutputMeshPointer;
		typedef typename OutputMeshType::CellTraits      OutputCellTraits;
		typedef typename OutputMeshType::CellIdentifier  OutputCellIdentifier;
		typedef typename OutputMeshType::CellType        OutputCellType;
		typedef typename OutputMeshType::CellAutoPointer OutputCellAutoPointer;
		typedef typename OutputMeshType::PointIdentifier OutputPointIdentifier;
		typedef typename OutputCellTraits::PointIdIterator     OutputPointIdIterator;

		typedef typename OutputMeshType::PointsContainerPointer
			OutputPointsContainerPointer;

		typedef typename OutputMeshType::PointsContainer
			OutputPointsContainer;

		typedef typename OutputMeshType::CellsContainer
			OutputCellsContainer;

		typedef typename OutputMeshType::CellsContainerPointer
			OutputCellsContainerPointer;

		typedef PolylineCell<OutputCellType>                      PolylineCellType;


		int GetNumberOfClusters()
		{
			return this->m_numberOfClusters;
		}

		void SetNumberOfClusters(int n)
		{
			this->m_numberOfClusters = n;
		}
		void SetHierarchyFilename(std::string filename)
		{
			this->m_hierarchyFilename = filename;
		}
		std::string GetHierarchyFilename()
		{
			return this->m_hierarchyFilename;
		}
		void SetClustersPath(std::string path)
		{
			this->m_clustersPath = path;
		}
		std::string GetClustersPath()
		{
			return this->m_clustersPath;
		}
	 	void SetExtension(FiberFormat f)
		{
			this->m_fiberFormat = f;
		}

		FiberFormat GetExtension()
		{
			return this->m_fiberFormat;
		}
		std::vector<vtkSmartPointer<vtkPolyData>> GetOutputBundles()
		{
			return this->m_outputBundles;
		}
		std::vector<long long> GetClustersIds()
		{
			return this->m_clustersIds;
		}
		void SetReferenceImage(typename ImageType::Pointer image)
		{
			this->m_referenceImage=image;
		}	
			
	protected:
	
		HierarchicalClusteringPruner();
		~HierarchicalClusteringPruner() {};

		virtual void GenerateData (void);

	private:
		int m_numberOfClusters;
	        std::string m_hierarchyFilename;
		std::string m_clustersPath;
		FiberFormat m_fiberFormat;
	        std::vector<vtkSmartPointer<vtkPolyData>> m_outputBundles;	
		std::vector<long long> m_clustersIds;
		typename ImageType::Pointer m_referenceImage;
};


#include "HierarchicalClusteringPruner.txx"

#endif
