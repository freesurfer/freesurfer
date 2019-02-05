#include <iostream>

#include "itkImage.h"
#include "itkVector.h"
#include "itkMesh.h"

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "TrkVTKPolyDataFilter.txx"
#include "vtkSplineFilter.h"
#include "HierarchicalClusteringPruner.h"
typedef std::vector<int>                  PointDataType;
typedef float PixelType;
const unsigned int PointDimension = 3;
const unsigned int MaxTopologicalDimension = 3;
typedef double CoordinateType;
typedef double InterpolationWeightType;
typedef itk::Mesh< PixelType, PointDimension > MeshType;
typedef MeshType::CellType              CellType;
typedef CellType::CellAutoPointer       CellAutoPointer;
typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;
typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
typedef  itk::Image< double,3> ImageType;
typedef ImageType::IndexType IndexType;


std::vector<MeshType::Pointer> FixSampleClusters(std::vector<vtkSmartPointer<vtkPolyData>> polydatas)
{

	std::vector<MeshType::Pointer> meshes;
	for (unsigned int i=0;i<polydatas.size(); i++)
	{
		vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
		spline->SetInput(polydatas[i]);
		spline->SetNumberOfSubdivisions(19);
		spline->Update();

		MeshConverterType::Pointer converter = MeshConverterType::New();
		converter->SetVTKPolyData ( spline->GetOutput() );
		converter->GenerateData2();

		MeshType::Pointer mesh =  converter->GetOutput();
		std::cout << mesh->GetNumberOfCells() << std::endl;
		meshes.push_back(mesh);	
	}
	return meshes;
}

int main(int narg, char*  arg[])
{
	GetPot cl(narg, const_cast<char**>(arg));
	GetPot cl2(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -i anatomicutsFolder -n numClusters -c correspondenceFile -m <numMeasures> <measure1Name> <measure1> ... <measureNName> <measureN> -o output"  << std::endl;   
		return -1;
	}
	int numClusters = cl.follow(0,"-n");
	const char *fileName =cl.follow("output.csv",2,"-o","-O"); 
	const char *fileCorr =cl.follow("output.csv",2,"-c","-C"); 
	std::string hierarchyFilename  = std::string(cl.follow("","-i"));
	
	std::vector<ImageType::Pointer> measures;
	std::vector<std::string> measuresNames;
	int numMeasures = cl.follow(0,"-m");
	for(int i=0;i<numMeasures;i++)
	{
		measuresNames.push_back(std::string( cl.next ("")));
		const char *imFile = cl.next ("");
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		ImageReaderType::Pointer readerS = ImageReaderType::New();
		readerS->SetFileName (imFile);
		readerS->Update();
		ImageType::Pointer image  = readerS->GetOutput();
		measures.push_back(image);
	}

	HierarchicalClusteringPruner<MeshType, ImageType>::Pointer hierarchyPruner =  HierarchicalClusteringPruner<MeshType,ImageType>::New();
	hierarchyPruner->SetNumberOfClusters(numClusters);
	hierarchyPruner->SetHierarchyFilename(hierarchyFilename+"/HierarchicalHistory.csv"); 
	hierarchyPruner->SetClustersPath(hierarchyFilename);
	hierarchyPruner->SetExtension(FiberFormat::TRK);
	hierarchyPruner->SetReferenceImage(measures[0]);
	hierarchyPruner->Update();
	std::vector<vtkSmartPointer<vtkPolyData>> polydatas= hierarchyPruner->GetOutputBundles();
	std::vector<MeshType::Pointer> bundles= FixSampleClusters(polydatas);
	std::vector<long long> meshesIds = hierarchyPruner->GetClustersIds(); 
	std::map<long long,int> bundlesIndeces;
	for(unsigned int i =0;i<meshesIds.size();i++)
	{
		bundlesIndeces[meshesIds[i]]=i;
		std::cout << meshesIds[i]<< " " << i <<  std::endl;
	}
	std::ifstream file( fileCorr); 
	std::string value;
	getline ( file, value, ',' ); 
	getline ( file, value, ',' ); 
	std::vector<long long> correspondences ;
	while ( file.good() )
	{
		getline ( file, value, ',' );
		// long long v1 = atoll(value.c_str());
		getline ( file, value, ',' ); 
		long long v2 = atoll(value.c_str());
		correspondences.push_back(v2);		
		std::cout << v2 << std::endl;
	} 	
	std::ofstream csv_file;
	csv_file.open (fileName);

	csv_file << "Cluster,  N" ;
	for (unsigned int i =0; i< measuresNames.size();i++)
	{
		csv_file <<" , mean"<< measuresNames[i]<< " ,   stde"<< measuresNames[i]<< " ,   meanCentroid"<< measuresNames[i]<< " ,   stdeCentroid"<< measuresNames[i]; 
	}
	csv_file << std::endl;

	for(unsigned int k=0;k<correspondences.size();k++)  
	{
		std::cout << correspondences[k] << " " << bundlesIndeces[correspondences[k]]<< std::endl;
		MeshType::Pointer mesh =  bundles[bundlesIndeces[correspondences[k]]];
		csv_file << correspondences[k] <<  " , " << mesh->GetNumberOfCells() ;

		for(unsigned int m=0;m<measures.size();m++)
		{
			itk::Vector<itk::Point<float,3>,20> avgPoints;
			std::vector<float> FAs;
			std::vector<float>  averageFAs;

			float meanFA=0, meanAvgFA=0;


			MeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();

			for(unsigned  int i=0;i<20;i++)
			{
				MeshType::PointType pt;
				pt.Fill(0.0);
				avgPoints.SetElement(i,pt);
			}
			for(;cells!= mesh->GetCells()->End(); cells++)
			{
				MeshType::CellTraits::PointIdIterator  pointIdIt;
				double dist=0.0;
				double dist_inv=0.0;
				int j=0;
				//std::cout << " num points " << cells.Value()->GetNumberOfPoints() << std::endl;
				for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
				{
					//std::cout<< j<< std::endl;	
					MeshType::PointType pt = 0;
					mesh->GetPoint (*pointIdIt, &pt);
					dist +=	avgPoints[j].EuclideanDistanceTo(pt);
					dist_inv +=avgPoints[avgPoints.Size()-j-1].EuclideanDistanceTo(pt);
					ImageType::IndexType index;
					if( measures[m]->TransformPhysicalPointToIndex(pt,index))
					{
						meanAvgFA += measures[m]->GetPixel(index);		
						averageFAs.push_back( measures[m]->GetPixel(index));
					}
				}
				j=0;
				for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
				{	
					MeshType::PointType pt = 0;
					mesh->GetPoint (*pointIdIt, &pt);
					for(int k=0;k<3;k++)
					{
						if (dist <=dist_inv)
							avgPoints[j][k]+=pt[k]/mesh->GetNumberOfCells();
						else
							avgPoints[avgPoints.Size()-j-1][k]+=pt[k]/mesh->GetNumberOfCells();
					}
				}

			}
			for(unsigned int i=0;i<avgPoints.Size(); i++)
			{	
				ImageType::IndexType index;
				if( measures[m]->TransformPhysicalPointToIndex(avgPoints[i],index))
				{
					float val = measures[m]->GetPixel(index);		
					meanFA += val;		
					//std::cout <<  val << std::endl;
					FAs.push_back( val);
				}
			}
			meanFA/=FAs.size();
			meanAvgFA/=averageFAs.size();
			float stde=0, stdeAvg=0;
			for(unsigned int i=0;i<FAs.size();i++)
			{
				stde+= pow(FAs[i]-meanFA,2);					
				stdeAvg+= pow(averageFAs[i]-meanAvgFA,2);					
			}	
			stde/=FAs.size();
			stdeAvg/=averageFAs.size();
			
			csv_file <<" , " <<  meanFA << " , " << stde;
			csv_file << " , "<< meanAvgFA << " , " <<  stdeAvg;

		}
		csv_file<< std::endl;
	}
	csv_file.close();
}
