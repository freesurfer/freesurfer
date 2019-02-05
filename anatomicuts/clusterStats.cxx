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
int main(int narg, char*  arg[])
{

	enum {Dimension =3};
	typedef double                                                        PixelType;
	typedef itk::Image< PixelType,  Dimension >  ImageType;


	GetPot cl(narg, const_cast<char**>(arg));
	GetPot cl2(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << "-f numberOfBundles <list of vtk bundles> "  << std::endl;   
		std::cout << "-o metricFile.csv -i faImage" << std::endl;
		return -1;
	}
	typedef std::vector<int>                  PointDataType;
	const unsigned int PointDimension = 3;
	const unsigned int MaxTopologicalDimension = 3;
	typedef double CoordinateType;
	typedef double InterpolationWeightType;
	typedef itk::DefaultStaticMeshTraits<
		PointDataType, PointDimension, MaxTopologicalDimension,
		CoordinateType, InterpolationWeightType, int> MeshTraits;
	typedef itk::Mesh< PixelType, PointDimension, MeshTraits > MeshType;
	typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;
	typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
	float faAvg=0, semAvg=0;
	if(cl.search(1,"-i"))
	{
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		const char *faFileName = cl.follow ("", "-i");
		ImageReaderType::Pointer faReader = ImageReaderType::New();
		faReader->SetFileName ( faFileName);
		faReader->Update();
		ImageType::Pointer fa = faReader->GetOutput();

		const char *fileName =cl.follow("output.csv",2,"-o","-O"); 
		std::cout << "FA variance  file " << fileName << std::endl;
		std::ofstream csv_file;
		csv_file.open (fileName);

		int k = cl.follow(0,2,"-f","-F");
		std::vector<MeshType::Pointer> meshes;
		csv_file << " fileName , meanFA , SE , SEM "<< std::endl; 
		int norm =0;
		for(;k>0;k--)  
		{
			vtkPolyDataReader *reader = vtkPolyDataReader::New();
			reader->SetFileName ( cl.next("") );

			reader->GetOutput()->Update();

			MeshConverterType::Pointer converter = MeshConverterType::New();
			converter->SetVTKPolyData ( reader->GetOutput() );
			reader->Delete();
			converter->Update();
			MeshType::Pointer mesh =  converter->GetOutput();
			typedef FixedVTKSamplingFilter<MeshType,MeshType> SamplingFilterType;
			SamplingFilterType::Pointer samplingFilter =  SamplingFilterType::New();
			samplingFilter->SetInput(mesh);
			samplingFilter->SetSampling(10);
			samplingFilter->Update();
			mesh= samplingFilter->GetOutput();


			std::vector<float> fas;
			float meanFA=0;
			MeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();
			std::vector<MeshType::PointType> avgPoints(10,0);
			for(;cells!= mesh->GetCells()->End(); cells++)
			{
				MeshType::CellTraits::PointIdIterator  pointIdIt;
				double dist=0.0;
				double dist_inv=0.0;
				int j=0;
				for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
				{	
					MeshType::PointType pt;
					mesh->GetPoint (*pointIdIt, &pt);
					dist +=	avgPoints[j].EuclideanDistanceTo(pt);
					dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
				}
				j=0;
				for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
				{	
					MeshType::PointType pt;
					mesh->GetPoint (*pointIdIt, &pt);
					for(int k=0;k<3;k++)
					{
						if (dist <dist_inv)
							avgPoints[j][k]+=pt[k]/mesh->GetNumberOfCells();
						else
							avgPoints[avgPoints.size()-j-1][k]+=pt[k]/mesh->GetNumberOfCells();
					}
				}

			}
			for(int i=0;i<avgPoints.size(); i++)
			{	
				ImageType::IndexType index;
				if( fa->TransformPhysicalPointToIndex(avgPoints[i],index))
				{
					meanFA += fa->GetPixel(index);		
					//std::cout << index << " " << pt <<" " <<fa->GetPixel(index)<< std::endl;
					fas.push_back( fa->GetPixel(index));
				}
			}
			meanFA/=fas.size();
			csv_file << reader->GetFileName() << " , " << meanFA;
			float stde=0;
			for( int i=0;i<fas.size();i++)
			{
				stde+= pow(fas[i]-meanFA,2);					
			}	

			csv_file << " , "<< sqrt(stde/fas.size()) << " , " << sqrt(stde)/fas.size() << std::endl;
			//			std::cout << reader->GetFileName()<< " "<< meanFA << std::endl;
			//if( meanFA > 0)
			{	faAvg+=meanFA;
				semAvg += sqrt(stde)/fas.size();
				norm++;
			}
		}
		csv_file.close();

		std::cout << " meanFA " << faAvg/norm << "SEM  "<< semAvg/norm << std::endl;
	}
}
