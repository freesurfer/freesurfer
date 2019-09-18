
#include <fstream>
#include <iostream>

#include "itkMesh.h"
#include "vtkSplineFilter.h"


#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "itkPolylineCell.h"
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <cmath>
#include "itkArray.h"
#include "itkPolylineCell.h"
#include "GetPot.h"
#include "TrkVTKPolyDataFilter.txx"
#include "itkImage.h"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "LabelPerPointVariableLengthVector.h"
#include "EuclideanMembershipFunction.h"
#include "ClusterTools.h"
#include "itkDefaultStaticMeshTraits.h"


int main(int argc, char *argv[])
{		
	enum {Dimension =3};
	typedef int                                                        PixelType;
	const unsigned int PointDimension = 3;
	typedef std::vector<int>                  PointDataType;
	const unsigned int MaxTopologicalDimension = 3;
	typedef double CoordinateType;
	typedef double InterpolationWeightType;
	typedef itk::DefaultStaticMeshTraits<
		PointDataType, PointDimension, MaxTopologicalDimension,
		CoordinateType, InterpolationWeightType, PointDataType > MeshTraits;
	typedef itk::Mesh< PixelType, PointDimension, MeshTraits > HistogramMeshType;

	typedef itk::Image<float, 3> ImageType;
	
	typedef itk::Mesh< PixelType, PointDimension > ColorMeshType;
	typedef ColorMeshType::PointType PointType;
	typedef ColorMeshType::CellType        CellType;
	typedef itk::PolylineCell<CellType>                      PolylineCellType;
	//typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
	// typedef MeshType::CellsContainer::ConstIterator CellIteratorType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;

	GetPot cl(argc, const_cast<char**>(argv));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout << argv[0] << " -i streamlines -o streamlines -l minLength  -m mask/segmentation -nu (filter ushape fibers) -s subsample -c cleaningThreshold (std times) -d outputDirectory" << std::endl;
		return EXIT_FAILURE;
	}
	const char* outputDirectory = cl.follow("",2,"-d","-D");
	int minLenght  = cl.follow(-1,2,"-min","-MIN");
	int maxLenght  = cl.follow(500,2,"-max","-MAX");
	int subSample = cl.follow(-1,2,"-s","-S");
	bool filterUShape =cl.search("-nu");
	bool maskFibers  =cl.search("-m");
	float cleaningThreshold = cl.follow(0.0,2,"-c","-C");
	std::cout << " clean " << cleaningThreshold << std::endl;
	ImageType::Pointer refImage;
	ImageType::Pointer mask;

	if(cl.search(2,"-m","-M"))
	{
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName ( cl.next(""));
		reader->Update();
		mask = reader->GetOutput();
	}

	std::vector<std::string> inputFiles;
	for(std::string inputName = std::string(cl.follow("",2,"-i","-I")); access(inputName.c_str(),0)==0; inputName = std::string(cl.next("")))
	{
		inputFiles.push_back(inputName);
	}

	std::vector<std::string> outputFiles;
	for(std::string outputName = std::string(cl.follow("",2,"-o","-O")); outputName.size()>4; outputName = std::string(cl.next("")))
	{
		outputFiles.push_back(outputName);
		std::cout << outputName << std::endl;
	}
	std::vector<ColorMeshType::Pointer>* meshes;
	std::vector<ColorMeshType::Pointer>* fixMeshes;
	std::vector<vtkSmartPointer<vtkPolyData>> polydatas;
	
	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();

	std::vector<HistogramMeshType::Pointer>* histoMeshes;
	clusterTools->GetPolyDatas(inputFiles, &polydatas, mask) ;
	bool clean =(cleaningThreshold>0)?true:false; 
	std::cout << " clean " << clean << std::endl;
	if( clean)
	{	
		meshes = clusterTools->PolydataToMesh(polydatas);
		fixMeshes = clusterTools->FixSampleClusters(polydatas, 20);
		histoMeshes = clusterTools->ColorMeshToHistogramMesh(*fixMeshes, mask, false);
	
		clusterTools->SetDirectionalNeighbors(histoMeshes,  mask, *clusterTools->GetDirections(DirectionsType::ALL), false);

	
		std::cout << polydatas.size() <<  " " << fixMeshes->size() << std::endl;
	}
	else
	{
		meshes = clusterTools->PolydataToMesh(polydatas);
	}
		
	for(int i=0; i<meshes->size(); i++)
	{ 
		ColorMeshType::Pointer input = (*meshes)[i];
		int offset =(subSample>0)? input->GetNumberOfCells()/subSample:1;

		int averageId =0;
		float stdCluster=0;
		if ( clean)
		{
			std::cout << "hola "<< std::endl; 
			int averageId = clusterTools->GetAverageStreamline((*fixMeshes)[i]);	
			stdCluster =  clusterTools->GetStandardDeviation((*histoMeshes)[i],averageId)*cleaningThreshold;
			std::cout << (*fixMeshes)[i]->GetNumberOfCells() << " " <<(*histoMeshes)[i]->GetNumberOfCells() << " " << stdCluster << " " <<cleaningThreshold << std::endl;
		}

		std::set<int> unfilteredIds;

		int pointIndices =0;
		int cellIndices = 0;
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
		for (int cellId=0 ; inputCellIt!=input->GetCells()->End(); ++inputCellIt, cellId++)
		{

			PointType firstPt;
			firstPt.Fill(0);
			float val1=0, val2=0;

			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
			input->GetPoint(*it,&firstPt);
			double lenghtSoFar = 0;
			for( ; it!=inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				pt.Fill(0);
				input->GetPoint(*it,&pt);
				lenghtSoFar+= firstPt.EuclideanDistanceTo(pt);
				input->GetPoint(*it,&firstPt);
				if(filterUShape || maskFibers)
				{
					ImageType::IndexType  index ;
					float value=0;
					if( mask->TransformPhysicalPointToIndex(firstPt,index))
					{
						value= mask->GetPixel(index);
						//std::cout <<  value << std::endl;
						if (val1 ==0 && value!=0)
							val1 =value;

						if( value!=0)
							val2 = value;
					}		
				}
			}

			float dist= clean?clusterTools->GetDistance((*histoMeshes)[i],averageId, cellId):0;	
//			std::cout << dist << " " << stdCluster << std::endl;	
			if(lenghtSoFar >= minLenght  && lenghtSoFar <= maxLenght && cellId % offset ==0 &&( val1!=val2 || !filterUShape) && ( val1!= 0 || !maskFibers) && (dist<=stdCluster))
			{	

				unfilteredIds.insert(cellId);
			}

		}
	
		ColorMeshType::Pointer om = ColorMeshType::New();
		om->SetCellsAllocationMethod( ColorMeshType::CellsAllocatedDynamicallyCellByCell );

		inputCellIt = input->GetCells()->Begin();
		for (int cellId=0 ; inputCellIt!=input->GetCells()->End(); ++inputCellIt, cellId++)
		{

			if (unfilteredIds.count(cellId)>0)
			{	
				CellAutoPointer line;
				line.TakeOwnership ( new PolylineCellType);
				int k=0;
				CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
				for( ; it!=inputCellIt.Value()->PointIdsEnd(); it++)
				{
					PointType pt;
					input->GetPoint (*it, &pt);

					om->SetPoint (pointIndices, pt);
					line->SetPointId (k,pointIndices);

					k++;
					pointIndices++;

				}
				om->SetCell (cellIndices, line);
				ColorMeshType::CellPixelType cellData;
				input->GetCellData(cellId, &cellData);
				om->SetCellData(cellIndices, cellData) ;
				cellIndices++;
			}

		}
		std::string outputName;
		if(outputFiles.size() == inputFiles.size())
		{		
			outputName = outputFiles[i];
		}
		else
		{	
			std::string filename = inputFiles[i].substr(inputFiles[i].find_last_of("/\\") + 1);	
			outputName= std::string(outputDirectory)+ "/"+filename;
		}
	
		std::cout << " mesh name  " << outputName << std::endl;

		clusterTools->SaveMesh(om, mask, outputName, inputFiles[i]); 
	}
	delete meshes;
	if(clean)
		delete histoMeshes;
}

