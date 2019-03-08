
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
#include "VTKPolyDataToPolylineMeshFilter.h"

int main(int argc, char *argv[])
{
	typedef itk::Mesh<double, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef MeshType::CellType        CellType;
	typedef itk::PolylineCell<CellType>                      PolylineCellType;
	typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
	// typedef MeshType::CellsContainer::ConstIterator CellIteratorType;
	typedef MeshType::CellAutoPointer CellAutoPointer;
	typedef itk::Image<float, 3> ImageType;
	typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;

	GetPot cl(argc, const_cast<char**>(argv));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout << argv[0] << " -i streamlines -o streamlines -l minLength  -m mask/segmentation -nu (filter ushape fibers) -s subsample" << std::endl;
		return EXIT_FAILURE;
	}
	const char* inputName= cl.follow("",2,"-i","-I");
	const char* outputName = cl.follow("",2,"-o","-O");
	int maxLenght  = cl.follow(-1,2,"-l","-L");
	int subSample = cl.follow(-1,2,"-s","-S");
	bool filterUShape =cl.search("-nu");
	bool maskFibers  =cl.search("-m");

	ImageType::Pointer refImage;
	ImageType::Pointer mask;

	/*if(cl.search(2,"-r","-R"))
	{
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName ( cl.next(""));
		reader->Update();
		refImage = reader->GetOutput();
	}
	*/
	if(cl.search(2,"-m","-M"))
	{
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName ( cl.next(""));
		reader->Update();
		mask = reader->GetOutput();
	}

	MeshConverterType::Pointer converter = MeshConverterType::New();
	if(  std::string(inputName).find(".trk") !=std::string::npos)
	{
		itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
		trkReader->SetTrkFileName(inputName);
		if( cl.search(2,"-m","-M"))
		{
			trkReader->SetReferenceImage(mask);
		}
		trkReader->TrkToVTK();
		converter->SetVTKPolyData ( trkReader->GetOutputPolyData() );
	}
	else
	{
		vtkSmartPointer<vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
		vtkReader->SetFileName ( inputName);
		vtkReader->Update();
		converter->SetVTKPolyData ( vtkReader->GetOutput() );
	}
	
	converter->Update();
	MeshType::Pointer input = converter->GetOutput();
	MeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
        int offset =(subSample>0)? input->GetNumberOfCells()/subSample:1;

	MeshType::Pointer om = MeshType::New();
	om->SetCellsAllocationMethod( MeshType::CellsAllocatedDynamicallyCellByCell );

	int pointIndices =0;
	int cellIndices = 0;
	for (int cellId=0 ; inputCellIt!=input->GetCells()->End(); ++inputCellIt, cellId++)
	{

		PointType firstPt;
		firstPt.Fill(0);
		float val1=0, val2=0;

		CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
		input->GetPoint(*it,&firstPt);
		double lenghtSoFar = 0;
		bool masked = false;
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
					else
						masked=true;
				}		
			}
		}
		
		if(lenghtSoFar >= maxLenght &&  cellId % offset ==0 &&( (val1!=val2 && val1!= 0) || !filterUShape) &&(( !masked )|| !maskFibers))
		{	
			CellAutoPointer line;
			line.TakeOwnership ( new PolylineCellType);
			int k=0;
			it = inputCellIt.Value()->PointIdsBegin();
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
			cellIndices++;
		}

	}
	std::cout <<" output file: "<<  outputName<< om->GetNumberOfCells() <<std::endl;
	char meshName2[100];
	sprintf(meshName2, "%s", outputName);

	std::cout << " mesh name  " << meshName2 << std::endl;
	VTKConverterType::Pointer vtkConverter =  VTKConverterType::New();
	vtkConverter->SetInput(om);
	vtkConverter->Update();
	if(  std::string(outputName).find(".trk") != std::string::npos)
	{
		itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
		trkReader->SetInput(vtkConverter->GetOutputPolyData());
		if( cl.search(2,"-m","-M"))
		{
			trkReader->SetReferenceImage(mask);
		}
		trkReader->SetReferenceTrack(	std::string(inputName));
		trkReader->VTKToTrk(std::string(outputName));

	}
	else
	{
		vtkPolyDataWriter *writerFixed = vtkPolyDataWriter::New();
		writerFixed->SetFileName ( meshName2);
#if VTK_MAJOR_VERSION > 5
		writerFixed->SetInputData(vtkConverter->GetOutputPolyData());
#else
		writerFixed->SetInput(vtkConverter->GetOutputPolyData());
#endif
		writerFixed->SetFileTypeToBinary();
		writerFixed->Update();
	}

}

