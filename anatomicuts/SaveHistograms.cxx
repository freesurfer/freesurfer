#include <iostream>
#include "itkImage.h"
#include "itkVector.h"
#include "itkMesh.h"

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include "TrkVTKPolyDataFilter.txx"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "OrientationPlanesFromParcellationFilter.h"
#include "vtkSplineFilter.h"


int main(int narg, char*  arg[])
{

	enum {Dimension =3};
	typedef double                                                        PixelType;
	typedef itk::Image< PixelType,  Dimension >  ImageType;

	typedef ImageType::IndexType IndexType;


	GetPot cl(narg, const_cast<char**>(arg));
	GetPot cl2(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -p parcellation -f numberOfBundles <list of vtk bundles> -o output.csv -bb"  << std::endl;   
		return -1;
	}
	const unsigned int PointDimension = 3;
	typedef itk::Mesh< PixelType, PointDimension > MeshType;
	typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;
	typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;

	const char* filename =cl.follow("histograms.csv",2,"-O","-o"); 
	std::cout << "Output histogram file " << filename << std::endl;
	std::ofstream csv_file;
	csv_file.open (filename);
	bool bb= cl.search("-bb");
	const char *segFile = cl2.follow ("", "-p");

	typedef itk::ImageFileReader<ImageType> ImageReaderType;
	ImageReaderType::Pointer readerS = ImageReaderType::New();
	readerS->SetFileName ( segFile);
	readerS->Update();
	ImageType::Pointer segmentation  = readerS->GetOutput();
	OrientationPlanesFromParcellationFilter<ImageType,ImageType>::Pointer orientationFilter = OrientationPlanesFromParcellationFilter<ImageType, ImageType>::New();
	orientationFilter->SetInput(segmentation);
	orientationFilter->SetBabyMode(bb);
	orientationFilter->Update();
	std::vector<itk::Vector<float>> orientations;
	orientations.push_back(orientationFilter->GetUpDown());
	orientations.push_back(orientationFilter->GetFrontBack());
	orientations.push_back(orientationFilter->GetLeftRight());


	std::vector<IndexType> indeces;
	std::vector<itk::Vector<float>> direcciones;
	int possibles[3] = {0,1,-1};
	int numFiles = cl.follow(0,2,"-f", "-F");
	const char* fiberFile= cl.next("");
	for(int i=0;i<3;i++)
	{
		for(int k=0;k<3;k++)
		{
			for(int j=0;j<3;j++)
			{
				IndexType index;
				index[0] = possibles[i];
				index[1] = possibles[j];
				index[2] = possibles[k];

				itk::Vector<float> dir;
				for(int w=0;w<3;w++)
				{

					dir[w] = possibles[i] * orientations[0][w]+possibles[j] * orientations[1][w]+possibles[k] * orientations[2][w];
				}
				dir.Normalize();
				int howManyZeros=0;
				if(i==0)
					howManyZeros++;
				if(j==0)
					howManyZeros++;
				if(k==0)
					howManyZeros++;

				if(howManyZeros!=3)
				{
					direcciones.push_back(dir);
					indeces.push_back(index);
				}

			}
		}
	}

	while(numFiles)  
	{
		MeshConverterType::Pointer converter = MeshConverterType::New();
		vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
		spline->SetNumberOfSubdivisions(10);

		if(  std::string(fiberFile).find(".trk") !=std::string::npos)
		{
			itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
			trkReader->SetTrkFileName(fiberFile);
			//if(segmentation.IsNotNull())
			trkReader->SetReferenceImage(segmentation);
			trkReader->TrkToVTK();
#if VTK_MAJOR_VERSION > 5
			spline->SetInputData(trkReader->GetOutputPolyData());
#else
			spline->SetInput(trkReader->GetOutputPolyData());
#endif
		}
		else
		{
			vtkSmartPointer<vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
			vtkReader->SetFileName ( fiberFile);
			vtkReader->Update();
#if VTK_MAJOR_VERSION > 5
			spline->SetInputConnection( vtkReader->GetOutputPort() );
#else
			spline->SetInput( vtkReader->GetOutput() );
#endif
		}
		spline->Update();
		converter->SetVTKPolyData ( spline->GetOutput() );
		converter->GenerateData2();

		MeshType::Pointer mesh =  converter->GetOutput();

		MeshType::CellsContainer::ConstIterator itCell = mesh->GetCells()->Begin();
		int cellId = 0;

		std::vector<std::map<int,int>> histograms;
		for( int i=0;i<direcciones.size();i++)
			histograms.push_back(std::map<int,int>());

		for( ;itCell != mesh->GetCells()->End();itCell++)
		{
			MeshType::CellTraits::PointIdIterator  pointIdIt  = itCell.Value()->PointIdsBegin();
			for(;pointIdIt != itCell.Value()->PointIdsEnd(); pointIdIt++)
			{
				MeshType::PointType pt;
				mesh->GetPoint (*pointIdIt, &pt);
				MeshType::PointType lpt;

				IndexType index;	
				segmentation->TransformPhysicalPointToIndex(pt,index);				
				PixelType label = segmentation->GetPixel(index);
				if (histograms[0].count(label)>0)
					histograms[0][label] ++;
				else
					histograms[0][label]=1;

				for(int i=0;i<direcciones.size();i++)
				{
					PixelType vecino = label;
					itk::ContinuousIndex<float,3> continuousIndex = index;
					IndexType roundedIndex;
					MeshType::PointType point = pt;
					while(vecino == label )
					{
						for(int j=0; j<3;j++)
							continuousIndex[j] += direcciones[i][j];
						roundedIndex.CopyWithRound(continuousIndex);
						if(!segmentation->GetLargestPossibleRegion().IsInside(roundedIndex))
							break;
						vecino = segmentation->GetPixel(roundedIndex);
					}
					if(vecino!=0)
					{
						if(histograms[i].count(vecino)>0)
						{
							histograms[i][vecino]++;
						}
						else
						{
							histograms[i][vecino]=1;
						}
					}
					else
					{
						histograms[i][label]++;
					}

				}
			}
		}
		for(int i=0;i<indeces.size();i++)
		{
			csv_file << fiberFile << ",D:"<<indeces[i][0]<< ":"<<indeces[i][1]<<":"<<indeces[i][2]<<",";
			for(std::map<int,int>::iterator it=histograms[i].begin(); it!= histograms[i].end(); ++it)
			{
				csv_file<< it->first<<":"<<it->second<<",";
			}	
			csv_file<<std::endl;
		}
		fiberFile= cl.next("");
		segFile = cl2.next("");
		numFiles--;
		std::cout <<  fiberFile << " " << numFiles << std::endl;
	}
	csv_file.close();
	return 0;
}

