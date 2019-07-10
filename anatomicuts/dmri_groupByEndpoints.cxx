/* Andrew Zhang
 * Professor Siless
 * dmri_groupingByEndpoints.cxx
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>

#include "itkMesh.h"
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "itkPolylineCell.h"
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <cmath>
#include "itkPolylineCell.h"
#include "GetPot.h"
#include "TrkVTKPolyDataFilter.txx"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "ClusterTools.h"
#include "itkDefaultStaticMeshTraits.h"

using namespace std;

int main(int narg, char* arg[]) 
{
	//Take in inputs
	GetPot c1(narg, const_cast<char**>(arg));

	//Usage error
	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s streamlines -i imageFile -d outputDirectory" << endl;
		return -1; 
	}

	
	//Take in information
//	const char *stream_lines = c1.follow("string_file.trk", "-s"); 
	const char *image_file = c1.follow("image_file.nii.gz", "-i"); 
	const char *output = c1.follow("output_directory", "-d"); 

	//Read in files
	vector<string> inputFiles; 
	for (string inputName = string(c1.follow("", 2, "-s", "-S")); access(inputName.c_str(), 0) == 0; inputName = string(c1.next("")))
	{
		inputFiles.push_back(inputName); 
	}

	//Variable definitions
	enum {Dimension =3};
	typedef float                                                        PixelType;
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
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;

	ImageType::Pointer refImage;
	ImageType::Pointer mask; 

	vector<ColorMeshType::Pointer>* meshes; 
	vector<ColorMeshType::Pointer>* fixMeshes; 
	vector<vtkSmartPointer<vtkPolyData>> polydatas; 

	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType; 
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New(); 

	vector<HistogramMeshType::Pointer>* histoMeshes; 
	clusterTools->GetPolyDatas(inputFiles, &polydatas, mask); 

	meshes = clusterTools->PolydataToMesh(polydatas); 

	for (int i = 0; i < meshes->size(); i++) 
	{
		ColorMeshType::Pointer input = (*meshes)[i]; 

		int pointIndices = 0; 
		int cellIndices = 0; 
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin(); 
		for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
		{
			PointType firstPt, start, end; 
			firstPt.Fill(0); 

			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin(); 
			input->GetPoint(*it, &firstPt); 
			double lengthSoFar = 0; 

			start = firstPt; 

			cout << "First point: " << firstPt << endl; 

			for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt; 
				pt.Fill(0);
				input->GetPoint(*it, &pt);
				lengthSoFar += firstPt.EuclideanDistanceTo(pt); 
				input->GetPoint(*it, &firstPt); 
				end = pt; 	
			} 
			
			//Value of coordinates based on image
			using PixelType = short;
			constexpr unsigned int Dimension = 3; 

			using ImageType = Image<PixelType, Dimension>; 
			using ReaderType = ImageFileReader<ImageType>; 

			ReaderType::Pointer reader = ReaderType::New(); 

			reader->SetFileName(arg[4]); 
			reader->Update();
			
			ImageType::Pointer inputImage = reader->GetOutput(); 

			ImageType::IndexType index1, index2; 

			index1[0] = start[0]; 
			index1[1] = start[1]; 
			index1[2] = start[2]; 

			//index2[0] = end[0]; 
			//index2[1] = end[1]; 
			//index2[2] = end[2]; 

			cout << "First point: " << index1 << endl; 
			//cout << "Last point: " << index2 << endl; 
			
		//	const PixelType value1 = inputImage->GetPixel(index1); 
			//const PixelType value2 = inputImage->GetPixel(index2); 

			//cout << index1 << " = " << value1 << endl; 
			//cout << index2 << " = " << value2 << endl; 
		}
	}

	//TO DELETE
//	cout << stream_lines << endl << image_file << endl << output << endl; 

	return 0; 	
}



/*
 * Ushaped fibers go into a TRK file
 * Check endpoints and put into same region TRK file
 * If they connect two different structures, ignore
 *
 * First open/read files
 * Then identify endpoints
 * Then check if in same region
 * if so, then output them to TRK file and into a folder
 *
 */

/*	//Some variabel names are repeated but have _img added
	using PixelType_img = unsigned char; 
	constexpr unsigned int Dimension_img = 2; 

	using ImageType = Image<PixelType_img, Dimension_img>;
	using PointSetType = PointSet<PixelType_img, Dimension_img>; 
	using ReaderType = ImageFileReader<ImageType>; 

	ReaderType::Pointer reader = ReaderType::New(); 

	const char* inputFilename = arg[1]; //Test if image is being read 
	reader->SetFileName(inputFilename); 

	try
	{
		reader->Update(); 
	}
	catch(ExceptionObject &err)
	{
		cout << "ExceptionObject caught !" << endl; 
		cout << err << endl; 
		return EXIT_FAILURE; 
	}
*/
	//Extract pixels from image
/*	using ImageType = Image<unsigned short, 3>; 

	ImageType::Pointer image = ImageType::New(); 

	const ImageType::SizeType size = {{200, 200, 200}}; 
	const ImageType::IndexType start = {{0, 0, 0}}; 

	ImageType::RegionType region; 
	region.SetSize(size); 
	region.SetIndex(start); 

	image->SetRegions(region); 
	image->Allocate(true); 

	const ImageType::IndexType pixelIndex = {{144, 168, 75}}; 

	ImageType::PixelType pixelValue = image->GetPixel(pixelIndex); 

	image->SetPixel(pixelIndex, pixelValue); 
*/

