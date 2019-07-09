/* Andrew Zhang
 * Professor Siless
 * dmri_groupingByEndpoints.cxx
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkPointSet.h>
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
	//Define variables
	GetPot c1(narg, const_cast<char**>(arg));
	GetPot c2(narg, const_cast<char**>(arg)); 

	//Usage error
	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s streamlines -i imageFile -d outputDirectory" << endl;
		return -1; 
	}

	//More variable definitions
	enum {Dimension = 3};
	typedef int                                                        PixelType;
	const unsigned int PointDimension = 3;
	typedef vector<int>                  PointDataType;
	
	typedef Mesh<PixelType, PointDimension> ColorMeshType; //Error with PixelType with no _
	typedef ColorMeshType::PointType PointType;
	typedef ColorMeshType::CellType        CellType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;

	//Image variables
	//Some variabel names are repeated but have _img added
	using PixelType_img = unsigned char; 
	constexpr unsigned int Dimension_img = 2; 

	using ImageType = Image<PixelType_img, Dimension_img>;
	using PointSetType = PointSet<PixelType_img, Dimension_img>; 
	using ReaderType = ImageFileReader<ImageType>; 

	ReaderType::Pointer reader = ReaderType::New(); 

	const char* inputFilename = arg[4]; 
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

	PointSetType::Pointer pointSet = PointSetType::New(); 

	using IteratorType = ImageRegionConstIterator<ImageType>; 

	const ImageType* image = reader->GetOutput(); 

	IteratorType it(image, image->GetBufferedRegion());

	it.GoToBegin(); 

	using PointType_img = PointSetType::PointType; 
	PointType_img point; 

	unsigned long pointId = 0; 

	//Take in information
	const char *stream_lines = c1.follow("string_file.trk", "-s"); 
	const char *image_file = c1.follow("image_file.nii.gz", "-i"); 
	const char *output = c1.follow("output_directory", "-o"); 

	//Image processing function
	while (!it.IsAtEnd()) 
	{
		//Convert the pixel position into a Point
		image->TransformIndexToPhysicalPoint(it.GetIndex(), point); 
		pointSet->SetPoint(pointId, point); 

		//Transfer the pixel data to the value associated with the point
		pointSet->SetPointData(pointId, it.Get()); 

		++it; 
		++pointId; 
	}

	cout << "Number of Points = " << pointSet->GetNumberOfPoints() << endl; 

	//Declaration of vectors
	vector<ColorMeshType::Pointer>* meshes; 
	vector<ColorMeshType::Pointer>* fixMeshes; 
	vector<vtkSmartPointer<vtkPolyData>> polydatas; 

	//TO DELETE
	//cout << stream_lines << endl << image_file << endl << output << endl; 

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
