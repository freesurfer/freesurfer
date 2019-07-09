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
/*	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s streamlines -i imageFile -d outputDirectory" << endl;
		return -1; 
	}
*/
	//More variable definitions
/*	enum {Dimension = 3};
	typedef int                                                        PixelType;
	const unsigned int PointDimension = 3;
	typedef vector<int>                  PointDataType;
	
	typedef Mesh<PixelType, PointDimension> ColorMeshType; //Error with PixelType with no _
	typedef ColorMeshType::PointType PointType;
	typedef ColorMeshType::CellType        CellType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;
*/
	//Image variables
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

	using PixelType = short;
	constexpr unsigned int Dimension = 3; 

	using ImageType = Image<PixelType, Dimension>; 
	using ReaderType = ImageFileReader<ImageType>; 

	ReaderType::Pointer reader = ReaderType::New(); 

	reader->SetFileName(arg[1]); 
	reader->Update();

	ImageType::Pointer inputImage = reader->GetOutput(); 

	ImageType::IndexType index; 

	index[0] = 144; 
	index[1] = 168; 
	index[2] = 75; 

	const PixelType value = inputImage->GetPixel(index); 

	cout << index << " = " << value << endl; 
	
/*
	using WriterType = ImageFileWriter<ImageType>; 
	WriterType::Pointer writer = WriterType::New(); 
	writer->SetInput(image); 
	writer->SetFileName(arg[1]); 

	try
	{
		writer->Update(); 
	}
	catch (ExceptionObject & error)
	{
		cerr << "Error: " << error << endl; 
		return EXIT_FAILURE; 
	}

	cout << pixelValue << endl; 
*/
	//Take in information
	const char *stream_lines = c1.follow("string_file.trk", "-s"); 
	const char *image_file = c1.follow("image_file.nii.gz", "-i"); 
	const char *output = c1.follow("output_directory", "-o"); 

	//Declaration of vectors
//	vector<ColorMeshType::Pointer>* meshes; 
//	vector<ColorMeshType::Pointer>* fixMeshes; 
//	vector<vtkSmartPointer<vtkPolyData>> polydatas; 

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
