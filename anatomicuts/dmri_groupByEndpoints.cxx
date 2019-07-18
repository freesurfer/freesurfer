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
#include <algorithm>
#include <map>

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
#include "LabelPerPointVariableLengthVector.h"
#include "EuclideanMembershipFunction.h"
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
	typedef float PixelType;
	const unsigned int PointDimension = 3;
	typedef std::vector<int> PointDataType;
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

	vector<ColorMeshType::Pointer>* meshes; 
	vector<vtkSmartPointer<vtkPolyData>> polydatas; 
	ImageType::Pointer inputImage;  

	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType; 
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New(); 

	clusterTools->GetPolyDatas(inputFiles, &polydatas, inputImage); 

/*	using ImageType = Image<PixelType, Dimension>; 
	using ReaderType = ImageFileReader<ImageType>; 

	ReaderType::Pointer reader = ReaderType::New(); 

	reader->SetFileName(image_file); 
	reader->Update();
*/
	//Variable to read in the image file
	typedef ImageFileReader<ImageType> ImageReaderType; 
	ImageReaderType::Pointer reader = ImageReaderType::New(); 
	reader->SetFileName(c1.next("")); 
	reader->Update(); 	
	inputImage = reader->GetOutput(); 

	//Take in input trk file
	meshes = clusterTools->PolydataToMesh(polydatas); 
	ColorMeshType::Pointer input = (*meshes)[0];
	ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin(); 

	//Variables to hold the start and end values of a streamline
	int val1, val2; 

	//Map of the region values and their corresponding meshes
	map<int, ColorMeshType::Pointer> sorted_meshes; 
	map<int, int> pointIndices; 
	map<int, int> cellIndices; 

	//Variables for testing
	int stream_count = 0; 
	ImageType::IndexType index1, index2; 

	//Cycles through each streamline
	for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
	{
		cerr << cellId << endl; 

		PointType start, end; 
		start.Fill(0); 
		val1 = 0, val2 = 0; 

		//Make a variable to iterate thru one stream at a time
		CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin(); 
		input->GetPoint(*it, &start); 

		//Goes through each point in a streamline
		for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
		{
			PointType pt; 
			pt.Fill(0);
			input->GetPoint(*it, &pt);
	
			ImageType::IndexType index; 
			int value = 0; 

			//Find the first and last nonzero values
			if (inputImage->TransformPhysicalPointToIndex(pt, index))
			{
				value = inputImage->GetPixel(index); 
				if (val1 == 0 and value != 0)
				{
					val1 = value; 
					index1 = index; 
				}
				if (value != 0)
				{
					val2 = value; 
					index2 = index;
				        end = pt; 	
				}
			}

		}  	
		
		cerr << "First point: " << start << endl; 
		cerr << "Last point: " << end << endl; 

		//If start and end values match, take in that cell Id
		if (val1 != 0 and val1 == val2)
		{
			cout << cellId << endl; 
			cout << "First point: " << start << endl; 
			cout << "End point: " << end << endl; 
			cout << "Start and ends match: " << endl; 
			cout << index1 << " = " << val1 << endl; 
			cout << index2 << " = " << val2 << endl;  
			stream_count++;

			//Obtain the mesh associated with the value
			if (sorted_meshes.count(val1) == 0)
			{
				ColorMeshType::Pointer om = ColorMeshType::New(); 
				om->SetCellsAllocationMethod(ColorMeshType::CellsAllocatedDynamicallyCellByCell); 

				sorted_meshes.insert(pair<int, ColorMeshType::Pointer> (val1, om)); 
			} 
			map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.find(val1); 
			ColorMeshType::Pointer target_mesh = iter->second; 

			if (pointIndices.count(val1) == 0)
			{
				pointIndices.insert(pair<int, int> (val1, 0)); 
			}

			//INSERT THE CELLID OF WHAT IS ADDED TO THE MESH?

			CellAutoPointer line;
			line.TakeOwnership (new PolylineCellType);
			int k = 0;
			it = inputCellIt.Value()->PointIdsBegin();

			//Copy over the points of the streamlines that are going to be outputted
			for( ; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				input->GetPoint (*it, &pt);

				//Shows correct original coordinates
				//if (it == inputCellIt.Value()->PointIdsBegin())
				//	cerr << pt << endl; 	

				target_mesh->SetPoint (pointIndices.at(val1), pt);

				//Sets the point indicated by pointIndices to index k
				line->SetPointId (k, pointIndices.at(val1));

				k++;
				pointIndices.at(val1)++;
			}

			//No inverted coordinates by this point?

			if (cellIndices.count(val1) == 0)
			{
				cellIndices.insert(pair<int, int> (val1, 0)); 
			}

			//Cell is inserted into the mesh
			target_mesh->SetCell(cellIndices.at(val1), line);
			ColorMeshType::CellPixelType cellData;
			input->GetCellData(cellId, &cellData);
			target_mesh->SetCellData(cellIndices.at(val1), cellData) ;
			cellIndices.at(val1)++;
		}
	}
			
	//The x and y coordinates are turning into their opposites???

	//Check that the starting points in the map are correct, which they are
	for (map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.begin(); iter != sorted_meshes.end(); iter++)
	{
		ColorMeshType::Pointer meshh = iter->second; 
		ColorMeshType::CellsContainer::Iterator test = meshh->GetCells()->Begin(); 
	
		for (; test != meshh->GetCells()->End(); test++)
		{
			cerr << iter->first << endl; 

			PointType new_start; 
			new_start.Fill(0);

			CellType::PointIdIterator it = test.Value()->PointIdsBegin(); 
			meshh->GetPoint(*it, &new_start); 
			cerr << "Output's first point: " << new_start << endl; 
		}
	}

	//Print out the output trk file for each mesh with a unique key
	for (map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.begin(); iter != sorted_meshes.end(); iter++)
	{
		string outputName; 
		
		//Naming the trk file based on the index value of the start and end points
		string number = to_string(iter->first); 
		
		string filename = number + ".trk"; 
		outputName = string(output) + "/" + filename; 

		cout << "Mesh name: " << outputName << endl; 

		clusterTools->SaveMesh(iter->second, inputImage, outputName, inputFiles[0]);  
	
		cerr << "trk file made" << endl; 
	}

	cerr << "Total of " << stream_count << " streamlines" << endl; 

	delete meshes;

	return 0; 	
}


/* Questions:
 * -The x and y coordinates become inverted? 
 * -Need multiple trk files as inputs?
 *
 */

//Check points before endpoint if endpoint gives value zero
//What if starts at 0?

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
