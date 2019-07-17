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

	typedef ImageFileReader<ImageType> ImageReaderType; 
	ImageReaderType::Pointer reader = ImageReaderType::New(); 
	reader->SetFileName(c1.next("")); 
	reader->Update(); 	
	inputImage = reader->GetOutput(); 

	//Take in input trk file
	meshes = clusterTools->PolydataToMesh(polydatas); 

	ColorMeshType::Pointer input = (*meshes)[0];
	ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin(); 

	set<int> filteredIds; 
	vector<int> regions;
        vector<int> repeat_check;
	bool not_repeat = true; 	

//	vector<int> pointIndices; 
//	vector<int> cellIndices; 
//	int previous_point_count = 0; 
//	int previous_cell_count = 0; 

	int pointIndices = 0; 
	int cellIndices = 0; 	
	int val1, val2; 

	//Map of the region values and their corresponding meshes
	map<int, ColorMeshType::Pointer> sorted_meshes; 

	//Variables for testing
	int stream_count = 0; 
	ImageType::IndexType index1, index2; 

	//Cycles through each streamline
	for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
	{
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

			if (sorted_meshes.count(val1) == 0)
			{
				ColorMeshType::Pointer om = ColorMeshType::New(); 
				om->SetCellsAllocationMethod(ColorMeshType::CellsAllocatedDynamicallyCellByCell); 

				sorted_meshes.insert(pair<int, ColorMeshType::Pointer> (val1, om)); 
			} 
			//map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.at(val1); 
			ColorMeshType::Pointer target_mesh = sorted_meshes.at(val1); 

			CellAutoPointer line;
			line.TakeOwnership (new PolylineCellType);
			int k = 0;
			it = inputCellIt.Value()->PointIdsBegin();
			for( ; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				input->GetPoint (*it, &pt);

				target_mesh->SetPoint (pointIndices, pt);
				//Sets the point indicated by pointIndices to index k
				line->SetPointId (k, pointIndices);

				k++;
				pointIndices++;
			}

			//Cell is inserted into the mesh
			//Assign lines of a specific region into that corresponding file?
			target_mesh->SetCell(cellIndices, line);
			ColorMeshType::CellPixelType cellData;
			input->GetCellData(cellId, &cellData);
			target_mesh->SetCellData(cellIndices, cellData) ;
			cellIndices++;
		}
	}
			
			
//		        filteredIds.insert(cellId); 	
//		}
		
//		regions.push_back(val1); 
		
//	}

	//The x and y coordinates are turning into their opposites???


	for (map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.begin(); iter != sorted_meshes.end(); iter++)
	{
		string outputName; 
		string number = to_string(iter->first); 
		string filename = number + ".trk"; 
		outputName = string(output) + "/" + filename; 

		cout << "Mesh name: " << outputName << endl; 

		//State for each unique trk file?
		clusterTools->SaveMesh(iter->second, inputImage, outputName, inputFiles[0]);  
	
		cerr << "trk file made" << endl; 
	}



	//Make a separate function with om as an input????
/*	inputCellIt = input->GetCells()->Begin(); 
	for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
	{
		//pointIndices.push_back(0); 
		//cellIndices.push_back(0); 

		//Checks for streamlines that start and end in a region already seen
		for (int i = 0; i < repeat_check.size(); i++)
		{
			if (regions[cellId] == repeat_check[i])
			{
//				cerr << "Reject" << endl; 
				not_repeat = false;
				break; 
			}
		}

		//Runs for every cell that applies, including repeats
		if (filteredIds.count(cellId) > 0)
		{
			cerr << "Output start" << endl; 

			//Should this be declared here? Or stay outside all the loops?
			//Output files
			
			ColorMeshType::Pointer om = ColorMeshType::New(); 
			om->SetCellsAllocationMethod(ColorMeshType::CellsAllocatedDynamicallyCellByCell); 

			CellAutoPointer line;
			line.TakeOwnership (new PolylineCellType);
			int k = 0;
			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
			for( ; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				input->GetPoint (*it, &pt);

				om->SetPoint (pointIndices, pt);
				//Sets the point indicated by pointIndices to index k
				line->SetPointId (k, pointIndices);

				k++;
				pointIndices++;
			}

			cerr << "Set points" << endl; 

			//Cell is inserted into the mesh
			//Assign lines of a specific region into that corresponding file?
			om->SetCell(cellIndices, line);
			ColorMeshType::CellPixelType cellData;
			input->GetCellData(cellId, &cellData);
			om->SetCellData(cellIndices, cellData) ;
			cellIndices++;

			string outputName; 
			string number = to_string(regions[cellId]); 
			string filename = number + ".trk"; 
			outputName = string(output) + "/" + filename; 

			cout << "Mesh name: " << outputName << endl; 

			//State for each unique trk file?
			clusterTools->SaveMesh(om, inputImage, outputName, inputFiles[0]);  

			repeat_check.push_back(regions[cellId]); 

			//Iterate thru the streamline again or take the .value() of the streamline
			//which can have getcell() to get streamline later

			//save id and label of the streamline to know what output it will go to

			cerr << "Repeat loop" << endl; 
		}	

		not_repeat = true; 
	}

	//Store streamlines that lie in the same structure, then place into one trk file

	//Find a way to change the order of the streamlines as well
	//Keep them unorganized, check previous values if they've been used?
*/

	cerr << "Total of " << stream_count << " streamlines" << endl; 

	delete meshes;

	return 0; 	
}


/* Questions:
 * -Need a new mesh type for each unique valued cell?
 * -Check if the output file already exists?
 * -For the output files, only all the streamlines of the first region are transferred over
 * -The x and y coordinates become inverted? 
 * -Sort so that the regions are ordered?
 * -.Value() returns a memory address, even when dereferenced
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
