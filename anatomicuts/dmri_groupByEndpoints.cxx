/* Andrew Zhang
 * Professor Siless
 * dmri_groupingByEndpoints.cxx
 * July 2019
 *
 * A trk file with many streamlines is filtered based on if the starting and ending points of each 
 * streamline lie in the same structure. They are placed into separate trk files based on an inputted
 * image that acts as a reference for the structures the streamlines start and end at. 
 */

#include <iostream>
#include <string>
#include <map>

#include <itkImage.h>
#include <itkImageFileReader.h>

#include "itkMesh.h"
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "itkPolylineCell.h"
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include "GetPot.h"
#include "TrkVTKPolyDataFilter.txx"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "ClusterTools.h"

using namespace std;

int main(int narg, char* arg[]) 
{
	//Receive inputs
	GetPot c1(narg, const_cast<char**>(arg));

	//Usage error
	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s streamlineFile -i imageFile -d outputDirectory" << endl;
		return -1; 
	}

	
	//Take in information
	const char *image_file = c1.follow("image_file.nii.gz", "-i"); 
	const char *output = c1.follow("output_directory", "-d"); 

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

	//Create a color table
	typedef struct {
		unsigned char r;
		unsigned char g;
		unsigned char b;
	} color_triplet2;

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

	//Map of the region values to their corresponding meshes and other information
	map<int, ColorMeshType::Pointer> sorted_meshes; 
	//Holds the number of points for each mesh
	map<int, int> pointIndices; 
	//Holds the number of streamlines for each mesh
	map<int, int> cellIndices; 

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

			//Find the first and last nonzero values based on the transformation of the point
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
		
		//cerr << "First point: " << start << endl; 
		//cerr << "Last point: " << end << endl; 

		//If start and end values match, take in that cell Id
		if (val1 != 0 and val1 == val2)
		{
			//Print text for testing purposes
			//cout << cellId << endl; 
			//cout << "First point: " << start << " "; 
			//cout << "End point: " << end << endl; 
			//cout << "Start and ends match: " << " "; 
			//cout << index1 << " = " << val1 << ", "; 
			//cout << index2 << " = " << val2 << endl;  
			//stream_count++;

			//Obtain the mesh and related information associated with the value
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

			//Allocates and can release the cell's memory 
			CellAutoPointer line;
			line.TakeOwnership (new PolylineCellType);
			
			//Holds an index number for each unique point per streamline
			int k = 0;
			it = inputCellIt.Value()->PointIdsBegin();

			//Copy over the points of the streamlines that are going to be outputted
			for( ; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				input->GetPoint (*it, &pt);

				target_mesh->SetPoint (pointIndices.at(val1), pt);
				line->SetPointId (k, pointIndices.at(val1));

				k++;
				pointIndices.at(val1)++;
			}

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
/*
	//Check that the starting points in the map are correct, which they are
	for (map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.begin(); iter != sorted_meshes.end(); iter++)
	{
		ColorMeshType::Pointer meshh = iter->second; 
		ColorMeshType::CellsContainer::Iterator test = meshh->GetCells()->Begin(); 
	
		for (; test != meshh->GetCells()->End(); test++)
		{
	//		cerr << iter->first << " "; 

			PointType new_start; 
			new_start.Fill(0);

			CellType::PointIdIterator it = test.Value()->PointIdsBegin(); 
			meshh->GetPoint(*it, &new_start); 
	//		cerr << "Output's first point: " << new_start << endl; 
		}
	}
*/
	//Initiate the color table to give output meshes unique colors
	int i;
	int red, green, blue; 
	color_triplet2 table[256] = {{0, 0, 0}}; 

	const color_triplet2 black = {0, 0, 0}; 
	const color_triplet2 white = {255, 255, 255}; 

	table[0] = white; 
	table[255] = black; 

	i = 20; 
	for (red = 0; red <= 255; red += 51) {
		for (green = 0; green <= 255; green += 51) {
			for (blue = 0; blue <= 255; blue += 51) {
				table[i].r = red; 
				table[i].g = green; 
				table[i].b = blue; 
				++i; 
			}
		}
	}
	table[0] = white; 
	table[255] = black; 

	i = 0; 

	//Print out the output trk file for each mesh with a unique key	
	for (map<int, ColorMeshType::Pointer>::iterator iter = sorted_meshes.begin(); iter != sorted_meshes.end(); iter++)
	{
		string outputName; 
		
		//Naming the trk file based on the index value of the start and end points' corresponding region
		string number = to_string(iter->first); 
		
		string filename = number + ".trk"; 
		outputName = string(output) + "/" + filename; 

		//cout << "Mesh name: " << outputName << endl; 

		//Assign the unique color to each mesh		
		int index = ((int)47.*((i % sorted_meshes.size()) % (150))) % 197 + 5; 
		index = (int)(13 * (i % sorted_meshes.size())) % 150 + 65; 

		unsigned char color[3] = {table[index].r, table[index].g, table[index].b};
	
		//Create the output trk
		typedef PolylineMeshToVTKPolyDataFilter<ColorMeshType> VTKConverterType;
		typename VTKConverterType::Pointer vtkConverter = VTKConverterType::New(); 
		vtkConverter->SetInput(iter->second);
		vtkConverter->Update(); 

		SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader = TrkVTKPolyDataFilter<ImageType>::New(); 
		trkReader->SetInput(vtkConverter->GetOutputPolyData()); 
		trkReader->SetReferenceTrack(inputFiles[0]); 
		trkReader->SetColor(color);
		trkReader->VTKToTrk(outputName); 

		i++;
	}

	//cerr << "Total of " << stream_count << " streamlines" << endl; 

	//Clear associated memory
	delete meshes;
	sorted_meshes.clear(); 
	pointIndices.clear(); 
	cellIndices.clear(); 

	return 0; 	
}


/* Questions:
 * -The x and y coordinates become inverted? 
 * -Need multiple trk files as inputs?
 */

/*
 * Ushaped fibers go into a TRK file
 * Check endpoints and put into same region TRK file
 * If they connect two different structures, ignore
 *
 * First open/read files
 * Then identify endpoints
 * Then check if in same region
 * if so, then output them to TRK file and into a folder
 */
