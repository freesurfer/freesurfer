/* Andrew Zhang
 * Professor Siless
 * dmri_coloredFA.cxx
 * July 2019
 *
 * Find the FA values of all the points in a streamline and assign colors to them, outputting versions of the inputted files with colored streamlines. 
 *
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

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

using namespace std; 

int main(int narg, char* arg[])
{
        //Receive inputs
        GetPot c1(narg, const_cast<char**>(arg));

        //Usage error
	//Want to have a directory for streamline inputs?
        if (c1.size() == 1 || c1.search(2, "--help", "-h"))
        {
                cout << "Usage: " << endl;
                cout << arg[0] << " -s streamlines -i imageFile -d outputDirectory" << endl;
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
	//vector<ColorMeshType::Pointer>* colored_meshes; 
        vector<vtkSmartPointer<vtkPolyData>> polydatas;
        //vector<vtkSmartPointer<vtkPolyData>> colored_polydatas;
	ImageType::Pointer inputImage;

	//Variable to read in the image file
        typedef ImageFileReader<ImageType> ImageReaderType;
        ImageReaderType::Pointer reader = ImageReaderType::New();
        reader->SetFileName(c1.next(""));
        reader->Update();
        inputImage = reader->GetOutput();

        typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
        ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();

        clusterTools->GetPolyDatas(inputFiles, &polydatas, inputImage);

        //Variable to take in input trk file
        meshes = clusterTools->PolydataToMesh(polydatas);

	//Values and color association
	vector<float> FA_value; 
	//vtkSmartPointer<vtkPolyData> pointsPolyData = vtkSmartPointer<vtkPolyData>::New(); 
	//vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New(); 

			//Test colors
			unsigned char red[3] = {255, 0, 0}; 
			unsigned char green[3] = {0, 255, 0}; 
			unsigned char blue[3] = {0, 0, 255}; 

			colors->SetNumberOfComponents(3); 
			colors->SetName("Colors"); 
	

	for (int i = 0; i < meshes->size(); i++)
	{
		ColorMeshType::Pointer input = (*meshes)[i];
      		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New(); 

		for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
	  	{
			PointType start, end;
                	start.Fill(0);

			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
               		input->GetPoint(*it, &start);

              		//Goes through each point in a streamline
               		for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
                	{
				PointType pt;
	                        pt.Fill(0);
        	                input->GetPoint(*it, &pt);

                	        ImageType::IndexType index;
                	        //int FAvalue = 0;

                	       	//Find the first and last nonzero values based on the transformation of the point
                        	if (inputImage->TransformPhysicalPointToIndex(pt, index))
				{
					FA_value.push_back(inputImage->GetPixel(index));
					points->InsertNextPoint(pt[0], pt[1], pt[2]); 
					
					//Test color
					colors->InsertNextTupleValue(red); 
				}	
			}
		}

		//pointsPolyData->SetPoints(points); 
		polydatas[i]->SetPoints(points); 

		polydatas[i]->GetPointData()->SetScalars(colors); 

		//delete points;  

		//colored_polydatas.push_back(polydata); 

		//<vtkFloatArray>?
	}

	//FIX OUTPUT SO THE MODIFIED MESHES ARE PRINTED
	//colored_meshes = clusterTools->PolydataToMesh(colored_polydatas);
	for (int i = 0; i < polydatas.size(); i++)
	{
		//Create an output file for each mesh
		string outputName = string(output) + "/"; 

		//string fileName; 
		int place; 
		for (int j = 0; j < inputFiles[i].length(); j++)
		{
			if (inputFiles[i][j] == '/')
			{
				place = j; 
			}
		}

		for (int j = place + 1; j < inputFiles[i].length(); j++)
		{
			outputName += inputFiles[i][j]; 
		}

		cerr << outputName << endl; 

		//typedef PolylineMeshToVTKPolyDataFilter<ColorMeshType> VTKConverterType;
                //typename VTKConverterType::Pointer vtkConverter =  VTKConverterType::New();
		//vtkConverter->SetInput(input);
		//vtkConverter->Update();

		SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader = TrkVTKPolyDataFilter<ImageType>::New();
		trkReader->SetInput(polydatas[i]); //(vtkConverter->GetOutputPolyData());
                trkReader->SetReferenceTrack(inputFiles[i]);
		trkReader->SetReferenceImage(inputImage); 
                trkReader->VTKToTrk(outputName);	
	}

	return 0; 
}
