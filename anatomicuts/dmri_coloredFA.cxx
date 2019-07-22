/* Andrew Zhang
 * Professor Siless
 * dmri_coloredFA.cxx
 * July 2019
 *
 * Find the FA values of all the points in a streamline and assign colors to them
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
        vector<vtkSmartPointer<vtkPolyData>> polydatas;
        ImageType::Pointer inputImage;

        typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
        ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();

        clusterTools->GetPolyDatas(inputFiles, &polydatas, inputImage);

	//Variable to read in the image file
        typedef ImageFileReader<ImageType> ImageReaderType;
        ImageReaderType::Pointer reader = ImageReaderType::New();
        reader->SetFileName(c1.next(""));
        reader->Update();
        inputImage = reader->GetOutput();

        //Variable to take in input trk file
        meshes = clusterTools->PolydataToMesh(polydatas);

	for (int i = 0; i < meshes.size(); i++)
	{
		ColorMeshType::Pointer input = (*meshes)[i];
      		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();


		for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
	  	{
			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
               		input->GetPoint(*it, &start);

              		//Goes through each point in a streamline
               		for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
                	{
				PointType pt;
	                        pt.Fill(0);
        	                input->GetPoint(*it, &pt);

                	        ImageType::IndexType index;
                	        int FAvalue = 0;

                	       	//Find the first and last nonzero values based on the transformation of the point
                        	if (inputImage->TransformPhysicalPointToIndex(pt, index))
				{
					FAvalue = inputImage->GetPixel(index);
				}
			}
		}
	}




	return 0; 
}
