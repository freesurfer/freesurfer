/* Alexander Zsikla
 * Professor Siless
 * August 2019
 * dmri_extractSurfaceMeasurements.cxx
 *
 * Takes in a cluster and outputs metrics, specifically thickness and curvature, based on the endpoints of the connections and outputs it to an external CSV file
 *
 */

//Libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

// Input Splicing
#include "GetPot.h"

// TRK Reading
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
#include "LabelPerPointVariableLengthVector.h"
#include "EuclideanMembershipFunction.h"
#include "ClusterTools.h"
#include "itkDefaultStaticMeshTraits.h"

// Surface Reading
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include <set>
#include "colortab.h"
#include "fsenv.h"
#include "mrisurf.h"

using namespace std;

int main(int narg, char* arg[])
{
	GetPot num1(narg, const_cast<char**>(arg));
	
	// Checking for correct parameters
	if ((num1.size() <= 6) or (num1.search(2, "--help", "-h")))
	{
		cout << "Usage: " << endl
		     << arg[0] << " -i streamlineFile.trk"
		     << " -s surfaceFile.orig -o output.csv"
		     << endl;

		return EXIT_FAILURE;
	}

	//
	// Input Parsing
	//

	// Looping for multiple file input of a certain kind of file
	vector<string> inputFiles;
	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);
	

	const char *surface = num1.follow("lh.orig", "-s"); 	// this is a surface
	const char *output  = num1.follow("output.csv", "-o");


	// TO DELETE
	// Testing that files are saved and can be outputted
	
	for (int i = 0; i < inputFiles.size(); i++)
		cout << inputFiles.at(i) << endl;

	cout << surface << endl << output << endl;

	//
	// Reading in TRK File
	//

	enum {Dimension =3};
	typedef int PixelType;
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
	typedef ColorMeshType::CellType CellType;
	typedef itk::PolylineCell<CellType> PolylineCellType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;

	ImageType::Pointer refImage;
	ImageType::Pointer mask;	

	vector<ColorMeshType::Pointer>* meshes;
	vector<ColorMeshType::Pointer>* fixMeshes;
	vector<vtkSmartPointer<vtkPolyData>> polydatas;

	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();

	std::vector<HistogramMeshType::Pointer>* histoMeshes;
	clusterTools->GetPolyDatas(inputFiles, &polydatas, mask);

	meshes = clusterTools->PolydataToMesh(polydatas);

	for(int i = 0; i < meshes->size(); i++)
	{ 
		ColorMeshType::Pointer input = (*meshes)[i];

		int averageId = 0;
		float stdCluster = 0;

		set<int> unfilteredIds;

		int pointIndices = 0;
		int cellIndices  = 0;
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
		for (int cellId = 0; inputCellIt != input->GetCells()->End(); ++inputCellIt, cellId++)
		{

			PointType firstPt;
			firstPt.Fill(0);
			float val1 = 0, val2 = 0;

			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
			input->GetPoint(*it,&firstPt);
			double lenghtSoFar = 0;
			// iterates through the points in one stream line and computes the distance between each pairs of points
			for(; it != inputCellIt.Value()->PointIdsEnd(); it++)
			{
				PointType pt;
				pt.Fill(0);
				input->GetPoint(*it,&pt);
				lenghtSoFar += firstPt.EuclideanDistanceTo(pt);
				input->GetPoint(*it,&firstPt);
				cout << "Length of Stream: " << lenghtSoFar << endl;
			}


		}
	}

	//
	// Reading in surface File
	//

	MRI_SURFACE *surf;
        surf = MRISread(surface);
	
	// TODO with VIV

	//
	// Outputing to an extneral file
	//

	ofstream oFile;
	oFile.open(output);

	if (not oFile.is_open()) {
		cerr << "Could not open output file" << endl;
		return -1;
	}

	oFile.close();

	return EXIT_SUCCESS;
}

/*
 * recieve a list of streamlines (TRK)
 * connecting same structure
 *
 * checking endpoints for cortical thickness (there is a file that contains the info, but it is surface based)
 *
 * table
 * streamline name, thickness of first point, thickness of second point, curvature of first point, curvature of second point (in same file as the thickness)
 *
 */


	/*for (unsigned i = 0; i < surf->nvertices; i++)
	{
		double x, y, z;
		float pdx, pdy, pdz;

		if (surf->vertices[j].x > 1)
		{
			MRISsurfaceRASTToVoxel(surf, images, surf->vertices[j].x, surf->vertices[j].y, surf->vertices[j].z, &x, &y, &z);
			float magnitud = MRIvoxelGradient(images, (float) x, (float) y, (float) z, &pdx, &pdy, &pdz);
	}*/
