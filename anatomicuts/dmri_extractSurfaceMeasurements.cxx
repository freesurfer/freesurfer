/* Alexander Zsikla
 * Professor Siless
 * August 2019
 * dmri_extractSurfaceMeasurements.cxx
 *
 * Takes in a cluster and outputs metrics, specifically thickness and curvature, based on the endpoints of the connections and outputs it to an external CSV file
 *
 * export FREESURFER_HOME /home/fsuser2/alex_zsikla/install
 * source /home/fsuser2/alex_zsikla/install/SetUpFreeSurfer.sh
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
#include "fsSurface.h"
#include "itkTriangleCell.h"
#include <set>
#include "colortab.h"
#include "fsenv.h"
#include "itkVTKPolyDataWriter.h"
#include "itkSmoothingQuadEdgeMeshFilter.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
	
#include "vtkFillHolesFilter.h" 
#include "vtkPolyDataNormals.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkDecimatePro.h"
#include "vtkCleanPolyData.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTriangleFilter.h"

#include "vtkDelaunay3D.h"
#include "macros.h"
#include "mrisurf.h"
#include "mri.h"
#include "vtkKdTreePointLocator.h"

//#include "vtkPCACurvatureEstimation.h"
#include "vtkCurvatures.h"

/*#if VTK_MAJOR_VERSION > 5	
	#include "vtkPCACurvatureEstimation.h"
#else
	#include "vtkCurvatures.h"*/

using namespace std;

// HELPER FUNCTIONS
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf);

int main(int narg, char* arg[])
{
	GetPot num1(narg, const_cast<char**>(arg));
	
	// Checking for correct parameters
	if ((num1.size() <= 6) or (num1.search(2, "--help", "-h")))
	{
		cout << "Usage: " << endl
		     << arg[0] << " -i streamlineFile.trk -s surfaceFile.orig -t overlayFile.thickness -c overlayFile.curv -o output.csv"
		     << endl;

		return EXIT_FAILURE;
	}

	// INPUT PARSING
	// Looping for multiple file input of a certain kind of file
	vector<string> inputFiles;
	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);
	

	const char *surfaceFile = num1.follow("surfacy", "-s");
	const char *thickFile   = num1.follow("thicky", "-t");
	const char *curvFile    = num1.follow("curvy", "-c");
	const char *outputFile  = num1.follow("outputy", "-o");

	// TO DELETE
	// Testing that files are saved and can be outputted
	
	for (int i = 0; i < inputFiles.size(); i++)
		cout << inputFiles.at(i) << endl;

	cout << "Surface: " << surfaceFile << endl << "Thickness: " << thickFile << endl << "Curvature: " << curvFile << endl << "Output: " << outputFile << endl;
	
	// Opening External File for Outputting
	ofstream oFile;
	oFile.open(outputFile);

	if (not oFile.is_open()) {
		cerr << "Could not open output file" << endl;
		return -1;
	}
	
	// Declaration of Variables for Program to Function
	// TRK file Definitions
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
	
	// Surface Defintions
	typedef float CoordType;
	typedef fs::Surface< CoordType, Dimension> SurfType;

	//Reading in surface from file
	MRI_SURFACE *surf;
        surf = MRISread(surfaceFile);
	
	SurfType::Pointer surface =  SurfType::New();
	surface->Load(&*surf);
	
	surf = surface->GetFSSurface(&*surf);
	
	// Curvature Defintions
	vtkSmartPointer<vtkCurvatures> curvature = vtkSmartPointer<vtkCurvatures>::New();
	curvature->SetInput(FSToVTK(surf));
	curvature->SetCurvatureTypeToGaussian();
	
	curvature->Update();
	vtkSmartPointer<vtkPolyData> polydata =curvature->GetOutput();	

	// The first and last point of a stream
	PointType firstPt, lastPt;
	firstPt.Fill(0);
	lastPt.Fill(0);

	// Finding the Thickness
	vtkSmartPointer<vtkPolyData> surfVTK = FSToVTK(surf);	
	
	vtkSmartPointer<vtkKdTreePointLocator> surfTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	surfTree->SetDataSet(surfVTK);
	surfTree->BuildLocator();	

	// Holding the coordinates of the points
	double firstPt_array[3];
	double lastPt_array[3];

	MRISreadCurvature(surf, curvFile);

	for(int i = 0; i < meshes->size(); i++)
	{ 
		ColorMeshType::Pointer input = (*meshes)[i];

		int averageId = 0;
		float stdCluster = 0;

		set<int> unfilteredIds;

		int pointIndices = 0;
		int cellIndices  = 0;
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
		for (; inputCellIt != input->GetCells()->End(); ++inputCellIt)
		{
			// Creating streamline variable and finding first point
			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
			input->GetPoint(*it,&firstPt);

			// Finding the last point
			for (; it != inputCellIt.Value()->PointIdsEnd(); it++) 
				input->GetPoint(*it, &lastPt);
			
			for (int j = 0; j < 3; j++) {
				firstPt_array[j] = firstPt[j];
				lastPt_array[j]	 = lastPt[j];
			}

			double dist1 = 0, dist2 = 0;
			vtkIdType ID1 = surfTree->FindClosestPointWithinRadius(20.0, firstPt_array, dist1);
			vtkIdType ID2 = surfTree->FindClosestPointWithinRadius(20.0, lastPt_array, dist2);

			/*cout << "ID1: " << ID1 << " ID2: " << ID2 << endl;
			cout << "First Point: " << firstPt << " | [" << surf->vertices[ID1].x << ", " << surf->vertices[ID1].y << ", " << surf->vertices[ID1].z << "]" << endl;
			cout << "Last Point:  " << lastPt << " | [" << surf->vertices[ID2].x << ", " << surf->vertices[ID2].y << ", " << surf->vertices[ID2].z << "]" << endl;
			cout << "Curvature of First Point: " << surf->vertices[ID1].curv << " Curvature of Last Point:  " << surf->vertices[ID2].curv << endl;*/
		}
	}
		
	oFile.close();

	return EXIT_SUCCESS;
}

//
// Converts a surface to a VTK
//
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	cout << points->GetNumberOfPoints() << endl;

	for(int i = 0; i < surf->nvertices; i++)
		points->InsertNextPoint(surf->vertices[i].x, surf->vertices[i].y, surf->vertices[i].z);

	for( int i = 0; i < surf->nfaces; i++ ) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		for(int j = 0;j < 3; j++)
			triangle->GetPointIds()->SetId(j, surf->faces[i].v[j]);

		triangles->InsertNextCell(triangle);
	}	

	vtkSmartPointer<vtkPolyData> vtkSurface = vtkSmartPointer<vtkPolyData>::New();
	vtkSurface->SetPoints(points);
	vtkSurface->SetPolys(triangles);

	return vtkSurface;
}

/*
 * Things I Have Changed
 * 1. Commented out some directories that were not recognized
 * 2. Copied directories from freesurfer/resurf/Code to freesurfer/anatomicuts/Code to compile
 * 3.
 *
 */


/*
 * Things to Do:
 *
 * recieve a list of streamlines (TRK)
 * connecting same structure
 *
 * checking endpoints for cortical thickness (there is a file that contains the info, but it is surface based)
 *
 * table
 * streamline name, thickness of first point, thickness of second point, curvature of first point, curvature of second point (in same file as the thickness)
 *
 */

//Code for cycling through vertices
	/*for (unsigned i = 0; i < surf->nvertices; i++)
	{
		double x, y, z;
		float pdx, pdy, pdz;

		if (surf->vertices[j].x > 1)
		{
			MRISsurfaceRASTToVoxel(surf, images, surf->vertices[j].x, surf->vertices[j].y, surf->vertices[j].z, &x, &y, &z);
			float magnitud = MRIvoxelGradient(images, (float) x, (float) y, (float) z, &pdx, &pdy, &pdz);
	}*/

//Code that did not do exactly what it had to
			// Gave the same First and Last point
			/*CellType::PointIdIterator it2 = inputCellIt.Value()->PointIdsEnd();
			input->GetPoint(*it2, &lastPt);

			outFile << "Last Point: " << lastPt << endl;*/

