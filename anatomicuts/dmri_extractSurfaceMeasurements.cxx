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
//#include "fsSurface.h"
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

#include "vtkPCACurvatureEstimation.h"
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

	//
	// Input Parsing
	//

	// Looping for multiple file input of a certain kind of file
	vector<string> inputFiles;
	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);
	

	const char *surface = num1.follow("lh.orig", "-s");
	const char *thick   = num1.follow("lh.thickness", "-t");
	const char *curv    = num1.follow("lh.curv", "-c");
	const char *output  = num1.follow("output.csv", "-o");


	// TO DELETE
	// Testing that files are saved and can be outputted
	
	for (int i = 0; i < inputFiles.size(); i++)
		cout << inputFiles.at(i) << endl;

	cout << surface << endl << thick << endl << curv << endl << output << endl;

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

	PointType first;
	PointType last;

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

			PointType firstPt, lastPt;
			firstPt.Fill(0);
			lastPt.Fill(0);

			// Creating streamline variable and finding first point
			CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
			input->GetPoint(*it,&firstPt);
			cout << "First Point: " << firstPt << endl;

			// Finding the last point
			for (; it != inputCellIt.Value()->PointIdsEnd();it++)
				input->GetPoint(*it, &lastPt);
			
			input->GetPoint(*it, &lastPt);	
			cout << "Last Point: " << lastPt << endl;
			
			// Gave the same First and Last point
			/*CellType::PointIdIterator it2 = inputCellIt.Value()->PointIdsEnd();
			input->GetPoint(*it2, &lastPt);

			outFile << "Last Point: " << lastPt << endl;*/
		}
	}

	//
	// Reading in Surface and find data values
	//

	//Reading in surface from file
	MRI_SURFACE *surf;
        surf = MRISread(surface);

	SurfType::Pointer surface =  SurfType::New();
	surface->Load(&*surf);

	surf = surface->GetFSSurface(&*surf);

	// Finding the Thickness

	vtkSmartPointer<vtkPolyData> surfVTK = FSToVTK(surf);	
	
	vtkSmartPointer<vtkKdTreePointLocator> surfTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	surfTree->SetDataSet(surfVTK);
	surfTree->BuildLocator();

	vtkPoints* points = vtkPoints::New();
	for (int i = 0; i < surfVTK->GetNumberOfPoints(); ++i)
	{

		double* point = surfVTK->GetPoint(i);
		vtkIdType iD = surfTree->FindClosestPoint(point);
		double* point2 = surfVTK->GetPoint( iD);
		float distance =  vtkMath::Distance2BetweenPoints(point,point2);
		cout << point [0] << " " << point2[0] << " " << distance << endl;
		
		/*if( distance > 0.01)
		{
			points->InsertPoint(i,point[0], point[1], point[2]);
		}*/
	}
	
	// Finding the Curvature
	//TODO



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
