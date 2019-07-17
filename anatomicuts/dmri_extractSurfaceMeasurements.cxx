/* Alexander Zsikla
 * Professor Siless
 * August 2019
 * dmri_extractSurfaceMeasurements.cxx
 *
 * Takes in a cluster and outputs metrics, specifically thickness and curvature, based on the endpoints of the connections and outputs it to an external CSV file
 *
 * export FREESURFER_HOME /home/fsuser2/alex_zsikla/install
 * source $FREESURFER_HOME/SetUpFreeSurfer.sh
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
#include "vtkCurvatures.h"

using namespace std;

// HELPER FUNCTIONS
string cleanFile(string file);
string makeCSV(string dir, string number);
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf);

int main(int narg, char* arg[])
{
	GetPot num1(narg, const_cast<char**>(arg));
	
	// Checking for correct parameters
	if ((num1.size() <= 6) or (num1.search(2, "--help", "-h")))
	{
		cerr << "Usage: " << endl
		     << arg[0] << " -i streamlineFile.trk -s surfaceFile.orig -t overlayFile.thickness -c overlayFile.curv -o outputDirectory"
		     << endl << "OPTION: -fa FA_file.nii.gz" << endl;

		return EXIT_FAILURE;
	}

	// INPUT PARSING
	// Looping for multiple file input of a certain kind of file
	vector<string> inputFiles;
	
	int counter = 1;
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);

	const char *surfaceFile = num1.follow("rh.orig", "-s");
	const char *thickFile   = num1.follow("rh.thickness", "-t");
	const char *curvFile    = num1.follow("rh.curv", "-c");
	const char *outputDir   = num1.follow("/home/fsuser2/alex_zsikla/freesurfer/anatomicuts/output_TRK", "-o");

	if (num1.search("-fa"))
		const char *FA_file = num1.follow("fa.nii.gz", "-fa");

	// TO DELETE
	// Testing that files are saved and can be outputted
	
	cerr << endl;
	for (int i = 0; i < inputFiles.size(); i++)
		cerr << "TRK File " << i + 1 << ": " << inputFiles.at(i) << endl;

	cerr << "Surface:    " << surfaceFile << endl << "Thickness:  " << thickFile << endl << "Curvature:  " << curvFile << endl << "Output:     " << outputDir << endl;


	// Opening External File for Outputting
	string outputFilename = cleanFile(inputFiles.at(0));
	
	ofstream oFile;
	string aux = makeCSV(outputDir, outputFilename);
	oFile.open(makeCSV(outputDir, outputFilename));

	cerr << aux << endl;

	if (not oFile.is_open()) {
		cerr << "Could not open output file" << endl;
		return -1;
	}

	// Adds option for finding FA values	
	if (num1.search("-fa"))
		oFile << "Streamline Name , Curv of Start Point , Curv of Last Point , Thickness of Start Point , Thickness of Last Point , meanFA , avgFA , stdFA" << endl;
	else 
		oFile << "Streamline Name , Curvature of Start Point , Curvature of Last Point , Thickness of Start Point , Thickness of Last Point" << endl;

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
	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
	
	// Surface Defintions
	typedef float CoordType;
	typedef fs::Surface< CoordType, Dimension> SurfType;

	// Variable Declaration
	ImageType::Pointer mask;	

	vector<ColorMeshType::Pointer>* meshes;
	vector<vtkSmartPointer<vtkPolyData>> polydatas;
	
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();
	clusterTools->GetPolyDatas(inputFiles, &polydatas, mask);
	meshes = clusterTools->PolydataToMesh(polydatas);
	
	//Reading in surface from file
	//For Curvature
	MRI_SURFACE *surf;
        surf = MRISread(surfaceFile);
	
	SurfType::Pointer surface =  SurfType::New();
	surface->Load(&*surf);
	
	surf = surface->GetFSSurface(&*surf);

	MRISreadCurvature(surf, curvFile);	

	//For Thickness
	MRI_SURFACE *surf_t;
	surf_t = MRISread(surfaceFile);

	SurfType::Pointer surface_t = SurfType::New();
	surface_t->Load(&*surf_t);

	surf_t = surface_t->GetFSSurface(&*surf_t);

	MRISreadCurvature(surf_t, thickFile);

	// Initializing the KdTree
	vtkSmartPointer<vtkPolyData> surfVTK = FSToVTK(surf);
	
	vtkSmartPointer<vtkKdTreePointLocator> surfTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	surfTree->SetDataSet(surfVTK);
	surfTree->BuildLocator();	

	// Working with FA file
	


	PointType firstPt, lastPt;	// The first and last point of a stream
	firstPt.Fill(0);
	lastPt.Fill(0);

	double firstPt_array[3];	// Holding the coordinates of the points
	double lastPt_array[3];

	counter = 1;
	for(int i = 0; i < meshes->size(); i++)
	{ 
		ColorMeshType::Pointer input = (*meshes)[i];
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
		
		for (; inputCellIt != input->GetCells()->End(); ++inputCellIt, ++counter)
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

			vtkIdType ID1 = surfTree->FindClosestPoint(firstPt_array);
			vtkIdType ID2 = surfTree->FindClosestPoint(lastPt_array);

			//oFile << "StreamLine " << counter << "," << surf->vertices[ID1].curv << "," << surf->vertices[ID2].curv << ","
			//      << surf_t->vertices[ID1].curv << "," << surf_t->vertices[ID2].curv << endl;
		}
	}
		
	//oFile.close();

	return EXIT_SUCCESS;
}

string cleanFile(string file)
{
	int front = file.find_last_of("/");
       	int back  = file.find_last_of(".");
	
	return file.substr(front + 1, back - front - 1);
}

string makeCSV(string dir, string number)
{
	dir.append("/");
	dir.append(number);
	dir.append(".csv");

	return dir;
}

//
// Converts a surface to a VTK
//
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

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
 */

/*
 * Questions about the FA File Reading Code
 * 1. Line 27: This definition has a conflict with a previous definition that is described above (MeshType vs ColorMeshType) --> I think they are just names so I don't think it would matter
 * 2. Line 32: This definition has a conflict with something that is already made inside of another .h file, should i just rename and continue?
 * 3. Line 74: Is this the file that I pass in (AKA fa.nii.gz)
 * 4. Line 110 and 115: Why do you have two get lines that save to the same value?
 * 5. Line 140: What is the purpose of Avg Points?
 * 6. For my table, all I am doing is adding another column for FA values right? The program that I run on my computer automatically computes mean and std
 */

/*
 * Things to Do
 * 1. Have a different CSV file for each TRK file
 * 2. Find meanFA, avgFA, and stdFA for every streamline and output that to the CSV file
 * 3. Find the mean, avg, and std for each of those columns (probably using python) and output it to the bottom
 * 4.
 */

