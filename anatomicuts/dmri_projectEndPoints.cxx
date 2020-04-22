/* Author: Alexander Zsikla
 * Advisor: Professor Siless
 * Time: Summer 2019
 * Name: dmri_projectEndPoints.cxx
 *
 * Description:
 * Changes the endpoints of a streamline to a different value and save it to a new overlay file
 *
 */

const int ENDPOINT_VALUE = 1;

// Libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

// Input Splicing
#include "GetPot.h"

// TRK Loading
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

// Surface Loading
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

// Helper Functions
vtkIdType which_ID(double n1, double n2, vtkIdType ID1, vtkIdType ID2);
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf);

int main(int narg, char* arg[])
{
	GetPot gp(narg, const_cast<char**>(arg));

	// Checking for correct parameters
	if ((gp.size() <= 8) or (gp.search(2, "--help", "-h")))
        {
                cerr << "Usage: " << endl
                     << arg[0] << " -i streamlineFile.trk -sl surfaceFile_lh.orig -sr surfaceFile_rh.orig" << endl
		     << "                       -ol left_overlayFile -or right_overlayFile" << endl 
		     << "		       -ri reference_Image" << endl;

                return EXIT_FAILURE;
        }

	// Declaration of Variables for Program to Function
	// TRK file Definition
	enum {Dimension = 3};
        typedef int PixelType;
        const unsigned int PointDimension = 3;
        typedef vector<int> PointDataType;
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

	// Surface file Definition
	typedef float CoordType;
        typedef fs::Surface< CoordType, Dimension> SurfType;	

	// Input Parsing
	vector<string> TRKFile;
	TRKFile.push_back(gp.follow("Could not find TRK file", "-i"));
	const char *surfaceFileL   = gp.follow("Could not find Surface File", "-sl");
	const char *surfaceFileR   = gp.follow("Could not find Surface File", "-sr");
	const char *overlayFileL   = gp.follow("Could not find Overlay File", "-ol");
	const char *overlayFileR   = gp.follow("Could not find Overlay File", "-or");
	const char *refImage       = gp.follow("Could not find Reference Image", "-ri");

	// Reading in the Image
	// ITK Version
	typedef itk::ImageFileReader<ImageType> ImageReaderType;
	ImageReaderType::Pointer readerS = ImageReaderType::New();
	readerS->SetFileName(refImage);
	readerS->Update();
	ImageType::Pointer volume = readerS->GetOutput();

	// FS Version
	MRI *image = MRIread(refImage);

	//Outputting the Files to Ensure the correct files were input
	cerr << endl 
	     << "TRK File:           " << TRKFile.at(0)  << endl 
	     << "Left Surface File:  " << surfaceFileL   << endl 
	     << "Left Overlay File:  " << overlayFileL   << endl
	     << "Right Surface File: " << surfaceFileR   << endl 
	     << "Right Overlay File: " << overlayFileR   << endl
	     << "Reference Image:    " << refImage << endl;

	// Loading the TRK files into a mesh
        vector<ColorMeshType::Pointer>* meshes;
        vector<vtkSmartPointer<vtkPolyData>> polydatas;

        ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();
        clusterTools->GetPolyDatas(TRKFile, &polydatas, volume);
        meshes = clusterTools->PolydataToMesh(polydatas);

	// Loading the Surface files and initialization of the KdTree
	// Left Hemisphere
	MRI_SURFACE *surfL;
        surfL = MRISread(surfaceFileL);
	
	vtkSmartPointer<vtkPolyData> surfVTK_L = FSToVTK(surfL);

        vtkSmartPointer<vtkKdTreePointLocator> surfTreeL = vtkSmartPointer<vtkKdTreePointLocator>::New();
        surfTreeL->SetDataSet(surfVTK_L);
        surfTreeL->BuildLocator();
	
	// Right Hemisphere
	MRI_SURFACE *surfR;
        surfR = MRISread(surfaceFileR);

	vtkSmartPointer<vtkPolyData> surfVTK_R = FSToVTK(surfR);

        vtkSmartPointer<vtkKdTreePointLocator> surfTreeR = vtkSmartPointer<vtkKdTreePointLocator>::New();
        surfTreeR->SetDataSet(surfVTK_R);
        surfTreeR->BuildLocator();

	// Variables able to hold a point in two different types
	PointType point;
        point.Fill(0);
        double point_array[3];

	// Initialization of a streamline
	ColorMeshType::Pointer input = (*meshes)[0];
        ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();

	// Cycles through the points of the streamlines
	for (; inputCellIt != input->GetCells()->End(); ++inputCellIt)
	{
		CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
                input->GetPoint(*it, &point);

		ImageType::IndexType index;
		volume->TransformPhysicalPointToIndex(point, index);
		
		MRIvoxelToSurfaceRAS(image, index[0], index[1], index[2], &point_array[0], &point_array[1], &point_array[2]);

		// Finds closest point and sets value equal to ENDPOINT_VALUE
		double distL, distR;
		vtkIdType Left_ID  = surfTreeL->FindClosestPointWithinRadius(100000, point_array, distL);
                vtkIdType Right_ID = surfTreeR->FindClosestPointWithinRadius(100000, point_array, distR);
                vtkIdType ID = which_ID(distL, distR, Left_ID, Right_ID);
	
		if (ID == Left_ID)
			surfL->vertices[ID].curv = ENDPOINT_VALUE;
		else
			surfR->vertices[ID].curv = ENDPOINT_VALUE;

		// Finding last point in the stream
		for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
                	input->GetPoint(*it, &point);

		volume->TransformPhysicalPointToIndex(point, index);

		MRIvoxelToSurfaceRAS(image, index[0], index[1], index[2], &point_array[0], &point_array[1], &point_array[2]);

		Left_ID  = surfTreeL->FindClosestPointWithinRadius(100000, point_array, distL);
                Right_ID = surfTreeR->FindClosestPointWithinRadius(100000, point_array, distR);
                ID = which_ID(distL, distR, Left_ID, Right_ID);
	
		if (ID == Left_ID)
			surfL->vertices[ID].curv = ENDPOINT_VALUE;
		else
			surfR->vertices[ID].curv = ENDPOINT_VALUE;

	}

	MRISwriteCurvature(surfL, overlayFileL);
	MRISwriteCurvature(surfR, overlayFileR);

	return 0;
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

/* Function: which_ID
 * Input: the two distances and the two vertice IDs
 * Return: whichever vertice is closer to the point
 * Does: Compares the two distances and returns the vertice of the shorter distance
 */
vtkIdType which_ID(double n1, double n2, vtkIdType ID1, vtkIdType ID2)
{
        if (n1 < n2)
                return ID1;

        return ID2;
}
