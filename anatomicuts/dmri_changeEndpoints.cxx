/* Alexander Zsikla
 * Professor Siless
 * Summer 2019
 * dmri_changeEndpoints.cxx
 *
 * Changing the endpoint values to a different value
 *
 */

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

int main(int narg, char* arg[])
{
	GetPot input(narg, const_cast<char**>(arg));

	// Checking for correct parameters
	if ((input.size() <= 3) or (input.search(2, "--help", "-h")))
        {
                cerr << "Usage: " << endl
                     << arg[0] << " -i streamlineFile.trk -s surfaceFile.orig -o outputDirectory" << endl;

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
	TRKFile.push_back(input.follow("Could not find TRK file", "-i"));
	const char *surface    = input.follow("Could not find Surface File", "-s");
	const char *outputDir  = input.follow("Could not find Output Directory", "-o");

	//Outputting the Files to Ensure the correct files were input
	cerr << endl 
	     << "TRK File:     " << TRKFile.at(0) << endl 
	     << "Surface File: " << surface << endl 
	     << "Output Dir:   " << outputDir << endl;

	// Loading the TRK files into a mesh
	ImageType::Pointer mask;

        vector<ColorMeshType::Pointer>* meshes;
        vector<vtkSmartPointer<vtkPolyData>> polydatas;

        ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();
        clusterTools->GetPolyDatas(TRKFile, &polydatas, mask);
        meshes = clusterTools->PolydataToMesh(polydatas);

	return 0;
}












