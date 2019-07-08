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
#include "ClusterTools.h"
#include "itkDefaultStaticMeshTraits.h"

using namespace std;

int main(int narg, char* arg[]) 
{
	//Define variables
	GetPot c1(narg, const_cast<char**>(arg));
	GetPot c2(narg, const_cast<char**>(arg)); 

	enum {Dimension =3};
	typedef int                                                        PixelType;
	const unsigned int PointDimension = 3;
	typedef std::vector<int>                  PointDataType;
	
	typedef itk::Mesh< PixelType, PointDimension > ColorMeshType;
	typedef ColorMeshType::PointType PointType;
	typedef ColorMeshType::CellType        CellType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;

	//Usage error
	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s streamlines -i imageFile -d outputDirectory" << endl;
		return -1; 
	}

	//Take in information
	const char *stream = c1.follow("string_file.trk", "-s"); 
	const char *image = c1.follow("image_file.nii.gz", "-i"); 
	const char *output = c1.follow("output_directory", "-o"); 

	//Declaration of vectors
	vector<ColorMeshType::Pointer>* meshes; 
	vector<ColorMeshType::Pointer>* fixMeshes; 
	vector<vtkSmartPointer<vtkPolyData>> polydatas; 

	//TO DELETE
	cout << stream << endl << image << endl << output << endl; 

	return 0; 	

}

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
