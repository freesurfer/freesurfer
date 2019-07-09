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

// to use access
/*#include <vtkPolyData.h>
#include "TrkVTKPolyDataFilter.txx"*/

// Input Splicing
#include "GetPot.h"

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
	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); inputName != "-s" and inputName != "-o"; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);

	// Could not figure out how access() worked so I tried hardcoding it in
	/*vector<string> inputFiles;
	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);
	*/

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
	



	//
	// Reading in surface File
	//

	MRI_SURFACE *surf;
        //surf = MRISread(surface);
	
	/*for (unsigned i = 0; i < surf->nvertices; i++)
	{
		double x, y, z;
		float pdx, pdy, pdz;

		if (surf->vertices[j].x > 1)
		{
			MRISsurfaceRASTToVoxel(surf, images, surf->vertices[j].x, surf->vertices[j].y, surf->vertices[j].z, &x, &y, &z);
			float magnitud = MRIvoxelGradient(images, (float) x, (float) y, (float) z, &pdx, &pdy, &pdz);
	}*/

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

