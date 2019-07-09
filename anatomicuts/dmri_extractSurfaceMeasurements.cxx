/* Alexander Zsikla
 * Professor Siless
 * August 2019
 * dmri_extractSurfaceMeasurements.cxx
 *
 * Takes in a cluster and outputs metrics, specifically thickness and curvature, based on the endpoints of the connections and outputs it to an external CSV file
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstlib>

#include <vtkPolyData.h>
#include "TrkVTKPolyDataFilter.txx"
#include "GetPot.h"
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

	// Input Parsing
	
	// Looping for multiple file input of a certain kind of file
	vector<string> inputFiles;

	for (string inputName = string(num1.follow("input.trk", "-i")); access(inputName.c_str(),0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);

	const char *surface = num1.follow("lh.orig", "-s"); 	// this is a surface
	const char *output  = num1.follow("output.csv", "-o");


	// TO DELETE
	cout << surface << endl << output << endl;

	// Reading in TRK File
	

	// Opening output File
	ofstream oFile;
	oFile.open(output);

	if (not oFile.is_open()) {
		cerr << "Could not open output file" << endl;
		return -1;
	}

	oFile.close();

	return 0;
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

