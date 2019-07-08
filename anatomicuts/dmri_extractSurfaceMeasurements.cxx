/* Alexander Zsikla
 * Professor Siless
 * August 2019
 * dmri_extractSurfaceMeasurements.cxx
 *
 * 
 *
 */

#include <iostream>
#include <fstream>
#include <string>

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
		     << arg[0] << " -s streamlineFile.trk"
		     << " -f surfaceFile.orig -o output.csv"
		     << endl;

		return -1;
	}

	// Input Parsing
	const char *stream  = num1.follow("input.trk", "-s");
	const char *surface = num1.follow("lh.orig", "-f");
	const char *output  = num1.follow("output.csv", "-o");

	// Looping for multiple file input of a certain kind of file
	/*vector<string> inputFiles;

	for (string inputName = string(num1.follow("files for whatever", "-s")); access(inputName.c_str(),0) == 0; inputName = string(num1.next("")))
		inputFiles.push_back(inputName);
	*/

	// TO DELETE
	cout << stream << endl << surface << endl << output << endl;

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
 * 1) Why do you have GetPot cl2 when it is only used once (line 65)
 * 2) 
 */

