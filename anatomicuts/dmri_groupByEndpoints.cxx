/* Andrew Zhang
 * Professor Siless
 * dmri_groupingByEndpoints.cxx
 *
 */

#include <iostream>
#include <string>
#include <string.h>

#include "GetPot.h"

using namespace std;

int main(int narg, char* arg[]) 
{
	GetPot c1(narg, const_cast<char**>(arg));
	GetPot c2(narg, const_cast<char**>(arg)); 

	if (c1.size() == 1 || c1.search(2, "--help", "-h"))
	{
		cout << "Usage: " << endl; 
		cout << arg[0] << " -s " 
		     << "streamline -i parcelled_image "
		     << "-o folder" 
		     << endl;
		return -1; 
	}

	const char *stream = c1.follow("string_file.trk", "-s"); 
	const char *image = c1.follow("image_file.nii.gz", "-i"); 
	const char *folder = c1.follow("output_folder", "-o"); 

	cout << stream << endl << image << endl << folder << endl; 

	return 0; 	

}
