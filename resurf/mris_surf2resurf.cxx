#include <iostream>
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "colortab.h"
#include "fsenv.h"
extern "C" {
#include "mrisurf.h"
}

#include "mrisurf.h"

int main(int narg, char*  arg[])
{

/*COLOR_TABLE *ct;

  FSENV *fsenv = FSENVgetenv();
  char tmpstr[2000];
  sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
  ct = CTABreadASCII(tmpstr);
  */

   MRI_SURFACE *surf;

  surf = MRISread(arg[1]);

	try{

		enum {Dimension =3};
		typedef double                                                        PixelType;


		GetPot cl(narg, const_cast<char**>(arg));
		if(cl.size()==1 || cl.search(2,"--help","-h"))
		{
			std::cout<<"Usage: " << std::endl;
			std::cout<< arg[0] << " -i surface -o surface"  << std::endl;   
			return -1;
		}
		int numClusters = cl.follow(0,"-c");
		const char *segFile = cl.follow ("", "-s1");

		segFile = cl.follow ("", "-s2");
	}catch(...)
	{
		std::cout << "Error --> ";
		for(int i=0;i<narg;i++)
		{
			std::cout << arg[i];
		}
		std::cout << std::endl;

	}
	return 0;
}
