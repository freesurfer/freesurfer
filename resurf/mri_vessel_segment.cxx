#include "mri.h"
#include "mrisurf.h"
#include "fsenv.h"
#include "cma.h"

#include "GetPot.h"

int main(int narg, char * arg[])
{
	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -t1 t1.mgz -t2 t2.mgz -aseg aseg.mgz -o output.mgz "  << std::endl;   
		return -1;
	}

	const char *imageNameT1= cl.follow ("", "-t1");
	const char *imageNameT2 = cl.follow ("", "-t2");
	const char *outputName = cl.follow ("", "-o");
	
	MRI* imageAseg = NULL;
	if(cl.search("-aseg"))
	{
		const char *imageNameAseg = cl.follow ("", "-aseg");
 		imageAseg =  MRIread(imageNameAseg) ;
	}
	
	
 	MRI *imageT1 =  MRIread(imageNameT1) ;
 	MRI *imageT2 =  MRIread(imageNameT2) ;
 	MRI *output =  MRIcopy(imageT1, NULL) ;


	for (int x = 0 ; x < imageT1->width ; x++)
	{
		for (int y = 0 ; y < imageT1->height ; y++)
		{
			for (int z = 0 ; z < imageT1->depth ; z++)
			{
			        int label =(imageAseg)? MRIgetVoxVal(imageAseg, x, y, z, 0) :0;
				float T1 = MRIgetVoxVal(imageT1, x, y, z, 0) ;
				float T2 = MRIgetVoxVal(imageT2, x, y, z, 0) ;
				if ( IS_CORTEX(label) &&  1.1*T1 > T2)
				{
					MRIsetVoxVal(output, x, y, z, 0, 1) ;
				}
				else
				{
					MRIsetVoxVal(output, x, y, z, 0, 0) ;
				}
			}
		}
	}

	MRIwrite(output,outputName) ;
 	MRIfree(&imageT1);	
 	MRIfree(&imageT2);	
 	if(imageAseg)
		MRIfree(&imageAseg);	
 	MRIfree(&output);	
}
