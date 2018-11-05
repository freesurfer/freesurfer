#include <iostream>
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "fsSurface.h"
#include "itkTriangleCell.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "colortab.h"
#include "fsenv.h"

extern "C" 
{
	#include "mrisurf.h"
}

int main(int narg, char*  arg[])
{
	try{

		constexpr unsigned int Dimension = 3;
		typedef double CoordType;
		typedef fs::Surface< CoordType> SurfType;

		GetPot cl(narg, const_cast<char**>(arg));
		if(cl.size()==1 || cl.search(2,"--help","-h"))
		{
			std::cout<<"Usage: " << std::endl;
			std::cout<< arg[0] << " -i surface -o surface"  << std::endl;   
			return -1;
		}
		const char *inSurfFilename= cl.follow ("", "-i");
		const char *outSurfFilename = cl.follow ("", "-o");


		MRI_SURFACE *surf;
		surf = MRISread(inSurfFilename);

		SurfType::Pointer surface =  SurfType::New();
		surface->Load(&*surf);
		surf = surface->GetFSSurface(&*surf);
		/*typedef itk::VTKPolyDataWriter<MeshType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(mesh);
		writer->SetFileName(outSurf);
		writer->Update();*/
		
		MRISwrite(surf,outSurfFilename);
 		MRISfree(&surf);	
	}catch(...)
	{
		std::cout << "Error --> ";
		for(int i=0;i<narg;i++)
		{
			std::cout << arg[i];
		}
		std::cout << std::endl;

		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
