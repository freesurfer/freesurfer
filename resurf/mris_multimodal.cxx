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
#include "fsSurfaceOptimizationFilter.h"
#include "itkVTKPolyDataWriter.h"
		

 

	#include "mrisurf.h"


int main(int narg, char*  arg[])
{
//	try{

		constexpr unsigned int Dimension = 3;
		typedef float CoordType;
		typedef fs::Surface< CoordType, Dimension> SurfType;
		typedef fs::SurfaceOptimizationFilter< SurfType, SurfType> SurfFilterType;

		GetPot cl(narg, const_cast<char**>(arg));
		if(cl.size()==1 || cl.search(2,"--help","-h"))
		{
			std::cout<<"Usage: " << std::endl;
			std::cout<< arg[0] << " -i surface -o surface -vtk "  << std::endl;   
			return -1;
		}
		const char *inSurfFilename= cl.follow ("", "-i");
		const char *outSurfFilename = cl.follow ("", "-o");


		MRI_SURFACE *surf;
		surf = MRISread(inSurfFilename);

		SurfType::Pointer surface =  SurfType::New();
		surface->Load(&*surf);

		SurfFilterType::Pointer filter =  SurfFilterType::New();
		filter->SetInput(surface);
		filter->Update();
		surface =  filter->GetOutput();
		surf = surface->GetFSSurface(&*surf);

		if (cl.search("-vtk") )
		{
			typedef itk::VTKPolyDataWriter<SurfType> WriterType;
			WriterType::Pointer writer = WriterType::New();
			writer->SetInput(surface);
			writer->SetFileName(outSurfFilename);
			writer->Update();
		}
		else
		{	
			MRISwrite(surf,outSurfFilename);
 			MRISfree(&surf);	
		}	
/*	}catch(...)
	{
		std::cout << "Error --> ";
		for(int i=0;i<narg;i++)
		{
			std::cout << arg[i];
		}
		std::cout << std::endl;

		return EXIT_FAILURE;
	}
*/	return EXIT_SUCCESS;
}
