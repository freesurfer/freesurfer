#include <iostream>
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "colortab.h"
#include "fsenv.h"
 
#include "mrisurf.h"
#include "itkVTKPolyDataWriter.h"

int main(int narg, char*  arg[])
{
	try{

		constexpr unsigned int Dimension = 3;
		typedef double CoordType;
		typedef itk::Mesh< CoordType, Dimension > MeshType;
		typedef MeshType::PointsContainer PointsContainer;
		typedef MeshType::PointType PointType;
		//typedef MeshType::PointIdentifier PointIdentifier;
		typedef MeshType::CellType CellType;
		typedef itk::TriangleCell< CellType > TriangleType;

		GetPot cl(narg, const_cast<char**>(arg));
		if(cl.size()==1 || cl.search(2,"--help","-h"))
		{
			std::cout<<"Usage: " << std::endl;
			std::cout<< arg[0] << " -i surface -o surface -l gradientDescentlambda -k iterations -m numberOfImages image1 image2 image3"  << std::endl;   
			return -1;
		}
		const char *inSurf= cl.follow ("", "-i");
		const char *outSurf = cl.follow ("", "-o");
		float lambda= cl.follow (0.001, "-l");
		int iterations = cl.follow (10, "-k");


		MRI_SURFACE *surf;
		surf = MRISread(inSurf);

		int imageNumber = cl.follow(0,"-m");
		
		std::vector<MRI*> images; 

		for(;imageNumber>0;imageNumber--)
		{
			
			MRI *imageFS =  MRIread(cl.next("")) ;
			images.push_back(imageFS);
		}

		for(int i=0;i<iterations;i++)
		{ 
			for (unsigned j=0;j<surf->nvertices;j++)
			{

				for(unsigned k=0;k<images.size();k++)
				{
					double x,y,z;
					float pdx, pdy, pdz;

					if(surf->vertices[j].x >1)
					{
						MRISsurfaceRASToVoxel(surf, images[k], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &x,&y,&z);
						float magnitud = MRIvoxelGradient(images[k], (float) x, (float) y,(float) z, &pdx,  &pdy, &pdz);	
						float norm= magnitud/images.size();
						if(norm>1)
						{
//						std::cout << magnitud << " "<< pdx << " "<< pdy << " " << pdz << std::endl;
						surf->vertices[j].x += lambda *  pdx /norm;
						surf->vertices[j].y += lambda *  pdy/norm;
						surf->vertices[j].z += lambda *  pdz/norm;
						}
					}


				}
			}

		}
				
		MRISwrite(surf,outSurf);
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
