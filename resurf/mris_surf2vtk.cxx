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
			std::cout<< arg[0] << " -i surface -o surface"  << std::endl;   
			return -1;
		}
		const char *inSurf= cl.follow ("", "-i");
		const char *outSurf = cl.follow ("", "-o");


		MRI_SURFACE *surf;
		surf = MRISread(inSurf);


		MeshType::Pointer mesh = MeshType::New();
		PointsContainer::Pointer points = PointsContainer::New();
		points->Reserve( surf->nvertices );

		for( int i = 0; i < surf->nvertices; i++ )
		{
			PointType p;
			p[0]=surf->vertices[i].x;
			p[1]=surf->vertices[i].y;
			p[2]=surf->vertices[i].z;
			points->SetElement( i, p );
		}

		mesh->SetPoints( points );
		CellType::CellAutoPointer cellpointer;

		for( int i = 0; i < surf->nfaces; i++ )
		{
			cellpointer.TakeOwnership( new TriangleType );
			cellpointer->SetPointId( 0, surf->faces[i].v[0]);
			cellpointer->SetPointId( 1, surf->faces[i].v[1]);
			cellpointer->SetPointId( 2, surf->faces[i].v[2]);
			mesh->SetCell( i, cellpointer );
		}

		typedef itk::VTKPolyDataWriter<MeshType> WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(mesh);
		writer->SetFileName(outSurf);
		writer->Update();
	
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
