#include <iostream>
#include <fstream>
#include <string>         
#include "itkImageFileReader.h"
#include "GetPot.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include "itkPolylineCell.h"
#include "TrkVTKPolyDataFilter.txx"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "itkImageDuplicator.h"
#include "itkNeighborhoodIterator.h"
#include <time.h>

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
#include "mris_multimodal_refinement.h"
#include "itkImage.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "itkImageFileReader.h"

#include <vtkPolyLine.h>
#include <vtkPolyData.h>

int main(int narg, char*  arg[])
{
//	try{

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
			std::cout<< arg[0] << " -i surface -v overlay -o csvfile -m numberOfImages image1 image2 imaga3"  << std::endl;   
			return -1;
		}
		const char *inSurf= cl.follow ("", "-i");
		const char *outCsv = cl.follow ("", "-o");
		const char *overlayFilename= cl.follow ("", "-v");
		MRI_SURFACE *surf;
		surf = MRISread(inSurf);
		MRIScomputeMetricProperties(surf);
		MRISstoreMetricProperties(surf);

		int imageNumber = cl.follow(0,"-m");
		
		std::vector<MRI*> images; 

		std::vector<std::string> fileNames;
		for(;imageNumber>0;imageNumber--)
		{
			fileNames.push_back(cl.next(""));
			MRI *imageFS =  MRIread(fileNames[fileNames.size()-1].c_str()) ;
			images.push_back(imageFS);
		}

		double x,y,z;
		MRISreadCurvature(surf,overlayFilename) ;
		
		std::ofstream csv;
		csv.open(outCsv);
		csv << " Intensity, PrevIntensity, NextIntensity, Overlay "<< std::endl;

		for (unsigned j=0;j<surf->nvertices;j++)
		{
			
			double xv, yv,zv, x,y,z, nx,ny,nz,val;
		        float point[3]={0,0,0};
			float normal[3]={0,0,0};
	
			nx= surf->vertices[j].nx;
			ny=surf->vertices[j].ny;
			nz=surf->vertices[j].nz;
			float dist = sqrt(nx*nx + ny*ny + nz*nz);
			if( dist>0)
			{
				nx /= dist;
				ny /= dist;
				nz /= dist;
			}	

			MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &xv,&yv,&zv);
			MRIsampleVolume(images[0], xv, yv, zv, &val);
			csv << val <<", ";
			x=surf->vertices[j].x -nx;
			y=surf->vertices[j].y -ny;
			z=surf->vertices[j].z -nz;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);

			MRIsampleVolume(images[0], xv, yv, zv, &val);
			csv << val <<", ";
			x=surf->vertices[j].x +nx;
			y=surf->vertices[j].y +ny;
			z=surf->vertices[j].z +nz;

			MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);

			MRIsampleVolume(images[0], xv, yv, zv, &val);
			csv << val <<", ";
			csv<< surf->vertices[j].curv << std::endl;
		}
		MRISfree(&surf);	
		csv.close();

	/*
	}catch(...)
	{
		std::cout << "Error --> ";
		for(int i=0;i<narg;i++)
		{
			std::cout << arg[i];
		}
		std::cout << std::endl;

		return EXIT_FAILURE;
	}*/
	return EXIT_SUCCESS;
}
