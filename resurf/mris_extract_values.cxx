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

#include "annotation.h"
#include "surfcluster.h"
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
			std::cout<< arg[0] << " -i surface -v overlay -a annotation -o csvfile -m numberOfImages image1 image2 imaga3"  << std::endl;   
			return -1;
		}
		const char *inSurf= cl.follow ("", "-i");
		const char *outCsv = cl.follow ("", "-o");
		const char *overlayFilename= cl.follow ("", "-v");
		const char *annotFilename= cl.follow ("", "-a");
		
		MRI_SURFACE *surf;
		surf = MRISread(inSurf);
		MRIScomputeMetricProperties(surf);
		MRISstoreMetricProperties(surf);

		int imageNumber = cl.follow(0,"-m");
		if( cl.search("-a") )
		{
			MRISreadCurvature(surf,overlayFilename) ;
			MRISreadAnnotation(surf,annotFilename) ;
			sclustCountClusters(surf);
			//read_named_annotation_table(annotFilename);
			std::vector<int> vector = readAnnotationIntoVector(annotFilename);
			std::map<int,int> myannot; 
			for(int i=0;i<vector.size();i++)
			{
				myannot[vector[i]]=i;
			}
		
			std::map<int,std::vector<float>> clusterVals;
			for (unsigned j=0;j<surf->nvertices;j++)
			{
				//std::cout <<  surf->vertices[j].undefval ;
				if( clusterVals.count(surf->vertices[j].annotation) ==0)
					clusterVals[surf->vertices[j].annotation]= std::vector<float>();
				clusterVals[surf->vertices[j].annotation].push_back(surf->vertices[j].curv);

			}

			MRISfree(&surf);	
			std::ofstream csv;
			csv.open(outCsv);
			csv << " cluster, mean,std  "<< std::endl;
			int i =0;
			for( auto& x:clusterVals)
			{	
				float mean = std::accumulate(x.second.begin(), x.second.end(), 0.0)/x.second.size();
				double accum = 0.0;
				std::for_each (std::begin(x.second), std::end(x.second), [&](const double d) {
					accum += (d - mean) * (d - mean);
				});

				double stdev = sqrt(accum / (x.second.size()-1));
				csv << i << "," <<   x.first <<  ","<< mean << ","<<stdev << std::endl;
				i++;
			}
			csv.close();
		}
		else
		{
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
		}
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
