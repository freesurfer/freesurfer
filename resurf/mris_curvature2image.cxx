#include "itkKdTree.h"

#include "vtkPolyData.h"
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
#include "itkVTKPolyDataWriter.h"
#include <vnl/vnl_cross.h>
#include <cmath>
#include "vtkIdList.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkKdTreePointLocator.h"
#include "mrisurf.h"
int main(int narg, char * arg[])
{
	constexpr unsigned int Dimension = 3;
	typedef float CoordType;
  	
	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -s surface -m mask -o outputOverlayImage -d outputDistanceImage -c overlay -l label -inv -r radius"  << std::endl;   
		return -1;
	}

	const char *surfFilename= cl.follow ("", "-s");
	const char *outputImageFilename = cl.follow ("", "-o");
	const char *displayImageFilename = cl.follow ("", "-d");
	const char *overlayFilename = cl.follow ("", "-c");
	const char *maskFilename = cl.follow("","-m");
	const int label = cl.follow(1,"-l");
	const float  radius= cl.follow(1.1,"-r");
	
	MRI_SURFACE *surf;
	surf = MRISread(surfFilename);

 	MRI *image =  MRIread(maskFilename) ;
	MRI *output = MRIalloc(image->width, image->height, image->depth, MRI_FLOAT);
	MRIcopyHeader(image, output);
	
	MRISreadCurvature(surf,overlayFilename) ;
	
	if(! cl.search("-inv"))
	{
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		for( int i = 0; i < surf->nvertices; i++ )
		{
			points->InsertNextPoint(surf->vertices[i].x, surf->vertices[i].y, surf->vertices[i].z);
		}
		vtkSmartPointer<vtkPolyData> pdPoints = vtkSmartPointer<vtkPolyData>::New();
		pdPoints->SetPoints(points);

		vtkSmartPointer<vtkKdTreePointLocator> kDTree =	vtkSmartPointer<vtkKdTreePointLocator>::New();
		kDTree->SetDataSet(pdPoints);
		kDTree->BuildLocator();


		for(int i=0;i<image->width;i++)
		{
			for(int j=0;j<image->height;j++)
			{
				for(int k=0;k<image->depth;k++)
				{
					float val = MRIgetVoxVal(image,i,j,k,0);
					if(val  == label)
					{
						double point[3];
						//MRIvoxelToWorld( image, i,j,k, &point[0],&point[1],&point[2]);
						MRIvoxelToSurfaceRAS( image, i,j,k, &point[0],&point[1],&point[2]);
						int id = kDTree->FindClosestPoint(point);
						//					std::cout << point[0] << " " << point[1] << " " << point[2] << " " << id <<  " " << surf->vertices[id].curv << std::endl;
						MRIsetVoxVal(output,i,j,k,0, (float)surf->vertices[id].curv);			
					}
					else
					{
						MRIsetVoxVal(output,i,j,k,0, (float)0.0);			
					}
				}
			}
		}

		//	MRISsurfaceRASToVoxel(surf, image, pt[0], pt[1],pt[2], &x,&y,&z);
		MRIwrite(output,outputImageFilename);
		MRISfree(&surf);	
	}
	else
	{	
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		int w =0;
		for(int i=0;i<image->width;i++)
		{
			for(int j=0;j<image->height;j++)
			{
				for(int k=0;k<image->depth;k++)
				{
					//double point[3];
					//MRIvoxelToSurfaceRAS( image, i,j,k, &point[0],&point[1],&point[2]);
					points->InsertNextPoint(i,j,k);
					w++;
				}
			}
		}
		
		vtkSmartPointer<vtkPolyData> pdPoints = vtkSmartPointer<vtkPolyData>::New();
		pdPoints->SetPoints(points);

		vtkSmartPointer<vtkKdTreePointLocator> kDTree =	vtkSmartPointer<vtkKdTreePointLocator>::New();
		kDTree->SetDataSet(pdPoints);
		kDTree->BuildLocator();

		for(int j=0;j<surf->nvertices;j++)
		{
			double pt[3];
			MRISsurfaceRASToVoxel(surf, image, surf->vertices[j].x, surf->vertices[j].y, surf->vertices[j].z,&pt[0],&pt[1],&pt[2]);
			vtkNew<vtkIdList> listIds;
			kDTree->FindPointsWithinRadius(radius, pt,listIds.GetPointer());
			double val=0;
			int total = 0;
			for(int i=0; i< listIds->GetNumberOfIds();i++)
			{	
				points->GetPoint(listIds->GetId(i), pt);
				double tmpval =  MRIgetVoxVal(image, pt[0],pt[1],pt[2],0);
				if (tmpval >0)
				{
					val+= tmpval ;
					total++;
				}
			}
			surf->vertices[j].curv = val/total;
		}

		//MRIwrite(output,outputImageFilename);
		
		MRISwriteCurvature(surf,overlayFilename) ;

		MRISfree(&surf);	
	}
}
