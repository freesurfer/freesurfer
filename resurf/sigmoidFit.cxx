#include <iostream>
//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interpolation.h"

#include "itkVector.h"

#include "vtkVersion.h"
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


#include "fsSurfaceOptimizationFilter.h"
#include "itkVTKPolyDataWriter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVariableLengthVector.h"
#include "mrisurf.h"
#include "colortab.h"
#include "fsenv.h"

#include "GetPot.h"
using namespace alglib;
typedef float CoordType;
typedef itk::Image< CoordType, 4 >         ProfileType;
typedef itk::ImageFileReader<ProfileType> ReaderType;
typedef fs::Surface< CoordType, 3> SurfaceType;


void function_cx_1_func(const real_1d_array &c, const real_1d_array &x, double &func, void *ptr) 
{
	// this callback calculates f(c,x)=exp(-c0*sqr(x0))
	// where x is a position on X-axis and c is adjustable parameter
	//func = exp(-c[0]*pow(x[0],2));
	/*float a1 = c[0];
	float x1 = c[1];
	float e1= c[2];
	float a2 = c[3];
	float x2 = c[4];
	float e2 = c[5];
	float offset = c[6];

	float s1=tanh(e1*(x[0]-x1));
	float s2=tanh(e2*(x[0]-x2));
	float y =  offset +  a1*s1  - a2*s2 ;
	*/
	func = c[6] + c[0]*tanh(c[1]*(x[0]-c[2])) - c[3]*tanh(c[4]*(x[0]-c[5])) ;
}

int main(int narg, char **arg)
{
  	
	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -s surface -i image -p profile 4D image -o output fitted image 4D "  << std::endl;   
		return -1;
	}

	const char *surfaceFile= cl.follow ("", "-s");
	const char *imageFile= cl.follow ("", "-i");
	const char *profileFile= cl.follow ("", "-p");
	const char *fittingFile = cl.follow ("", "-o");
	

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(profileFile);
	reader->Update();

	ProfileType::Pointer image =reader->GetOutput();

	MRI_SURFACE *surf;
	surf = MRISread(surfaceFile);

	SurfaceType::Pointer surface =  SurfaceType::New();
	surface->Load(&*surf);

 	MRI *imageFS =  MRIread(imageFile) ;
	
	std::vector<SurfaceType::PointType>  points;
	int halfSize= image->GetLargestPossibleRegion().GetSize()[3]/2;
	for(SurfaceType::CellsContainerConstIterator itCells = surface->GetCells()->Begin();itCells < surface->GetCells()->End();itCells+=10)
	{
		SurfaceType::CellType::PointIdIterator pointsIt = itCells.Value()->PointIdsBegin();
		SurfaceType::PointType p1 = surface->GetPoint(*pointsIt);
		SurfaceType::PointType p2 = surface->GetPoint(*(++pointsIt));
		SurfaceType::PointType p3 = surface->GetPoint(*(++pointsIt));
		SurfaceType::PointType edge1 = p1-p2;
		SurfaceType::PointType edge2 = p1-p3;
		
		SurfaceType::PointType mean ;
		mean.Fill(0.0);
		for(int i=0;i<3;i++)
			mean[i]= (p1[i]+p2[i]+p3[i])/3.0;
	
		ProfileType::IndexType index; 
	
		double x,y,z;
		MRISsurfaceRASToVoxel(surf, imageFS, mean[0], mean[1],mean[2], &x,&y,&z);
		index[0]=x;
		index[1]=y;
		index[2]=z;
		index[3]=0;
		std::string y_fit = std::string("[") + std::to_string(image->GetPixel(index)) ;
		std::string x_fit = std::string("[[")+ std::to_string(-halfSize)+std::string("]");
		for(int  t=1; t<image->GetLargestPossibleRegion().GetSize()[3];t++)
		{
			index[3]=t;
			image->GetPixel(index);
			y_fit+= std::string(",") + std::to_string(image->GetPixel(index)) ;
			x_fit+= std::string(",[")+std::to_string(t-halfSize)+std::string("]");
		}

		y_fit+= std::string("]");
		x_fit+= std::string("]");
		real_1d_array yval( y_fit.c_str());
		real_2d_array xval( x_fit.c_str());

		real_1d_array c = "[1, .1, 2,1, .1,2, 120]";
		double epsx = 0.01;
		ae_int_t maxits = 1000;
		ae_int_t info;
		lsfitstate state;
		lsfitreport rep;
		double diffstep = 0.01;

		lsfitcreatef(xval, yval, c, diffstep, state);
		lsfitsetcond(state, epsx, maxits);
		alglib::lsfitfit(state, function_cx_1_func);
		lsfitresults(state, info, c, rep);
		//printf("%d\n", int(info)); // EXPECTED: 2
		//printf("%s\n", c.tostring(1).c_str()); // EXPECTED: [1.5]
	
		for(int  t=0; t<image->GetLargestPossibleRegion().GetSize()[3];t++)
		{
			index[3]=t;
			float value = c[6] + c[0]*tanh(c[1]*(t-halfSize-c[2])) - c[3]*tanh(c[4]*(t-halfSize-c[5])) ;
			image->SetPixel(index, value);
		}

	}
	typedef  itk::ImageFileWriter< ProfileType  > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fittingFile);
	writer->SetInput(image);
	writer->Update();

	return 0;
}

