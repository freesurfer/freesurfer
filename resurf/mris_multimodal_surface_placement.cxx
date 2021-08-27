#include <iostream>
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
			std::cout<< arg[0] << " -i surface -o surface -b spheresurf -n normals.vtk -v values.vtk -d debugVertex -s step_size -k numberOfSteps  -g gradientSigma -a aseg.aparc  -w whitesurface -p pOfCSF  -min/max -t1 image -t2 image -flair image"  << std::endl;   
			return -1;
		}
		const char* inSurf= cl.follow ("", "-i");
		const char* inSphere= cl.follow ("", "-b");
		const char* whSurf= cl.follow ("", "-w");
		const char* outSurf = cl.follow ("", "-o");
		const char* outNormals = cl.follow ("", "-n");
		const char* outValues= cl.follow ("", "-v");
		const char* overlayFilename= cl.follow ("", "-p");
		int debugVertex= cl.follow (-1, "-d");
		float step_size= cl.follow (.4, "-s");
		int numberOfSteps= cl.follow (20, "-k");
		const char* asegFile= cl.follow ("", "-a");
		float gradientSigma= cl.follow (.20, "-g");
		bool maxGradient = !cl.search("-min");
		std::cout << maxGradient << std::endl;

		MRI_SURFACE* whiteSurf;
		whiteSurf = MRISread(whSurf);
		MRI_SURFACE* surf;
		surf = MRISread(inSurf);
		MRIScomputeMetricProperties(surf);
		MRISstoreMetricProperties(surf);

		MRIS* sph = MRISread(inSphere);
		//	MRIScomputeMetricProperties(sphere);
		MRI_SP* sphere =  MRIStoParameterization(sph,NULL, 1,0);

		std::cout << " sphere " << sph->vertices[9].phi << " " << sph->vertices[9].theta << " " << MRISaverageRadius(sph) << std::endl;
/*		MRISedges(surf);
		MRIScorners(surf);
		MRISfaceMetric(surf,0);
		MRISedgeMetric(surf,0);
		MRIScornerMetric(surf,0);	
*/
		std::cout << "debug vertex " << debugVertex << " " <<surf->vertices[debugVertex].x << " "  << surf->vertices[debugVertex].y<< " " <<surf->vertices[debugVertex].z << std::endl; 


		MRIS_MultimodalRefinement* t2refinement = new MRIS_MultimodalRefinement();	
		std::vector<MRI*> images; 
		int modality=0;
		if( cl.search("-t1"))
		{
			MRI *t1 = MRIread(cl.follow("","-t1"));
			images.push_back(t1);
			t2refinement->addImage(t1);
		}
		if ( cl.search("-t2"))	
		{
			MRI *t2 = MRIread(cl.follow("","-t2"));
			images.push_back(t2);
			t2refinement->addImage(t2);
			modality=1;
		}
		if( cl.search("-flair"))
		{
			MRI *flair = MRIread(cl.follow("","-flair"));
			images.push_back(flair);
			t2refinement->addImage(flair);
			modality=2;
		}

		MRI* whiteMR= MRIcopy(images[0], NULL);
		MRI* vesselMR= MRIcopy(images[0], NULL);
		t2refinement->SegmentVessel(images[0], images[1],vesselMR, modality);
		t2refinement->SegmentWM(images[0], images[1],whiteMR,modality);
		t2refinement->SetWhiteMR(whiteMR);
		t2refinement->SetVesselMR(vesselMR);

		t2refinement->SetSphere (sph);
		for (unsigned j=0;j<surf->nvertices;j++)
		{
			
			surf->vertices[j].whitex = whiteSurf->vertices[j].x;
			surf->vertices[j].whitey = whiteSurf->vertices[j].y;
			surf->vertices[j].whitez = whiteSurf->vertices[j].z;
		}

		MRI *aseg= MRIread(asegFile);
	        t2refinement->SetSegmentation(aseg);
		
		t2refinement->SetVertexDebug(debugVertex);
		t2refinement->SetStep(step_size);
		t2refinement->SetNumberOfSteps(numberOfSteps);
		t2refinement->SetGradientSigma(gradientSigma);
		t2refinement->FindMaximumGradient(maxGradient);
		t2refinement->SetWhite(whiteSurf); //, debugVertex);
		t2refinement->getTarget(surf); //, debugVertex);
		double x,y,z;
		vtkSmartPointer<vtkPoints> points = vtkPoints::New();
		int totalPoints=0;
  		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
/*		typedef itk::Image<float,3> ImageType;
		typedef itk::ImageFileReader<itk::Image<float,3>> ImageReaderType;
		ImageReaderType::Pointer reader = ImageReaderType::New();
		reader->SetFileName ( fileNames[0]);
		reader->Update();
	
		ImageType::Pointer image = reader->GetOutput();
*/
		for (unsigned j=0;j<surf->nvertices;j++)
		{
			
			surf->vertices[j].x = surf->vertices[j].targx;
			surf->vertices[j].y = surf->vertices[j].targy;
			surf->vertices[j].z = surf->vertices[j].targz;

		
			double x,y,z, nx,ny,nz;
		        float point[3]={0,0,0};
			float normal[3]={0,0,0};
			/*
			MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].x,surf->vertices[j].y,surf->vertices[j].z, &x,&y,&z);
			MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].nx,surf->vertices[j].ny,surf->vertices[j].nz, &nx,&ny,&nz);

			ImageType::IndexType index;
			index[0]=x;
			index[1]=y;
			index[2]=z;
	
			image->TransformIndexToPhysicalPoint(index, point);
			index[0]=nx;
			index[1]=ny;
			index[2]=nz;
	
			image->TransformIndexToPhysicalPoint(index, normal);
			*/

			point[0]= surf->vertices[j].x;
			point[1]= surf->vertices[j].y;
			point[2]= surf->vertices[j].z;

			normal[0]= surf->vertices[j].nx;
			normal[1]= surf->vertices[j].ny;
			normal[2]= surf->vertices[j].nz;

			vtkIdType *ids = new vtkIdType [2];
			points->InsertPoint (totalPoints,point[0],point[1],point[2]);
			ids[0] = totalPoints;
			totalPoints++;
			points->InsertPoint (totalPoints, point[0]+normal[0],point[1]+normal[1],point[2]+normal[2]);
			ids[1] = totalPoints;
			totalPoints++;
			
			vtkSmartPointer<vtkPolyLine> polyLine =    vtkSmartPointer<vtkPolyLine>::New();
			polyLine->GetPointIds()->SetNumberOfIds(2);
			polyLine->GetPointIds()->SetId(0,ids[0]);
			polyLine->GetPointIds()->SetId(1,ids[1]);
		//	vtk->InsertNextCell (VTK_POLY_LINE, 2, ids);
			cells->InsertNextCell(polyLine);
		}

                 vtkSmartPointer<vtkPolyData> polyData  = vtkSmartPointer<vtkPolyData>::New();
                 polyData->SetPoints(points);
                polyData->SetLines(cells);

		vtkSmartPointer<vtkPolyDataWriter> pdWriter =  vtkSmartPointer<vtkPolyDataWriter>::New();
		
		#if VTK_MAJOR_VERSION <= 5	
			pdWriter->SetInput(polyData);
		#else
			pdWriter->SetInputData(polyData);
		#endif
		pdWriter->SetFileName(outNormals);
		pdWriter->Update();


		vtkSmartPointer<vtkPoints> pointsValues= vtkPoints::New();
  		vtkSmartPointer<vtkCellArray> cellsValues = vtkSmartPointer<vtkCellArray>::New();
		totalPoints = 0;

		numberOfSteps/=3;
		step_size*=3;
		for (unsigned j=0;j<surf->nvertices;j++)
		{
			double x,y,z, nx,ny,nz;
			double xv, yv, zv, xvp, yvp, zvp, xvn, yvn, zvn;
		        float point[3]={0,0,0};
			float normal[3]={0,0,0};
	
			vtkSmartPointer<vtkPolyLine> polyLine =    vtkSmartPointer<vtkPolyLine>::New();
			vtkIdType *ids = new vtkIdType [numberOfSteps*2+1];
			polyLine->GetPointIds()->SetNumberOfIds(numberOfSteps*2+1);

			MRISsurfaceRASToVoxel(surf, images[0], surf->vertices[j].nx,surf->vertices[j].ny,surf->vertices[j].nz, &nx,&ny,&nz);
			float dist = sqrt(nx*nx + ny*ny + nz*nz);
			if( dist>0)
			{
				nx /= dist;
				ny /= dist;
				nz /= dist;
			}	

			double max_mag=0;
			double max_val=0;
			for (int d=-1; d<2;d+=2)
			{
				for(int t=-numberOfSteps, i=0; t<=numberOfSteps;t++,i++)
				{
					x=surf->vertices[j].x +nx*t*step_size;
					y=surf->vertices[j].y +ny*t*step_size;
					z=surf->vertices[j].z +nz*t*step_size;

					MRISsurfaceRASToVoxel(surf, images[0], x,y,z, &xv,&yv,&zv);

					double mag, val;

					MRIsampleVolume(images[0], xv, yv, zv, &val);
					MRIsampleVolumeDerivativeScale(images[0], xv, yv, zv, -nx, -ny, -nz, &mag, 1.0);   // expensive
				//	std::cout << ((float)t)*step << " " << val << " " <<  mag << std::endl;
					pointsValues->InsertPoint (totalPoints,((float)t)*step_size, val, mag);
					ids[i] = totalPoints;
					totalPoints++;

					polyLine->GetPointIds()->SetId(i,ids[i]);
					polyLine->GetPointIds()->SetId(i,ids[i]);
					cellsValues->InsertNextCell(polyLine);
				}
			}
		}
                 vtkSmartPointer<vtkPolyData> polyDataValues  = vtkSmartPointer<vtkPolyData>::New();
                 polyDataValues->SetPoints(pointsValues);
                polyDataValues->SetLines(cellsValues);

		vtkSmartPointer<vtkPolyDataWriter> pdWriterValues =  vtkSmartPointer<vtkPolyDataWriter>::New();
		
		#if VTK_MAJOR_VERSION <= 5	
		pdWriterValues->SetInput(polyDataValues);
		#else
		pdWriterValues->SetInputData(polyDataValues);
		#endif	
		pdWriterValues->SetFileName(outValues);
		pdWriterValues->Update();


/*			itk::SmartPointer<TrkVTKPolPyDataFilter<ImageType>> vtk2trk  = TrkVTKPolyDataFilter<ImageType>::New();
		vtk2trk->SetReferenceImage(reader->GetOutput());		
		vtk2trk->SetInput( polyData);
		vtk2trk->VTKToTrk(outNormals);
*/
		MRISwrite(surf,outSurf);		
		MRISwriteCurvature(surf,overlayFilename) ;
		MRISfree(&surf);	


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
