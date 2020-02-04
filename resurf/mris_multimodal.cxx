#define COMPILING_MRISURF_TOPOLOGY
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
#include "itkSmoothingQuadEdgeMeshFilter.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
	
#include "vtkFillHolesFilter.h" 
#include "vtkPolyDataNormals.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkDecimatePro.h"
#include "vtkCleanPolyData.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTriangleFilter.h"

#include "vtkDelaunay3D.h"
#include "macros.h"
#include "mrisurf.h"
#include "mri.h"
#include "vtkKdTreePointLocator.h"
#if VTK_MAJOR_VERSION > 5	
	#include "vtkPCACurvatureEstimation.h"
#else
	#include "vtkCurvatures.h"
#endif

MRIS* VTKToSurf(vtkSmartPointer<vtkPolyData> vtkSurface) 
{
	MRIS* surf = MRISalloc( vtkSurface->GetNumberOfPoints(), vtkSurface->GetNumberOfPolys());
	surf->type = MRIS_TRIANGULAR_SURFACE;

	for(int i=0; i<vtkSurface->GetNumberOfPoints();i++)
	{	
		double* point = vtkSurface->GetPoint( i);
		double* point2 = vtkSurface->GetPoint( i);

		surf->vertices[i].x = point2[0];
		surf->vertices[i].y = point2[1];
		surf->vertices[i].z =point2[2];
		//face = &surf->faces[i];	
	}

	// Copy in the faces.
	vtkIdType cPointIDs = 0;
	vtkIdType* pPointIDs = NULL;
	vtkCellArray* polys = vtkSurface->GetPolys();
	assert( polys );
	vtkIdType nFace = 0;
	for( polys->InitTraversal();polys->GetNextCell( cPointIDs, pPointIDs ); nFace++ ) 
	{
		if( cPointIDs == 3 ) 
		{
			for( int nPointID = 0; nPointID < 3; nPointID++ )
			{	
				surf->faces[nFace].v[nPointID] = pPointIDs[nPointID];
			}
		}
	}

	return surf;
}

vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	//points->SetNumberOfPoints(surf->nvertices);
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	//triangles->SetNumberOfCells(surf->nfaces);
	std::cout << points->GetNumberOfPoints() << std::endl;
	for( int i = 0; i < surf->nvertices; i++ )
	{
		points->InsertNextPoint(surf->vertices[i].x, surf->vertices[i].y, surf->vertices[i].z);
	}

	for( int i = 0; i < surf->nfaces; i++ )
	{
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		for(int j=0;j<3;j++)
		{
			triangle->GetPointIds()->SetId(j, surf->faces[i].v[j]);
		}
		triangles->InsertNextCell(triangle);
	}	

	vtkSmartPointer<vtkPolyData> vtkSurface = vtkSmartPointer<vtkPolyData>::New();
	vtkSurface->SetPoints(points);
	vtkSurface->SetPolys(triangles);
	return vtkSurface;
}
int main(int narg, char*  arg[])
{
	constexpr unsigned int Dimension = 3;
	typedef float CoordType;
	typedef fs::Surface< CoordType, Dimension> SurfType;
	//	typedef fs::SurfaceOptimizationFilter< SurfType, SurfType> SurfFilterType;

	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -i surface -t surface -o surface -fillHoles --curvature --thickness -a anotationOutput -v overlayOutput -c output.csv -vtk "  << std::endl;   
		return -1;
	}
	const char *inSurfFilename= cl.follow ("", "-i");
	const char *targSurfFilename= cl.follow ("", "-t");
	const char *outSurfFilename = cl.follow ("", "-o");
	const char *annotationFilename = cl.follow("","-a");
	const char *overlayFilename = cl.follow("","-v");
	const char *csvFilename = cl.follow("","-c");

	MRI_SURFACE *surf;
	surf = MRISread(inSurfFilename);

	MRI_SURFACE *targetSurf;
	targetSurf = MRISread(targSurfFilename);


	SurfType::Pointer surface =  SurfType::New();
	surface->Load(&*surf);

	surf = surface->GetFSSurface(&*surf);
	if(cl.search("--fillHoles"))
	{
		vtkSmartPointer<vtkPolyData> hola = FSToVTK(surf);
		vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();

		#if VTK_MAJOR_VERSION <= 5	
		decimate->SetInput(hola);
		#else
		decimate->SetInputData(hola);
		#endif

		decimate->SetPreserveTopology(true);
		decimate->SplittingOff();
		decimate->BoundaryVertexDeletionOn();
		decimate->SetTargetReduction(.5); //99% reduction (if there was 100 triangles, now there will be 1)

		vtkSmartPointer<vtkTriangleFilter> stripper =		vtkSmartPointer<vtkTriangleFilter>::New();
		stripper->SetInputConnection( decimate->GetOutputPort() );
		stripper->PassVertsOff();
		stripper->PassLinesOff();

		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->SetInputConnection(stripper->GetOutputPort());
		cleaner->Update();

		vtkSmartPointer<vtkDelaunay3D> delaunay3D =
		vtkSmartPointer<vtkDelaunay3D>::New();
		delaunay3D->SetInputConnection (cleaner->GetOutputPort());
	
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =     vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoother->SetInputConnection( delaunay3D->GetOutputPort() );
		smoother->SetNumberOfIterations(1500);
	//	smoother->SetFeatureAngle(35);
	//	smoother->SetEdgeAngle(35);
		smoother->SetRelaxationFactor(.7);
		smoother->FeatureEdgeSmoothingOff();
		smoother->BoundarySmoothingOn();
		smoother->SetConvergence(0);
		smoother->Update();
		hola = smoother->GetOutput();
		vtkSmartPointer<vtkFillHolesFilter> fillHoles =	vtkSmartPointer<vtkFillHolesFilter>::New();
		#if VTK_MAJOR_VERSION <= 5	
		fillHoles->SetInput(hola);
		#else
		fillHoles->SetInputData(hola);
		#endif
		fillHoles->SetHoleSize(100000000.0);
		fillHoles->Update();
		surf = VTKToSurf(fillHoles->GetOutput());
	}
	if( cl.search("--smooth"))
	{

/*		itk::SmoothingQuadEdgeMeshFilter<SurfType, SurfType> smoothing =  itk::SmoothingQuadEdgeMeshFilter<SurfType, SurfType>::New();
		smoothing->SetInput(surface);
		smoothing->Update();
		surface = smoothing->GetOutput();
		surf = surface->GetFSSurface(surf);
*/
	}
	if( cl.search("--thickness"))
	{

		vtkSmartPointer<vtkPolyData> surfVTK =  FSToVTK(surf);
		vtkSmartPointer<vtkPolyData> targetVTK =  FSToVTK(targetSurf);

		vtkSmartPointer<vtkKdTreePointLocator> surfTree =	vtkSmartPointer<vtkKdTreePointLocator>::New();
		surfTree->SetDataSet(surfVTK);
		surfTree->BuildLocator();

		vtkPoints* points = vtkPoints::New();
		for (int i=0; i<targetVTK->GetNumberOfPoints(); ++i)
		{

			double* point = targetVTK->GetPoint( i);
			vtkIdType iD = surfTree->FindClosestPoint(point);
			double* point2 = surfVTK->GetPoint( iD);
			float distance =  vtkMath::Distance2BetweenPoints(point,point2);
		//	std::cout << point [0] << " " << point2[0]<< " " <<distance << std::endl;
			//if( distance > 0.01)
			{
				points->InsertPoint(i,point[0], point[1], point[2]);
			}
		}


  		vtkPolyData* polydata = vtkPolyData::New();
		polydata->SetPoints(points);

		vtkSmartPointer<vtkKdTreePointLocator> kDTree =	vtkSmartPointer<vtkKdTreePointLocator>::New();
		kDTree->SetDataSet(polydata);
//		kDTree->SetPoints(points);
		kDTree->BuildLocator();

		//		polydata->GetCellData()->GetScalars();
		//surf = VTKToSurf(polydata);
		for(int i=0; i<surfVTK->GetNumberOfPoints();i++)
		{	
			double* point = surfVTK->GetPoint( i);
			vtkIdType iD = kDTree->FindClosestPoint(point);
			
			double* point2 = targetVTK->GetPoint( iD);
			float distance =  vtkMath::Distance2BetweenPoints(point,point2);
			//kDTree->GetDataSet()->GetPoint(iD, closestPoint);

		//	std::cout << point[0] <<  " " << point2[0] << distance << std::endl;	
			//if (distance <2)
			{
				surf->vertices[i].curv=distance;
			}
	
		}
		fstream fout; 
		fout.open(csvFilename, ios::out | ios::app); 


		MRISwriteCurvature(surf,overlayFilename) ;
		//COLOR_TABLE *ct;
		//int annot;

		//ct = CTABalloc(100);
		//surf->ct = ct;
		for(int i=0; i<surfVTK->GetNumberOfPoints();i++)
		{	
			double* point = surfVTK->GetPoint( i);
			vtkIdType iD = kDTree->FindClosestPoint(point);
			
			double* point2 = targetVTK->GetPoint( iD);
			float distance =  vtkMath::Distance2BetweenPoints(point,point2);
			//kDTree->GetDataSet()->GetPoint(iD, closestPoint);

		//	std::cout << point[0] <<  " " << point2[0] << distance << std::endl;	
			//if (distance <2)
			{
				//CTABannotationAtIndex(surf->ct, int(distance*10),  &annot);
				//surf->vertices[i].annotation=annot;

				fout << i <<  ", "<< distance <<  "\n"; 

			}
	
		}
		//MRISwriteAnnotation(surf,annotationFilename) ;
		fout.close();
	}
	if( cl.search("--curvature"))
	{
		#if VTK_MAJOR_VERSION <= 5	
		vtkSmartPointer<vtkCurvatures> curvature=   vtkSmartPointer<vtkCurvatures>::New();
		curvature->SetInput(FSToVTK(surf));
		curvature->SetCurvatureTypeToGaussian();
		#else
		vtkSmartPointer<vtkPCACurvatureEstimation> curvature =   vtkPCACurvatureEstimation::New();
		//curvature->SetSampleSize(100);
		curvature->SetInputData(FSToVTK(surf));
		#endif

		//curvature->SetCurvatureTypeToMinimum();
		//curvature->SetCurvatureTypeToMaximum();
		//curvature->SetCurvatureTypeToMean();
		curvature->Update();
		vtkSmartPointer<vtkPolyData> polydata =curvature->GetOutput();

		//		polydata->GetCellData()->GetScalars();
		//surf = VTKToSurf(polydata);

		for(int i=0;i<polydata->GetNumberOfPoints();i++)
		{	
			double curv =0;
		#if VTK_MAJOR_VERSION <= 5	
			curv= 	polydata->GetPointData()->GetScalars()->GetTuple(i)[0];
		#else  			
			double* curvs = dynamic_cast<vtkDataArray*>(polydata->GetPointData()->GetArray("PCACurvature"))->GetTuple3(i);
			curv = curvs[1]/( curvs[0] + curvs[1] +curvs[2]); //,  curvs[2])*100;
		#endif
			surf->vertices[i].curv= curv;
		}	
		MRISwriteCurvature(surf,overlayFilename) ;
	
		fstream fout; 
		fout.open(csvFilename, ios::out | ios::app); 


		
		//COLOR_TABLE *ct;
		//int annot;

		//double scalarRange[2];
    		//polydata->GetScalarRange(scalarRange);
		//std::cout << scalarRange[1] << " ," <<scalarRange[0] <<std::endl;
		//ct = CTABalloc(100);
		//surf->ct = ct;
		

		for(int i=0;i<polydata->GetNumberOfPoints();i++)
		{	
				
  			double* curvs = dynamic_cast<vtkDataArray*>(polydata->GetPointData()->GetArray("PCACurvature"))->GetTuple3(i);
			//double curv = std::max(std::max( curvs[0] , curvs[1]),  curvs[2])*100;
			double curv =  curvs[0]*.5+ curvs[1]*.5; //,  curvs[2])*100;
			//CTABannotationAtIndex(surf->ct,curv,  &annot);
			//surf->vertices[i].annotation=annot;
			fout << i <<  ", "<< curv <<  "\n"; 
		}	
		//MRISwriteAnnotation(surf,annotationFilename) ;

		fout.close();
	}
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
	return EXIT_SUCCESS;
}


