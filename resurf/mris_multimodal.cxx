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


		//MRISwriteCurvature(surf,overlayFilename) ;
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
		curvature->SetSampleSize(500);
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
			curv =  curvs[0]*50+ curvs[1]*50; //,  curvs[2])*100;
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

/*
 *
 *vtkSmartPointer<vtkPolyDataNormals> normals =   vtkSmartPointer<vtkPolyDataNormals>::New();
		normals->SetInput( FSToVTK(surf) );
		//	normals->SetInput( polydata );
		normals->SetFeatureAngle( 60.0 );
		normals->ComputePointNormalsOff();
		normals->ComputeCellNormalsOn();
		normals->ConsistencyOn();
		normals->SplittingOff();
		normals->Update();


		vtkSmartPointer<vtkPolyData> polydata = normals->GetOutput();

		vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works
		std::cout << "There are " << normalsGeneric->GetNumberOfTuples() << std::endl;
		std::cout << "There are " <<polydata->GetNumberOfCells() << std::endl;
		int num = normalsGeneric->GetNumberOfTuples() ;
		vtkSmartPointer<vtkFloatArray> colors= vtkSmartPointer<vtkFloatArray>::New();
		std::cout << image->GetLargestPossibleRegion()<< std::endl;
		for(int i=0;i<num;i++)
		{
			double normal[3];
			normalsGeneric->GetTuple(i, normal);

			double point[3];
			//polydata->GetPoint(i, point);
			polydata->GetCell(i)->GetPoints()->GetPoint(1, point);

			ImageType::PointType point1,point2, point0;
			ImageType::IndexType index;
			int dir=1;
			for(int j=0;j<3;j++)
			{
				normal[j] =normal[j]*spacing[0]/2;
				point1[j] = -point[j];			
				point2[j] = -point[j];			
				point0[j] = -point[j];			

			}
			point1[2]=-point1[2];
			point2[2]=-point2[2];
			point0[2]=-point0[2];
			image->TransformPhysicalPointToIndex(point1, index);
			float thickness = 0;

			float thicknessPos = 0, thicknessNeg =0;
			int currentLabel = image->GetPixel(index);
			while(currentLabel != label && currentLabel != 0 && currentLabel!= nolabel)
			{
				for(int j=0;j<3;j++)
				{
					point1[j] = point1[j] +  normal[j];			
				}		
				if (image->TransformPhysicalPointToIndex(point1, index))
				{ //std::cout << index << std::endl;

					currentLabel = image->GetPixel(index) ;
					thicknessPos ++;
				}
				else
				{
					break;
				}
			}
			thicknessPos = point1.EuclideanDistanceTo(point0);
			if( currentLabel == label )
				thickness = thicknessPos;
			std::cout << currentLabel <<  " " << label  << " " << thicknessPos ;
			image->TransformPhysicalPointToIndex(point2, index);
			currentLabel = image->GetPixel(index);
			while(currentLabel != label && currentLabel!=0 && currentLabel != nolabel)
			{
				for(int j=0;j<3;j++)
				{
					point2[j] = point2[j] - normal[j];			
				}		
				if(image->TransformPhysicalPointToIndex(point2, index))
				{
					currentLabel = image->GetPixel(index) ;
					thicknessNeg ++;
				}
				else
				{
					break;
				}
			}
			thicknessNeg = point2.EuclideanDistanceTo(point0);
			std::cout << " " <<thicknessNeg << " " << currentLabel  << std::endl;
			if( currentLabel == label )
				thickness = thicknessNeg;
			//if(thickness == 0)
			{
				if( thicknessNeg < 2.5 && thicknessNeg>.1)
					thickness =  thicknessNeg;
				else if( thicknessPos < 2.5 && thicknessPos>.1)
					thickness = thicknessPos;

			}	
			colors->InsertNextValue(thickness);

		}
		polydata->GetCellData()->SetScalars(colors);


 */	/*if( cl.search("--turvatures"))
	{

		vtkSmartPointer<vtkPolyData> surfVTK =  FSToVTK(surf);

		vtkSmartPointer<vtkKdTreePointLocator> surfTree =	vtkSmartPointer<vtkKdTreePointLocator>::New();
		surfTree->SetDataSet(surfVTK);
		surfTree->BuildLocator();


		double *a[3], a0[3], a1[3], a2[3], xp[3];
		a[0] = a0; a[1] = a1; a[2] = a2;

		for(int i=0; i<surfVTK->GetNumberOfPoints();i++)
		{	
			double* point = surfVTK->GetPoint( i);
	
      			vtkIdList*& pIds;
			surfTree->FindClosestNPoints(20, point, pIds);
			int numPts = pIds->GetNumberOfIds();

			// First step: compute the mean position of the neighborhood.
			double mean[3];
			mean[0] = mean[1] = mean[2] = 0.0;
			for (int sample=0; sample<numPts; ++sample)
			{
				double* point2 = surfVTK->GetPoint( sample);
				for(int j=0;j<3;j++)
					mean[j] += sample[j];
			}
			for(int j=0;j<3;j++)
				mean[j] /= numPts;

			// Now compute the covariance matrix
			a0[0] = a1[0] = a2[0] = 0.0;
			a0[1] = a1[1] = a2[1] = 0.0;
			a0[2] = a1[2] = a2[2] = 0.0;
			for (int sample=0; sample < numPts; ++sample )
			{
				double* point2 = surfVTK->GetPoint( sample);
				for(int j=0;j<3;j++)
					xp[j] = sample[j] - mean[j];
				for (i=0; i < 3; i++)
				{
					a0[i] += xp[0] * xp[i];
					a1[i] += xp[1] * xp[i];
					a2[i] += xp[2] * xp[i];
				}
			}
			for (i=0; i < 3; i++)
			{
				a0[i] /= numPts;
				a1[i] /= numPts;
				a2[i] /= numPts;
			}

			// Next extract the eigenvectors and values
			vtkMath::Jacobi(a,eVal,v);

			// Finally compute the curvatures
			double den = eVal[0] + eVal[1] + eVal[2];
			c[0] = (eVal[0] - eVal[1]) / den;
			c[1] = 2.0*(eVal[1] - eVal[2]) / den;
			c[3] = 3.0*eVal[2] / den;

			

	}*/

