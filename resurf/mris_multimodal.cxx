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
#include "fsSurfaceOptimizationFilter.h"
#include "itkVTKPolyDataWriter.h"
#include "itkSmoothingQuadEdgeMeshFilter.h"
	
#include "vtkFillHolesFilter.h" 
#include "vtkPolyDataNormals.h"
#include "vtkCurvatures.h"
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
	typedef fs::SurfaceOptimizationFilter< SurfType, SurfType> SurfFilterType;

	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -i surface -o surface -fillHoles -vtk "  << std::endl;   
		return -1;
	}
	const char *inSurfFilename= cl.follow ("", "-i");
	const char *outSurfFilename = cl.follow ("", "-o");


	MRI_SURFACE *surf;
	surf = MRISread(inSurfFilename);

	SurfType::Pointer surface =  SurfType::New();
	surface->Load(&*surf);

	surf = surface->GetFSSurface(&*surf);
	if(cl.search("--fillHoles"))
	{
		vtkSmartPointer<vtkPolyData> hola = FSToVTK(surf);
		vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
		decimate->SetInput(hola);
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
		fillHoles->SetInput(hola);
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
/*		vtkSmartPointer<vtkPolyDataNormals> normals =   vtkSmartPointer<vtkPolyDataNormals>::New();
		normals->SetInputData( surface );
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
			input->GetCell(i)->GetPoints()->GetPoint(1, point);

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
		input->GetCellData()->SetScalars(colors);

		vtkSmartPointer<vtkPolyDataWriter> pdWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		pdWriter->SetFileName(outputFileName);
		pdWriter->SetInputData(input);
		pdWriter->Update();

*/	}
	if( cl.search("--curvature"))
	{
		vtkSmartPointer<vtkCurvatures> curvature=   vtkSmartPointer<vtkCurvatures>::New();
		curvature->SetInput(FSToVTK(surf));
		//curvature->SetCurvatureTypeToMinimum();
		//curvature->SetCurvatureTypeToMaximum();
		//curvature->SetCurvatureTypeToMean();
		curvature->SetCurvatureTypeToGaussian();
		curvature->Update();
		vtkSmartPointer<vtkPolyData> polydata =curvature->GetOutput();

		//input->GetCellData()->SetScalars(colors);

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
