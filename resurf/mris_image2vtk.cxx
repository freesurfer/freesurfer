#define COMPILING_MRISURF_TOPOLOGY
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkSmartPointer.h"

#include "itkMesh.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMeshFileWriter.h"

#include "vtkContourFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkImageExport.h"
#include "itkMesh.h"
#include "itkImageFileWriter.h"
#include "vtkMarchingCubes.h"
#include "vtkImageThreshold.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkPolyDataWriter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "vtkImageResample.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkFillHolesFilter.h"
//#include "vtkNIFTIImageReader.h"
#include "vtkDecimatePro.h"

#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "macros.h"
#include "vtkImageData.h"
#include "mri.h"
#include "mrisurf.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"


#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "cma.h"
#include "version.h"
#include "vtkDelaunay3D.h"
#include "vtkFillHolesFilter.h"
#include "vtkImageDilateErode3D.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageContinuousErode3D.h"
#include "vtkImageContinuousDilate3D.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkVoxelContoursToSurfaceFilter.h"
#include "vtkQuadricDecimation.h"
vtkSmartPointer<vtkImageData> CreateImage( MRI* rasMRI )
{
	// first copy mri data to image
	vtkDataArray *scalars = NULL;
	vtkUnsignedCharArray  *ucharScalars = NULL;
	vtkIntArray           *intScalars = NULL;
	vtkShortArray         *shortScalars = NULL;
	vtkLongArray          *longScalars = NULL;
	vtkFloatArray         *floatScalars = NULL;
	vtkIdType cValues;

	vtkIdType zX = rasMRI->width;
	vtkIdType zY = rasMRI->height;
	vtkIdType zZ = rasMRI->depth;
	vtkIdType zFrames = rasMRI->nframes;

	vtkSmartPointer<vtkImageData> imageData =vtkSmartPointer<vtkImageData>::New();

	// This object's output space is in voxel coordinates.
	std::cout << zX<< zY<<zZ << std::endl;
	imageData->SetDimensions( zX, zY, zZ );

	double origin[3] = {0,0,0};

	imageData->SetSpacing( rasMRI->xsize, rasMRI->ysize, rasMRI->zsize);
	imageData->SetOrigin( rasMRI->xstart, rasMRI->ystart, rasMRI->zstart );
	//  imageData->SetWholeExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );
	imageData->SetDimensions(zX, zY, zZ);
	if (rasMRI->type == MRI_RGB)
		zFrames = 4;

	// create the scalars for all of the images. set the element size
	// for the data we will read.
#if VTK_MAJOR_VERSION > 5
	switch ( rasMRI->type )
	{
		case MRI_UCHAR:
		case MRI_RGB:
			imageData->AllocateScalars(VTK_UNSIGNED_CHAR, zFrames);
			break;
		case MRI_INT:
			imageData->AllocateScalars(VTK_INT, zFrames);
			break;
		case MRI_LONG:
			imageData->AllocateScalars(VTK_LONG, zFrames);
			break;
		case MRI_FLOAT:
			imageData->AllocateScalars(VTK_FLOAT, zFrames);
			break;
		case MRI_SHORT:
			imageData->AllocateScalars(VTK_SHORT, zFrames);
			break;
		default:
			;  
	}
#else
	imageData->SetNumberOfScalarComponents(zFrames);
	switch ( rasMRI->type )
	{
		case MRI_UCHAR:
		case MRI_RGB:
			imageData->SetScalarTypeToUnsignedChar();
			break;
		case MRI_INT:
			imageData->SetScalarTypeToInt();
			break;
		case MRI_LONG:
			imageData->SetScalarTypeToLong();
			break;
		case MRI_FLOAT:
			imageData->SetScalarTypeToFloat();
			break;
		case MRI_SHORT:
			imageData->SetScalarTypeToShort();
			break;
		default:
			;  
	}
	imageData->AllocateScalars();
#endif

	return imageData;
}
void CopyMRIDataToImage( MRI* mri, vtkImageData* image )
{
	// Copy the slice data into the scalars.
	int zX = mri->width;
	int zY = mri->height;
	int zZ = mri->depth;
	int zFrames = mri->nframes;

	vtkIdType nTuple = 0;
	vtkDataArray *scalars = image->GetPointData()->GetScalars();
	for ( int nZ = 0; nZ < zZ; nZ++ )
	{
		for ( int nY = 0; nY < zY; nY++ )
		{
			for ( int nX = 0; nX < zX; nX++ )
			{
				for ( int nFrame = 0; nFrame < zFrames; nFrame++ )
				{
					//int hola = mri->slices[nZ+(nFrame)*mri->depth][nY][nX];
					//scalars->SetComponent( nTuple, nFrame,hola);

					switch ( mri->type )
					{
						case MRI_UCHAR:
							scalars->SetComponent( nTuple, nFrame,
									MRIseq_vox( mri, nX, nY, nZ, nFrame ) );
							break;
						case MRI_INT:
							scalars->SetComponent( nTuple, nFrame,
									MRIIseq_vox( mri, nX, nY, nZ, nFrame ) );
							break;
						case MRI_LONG:
							scalars->SetComponent( nTuple, nFrame,
									MRILseq_vox( mri, nX, nY, nZ, nFrame ) );
							break;
						case MRI_FLOAT:
							scalars->SetComponent( nTuple, nFrame,
									MRIFseq_vox( mri, nX, nY, nZ, nFrame ) );
							break;
						case MRI_SHORT:
							scalars->SetComponent( nTuple, nFrame,
									MRISseq_vox( mri, nX, nY, nZ, nFrame ) );
							break;
						default:
							break;
					}

				}
				if (mri->type == MRI_RGB)
				{
					int val = MRIIseq_vox(mri, nX, nY, nZ, 0);
					scalars->SetComponent( nTuple, 0, val & 0x00ff);
					scalars->SetComponent( nTuple, 1, (val >> 8) & 0x00ff);
					scalars->SetComponent( nTuple, 2, (val >> 16) & 0x00ff);
					scalars->SetComponent( nTuple, 3, 255);
				}
				nTuple++;
			}
		}

	}
}
int main( int argc, char* argv[] )
{
	if( argc != 8 )
	{
		std::cerr << "Usage: "<< std::endl;
		std::cerr << argv[0];
		std::cerr << " <InputFileName> <OutputFileName> <Lower Threshold> <Upper Threshold> <vtk smoothing iterations> <image smoothing size> <reduction percentage>";
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	const char * inputFileName = argv[1];
	const char * outputFileName = argv[2];

	MRI *imageFS =  MRIread(inputFileName) ;
	int dirx = imageFS->xstart - imageFS->xend;
	int diry = imageFS->ystart - imageFS->yend;
	int dirz = imageFS->zstart - imageFS->zend;

	dirx/=abs(dirx);
	diry/=abs(diry);
	dirz/=abs(dirz);
	std::cout <<dirx << " " << dirz << " "<< diry << std::endl;

	vtkSmartPointer<vtkImageData> vtkImage = CreateImage(imageFS);
	CopyMRIDataToImage(imageFS, vtkImage);	

	using PixelType = unsigned char;

	auto label = static_cast< PixelType >( std::stof( argv[3] ) );
	auto upperThreshold = static_cast< PixelType >( std::stof( argv[4] ) );
	auto vtkSmoothingIterations= static_cast< PixelType >( atoi( argv[5] ) );
	auto imageSmoothingIterations= static_cast< PixelType >( atoi( argv[6] ) );
	float redPercentage =  std::stof( argv[7] );
	std::cout << redPercentage << std::endl;	
	//vtkSmartPointer<vtkNIFTIImageReader> reader =  vtkSmartPointer<vtkNIFTIImageReader>::New();
	//reader->SetFileName(inputFileName);
	//reader->Update();

	vtkSmartPointer<vtkImageDilateErode3D> dilateErode =
		vtkSmartPointer<vtkImageDilateErode3D>::New();
#if VTK_MAJOR_VERSION <= 5	
	dilateErode->SetInput(vtkImage);
#else
	dilateErode->SetInputData(vtkImage);
	#endif
	dilateErode->SetDilateValue(label);
	dilateErode->SetErodeValue(label);
	dilateErode->SetKernelSize(5, 5, 5);
	dilateErode->ReleaseDataFlagOff();
	dilateErode->Update();

	vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
	#if VTK_MAJOR_VERSION <= 5	
	threshold->SetInput( vtkImage); //dilateErode->GetOutput());
	#else
	threshold->SetInputData( vtkImage) ;// dilateErode->GetOutput());
	#endif
	threshold->ThresholdBetween( label, upperThreshold );
	threshold->SetReplaceIn(true);
	threshold->SetInValue( 1);
	threshold->SetReplaceOut(true);
	threshold->SetOutValue( 0 );
	threshold->SetOutputScalarTypeToInt();
	threshold->Update();
	vtkImage= threshold->GetOutput();
	for(int i=0;i<imageSmoothingIterations;i++)
	{
		vtkSmartPointer<vtkImageContinuousErode3D> erode3D = vtkSmartPointer<vtkImageContinuousErode3D>::New();
		#if VTK_MAJOR_VERSION <= 5	
	erode3D->SetInput(vtkImage);
	#else
	erode3D->SetInputData(vtkImage);
	#endif

		erode3D->SetKernelSize(2, 2, 2);

		vtkSmartPointer<vtkImageContinuousDilate3D> dilate3D = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
		dilate3D->SetInputConnection(erode3D->GetOutputPort());
		dilate3D->SetKernelSize(2, 2, 2);
		dilate3D->Update();
		vtkImage = dilate3D->GetOutput();
		
		vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter =    vtkSmartPointer<vtkImageGaussianSmooth>::New();
		#if VTK_MAJOR_VERSION <= 5	
		gaussianSmoothFilter->SetInput(vtkImage);
	#else
		gaussianSmoothFilter->SetInputData(vtkImage);
	#endif
		gaussianSmoothFilter->SetRadiusFactors(2*imageFS->xsize,2*imageFS->ysize,2*imageFS->zsize);
		gaussianSmoothFilter->Update();
		vtkImage = gaussianSmoothFilter->GetOutput();
	}



	vtkSmartPointer<vtkPolyData> vtkSurface; 

//	vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
	vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();

//	contour->ComputeGradientsOn();
//	contour->ComputeScalarsOn();
#if VTK_MAJOR_VERSION <= 5	
	contour->SetInput( vtkImage);
#else
	contour->SetInputData( vtkImage);
#endif
	contour->SetValue(0, 1);	
	contour->SetNumberOfContours(1);
	contour->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> conn =	vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	conn->SetInputConnection( contour->GetOutputPort() );
	//conn->SetExtractionModeToLargestRegion();
	conn->SetExtractionModeToAllRegions();

	vtkSurface = conn->GetOutput();

	std::cout << " num points " << vtkSurface->GetNumberOfPoints()<< std::endl;

//	for(int i =0;i<vtkSmoothingIterations;i++)
	{
		vtkSmartPointer<vtkFillHolesFilter> fillHoles =	vtkSmartPointer<vtkFillHolesFilter>::New();
	#if VTK_MAJOR_VERSION <= 5	
		fillHoles->SetInput(vtkSurface);
	#else
		fillHoles->SetInputData(vtkSurface);
	#endif
		fillHoles->SetHoleSize(10000000000.0);
		fillHoles->Update();
		vtkSurface = fillHoles->GetOutput();

		vtkSmartPointer<vtkTriangleFilter> stripper2 =  vtkSmartPointer<vtkTriangleFilter>::New();
	#if VTK_MAJOR_VERSION <= 5	
		stripper2->SetInput( vtkSurface );
	#else
		stripper2->SetInputData( vtkSurface );
	#endif
		stripper2->Update();   	
		vtkSurface = stripper2->GetOutput();

		float reduction = 1.0 - redPercentage /vtkSurface->GetNumberOfPoints();
		std::cout << " reduction " << reduction << " num points " << vtkSurface->GetNumberOfPoints()<< std::endl;
		vtkSmartPointer<vtkQuadricDecimation> decimate = vtkSmartPointer<vtkQuadricDecimation>::New();
		//vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	#if VTK_MAJOR_VERSION <= 5	
		decimate->SetInput(vtkSurface);
	#else
		decimate->SetInputData(vtkSurface);
	#endif
	//	decimate->SetVolumePreservation(true);
	//	decimate->SetPreserveTopology(true);
	//	decimate->SplittingOff();
	//	decimate->BoundaryVertexDeletionOn();
		decimate->SetTargetReduction(reduction); //99% reduction (if there was 100 triangles, now there will be 1)
		decimate->Update();
		vtkSurface = decimate->GetOutput();

		std::cout << "num points " << vtkSurface->GetNumberOfPoints()<< std::endl;
		vtkSmartPointer<vtkTriangleFilter> stripper =		vtkSmartPointer<vtkTriangleFilter>::New();
	#if VTK_MAJOR_VERSION <= 5	
		stripper->SetInput( vtkSurface );
	#else
		stripper->SetInputData( vtkSurface );
	#endif
	//	stripper->PassVertsOff();
	//	stripper->PassLinesOff();

		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->SetAbsoluteTolerance(0.1);
		cleaner->SetInputConnection(stripper->GetOutputPort());
		cleaner->Update();

		vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =     vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		//smoother->SetInputConnection( contour->GetOutputPort() );
		smoother->SetInputConnection( cleaner->GetOutputPort() );
		smoother->SetNumberOfIterations(vtkSmoothingIterations);
		smoother->SetRelaxationFactor(0.3);
		smoother->FeatureEdgeSmoothingOff();
		smoother->BoundarySmoothingOn();
		smoother->Update();
		vtkSurface = smoother->GetOutput();

	}


	
	MRIS* surf = MRISalloc( vtkSurface->GetNumberOfPoints(), vtkSurface->GetNumberOfPolys());
	surf->type = MRIS_TRIANGULAR_SURFACE;
	for(int i=0; i<vtkSurface->GetNumberOfPoints();i++)
	{	
		double* point = vtkSurface->GetPoint( i);
		double* point2 = vtkSurface->GetPoint( i);
		//::MRIvoxelToWorld( imageFS, point[0]/0.013, point[1]/0.1, point[2]/0.013, &point2[0], &point2[1], &point2[2] );
		//::MRIworldToVoxel( m_volumeRef->m_MRITarget, ras[0], ras[1], ras[2], &cindex[0], &cindex[1], &cindex[2] );

		//	MRISsurfaceRASToVoxel(surf, imageFS,point[0], point[1], point[2], &point2[0], &point2[1], &point2[2] );
		
		// MRIvoxelToSurfaceRAS( imageFS,point[0], point[1], point[2], &point2[0], &point2[1], &point2[2] );

		surf->vertices[i].x = -  point2[0];
		surf->vertices[i].z = -  point2[1];
		surf->vertices[i].y =  point2[2];
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
			//surf->faces[nFace].v = pPointIDs;
			for( int nPointID = 0; nPointID < 3; nPointID++ )
			{	
				surf->faces[nFace].v[nPointID] = pPointIDs[nPointID];
				///MRIS::faces face = surf->faces[1];
				//face.v[1] = pPointIDs[1];
			}
		}
	}
	//surf->SetMatrix(MRIvoxelXformToRasXform(imageFS));
	
	//surf->vg= imageFS->vg; //=  MRIgetVoxelToRasXform(imageFS) ; //MRIvoxelXformToRasXform(imageFS); //)imageFS->GetMatrix();
	// Write the data.
	MRISwrite( surf, outputFileName);

return EXIT_SUCCESS;
}



/*	vtkSmartPointer<vtkImageResample> resampler = vtkSmartPointer<vtkImageResample>::New();
	resampler->SetAxisMagnificationFactor(0, 2.0);
	resampler->SetAxisMagnificationFactor(1, 2.0);
	resampler->SetAxisMagnificationFactor(2, 2.0);
	resampler->SetInputConnection(threshold->GetOutputPort());
*/	

/*




		vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter =    vtkSmartPointer<vtkImageGaussianSmooth>::New();
		gaussianSmoothFilter->SetInput(vtkImage);
		gaussianSmoothFilter->SetRadiusFactors(.2,.2,.2);
		gaussianSmoothFilter->Update();
		vtkImage = gaussianSmoothFilter->GetOutput();
		//vtkSmartPointer<vtkMarchingCubes> contour =vtkSmartPointer<vtkMarchingCubes>::New(); 
		
		
		   vtkSmartPointer<vtkFlyingEdges3D> isosurface  = vtk.vtkFlyingEdges3D();
		   isosurface->SetInputConnection(conn->GetOutputPort());
		   isosurface->SetValue(label-0.5, label+0.5);
		   isosurface->ComputeNormalsOn();
		   isosurface->ComputeGradientsOn();
		   isosurface->ComputeScalarsOn();
		   isosurface->InterpolateAttributesOff();
		   isosurface->Update();
		   vtkSurface = isosurface->GetOutput();


	vtkSmartPointer<vtkImageThreshold> threshold2 = vtkSmartPointer<vtkImageThreshold>::New();
	threshold2->SetInput( vtkImage);
	threshold2->ThresholdBetween(lowerThreshold,upperThreshold);
	threshold2->Update();

	vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();
	contour->SetInput(vtkImage);
	//contour->SetInputConnection(threshold2->GetOutputPort());
	contour->SetValue(0, lowerThreshold+1);

	vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	conn->SetInputConnection( contour->GetOutputPort() );
	conn->SetExtractionModeToLargestRegion();
	conn->Update();

	vtkSmartPointer<vtkPolyDataNormals> normals =   vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection( conn->GetOutputPort() );
	normals->ConsistencyOn();
	normals->Update();

	vtkSmartPointer<vtkTriangleFilter> stripper =  vtkSmartPointer<vtkTriangleFilter>::New();
	stripper->SetInputConnection( normals->GetOutputPort() );
	stripper->Update();   

	vtkSmartPointer<vtkCleanPolyData> clean=  vtkSmartPointer<vtkCleanPolyData>::New();
	clean->SetInputConnection( stripper->GetOutputPort() );
	clean->Update();   

	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	decimate->SetInput(contour->GetOutput());
	//decimate->SetPreserveTopology(true);
	//decimate->SplittingOff();
	decimate->SetTargetReduction(.95); //99% reduction (if there was 100 triangles, now there will be 1)
	decimate->Update();
  //
  //
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother2 = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother2->SetInputConnection( decimate->GetOutputPort() );
	smoother2->SetFeatureEdgeSmoothing(true);
	smoother2->SetFeatureAngle(.30);
	smoother2->SetNumberOfIterations(100);//smoothingIterations);
	smoother2->Update();
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =     vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputConnection(decimate->GetOutputPort());
//	smoothFilter->SetNumberOfIterations(10);
//	smoothFilter->SetFeatureAngle(15);
//	smoothFilter->SetEdgeAngle(15);
//	smoothFilter->SetRelaxationFactor(.9);
//	smoothFilter->FeatureEdgeSmoothingOn();
//	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();

	// Update normals on newly smoothed polydata
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
	normalGenerator->ComputePointNormalsOn();
	normalGenerator->ComputeCellNormalsOn();
	normalGenerator->Update();

	vtkSmartPointer<vtkPolyData> vtkSurface = smoothFilter->GetOutput();
	//vtkSmartPointer<vtkPolyData> vtkSurface = normalGenerator->GetOutput();
	//vtkSmartPointer<vtkPolyData> vtkSurface = smoother2->GetOutput();
	*/
/*vtkSmartPointer<vtkVoxelContoursToSurfaceFilter> contoursToSurface = vtkSmartPointer<vtkVoxelContoursToSurfaceFilter>::New();
#if VTK_MAJOR_VERSION <= 5
		contoursToSurface->SetInput(vtkSurface );
#else
		contoursToSurface->SetInputData( vtkSurface );
#endif
		contoursToSurface->Update();

		vtkSurface =  contoursToSurface->GetOutput();
		*/

		// polygonal data --> image stencil:
/*		vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =	vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	#if VTK_MAJOR_VERSION <= 5	
		pol2stenc->SetInput(vtkSurface);
	#else
		pol2stenc->SetInputData(vtkSurface);
	#endif
		pol2stenc->SetOutputOrigin(vtkImage->GetOrigin());
		pol2stenc->SetOutputSpacing(vtkImage->GetSpacing());
		pol2stenc->SetOutputWholeExtent(vtkImage->GetExtent());
		pol2stenc->Update();

		// cut the corresponding white image and set the background:
		vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
	#if VTK_MAJOR_VERSION <= 5	
		imgstenc->SetInput(vtkImage);
	imgstenc->SetStencil(pol2stenc->GetOutput());
	#else
		imgstenc->SetInputData(vtkImage);
	imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
	#endif
		imgstenc->ReverseStencilOff();
		imgstenc->SetBackgroundValue(0);
		imgstenc->Update();
		vtkImage =  imgstenc->GetOutput();
		
		vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter =    vtkSmartPointer<vtkImageGaussianSmooth>::New();
	#if VTK_MAJOR_VERSION <= 5	
		gaussianSmoothFilter->SetInput(vtkImage);
	#else
		gaussianSmoothFilter->SetInputData(vtkImage);
	#endif
		gaussianSmoothFilter->SetRadiusFactors(.2,.2,.2);
		gaussianSmoothFilter->Update();
		vtkImage = gaussianSmoothFilter->GetOutput();
*/
