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
//#include "itkVTKImageImport.h"
#include "itkMesh.h"
#include "itkImageFileWriter.h"
#include "vtkMarchingCubes.h"
#include "vtkImageThreshold.h"
//#include "itkImageToVTKImageFilter.h"
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
#include "vtkNIFTIImageReader.h"
#include "vtkDecimatePro.h"

#include "vtkCellArray.h"
extern "C" 
{
#include "mrisurf.h"
}

int main( int argc, char* argv[] )
{
	if( argc != 7 )
	{
		std::cerr << "Usage: "<< std::endl;
		std::cerr << argv[0];
		std::cerr << " <InputFileName> <OutputFileName> <Lower Threshold> <Upper Threshold> <vtk smoothing iterations> <image smoothing size>";
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	const char * inputFileName = argv[1];
	const char * outputFileName = argv[2];


	using PixelType = unsigned char;

	auto lowerThreshold = static_cast< PixelType >( std::stof( argv[3] ) );
	auto upperThreshold = static_cast< PixelType >( std::stof( argv[4] ) );
	auto smoothingIterations= static_cast< PixelType >( atoi( argv[5] ) );
	vtkSmartPointer<vtkNIFTIImageReader> reader =  vtkSmartPointer<vtkNIFTIImageReader>::New();
	reader->SetFileName(inputFileName);
	reader->Update();

	vtkSmartPointer<vtkImageThreshold> threshold2 = vtkSmartPointer<vtkImageThreshold>::New();
	threshold2->SetInputConnection( reader->GetOutputPort());
	threshold2->ThresholdBetween(lowerThreshold,upperThreshold);
	threshold2->Update();

	vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();
	contour->SetInputConnection(reader->GetOutputPort());
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

  vtkSmartPointer<vtkDecimatePro> decimate =
    vtkSmartPointer<vtkDecimatePro>::New();
  decimate->SetInput(clean->GetOutput());
  decimate->SetTargetReduction(.90); //99% reduction (if there was 100 triangles, now there will be 1)
   decimate->Update();
  //
  //
/*	vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
	smoother->SetInputConnection( stripper->GetOutputPort() );
	smoother->SetNumberOfIterations(500) ; 
//	smoother->FeatureEdgeSmoothingOn();
//	smoother->SetEdgeAngle(30);
	smoother->Update();
*/
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother2 = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother2->SetInputConnection( decimate->GetOutputPort() );
	smoother2->SetNumberOfIterations(smoothingIterations);
	smoother2->Update();


	// Write file
/*	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outputFileName);
	writer->SetInputConnection(conn->GetOutputPort());
	writer->Write();
*/
	vtkSmartPointer<vtkPolyData> vtkSurface = smoother2->GetOutput();
	MRIS* surf = MRISalloc( vtkSurface->GetNumberOfPoints(), vtkSurface->GetNumberOfPolys());
	surf->type = MRIS_TRIANGULAR_SURFACE;

	for(int i=0; i<vtkSurface->GetNumberOfPoints();i++)
	{	
		double* point = vtkSurface->GetPoint( i);
		surf->vertices[i].x = point[0];
		surf->vertices[i].y = point[1];
		surf->vertices[i].z = point[2];
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
			//surf->faces[nFace].v = pPointIDs;
			for( int nPointID = 0; nPointID < 3; nPointID++ )
			{	
				surf->faces[nFace].v[nPointID] = pPointIDs[nPointID];
				///MRIS::faces face = surf->faces[1];
				//face.v[1] = pPointIDs[1];
			}
		}
	}

	// Write the data.
	MRISwrite( surf, outputFileName);

return EXIT_SUCCESS;
}
//using ReaderType = itk::ImageFileReader< ImageType >;
//ReaderType::Pointer reader = ReaderType::New();
//reader->SetFileName( inputFileName );

/*	typedef itk::SmoothingRecursiveGaussianImageFilter<
	ImageType, ImageType >  SmoothingFilterType;

	SmoothingFilterType::Pointer smoothingRecursiveGaussianImageFilter = SmoothingFilterType::New();
	smoothingRecursiveGaussianImageFilter->SetInput(reader->GetOutput());
	smoothingRecursiveGaussianImageFilter->SetSigma(atoi(argv[6]));
	smoothingRecursiveGaussianImageFilter->Update();

*/
/*using FilterType2 = itk::ImageToVTKImageFilter< ImageType >;
  FilterType2::Pointer filter2 = FilterType2::New();
  filter2->SetInput( smoothingRecursiveGaussianImageFilter->GetOutput() );
//filter2->SetInput( reader->GetOutput() );

filter2->Update();
vtkSmartPointer<vtkImageResample> resampler = vtkSmartPointer<vtkImageResample>::New();
resampler->SetAxisMagnificationFactor(0, 2.0);
resampler->SetAxisMagnificationFactor(1, 2.0);
resampler->SetAxisMagnificationFactor(2, 2.0);
resampler->SetInput(filter2->GetOutput());
resampler->Update();

std::cout << atoi(argv[6]) << std::endl;


vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter = 
vtkSmartPointer<vtkImageGaussianSmooth>::New();
gaussianSmoothFilter->SetInputConnection(resampler->GetOutputPort());
gaussianSmoothFilter->SetRadiusFactors(std::stof(argv[6]), std::stof(argv[6]));
gaussianSmoothFilter->Update();
*/

/*
//Get the contour of the PolyData
vtkSmartPointer<vtkContourFilter> contour =
vtkSmartPointer<vtkContourFilter>::New();
contour->SetInput(threshold->GetOutput());
contour->ComputeScalarsOff();
contour->SetValue(0,0.5);
//    
//    //Mesh the largest contour
vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter =
vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
connectivityFilter->SetInput(contour->GetOutput());
connectivityFilter->SetExtractionModeToLargestRegion();
connectivityFilter->Update();
//        
//Export vtk Data to ITK
vtkImageExport* exportdata = vtkImageExport::New();
exportdata->SetInput(connectivityFilter->GetOutput());
//            
const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Mesh< PixelType, Dimension > MeshType;
typedef itk::ImageFileWriter<MeshType> WriteMeshType;
typedef itk::VTKImageImport<MeshType> importMeshType;
//                
WriteMeshType::Pointer WMeshType = WriteMeshType::New();
importMeshType::Pointer IMeshType = importMeshType::New();
exportdata->SetOutput(IMeshType);
IMeshType->Update();
WMeshType->SetInput(IMeshType->GetOutput());
*/
//typename Surface::Pointer mesh = Self::New();
/*	   vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
	   vtkSmartPointer<vtkFillHolesFilter>::New();
	   fillHolesFilter->SetInputData(input);
	   fillHolesFilter->SetHoleSize(1000.0);
	   */

