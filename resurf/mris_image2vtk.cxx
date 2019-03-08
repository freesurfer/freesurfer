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
#include "itkImageToVTKImageFilter.h"
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

	const unsigned int Dimension = 3;

	using PixelType = unsigned char;
	using ImageType = itk::Image< PixelType, Dimension >;

	using ReaderType = itk::ImageFileReader< ImageType >;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputFileName );

	typedef itk::SmoothingRecursiveGaussianImageFilter<
		ImageType, ImageType >  SmoothingFilterType;

	SmoothingFilterType::Pointer smoothingRecursiveGaussianImageFilter = SmoothingFilterType::New();
	smoothingRecursiveGaussianImageFilter->SetInput(reader->GetOutput());
	smoothingRecursiveGaussianImageFilter->SetSigma(atoi(argv[6]));
	smoothingRecursiveGaussianImageFilter->Update();

	auto lowerThreshold = static_cast< PixelType >( std::stof( argv[3] ) );
	auto upperThreshold = static_cast< PixelType >( std::stof( argv[4] ) );
	auto smoothingIterations= static_cast< PixelType >( atoi( argv[5] ) );

	using FilterType2 = itk::ImageToVTKImageFilter< ImageType >;
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




	vtkSmartPointer<vtkImageThreshold> threshold2 = 
		vtkSmartPointer<vtkImageThreshold>::New();
	threshold2->SetInput( filter2->GetOutput());
	//threshold2->SetInput( gaussianSmoothFilter->GetOutput());
	threshold2->ThresholdBetween(lowerThreshold,upperThreshold);
	//threshold2->ReplaceOutOn();
	//   threshold2->SetOutValue( 0 );
	threshold2->Update();
	vtkSmartPointer<vtkMarchingCubes> contour = 
		vtkSmartPointer<vtkMarchingCubes>::New();
	contour->SetInputConnection(threshold2->GetOutputPort());
	//   contour->SetInputConnection(resampler->GetOutputPort());
	contour->SetValue(0, lowerThreshold+1);

	vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = 
		vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	conn->SetInputConnection( contour->GetOutputPort() );
	conn->SetExtractionModeToLargestRegion();
	conn->Update();
/*	 vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
	   vtkSmartPointer<vtkFillHolesFilter>::New();
	   fillHolesFilter->SetInputConnection(conn->GetOutputPort());
	   //fillHolesFilter->SetHoleSize(1.0);
	   fillHolesFilter->Update();
*/	
	vtkSmartPointer<vtkPolyDataNormals> normals =   vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection( conn->GetOutputPort() );
//	normals->SetInput( polydata );
	//normals->SetFeatureAngle( 1 );
	normals->ConsistencyOn();
	//normals->SplittingOff();
	normals->Update();

	vtkSmartPointer<vtkTriangleFilter> stripper =  vtkSmartPointer<vtkTriangleFilter>::New();
	stripper->SetInputConnection( normals->GetOutputPort() );
	//stripper->PassVertsOn();   
	stripper->Update();   
	//stripper->SetInputConnection( smoother->GetOutputPort() );
	vtkSmartPointer<vtkCleanPolyData> cleaner =  vtkSmartPointer<vtkCleanPolyData>::New();
//	cleaner->SetInputConnection(stripper->GetOutputPort());
//	cleaner->Update();

	vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother = 
		vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();

//	smoother->SetInputConnection( fillHolesFilter->GetOutputPort() );

	smoother->SetInputConnection( stripper->GetOutputPort() );
	//smoother->SetInputConnection( contour->GetOutputPort() );
	//   }
	smoother->SetNumberOfIterations(500) ; // smoothingIterations );
	smoother->FeatureEdgeSmoothingOn();
	smoother->SetEdgeAngle(30);
	smoother->Update();
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother2 = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoother2->SetInputConnection( stripper->GetOutputPort() );
	smoother2->SetNumberOfIterations(1);
	smoother2->Update();

	
/*	   vtkSmartPointer<vtkFillHolesFilter> fillHolesFilter =
	   vtkSmartPointer<vtkFillHolesFilter>::New();
	   fillHolesFilter->SetInputData(input);
	   fillHolesFilter->SetHoleSize(1000.0);
	   */
	// Write file
	vtkSmartPointer<vtkPolyDataWriter> writer =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outputFileName);
	//writer->SetInputConnection(cleaner->GetOutputPort());
	writer->SetInputConnection(conn->GetOutputPort());
	writer->Write();

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

	return EXIT_SUCCESS;
}
