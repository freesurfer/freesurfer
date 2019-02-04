#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include <iostream>
#include <string>         
#include "itkImageFileReader.h"
#include "GetPot.h"
#include "itkImage.h"
#include "itkMesh.h"
#include "TrkVTKPolyDataFilter.txx"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"
#include "vtkPolyDataWriter.h"
#include "itkImageFileWriter.h"
#include "MeshToImageFilter.h"
#include "mri.h"


int main(int narg, char* arg[])
{
	enum {Dimension =3};
	typedef int                                                        PixelType;
	typedef itk::Image< PixelType,Dimension> ImageType;
	// typedef ImageType::IndexType 			IndexType;
	typedef itk::Mesh< PixelType, Dimension > MeshType;
	typedef itk::ImageFileWriter<ImageType>                               WriterType;
	typedef itk::ImageFileReader<ImageType>                               ReaderType;
	typedef MeshToImageFilter<MeshType,ImageType> MeshToImageType;
	typedef VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
	typedef itk::ImageFileReader<ImageType> ImageReaderType;
	typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;


	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -i referenceImage -f input.trk[options]"  << std::endl;
		std::cout<< " -o output.trk: output trk"  << std::endl;
		std::cout<< " -e output.nii.gz: export trk into an image"  << std::endl;
		std::cout<< " -u :update trk header with reference image"  << std::endl;
		std::cout<< " -v output.vtk: outputs streamlines in vtk format "  << std::endl;
		return -1;

	}
	const char *imageFile = cl.follow ("", "-i");
	const char *fiberFile = cl.follow ("", "-f");
	const char *output = cl.follow ("", "-o");
	const char *outputImage = cl.follow ("", "-e");
	//MRI *outref = 0;
	//MATRIX *outv2r;

	//outref = MRIread(imageFile);

	// Output space orientation information
	//outv2r = MRIgetVoxelToRasXform(outref);



	MeshConverterType::Pointer converter = MeshConverterType::New();

	ImageReaderType::Pointer reader = ImageReaderType::New();
	reader->SetFileName ( imageFile);
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();

	itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
	trkReader->SetTrkFileName(fiberFile);
	trkReader->SetReferenceImage(image);
	trkReader->TrkToVTK();
	converter->SetVTKPolyData ( trkReader->GetOutputPolyData() );

	converter->Update();

	MeshType::Pointer mesh = converter->GetOutput();

	if (cl.search("-u"))
	{

		VTKConverterType::Pointer vtkConverter =  VTKConverterType::New();
		vtkConverter->SetInput(mesh);
		vtkConverter->Update();

		if(cl.search("-v"))
		{
			vtkSmartPointer<vtkPolyDataWriter> writerFixed = vtkPolyDataWriter::New();
			writerFixed->SetFileName ( cl.follow("","-v"));
			writerFixed->SetInput(vtkConverter->GetOutputPolyData());
			writerFixed->SetFileTypeToBinary();
			writerFixed->Update();
		}

		trkReader->SetInput(vtkConverter->GetOutputPolyData());
		//	trkReader->SetReferenceTrack(fiberFile);
		trkReader->SetReferenceImage(image);
		trkReader->VTKToTrk(output);
		std::cout << " Saving trk with updated header  " << output << std::endl;

	}
	else if(cl.search("-e"))
	{

		MeshToImageType::Pointer filter = MeshToImageType::New();
		filter->SetOutputParametersFromImage(image); 
		filter->SetOutputSpacing(image->GetSpacing());
		filter->SetOutputDirection(image->GetDirection());
		filter->SetOutputOrigin(image->GetOrigin());
		filter->SetInput(mesh);
		filter->UpdateLargestPossibleRegion();

		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(outputImage);
		writer->SetInput(filter->GetOutput());
		writer->Update();

	}
}
