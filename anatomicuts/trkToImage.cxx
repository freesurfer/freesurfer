
#include <fstream>
#include <iostream>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <cmath>
#include "itkVTKPolyDataToPolylineMeshFilter.h"
#include "itkImage.h"
#include "itkMeshToImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMesh.h"
#include "interpolation.h"

int main(int argc, char *argv[])
{
  if(argc <2)
  {
    std::cout << " parameter mesh.vtk infoImage.mhd outputImage.mhd " << std::endl;
    return EXIT_FAILURE;
  }

  enum {Dimension =3};
  typedef double                                                        PixelType;
  typedef itk::Image< PixelType, Dimension >                            ImageType;

  typedef itk::ImageFileWriter<ImageType>                               WriterType;
  typedef itk::ImageFileReader<ImageType>                               ReaderType;

  typedef itk::Mesh<double, 3> MeshType;
  typedef MeshType::Pointer MeshTypePointer;
  typedef itk::VTKPolyDataToPolylineMeshFilter<MeshType> MeshConverterType;
  typedef itk::MeshToImageFilter<MeshType,ImageType> MeshToImageType;
  MeshType::Pointer mesh;;

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    std::cout <<" reading file: "<<  argv[1] << std::endl;
    reader->SetFileName ( argv[1] );
    reader->GetOutput()->Update();

    ReaderType::Pointer readerImage = ReaderType::New();
    readerImage->SetFileName(argv[2]);
    readerImage->Update();


    MeshConverterType::Pointer converter = MeshConverterType::New();
    converter->SetVTKPolyData ( reader->GetOutput());
    reader->Delete();
    converter->Update();
    mesh = converter->GetOutput();
  
    MeshToImageType::Pointer filter = MeshToImageType::New();
    filter->SetOutputParametersFromImage(readerImage->GetOutput()); 
    filter->SetOutputSpacing(readerImage->GetOutput()->GetSpacing());
    filter->SetOutputDirection(readerImage->GetOutput()->GetDirection());
    filter->SetOutputOrigin(readerImage->GetOutput()->GetOrigin());
    filter->SetInput(mesh);
    filter->UpdateLargestPossibleRegion();

  //  std::cout << " region before write " << filter->GetOutput()->GetLargestPossibleRegion() << std::endl;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName((argv[3]));
      writer->SetInput(filter->GetOutput());
      writer->Update();
 


}
